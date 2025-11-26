#!/usr/bin/env Rscript
# Title: Differential-expression workflow – Accelerated vs Decelerated (Horvath clock)
# Description: Automates Horvath-style differential expression by pairing cohort RNA-seq
#   counts with clock predictions, enforcing instructor-provided Accelerated/Decelerated
#   labels when available, and fitting DESeq2 models for each cohort. Writes annotated
#   result tables plus PCA and volcano plots per tissue. Designed to run standalone for
#   BRCA/KIRC/LIHC by default but can accept custom count/prediction inputs via CLI.

suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(ggplot2)
  library(SummarizedExperiment)
})

#======================================================================
# Helper functions
#======================================================================

ensure_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing ", pkg, " ...")
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
  library(pkg, character.only = TRUE)
}

invisible(lapply(c("data.table", "DESeq2", "ggplot2", "SummarizedExperiment"), ensure_pkg))

shorten_barcode <- function(x) substr(x, 1, 16)

sanitize_sample <- function(x) {
  x <- gsub("-", ".", x, fixed = TRUE)
  gsub("N\\.DNA", "N.RNA", x)
}

load_labels <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) return(NULL)
  dt <- tryCatch(fread(path), error = function(e) {
    warning(sprintf("Failed to read classification sheet %s: %s", path, e$message))
    return(NULL)
  })
  if (is.null(dt)) return(NULL)
  required <- c("sample_id", "classification")
  missing <- setdiff(required, names(dt))
  if (length(missing)) {
    warning(sprintf("Classification sheet %s missing columns: %s", path, paste(missing, collapse = ", ")))
    return(NULL)
  }
  dt[, sample_id := shorten_barcode(sample_id)]
  dt[, classification := factor(classification, levels = c("Accelerated", "Decelerated"))]
  if (all(is.na(dt$classification))) return(NULL)
  dt
}

#======================================================================
# Default parameters for the three main cohorts
#======================================================================

default_params <- list(
  BRCA = list(
    project     = "TCGA-BRCA",
    counts_path = "results/differential_analysis/TCGA-BRCA_tissue_quartiles_extreme_rnaseq_counts.tsv",
    clock_path  = "results/clock/gdc_pan/TCGA-BRCA_predictions.rds",
    clock_col   = "ageAcc2.Horvath",
    labels_path = "results/clock/derived/TCGA-BRCA_tissue_quartiles_extreme_samples.txt",
    out_dir     = "results/differential_analysis",
    padj_cut    = 0.05,
    lfc_cut     = 1.0
  ),
  KIRC = list(
    project     = "TCGA-KIRC",
    counts_path = "results/differential_analysis/TCGA-KIRC_tissue_quartiles_extreme_rnaseq_counts.tsv",
    clock_path  = "results/clock/gdc_pan/TCGA-KIRC_predictions.rds",
    clock_col   = "ageAcc2.Horvath",
    labels_path = "results/clock/derived/TCGA-KIRC_tissue_quartiles_extreme_samples.txt",
    out_dir     = "results/differential_analysis",
    padj_cut    = 0.05,
    lfc_cut     = 1.0
  ),
  LIHC = list(
    project     = "TCGA-LIHC",
    counts_path = "results/differential_analysis/TCGA-LIHC_balanced_global_extreme_rnaseq_counts.tsv",
    clock_path  = "results/clock/gdc_pan/TCGA-LIHC_predictions.rds",
    clock_col   = "ageAcc2.Horvath",
    labels_path = "results/clock/derived/TCGA-LIHC_balanced_global_extreme_samples.txt",
    out_dir     = "results/differential_analysis",
    padj_cut    = 0.05,
    lfc_cut     = 1.0
  )
)

#======================================================================
# Argument parsing – very flexible
#======================================================================

parse_args <- function(args) {
  if (length(args) == 0) {
    return(default_params[c("BRCA", "KIRC", "LIHC")])
  }

  # If first argument is an existing file → legacy single-run mode
  if (file.exists(args[1])) {
    return(list(custom = list(
      project     = args[3] %||% "custom",
      counts_path = args[1],
      clock_path  = args[2],
      clock_col   = args[4] %||% "ageAcc2.Horvath",
      out_dir     = args[5] %||% "results/differential_analysis",
      padj_cut    = as.numeric(args[6]) %||% 0.05,
      lfc_cut     = as.numeric(args[7]) %||% 1.0
    )))
  }

  # Otherwise treat arguments as cohort names (comma/space separated)
  wanted <- unique(trimws(unlist(strsplit(paste(args, collapse = " "), "[, ]+"))))
  wanted <- wanted[nzchar(wanted)]
  if (!length(wanted)) wanted <- names(default_params)

  unknown <- setdiff(wanted, names(default_params))
  if (length(unknown)) stop("Unknown cohort(s): ", paste(unknown, collapse = ", "))

  default_params[wanted]
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

#======================================================================
# Core workflow (one cohort)
#======================================================================

run_one <- function(p) {
  message("\n=== ", p$project, " ===\n")

  # 1. Load counts
  if (!file.exists(p$counts_path)) stop("Counts file not found: ", p$counts_path)
  counts_raw <- fread(p$counts_path)
  if (!"gene_id" %in% names(counts_raw)) stop("No gene_id column")
  gene_id <- counts_raw$gene_id
  count_mat <- as.matrix(counts_raw[, -"gene_id", with = FALSE])
  rownames(count_mat) <- gene_id
  storage.mode(count_mat) <- "integer"
  colnames(count_mat) <- sub("^X", "", colnames(count_mat))

  # 2. Load clock predictions
  if (!file.exists(p$clock_path)) stop("Clock file not found: ", p$clock_path)
  clock_df <- readRDS(p$clock_path)
  if (!"id" %in% names(clock_df)) stop("Clock RDS must contain column 'id'")
  clock_df$clock_val <- as.numeric(clock_df[[p$clock_col]])
  clock_df <- clock_df[!is.na(clock_df$clock_val), ]
  clock_df$short   <- shorten_barcode(clock_df$id)
  clock_df$clean   <- sanitize_sample(clock_df$id)

  # 3. Build sample annotation that matches the count matrix
  sample_info <- data.frame(
    raw_name   = colnames(count_mat),
    short      = shorten_barcode(colnames(count_mat)),
    clean      = sanitize_sample(colnames(count_mat)),
    stringsAsFactors = FALSE
  )

  label_map <- load_labels(p$labels_path)

  # Try to match first by short barcode, then by cleaned full barcode
  merged <- merge(sample_info, clock_df[, c("short", "clean", "clock_val")],
                  by = "short", all.x = FALSE)
  if (nrow(merged) == 0) {
    merged <- merge(sample_info, clock_df[, c("short", "clean", "clock_val")],
                    by.x = "clean", by.y = "clean", all.x = FALSE)
  }
  if (nrow(merged) < 10) stop("Fewer than 10 overlapping samples – something is wrong with IDs")

  if (!is.null(label_map)) {
    merged <- merge(merged, label_map, by.x = "short", by.y = "sample_id", all.x = TRUE, sort = FALSE)
    matched_labels <- sum(!is.na(merged$classification))
    if (matched_labels > 0) {
      message(sprintf("Applied %d provided classifications (Accelerated lower tail)", matched_labels))
      merged$Type <- merged$classification
    } else {
      warning("Classification sheet provided but no sample IDs overlapped; falling back to clock sign")
    }
    merged$classification <- NULL
  }

  if (!"Type" %in% names(merged)) {
    merged$Type <- NA_character_
  }

  if (anyNA(merged$Type)) {
    fallback_idx <- which(is.na(merged$Type))
    if (length(fallback_idx)) {
      merged$Type[fallback_idx] <- ifelse(merged$clock_val[fallback_idx] < 0, "Accelerated", "Decelerated")
      message(sprintf("Assigned %d samples from clock sign (Accelerated = lower residual)", length(fallback_idx)))
    }
  }
  merged$Type <- factor(merged$Type, levels = c("Accelerated", "Decelerated"))
  merged$Type <- factor(merged$Type, levels = c("Accelerated", "Decelerated"))

  message("Samples retained: ", nrow(merged),
          " (", sum(merged$Type == "Accelerated"), " Accel, ",
                sum(merged$Type == "Decelerated"), " Decel)")

  # Subset & order count matrix exactly to the merged table
  count_mat <- count_mat[, merged$raw_name, drop = FALSE]

  # 4. Filter low-count genes
  keep_genes <- rowSums(count_mat) >= 10
  message("Genes before/after filtering: ", nrow(count_mat), " → ", sum(keep_genes))
  count_mat <- count_mat[keep_genes, , drop = FALSE]

  # 5. DESeq2
  colData <- DataFrame(Type = merged$Type, row.names = merged$raw_name)
  dds <- DESeqDataSetFromMatrix(countData = count_mat,
                                colData = colData,
                                design = ~ Type)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("Type", "Accelerated", "Decelerated"))
  vst_mat <- vst(dds, blind = FALSE)

  # 6. Results table
  res_dt <- as.data.table(as.data.frame(res), keep.rownames = "gene_id")
  res_dt[, significant := padj < p$padj_cut & abs(log2FoldChange) >= p$lfc_cut & !is.na(padj)]

  message("Significant genes: ", sum(res_dt$significant))

  # 7. Plots
  pca_plot <- plotPCA(vst_mat, intgroup = "Type", returnData = FALSE) +
    ggtitle(paste(p$project, "– VST PCA")) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom")

  volcano_df <- as.data.frame(res)
  volcano_df$gene <- rownames(volcano_df)
  volcano_df$`-log10(padj)` <- -log10(pmax(volcano_df$padj, .Machine$double.xmin))
  volcano_df$sig <- "Not sig."
  volcano_df$sig[volcano_df$padj < p$padj_cut &
                 abs(volcano_df$log2FoldChange) >= p$lfc_cut] <- "Significant"

  volcano_plot <- ggplot(volcano_df,
                         aes(x = log2FoldChange, y = `-log10(padj)`, color = sig)) +
    geom_point(alpha = 0.7, size = 1.2) +
    scale_color_manual(values = c("grey60", "#d62728")) +
    geom_vline(xintercept = c(-p$lfc_cut, p$lfc_cut), linetype = "dashed", colour = "grey40") +
    geom_hline(yintercept = -log10(p$padj_cut), linetype = "dashed", colour = "grey40") +
    labs(title = paste(p$project, "– Volcano Plot"),
         x = expression(Log[2]~Fold~Change),
         y = expression(-Log[10]~adjusted~p)) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none")

  # 8. Write everything
  dir.create(p$out_dir, showWarnings = FALSE, recursive = TRUE)
  out_prefix <- file.path(p$out_dir, paste0(p$project, "_accel_vs_decel"))

  fwrite(res_dt, paste0(out_prefix, "_full_results.tsv"), sep = "\t")
  fwrite(res_dt[significant == TRUE], paste0(out_prefix, "_significant.tsv"), sep = "\t")
  ggsave(paste0(out_prefix, "_PCA.png"),    pca_plot,     width = 7, height = 6, dpi = 300)
  ggsave(paste0(out_prefix, "_volcano.png"), volcano_plot, width = 7, height = 6, dpi = 300)

  message("→ Done! Files written to ", p$out_dir)
}

#======================================================================
# Main
#======================================================================

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  cohorts <- parse_args(args)

  for (i in seq_along(cohorts)) {
    tryCatch({
      run_one(cohorts[[i]])
    }, error = function(e) {
      message("\n!!! Error in ", cohorts[[i]]$project, ": ", e$message, "\n")
    })
  }
}

if (identical(environment(), globalenv())) {
  main()
}
