#!/usr/bin/env Rscript
# Title: Differential-expression workflow for acceleration status contrasts
# Description: Reusable CLI script that ingests TCGA-style count matrices and
#   Horvath residuals/classifications, aligns samples by 16-character barcodes,
#   runs DESeq2, and emits harmonized QC plots plus DE tables.

ensure_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

required_pkgs <- c("data.table", "DESeq2", "ggplot2", "rlang")
invisible(lapply(required_pkgs, ensure_pkg))

defaults <- list(
  counts_path = "results/differential_analysis/TCGA-BRCA_tissue_quartiles_extreme_rnaseq_counts.tsv",
  clock_path = "results/clock/gdc_pan/TCGA-BRCA_predictions.rds",
  project_code = "TCGA-BRCA",
  clock_column = "ageAcc2.Horvath",
  output_dir = "results/differential_analysis",
  padj_cutoff = 0.05,
  lfc_cutoff = 1.0
)

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "sample_short", "clock_value", "column_name", "column_rank", "significant",
    "id", "padj", "log2FoldChange", "PC1", "PC2", "Type", "neg_log10_padj"
  ))
}

parse_args <- function(args) {
  params <- defaults
  if (length(args) >= 1) params$counts_path <- args[[1]]
  if (length(args) >= 2) params$clock_path <- args[[2]]
  if (length(args) >= 3) params$project_code <- args[[3]]
  if (length(args) >= 4) params$clock_column <- args[[4]]
  if (length(args) >= 5) params$output_dir <- args[[5]]
  if (length(args) >= 6) params$padj_cutoff <- as.numeric(args[[6]])
  if (length(args) >= 7) params$lfc_cutoff <- as.numeric(args[[7]])
  params
}

shorten_barcode <- function(x) substr(x, 1, 16)

derive_status <- function(values) {
  vals <- values
  if (is.factor(vals)) vals <- as.character(vals)
  if (is.character(vals)) {
    trimmed <- trimws(vals)
    uniq <- unique(trimmed)
    if (length(uniq) && all(uniq %in% c("Accelerated", "Decelerated"))) {
      return(factor(trimmed, levels = c("Decelerated", "Accelerated")))
    }
  }
  if (is.numeric(vals)) {
    uniq <- unique(na.omit(vals))
    if (length(uniq) && all(uniq %in% c(-1, 1))) {
      return(factor(ifelse(vals > 0, "Accelerated", "Decelerated"), levels = c("Decelerated", "Accelerated")))
    }
  }
  factor(ifelse(vals >= 0, "Accelerated", "Decelerated"), levels = c("Decelerated", "Accelerated"))
}

load_counts <- function(path) {
  if (!file.exists(path)) stop(sprintf("Counts file not found: %s", path))
  dt <- data.table::fread(path)
  if (!"gene_id" %in% names(dt)) stop("Counts file must contain a 'gene_id' column.")
  gene_ids <- dt$gene_id
  mat <- as.matrix(dt[, -"gene_id", with = FALSE])
  rownames(mat) <- gene_ids
  mode(mat) <- "numeric"
  mat
}

load_clock <- function(path, clock_column) {
  if (!file.exists(path)) stop(sprintf("Clock file not found: %s", path))
  df <- as.data.frame(readRDS(path))
  if (!"id" %in% names(df)) stop("Clock file must contain an 'id' column with TCGA barcodes.")
  if (!clock_column %in% names(df)) stop(sprintf("Column '%s' missing from %s", clock_column, path))
  df$sample_short <- shorten_barcode(df$id)
  df <- df[!is.na(df$sample_short) & nzchar(df$sample_short), ]
  df$clock_value <- df[[clock_column]]
  df <- df[!is.na(df$clock_value), c("sample_short", "clock_value"), drop = FALSE]
  df[!duplicated(df$sample_short), , drop = FALSE]
}

align_samples <- function(count_mat, clock_df) {
  sample_map <- data.frame(
    column_name = colnames(count_mat),
    sample_short = shorten_barcode(colnames(count_mat)),
    stringsAsFactors = FALSE
  )
  sample_map <- sample_map[!duplicated(sample_map$column_name), ]
  merged <- merge(sample_map, clock_df, by = "sample_short")
  if (!nrow(merged)) stop("No overlapping samples found between counts and clock predictions.")
  merged$column_rank <- match(merged$column_name, colnames(count_mat))
  merged <- merged[order(merged$column_rank), ]
  aligned_counts <- count_mat[, merged$column_name, drop = FALSE]
  col_data <- data.frame(
    row.names = merged$column_name,
    sample_short = merged$sample_short,
    clock_value = merged$clock_value,
    Type = derive_status(merged$clock_value)
  )
  col_data <- col_data[!is.na(col_data$Type), , drop = FALSE]
  aligned_counts <- aligned_counts[, rownames(col_data), drop = FALSE]
  list(counts = aligned_counts, col_data = col_data)
}

filter_counts <- function(count_mat, min_total = 10L) {
  keep <- rowSums(count_mat) >= min_total
  if (sum(keep) < 10) stop("Too few genes remain after filtering; adjust cohort or threshold.")
  message(sprintf("Filtered genes: %d retained of %d", sum(keep), length(keep)))
  count_mat[keep, , drop = FALSE]
}

run_deseq <- function(count_mat, col_data) {
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(count_mat), colData = col_data, design = ~ Type)
  message("Running DESeq2 …")
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds, contrast = c("Type", "Accelerated", "Decelerated"))
  list(dds = dds, res = res)
}

tidy_results <- function(res, padj_cutoff, lfc_cutoff) {
  df <- as.data.frame(res)
  df$gene_id <- rownames(df)
  df <- df[order(df$padj), ]
  df$significant <- !is.na(df$padj) & df$padj < padj_cutoff & abs(df$log2FoldChange) >= lfc_cutoff
  df
}

build_pca_plot <- function(dds, project_code) {
  vst <- DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
  pca <- stats::prcomp(t(SummarizedExperiment::assay(vst)))
  pca_df <- data.frame(
    row.names = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    Type = SummarizedExperiment::colData(vst)$Type
  )
  percent_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2, color = Type)) +
    ggplot2::geom_point(size = 2, alpha = 0.9) +
    ggplot2::scale_color_manual(values = c("Decelerated" = "#1b9e77", "Accelerated" = "#d95f02")) +
    ggplot2::labs(
      title = sprintf("%s: VST PCA", project_code),
      x = sprintf("PC1 (%.1f%%)", percent_var[1]),
      y = sprintf("PC2 (%.1f%%)", percent_var[2])
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom", plot.title = ggplot2::element_text(face = "bold"))
}

build_volcano_plot <- function(res_dt, padj_cutoff, lfc_cutoff, project_code) {
  plot_df <- res_dt[!is.na(res_dt$padj), , drop = FALSE]
  plot_df$neg_log10_padj <- -log10(pmax(plot_df$padj, .Machine$double.xmin))
  ggplot2::ggplot(plot_df, ggplot2::aes(x = log2FoldChange, y = neg_log10_padj)) +
    ggplot2::geom_point(ggplot2::aes(color = significant), alpha = 0.75, size = 1.4) +
    ggplot2::scale_color_manual(values = c("FALSE" = "grey80", "TRUE" = "firebrick"), guide = "none") +
    ggplot2::geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "#4d4d4d") +
    ggplot2::geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "#4d4d4d") +
    ggplot2::labs(
      title = sprintf("%s Accelerated vs Decelerated", project_code),
      x = "log2 Fold Change (Accelerated / Decelerated)",
      y = expression(-log[10](adj.~p~value))
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
}

write_outputs <- function(res_dt, output_dir, project_code, pca_plot, volcano_plot) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  results_path <- file.path(output_dir, sprintf("%s_accel_vs_decel_deseq2_results.tsv", project_code))
  sig_path <- file.path(output_dir, sprintf("%s_accel_vs_decel_sig.tsv", project_code))
  pca_path <- file.path(output_dir, sprintf("%s_accel_vs_decel_pca.png", project_code))
  volcano_path <- file.path(output_dir, sprintf("%s_accel_vs_decel_volcano.png", project_code))

  data.table::fwrite(data.table::as.data.table(res_dt), results_path, sep = "\t")
  data.table::fwrite(data.table::as.data.table(res_dt[res_dt$significant, , drop = FALSE]), sig_path, sep = "\t")
  ggplot2::ggsave(pca_path, plot = pca_plot, width = 7, height = 5, dpi = 300)
  ggplot2::ggsave(volcano_path, plot = volcano_plot, width = 7, height = 5, dpi = 300)

  message("Differential analysis complete.")
  message(sprintf("Results: %s", results_path))
  message(sprintf("Significant subset: %s", sig_path))
  message(sprintf("PCA plot: %s", pca_path))
  message(sprintf("Volcano plot: %s", volcano_path))
}

main <- function() {
  params <- parse_args(commandArgs(trailingOnly = TRUE))
  count_mat <- load_counts(params$counts_path)
  clock_dt <- load_clock(params$clock_path, params$clock_column)
  aligned <- align_samples(count_mat, clock_dt)
  message(sprintf("Clock entries available: %d", nrow(aligned$col_data)))
  if (length(unique(aligned$col_data$Type)) < 2) {
    stop("Need at least one sample in both Accelerated and Decelerated groups.")
  }
  filtered_counts <- filter_counts(aligned$counts)
  dds_res <- run_deseq(filtered_counts, aligned$col_data)
  res_dt <- tidy_results(dds_res$res, params$padj_cutoff, params$lfc_cutoff)
  message(sprintf("Significant genes (padj < %.3f): %d", params$padj_cutoff, sum(res_dt$significant, na.rm = TRUE)))
  pca_plot <- build_pca_plot(dds_res$dds, params$project_code)
  volcano_plot <- build_volcano_plot(res_dt, params$padj_cutoff, params$lfc_cutoff, params$project_code)
  write_outputs(res_dt, params$output_dir, params$project_code, pca_plot, volcano_plot)
}

if (identical(environment(), globalenv())) {
  main()
}
