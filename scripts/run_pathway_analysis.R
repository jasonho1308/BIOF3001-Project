#!/usr/bin/env Rscript
# Title: Hallmark ORA/GSEA workflow for DESeq2 tables
# Description: Ingests a DESeq2 result TSV, converts Ensembl identifiers to Entrez
#   IDs, and conducts both over-representation analysis (ORA) on significant genes
#   and ranked GSEA using Wald statistics or log2 fold-changes. Outputs include
#   tab-delimited Hallmark enrichment tables plus publication-ready plots (NES
#   bar chart and leading-edge trace) when discoveries surpass the significance
#   threshold. Designed to mirror the instructor’s preferred msigdbr/clusterProfiler
#   pipeline while remaining parameterizable via the command line.
# Usage:
#   Rscript run_pathway_analysis.R <results_tsv> <label>
#           [output_dir] [padj_cutoff] [lfc_cutoff]

ensure_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE)
  )
}

required_pkgs <- c(
  "data.table",
  "stringr",
  "dplyr",
  "tibble",
  "ggplot2",
  "msigdbr",
  "clusterProfiler",
  "enrichplot",
  "org.Hs.eg.db"
)
invisible(lapply(required_pkgs, ensure_pkg))

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("Description", "NES"))
}

# 1) Parse arguments describing DESeq2 results input and filtering thresholds
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) {
  stop(paste(
    "Usage:",
    "  run_pathway_analysis.R <results_tsv> <label> [output_dir] [padj_cutoff] [lfc_cutoff]",
    "  run_pathway_analysis.R <results_dir> [output_dir] [padj_cutoff] [lfc_cutoff]",
    "When a directory is supplied, all *_accel_vs_decel_deseq2_results.tsv files",
    "inside are processed sequentially, and labels are derived from filenames.",
    sep = "\n"
  ))
}

results_input <- args[[1]]
output_dir <- NULL
padj_cutoff <- 0.05
lfc_cutoff <- 1.0

args_consumed <- 1L
maybe_label <- NULL

if (!dir.exists(results_input)) {
  maybe_label <- if (length(args) >= 2) args[[2]] else NULL
  args_consumed <- if (is.null(maybe_label)) 1L else 2L
}

if (length(args) >= args_consumed + 1) {
  output_dir <- args[[args_consumed + 1]]
}
if (length(args) >= args_consumed + 2) {
  padj_cutoff <- as.numeric(args[[args_consumed + 2]])
}
if (length(args) >= args_consumed + 3) {
  lfc_cutoff <- as.numeric(args[[args_consumed + 3]])
}

if (is.null(output_dir)) {
  output_dir <- file.path("results", "pathway_analysis")
}

gene_id_col <- "gene_id"
logfc_col <- "log2FoldChange"
pval_col <- "pvalue"
padj_col <- "padj"

derive_label_from_path <- function(path) {
  base <- basename(path)
  label <- sub("_accel_vs_decel_deseq2_results\\.tsv$", "", base)
  label <- sub("_deseq2_results\\.tsv$", "", label)
  if (!nzchar(label)) {
    label <- tools::file_path_sans_ext(base)
  }
  label
}

collect_targets <- function(input_path, provided_label = NULL) {
  if (dir.exists(input_path)) {
    tsv_files <- sort(list.files(input_path, pattern = "_accel_vs_decel_deseq2_results\\.tsv$", full.names = TRUE))
    if (!length(tsv_files)) {
      stop(sprintf("No *_accel_vs_decel_deseq2_results.tsv files found in %s", input_path))
    }
    labels <- vapply(tsv_files, derive_label_from_path, character(1))
    return(stats::setNames(tsv_files, labels))
  }

  if (!file.exists(input_path)) {
    stop(sprintf("Results file not found: %s", input_path))
  }
  label <- if (!is.null(provided_label) && nzchar(provided_label)) provided_label else derive_label_from_path(input_path)
  stats::setNames(input_path, label)
}

run_pathway_for_dataset <- function(results_path, label, output_dir, padj_cutoff, lfc_cutoff) {
  message(sprintf("\n==== Processing %s ====", label))

  if (!file.exists(results_path)) {
    stop(sprintf("Results file not found: %s", results_path))
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE) # ensure results folder exists before writing artifacts

  message(sprintf("Reading differential results from %s", results_path))
  dt <- data.table::fread(results_path)
  needed_cols <- c(gene_id_col, logfc_col, pval_col, padj_col)
  missing_cols <- setdiff(needed_cols, names(dt))
  if (length(missing_cols)) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }

  dt$gene_base <- stringr::str_replace(dt[[gene_id_col]], "\\..*$", "")
  dt <- dt[!is.na(dt$gene_base) & nzchar(dt$gene_base)]
  dt <- dt[order(dt[[padj_col]], na.last = TRUE)]
  dt <- dt[!duplicated(dt$gene_base)] # keep best padj record per Ensembl gene before mapping

  # 3) Map Ensembl identifiers to Entrez IDs for Hallmark compatibility
  message("Mapping Ensembl IDs to Entrez IDs …")
  map_df <- clusterProfiler::bitr(
    unique(dt$gene_base),
    fromType = "ENSEMBL",
    toType = c("ENTREZID", "SYMBOL"),
    OrgDb = org.Hs.eg.db::org.Hs.eg.db
  )
  if (!nrow(map_df)) {
    stop("Failed to map any genes to Entrez IDs; cannot proceed with pathway analysis.")
  }
  data.table::setDT(map_df)
  merged_dt <- merge(dt, map_df, by.x = "gene_base", by.y = "ENSEMBL")
  if (!nrow(merged_dt)) {
    stop("No genes remain after mapping.")
  }

  sig_dt <- merged_dt[
    !is.na(merged_dt[[padj_col]]) & merged_dt[[padj_col]] <= padj_cutoff &
      !is.na(merged_dt[[logfc_col]]) & abs(merged_dt[[logfc_col]]) >= lfc_cutoff &
      !is.na(merged_dt$ENTREZID)
  ]

  if (!nrow(sig_dt)) {
    warning(sprintf("No significant genes detected for %s with padj <= %.3f and |log2FC| >= %.2f", label, padj_cutoff, lfc_cutoff))
  }

  universe_genes <- unique(merged_dt$ENTREZID[!is.na(merged_dt$ENTREZID)]) # reference set for ORA background
  sig_genes <- unique(sig_dt$ENTREZID) # foreground list that passed padj/LFC thresholds

  # 4) Get Hallmark gene sets (H collection) for Homo sapiens
  message("Preparing MSigDB Hallmark gene sets …")
  hallmark_df <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
  gene_col <- intersect(c("entrez_gene", "ncbi_gene"), names(hallmark_df))
  if (!length(gene_col)) {
    stop("Could not find an Entrez/Ncbi gene column in msigdbr output.")
  }
  term2gene <- hallmark_df[, c("gs_name", gene_col[1])] # build TERM2GENE mapping table once per run
  names(term2gene) <- c("gs_name", "entrez_gene")
  term2gene <- dplyr::distinct(term2gene)

  ora_path <- file.path(output_dir, sprintf("%s_hallmark_ora.tsv", label))
  gsea_path <- file.path(output_dir, sprintf("%s_hallmark_gsea.tsv", label))
  barplot_path <- file.path(output_dir, sprintf("%s_hallmark_gsea_barplot.png", label))
  enrichment_path <- file.path(output_dir, sprintf("%s_hallmark_gsea_top_enrichment.png", label))

  # 5) Run Over-Representation Analysis (ORA) on the significant subset
  if (length(sig_genes) >= 10) {
    message(sprintf("Running ORA for %d significant genes", length(sig_genes)))
    ora_res <- tryCatch({
      clusterProfiler::enricher(
        gene = sig_genes,
        TERM2GENE = term2gene,
        universe = universe_genes,
        pAdjustMethod = "BH",
        qvalueCutoff = 0.25
      )
    }, error = function(e) {
      warning(sprintf("ORA failed: %s", e$message))
      NULL
    })
    if (!is.null(ora_res) && nrow(ora_res@result)) {
      data.table::fwrite(data.table::as.data.table(ora_res@result), ora_path, sep = "\t")
      message(sprintf("ORA results written to %s", ora_path))
    } else {
      warning("ORA returned no pathways.")
    }
  } else {
    warning("Insufficient significant genes for ORA (need >= 10).")
  }

  # 6) Run GSEA on ranked statistics (full gene universe)
  message("Running GSEA …")
  merged_dt$stat_value <- ifelse(!is.na(merged_dt$stat), merged_dt$stat, merged_dt[[logfc_col]])
  rank_df <- merged_dt[!is.na(merged_dt$stat_value) & !is.na(merged_dt$ENTREZID)]
  rank_df <- rank_df[order(-abs(rank_df$stat_value))]
  rank_df <- rank_df[!duplicated(rank_df$ENTREZID)]
  rank_vector <- rank_df$stat_value
  names(rank_vector) <- rank_df$ENTREZID
  rank_vector <- sort(rank_vector, decreasing = TRUE)

  if (!length(rank_vector)) {
    warning("Ranked gene list is empty; skipping GSEA.")
    gsea_res <- NULL
  } else {
    set.seed(123)
    gsea_res <- tryCatch({
      clusterProfiler::GSEA(
        geneList = rank_vector,
        TERM2GENE = term2gene,
        minGSSize = 10,
        maxGSSize = 500,
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        seed = TRUE
      )
    }, error = function(e) {
      warning(sprintf("GSEA failed: %s", e$message))
      NULL
    })
  }

  if (!is.null(gsea_res) && nrow(gsea_res@result)) {
    gsea_dt <- data.table::as.data.table(gsea_res@result)
    data.table::setorder(gsea_dt, p.adjust)
    data.table::fwrite(gsea_dt, gsea_path, sep = "\t")
    message(sprintf("GSEA results written to %s", gsea_path))

    sig_gsea <- gsea_dt[!is.na(p.adjust) & p.adjust < 0.05]
    if (nrow(sig_gsea)) {
      plot_dt <- data.table::copy(sig_gsea)
      plot_dt$Description <- gsub("^HALLMARK_", "", plot_dt$ID)
      plot_dt <- plot_dt[order(p.adjust)]
      plot_dt$Description <- factor(plot_dt$Description, levels = plot_dt$Description[order(plot_dt$NES)])
      p_bar <- ggplot2::ggplot(plot_dt, ggplot2::aes_string(x = "Description", y = "NES")) +
        ggplot2::geom_col(ggplot2::aes(fill = -log10(p.adjust))) +
        ggplot2::coord_flip() +
        ggplot2::scale_fill_gradient(name = "-log10(padj)", low = "skyblue", high = "firebrick") +
        ggplot2::labs(
          title = sprintf("Significant Hallmark Pathways (n = %d)", nrow(plot_dt)),
          x = NULL,
          y = "Normalized Enrichment Score (NES)"
        ) +
        ggplot2::theme_minimal(base_size = 13)
      ggplot2::ggsave(barplot_path, p_bar, width = 8, height = 6, dpi = 300)
      message(sprintf("Saved GSEA barplot to %s", barplot_path))

      top_term <- gsea_dt$ID[1]
      ep <- enrichplot::gseaplot2(gsea_res, geneSetID = top_term, title = top_term)
      ggplot2::ggsave(enrichment_path, ep, width = 8, height = 6, dpi = 300)
      message(sprintf("Saved enrichment plot to %s", enrichment_path))
    } else {
      message("No significant hallmark pathways (padj < 0.05) for barplot/enrichment plot.")
    }
  } else {
    warning("GSEA returned no pathways.")
  }

  message(sprintf("Finished pathway analysis for %s.", label))
}

targets <- collect_targets(results_input, maybe_label)

for (label in names(targets)) {
  target_path <- targets[[label]]
  tryCatch(
    run_pathway_for_dataset(target_path, label, output_dir, padj_cutoff, lfc_cutoff),
    error = function(e) {
      warning(sprintf("Failed pathway analysis for %s (%s): %s", label, target_path, e$message))
    }
  )
}
