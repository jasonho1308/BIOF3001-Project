#!/usr/bin/env Rscript
# Title: Differential-expression workflow for acceleration status contrasts
# Description: Consumes a count matrix and a Horvath-classification RDS to compare
#   Accelerated vs Decelerated samples with DESeq2. The script harmonizes TCGA
#   barcodes, derives binary status labels (Accelerated/Decelerated) from either
#   categorical annotations or ±1 scores, and produces quality-control plots,
#   complete DESeq2 tables, and filtered significant gene lists.
# Usage:
#   Rscript run_accel_de_analysis.R <counts_tsv> <clock_rds>
#           [project_code] [clock_column] [output_dir] [padj_cutoff] [lfc_cutoff]
# Inputs:
#   counts_tsv    — tab-delimited matrix with Ensembl IDs plus TCGA columns.
#   clock_rds     — RDS containing a column with Horvath residuals/classification.
# Outputs (written to output_dir):
#   *_deseq2_results.tsv, *_sig.tsv, PCA PNG, volcano PNG.
# Notes:
#   • Samples are matched via 16-character barcodes.
#   • Genes with total counts < 10 are removed to avoid unstable dispersion.
#   • Wald statistics and log2FCs correspond to Accelerated/Decelerated contrast.

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

required_pkgs <- c("data.table", "DESeq2", "ggplot2")
invisible(lapply(required_pkgs, ensure_pkg)) # install/load dependencies if they are missing

# 1) Parse command-line arguments or fall back to project defaults
args <- commandArgs(trailingOnly = TRUE)
counts_path <- if (length(args) >= 1) args[[1]] else "results/differential_analysis/TCGA-BRCA_dec_rnaseq_counts.tsv"
clock_path <- if (length(args) >= 2) args[[2]] else "results/clock/gdc_pan/TCGA-BRCA_predictions.rds"
project_code <- if (length(args) >= 3) args[[3]] else "TCGA-BRCA"
clock_column <- if (length(args) >= 4) args[[4]] else "ageAcc2.Horvath"
output_dir <- if (length(args) >= 5) args[[5]] else "results/differential_analysis"
padj_cutoff <- if (length(args) >= 6) as.numeric(args[[6]]) else 0.05
lfc_cutoff <- if (length(args) >= 7) as.numeric(args[[7]]) else 1.0

# 2) Helper utilities for barcode alignment and status derivation
#    TCGA aliquot barcodes share participant+sample information in first 16 chars.
shorten_barcode <- function(x) substr(x, 1, 16)

derive_status <- function(values) {
  # Harmonize string labels and ±1 numeric residual encodings into a single factor
  vals <- values
  if (is.factor(vals)) vals <- as.character(vals)
  if (is.character(vals)) {
    trimmed <- trimws(vals)
    unique_vals <- unique(trimmed)
    if (length(unique_vals) && all(unique_vals %in% c("Accelerated", "Decelerated"))) {
      status <- trimmed
      return(factor(status, levels = c("Decelerated", "Accelerated")))
    }
  }
  if (is.numeric(vals)) {
    unique_vals <- unique(na.omit(vals))
    if (length(unique_vals) && all(unique_vals %in% c(-1, 1))) {
      status <- ifelse(vals > 0, "Accelerated", "Decelerated")
      return(factor(status, levels = c("Decelerated", "Accelerated")))
    }
  }
  status <- ifelse(vals >= 0, "Accelerated", "Decelerated")
  factor(status, levels = c("Decelerated", "Accelerated"))
}

# 3) Load expression counts and Horvath residuals/classification
if (!file.exists(counts_path)) stop(sprintf("Counts file not found: %s", counts_path))
if (!file.exists(clock_path)) stop(sprintf("Clock file not found: %s", clock_path))
message("Reading counts matrix …")
counts_dt <- data.table::fread(counts_path) # Expect gene_id + sample columns
if (!"gene_id" %in% names(counts_dt)) stop("Counts file must contain a 'gene_id' column.")

gene_ids <- counts_dt[["gene_id"]]
count_mat <- as.matrix(counts_dt[, -"gene_id", with = FALSE]) # convert data.table to numeric matrix
rownames(count_mat) <- gene_ids

message("Loading Horvath predictions …")
clock_df <- readRDS(clock_path) # Contains sample IDs plus clock residuals
if (!clock_column %in% colnames(clock_df)) {
  stop(sprintf("Column '%s' is missing from %s", clock_column, clock_path))
}

clock_dt <- data.table::as.data.table(clock_df)
clock_dt[, sample_short := shorten_barcode(id)]
clock_dt <- clock_dt[!is.na(sample_short) & nzchar(sample_short)]
clock_dt[, clock_value := get(clock_column)]
clock_dt <- clock_dt[!is.na(clock_value)]
clock_dt <- clock_dt[!duplicated(sample_short), .(sample_short, clock_value)]

# 4) Align samples between counts and residual table; drop any that cannot be matched
# Build lookup to align count matrix columns with residual entries
sample_map <- data.table::data.table(column_name = colnames(count_mat))
sample_map[, sample_short := shorten_barcode(column_name)]
sample_map <- sample_map[!duplicated(column_name)]
merged_map <- merge(sample_map, clock_dt, by = "sample_short")

if (!nrow(merged_map)) {
  stop("No overlapping samples found between counts and clock predictions.")
}

merged_map <- merged_map[order(column_name)]
count_mat <- count_mat[, merged_map$column_name, drop = FALSE] # ensure counts follow col_data row order

col_data <- data.frame(
  row.names = merged_map$column_name,
  sample_short = merged_map$sample_short,
  clock_value = merged_map$clock_value
)
col_data$Type <- derive_status(col_data$clock_value) # convert numeric/string classifications into factors

col_data <- col_data[!is.na(col_data$Type), , drop = FALSE] # drop samples lacking usable labels
count_mat <- count_mat[, rownames(col_data), drop = FALSE] # keep matrix columns synchronized

# 5) Filter to samples with usable labels and apply basic gene-count threshold
message(sprintf("Clock entries available: %d", nrow(col_data)))
if (length(unique(col_data$Type)) < 2) {
  stop("Need at least one sample in both Accelerated and Decelerated groups.")
}

keep_genes <- rowSums(count_mat) >= 10 # Basic low-count filter for DESeq2 stability
message(sprintf("Filtered genes: %d retained of %d", sum(keep_genes), length(keep_genes)))
if (sum(keep_genes) < 10) stop("Too few genes remain after filtering.")
count_mat <- count_mat[keep_genes, , drop = FALSE]

# 6) Build DESeq2 dataset and model Accelerated vs Decelerated contrast
dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(count_mat), colData = col_data, design = ~ Type) # Accelerated vs Decelerated contrast
message("Running DESeq2 …")
dds <- DESeq2::DESeq(dds)

# 7) QC step: visualize VST-transformed PCA by acceleration status
message("Generating PCA plot …")
vsd <- DESeq2::vst(dds, blind = FALSE)
pca_plot <- DESeq2::plotPCA(vsd, intgroup = "Type") + # visualize separation between acceleration states
  ggplot2::labs(title = sprintf("%s Accelerated vs Decelerated PCA", project_code)) +
  ggplot2::theme_minimal()

pca_path <- file.path(output_dir, sprintf("%s_accel_vs_decel_pca.png", project_code))
print(pca_plot)
ggplot2::ggsave(pca_path, plot = pca_plot, width = 7, height = 5, dpi = 300)

# 8) Extract Wald statistics and annotate significant genes
message("Extracting DESeq2 results …")
res <- DESeq2::results(dds, contrast = c("Type", "Accelerated", "Decelerated")) # log2FC positive => Accelerated up
res_dt <- data.table::as.data.table(as.data.frame(res), keep.rownames = "gene_id")
res_dt <- res_dt[order(padj)]
res_dt[, significant := !is.na(padj) & padj < padj_cutoff & abs(log2FoldChange) >= lfc_cutoff]
message(sprintf("Significant genes (padj < %.3f): %d", padj_cutoff, sum(res_dt$significant, na.rm = TRUE)))

results_path <- file.path(output_dir, sprintf("%s_accel_vs_decel_deseq2_results.tsv", project_code))
sig_path <- file.path(output_dir, sprintf("%s_accel_vs_decel_sig.tsv", project_code))

data.table::fwrite(res_dt, results_path, sep = "\t") # complete table for downstream stats
data.table::fwrite(res_dt[significant == TRUE], sig_path, sep = "\t") # compact list of hits for sharing

# 9) Plot Accelerated vs Decelerated contrast and persist outputs
message("Creating volcano plot …")
volcano_dt <- res_dt[!is.na(padj)] # drop NA padj rows to avoid plotting artifacts
volcano_dt[, neg_log10_padj := -log10(pmax(padj, .Machine$double.xmin))]
volcano_plot <- ggplot2::ggplot(volcano_dt, ggplot2::aes(x = log2FoldChange, y = neg_log10_padj)) +
  ggplot2::geom_point(ggplot2::aes(color = significant), alpha = 0.7, size = 1.2) +
  ggplot2::scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "firebrick"), guide = "none") +
  ggplot2::geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "black") +
  ggplot2::geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "black") +
  ggplot2::labs(
    title = sprintf("%s Accelerated vs Decelerated", project_code),
    x = "log2 Fold Change (Accelerated / Decelerated)",
    y = expression(-log[10](adj.~p~value))
  ) +
  ggplot2::theme_minimal()

volcano_path <- file.path(output_dir, sprintf("%s_accel_vs_decel_volcano.png", project_code))
ggplot2::ggsave(volcano_path, plot = volcano_plot, width = 7, height = 5, dpi = 300)

message("Differential analysis complete.")
message(sprintf("Results: %s", results_path))
message(sprintf("Significant subset: %s", sig_path))
message(sprintf("PCA plot: %s", pca_path))
message(sprintf("Volcano plot: %s", volcano_path))
