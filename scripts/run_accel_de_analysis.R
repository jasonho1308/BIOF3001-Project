#!/usr/bin/env Rscript
# Differential-expression workflow comparing epigenetic acceleration states in TCGA BRCA RNA-seq data

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
invisible(lapply(required_pkgs, ensure_pkg))

args <- commandArgs(trailingOnly = TRUE)
counts_path <- if (length(args) >= 1) args[[1]] else "results/differential_analysis/TCGA-BRCA_dec_rnaseq_counts.tsv"
clock_path <- if (length(args) >= 2) args[[2]] else "results/clock/gdc_pan/TCGA-BRCA_predictions.rds"
project_code <- if (length(args) >= 3) args[[3]] else "TCGA-BRCA"
clock_column <- if (length(args) >= 4) args[[4]] else "ageAcc2.Horvath"
output_dir <- if (length(args) >= 5) args[[5]] else "results/differential_analysis"
padj_cutoff <- if (length(args) >= 6) as.numeric(args[[6]]) else 0.05
lfc_cutoff <- if (length(args) >= 7) as.numeric(args[[7]]) else 1.0

shorten_barcode <- function(x) substr(x, 1, 16)

if (!file.exists(counts_path)) stop(sprintf("Counts file not found: %s", counts_path))
if (!file.exists(clock_path)) stop(sprintf("Clock file not found: %s", clock_path))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading counts matrix …")
counts_dt <- data.table::fread(counts_path)
if (!"gene_id" %in% names(counts_dt)) stop("Counts file must contain a 'gene_id' column.")

gene_ids <- counts_dt[["gene_id"]]
count_mat <- as.matrix(counts_dt[, -"gene_id", with = FALSE])
rownames(count_mat) <- gene_ids

message("Loading Horvath predictions …")
clock_df <- readRDS(clock_path)
if (!clock_column %in% colnames(clock_df)) {
  stop(sprintf("Column '%s' is missing from %s", clock_column, clock_path))
}

clock_dt <- data.table::as.data.table(clock_df)
clock_dt[, sample_short := shorten_barcode(id)]
clock_dt <- clock_dt[!is.na(sample_short) & nzchar(sample_short)]
clock_dt[, clock_value := get(clock_column)]
clock_dt <- clock_dt[!is.na(clock_value)]
clock_dt <- clock_dt[!duplicated(sample_short), .(sample_short, clock_value)]

sample_map <- data.table::data.table(column_name = colnames(count_mat))
sample_map[, sample_short := shorten_barcode(column_name)]
sample_map <- sample_map[!duplicated(column_name)]
merged_map <- merge(sample_map, clock_dt, by = "sample_short")

if (!nrow(merged_map)) {
  stop("No overlapping samples found between counts and clock predictions.")
}

merged_map <- merged_map[order(column_name)]
count_mat <- count_mat[, merged_map$column_name, drop = FALSE]

col_data <- data.frame(
  row.names = merged_map$column_name,
  sample_short = merged_map$sample_short,
  clock_value = merged_map$clock_value
)
col_data$Type <- ifelse(col_data$clock_value >= 0, "Accelerated", "Decelerated")
col_data$Type <- factor(col_data$Type, levels = c("Decelerated", "Accelerated"))

col_data <- col_data[!is.na(col_data$Type), , drop = FALSE]
count_mat <- count_mat[, rownames(col_data), drop = FALSE]

message(sprintf("Clock entries available: %d", nrow(col_data)))
if (length(unique(col_data$Type)) < 2) {
  stop("Need at least one sample in both Accelerated and Decelerated groups.")
}

keep_genes <- rowSums(count_mat) >= 10
message(sprintf("Filtered genes: %d retained of %d", sum(keep_genes), length(keep_genes)))
if (sum(keep_genes) < 10) stop("Too few genes remain after filtering.")
count_mat <- count_mat[keep_genes, , drop = FALSE]

dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(count_mat), colData = col_data, design = ~ Type)
message("Running DESeq2 …")
dds <- DESeq2::DESeq(dds)

message("Generating PCA plot …")
vsd <- DESeq2::vst(dds, blind = FALSE)
pca_plot <- DESeq2::plotPCA(vsd, intgroup = "Type") +
  ggplot2::labs(title = sprintf("%s Accelerated vs Decelerated PCA", project_code)) +
  ggplot2::theme_minimal()

pca_path <- file.path(output_dir, sprintf("%s_accel_vs_decel_pca.png", project_code))
print(pca_plot)
ggplot2::ggsave(pca_path, plot = pca_plot, width = 7, height = 5, dpi = 300)

message("Extracting DESeq2 results …")
res <- DESeq2::results(dds, contrast = c("Type", "Accelerated", "Decelerated"))
res_dt <- data.table::as.data.table(as.data.frame(res), keep.rownames = "gene_id")
res_dt <- res_dt[order(padj)]
res_dt[, significant := !is.na(padj) & padj < padj_cutoff & abs(log2FoldChange) >= lfc_cutoff]
message(sprintf("Significant genes (padj < %.3f): %d", padj_cutoff, sum(res_dt$significant, na.rm = TRUE)))

results_path <- file.path(output_dir, sprintf("%s_accel_vs_decel_deseq2_results.tsv", project_code))
sig_path <- file.path(output_dir, sprintf("%s_accel_vs_decel_sig.tsv", project_code))

data.table::fwrite(res_dt, results_path, sep = "\t")
data.table::fwrite(res_dt[significant == TRUE], sig_path, sep = "\t")

message("Creating volcano plot …")
volcano_dt <- res_dt[!is.na(padj)]
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
