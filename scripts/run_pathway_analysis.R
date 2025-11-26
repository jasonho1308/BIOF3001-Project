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

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(msigdbr)
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(org.Hs.eg.db)
})

# ===================================================================
# Helper: install missing packages
# ===================================================================
ensure_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing ", pkg, "...")
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
    library(pkg, character.only = TRUE)
  }
}
invisible(lapply(c("dplyr","tibble","readr","msigdbr","clusterProfiler","enrichplot","org.Hs.eg.db"), ensure_pkg))

set.seed(123)

# ===================================================================
# Parse arguments
# ===================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript run_hallmark_analysis.R <input> <label> [out_dir] [padj_cut=0.05] [lfc_cut=1.0]\n",
       "       Or provide a directory containing *_deseq2_results.tsv files.")
}

input_path   <- args[1]
label_in     <- if (length(args) >= 2) args[2] else NULL
out_dir      <- if (length(args) >= 3) args[3] else "results/pathway_analysis"
padj_cut     <- if (length(args) >= 4) as.numeric(args[4]) else 0.05
lfc_cut      <- if (length(args) >= 5) as.numeric(args[5]) else 1.0

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ===================================================================
# Load Hallmark gene sets (once)
# ===================================================================
message("Loading MSigDB Hallmark gene sets...")
hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene) %>%
  distinct()

term2gene <- hallmark %>%
  dplyr::select(term = gs_name, gene = entrez_gene)

# ===================================================================
# Process one DESeq2 results file
# ===================================================================
run_one <- function(tsv_path, label) {
  message("\n=== ", label, " ===\nReading: ", tsv_path)
  
  dt <- read_tsv(tsv_path, col_types = cols()) %>%
    dplyr::rename_with(~"gene_id", any_of(c("gene_id", "rowname", "...1"))) %>%
    mutate(gene_base = str_remove(gene_id, "\\..*$")) %>%
    filter(!is.na(gene_base), gene_base != "")
  
  if (!all(c("log2FoldChange", "padj") %in% names(dt))) {
    stop("Required columns missing: need log2FoldChange and padj")
  }
  
  # Prefer Wald statistic for ranking, fall back to log2FC
  dt <- dt %>%
    mutate(rank_stat = ifelse(is.na(stat) | is.nan(stat), log2FoldChange, stat))
  
  # Map Ensembl → Entrez (keep best per gene)
  message("Mapping Ensembl → Entrez...")
  gene_list <- unique(dt$gene_base)
  map <- bitr(gene_list, fromType = "ENSEMBL", toType = c("ENTREZID", "SYMBOL"), OrgDb = org.Hs.eg.db)
  
  dt_mapped <- dt %>%
    inner_join(map, by = c("gene_base" = "ENSEMBL")) %>%
    group_by(ENTREZID) %>%
    slice_max(order_by = abs(rank_stat), n = 1, with_ties = FALSE) %>%
    ungroup()
  
  if (nrow(dt_mapped) == 0) stop("No genes mapped to Entrez IDs")
  
  # ---------------------------------------------------------------
  # 1. ORA on significant genes
  # ---------------------------------------------------------------
  sig_genes <- dt_mapped %>%
    filter(padj <= padj_cut, abs(log2FoldChange) >= lfc_cut, !is.na(ENTREZID)) %>%
    pull(ENTREZID) %>% unique()
  
  universe <- dt_mapped %>% pull(ENTREZID) %>% unique()
  
  ora_file <- file.path(out_dir, paste0(label, "_hallmark_ora.tsv"))
  if (length(sig_genes) >= 10) {
    message("Running ORA (", length(sig_genes), " significant genes)...")
    ora <- enricher(
      gene = sig_genes,
      universe = universe,
      TERM2GENE = term2gene,
      pvalueCutoff = 0.25,
      pAdjustMethod = "BH",
      qvalueCutoff = 1
    )
    if (!is.null(ora) && nrow(ora@result) > 0) {
      write_tsv(as_tibble(ora@result), ora_file)
      message("→ ORA results: ", ora_file)
    }
  } else {
    message("Skipping ORA (<10 significant genes)")
  }
  
  # ---------------------------------------------------------------
  # 2. GSEA on full ranked list
  # ---------------------------------------------------------------
  message("Preparing ranked list for GSEA...")
  rank_vec <- dt_mapped$rank_stat
  names(rank_vec) <- dt_mapped$ENTREZID
  rank_vec <- sort(rank_vec, decreasing = TRUE)
  
  gsea_file <- file.path(out_dir, paste0(label, "_hallmark_gsea.tsv"))
  gsea_res <- GSEA(
    geneList = rank_vec,
    TERM2GENE = term2gene,
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    seed = TRUE,
    verbose = FALSE
  )
  
  gsea_tbl <- as_tibble(gsea_res@result) %>%
    arrange(p.adjust)
  
  write_tsv(gsea_tbl, gsea_file)
  message("→ GSEA results: ", gsea_file)
  
  # ---------------------------------------------------------------
  # 3. Plots (only if significant pathways exist)
  # ---------------------------------------------------------------
  sig_pathways <- gsea_tbl %>% filter(p.adjust < 0.05)
  
  if (nrow(sig_pathways) > 0) {
    message("Found ", nrow(sig_pathways), " significant Hallmark pathways → making plots")
    
    # Barplot
    plot_df <- sig_pathways %>%
      mutate(Description = str_remove(ID, "^HALLMARK_")) %>%
      arrange(NES)
    
    p_bar <- ggplot(plot_df, aes(x = reorder(Description, NES), y = NES, fill = -log10(p.adjust))) +
      geom_col(width = 0.8) +
      coord_flip() +
      scale_fill_gradient(name = "-log10(padj)", low = "steelblue", high = "firebrick") +
      labs(
        title = paste0("Significant Hallmark Pathways (n = ", nrow(plot_df), ") – ", label),
        x = NULL,
        y = "Normalized Enrichment Score (NES)"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        plot.title = element_text(face = "bold"),
        axis.text.y = element_text(size = 11)
      )
    
    bar_file <- file.path(out_dir, paste0(label, "_hallmark_gsea_barplot.png"))
    ggsave(bar_file, p_bar, width = 9, height = max(5, 0.25 * nrow(plot_df) + 2), dpi = 300, bg = "white")
    message("→ Barplot: ", bar_file)
    
    # Top enrichment plot
    top_id <- gsea_tbl$ID[1]
    p_enrich <- gseaplot2(gsea_res, geneSetID = top_id, title = str_remove(top_id, "^HALLMARK_"))
    
    enrich_file <- file.path(out_dir, paste0(label, "_hallmark_gsea_top_enrichment.png"))
    ggsave(enrich_file, p_enrich, width = 8, height = 6, dpi = 300, bg = "white")
    message("→ Top enrichment plot: ", enrich_file)
    
  } else {
    message("No significant Hallmark pathways (padj < 0.05)")
  }
}

# ===================================================================
# Main: collect files and run
# ===================================================================
if (dir.exists(input_path)) {
  files <- list.files(input_path, pattern = "_deseq2_results\\.tsv$|_full_results\\.tsv$", full.names = TRUE)
  if (length(files) == 0) stop("No DESeq2 results found in directory")
  
  labels <- basename(files) %>%
    str_remove("_accel_vs_decel.*|_full_results\\.tsv$|_deseq2_results\\.tsv$") %>%
    str_remove("\\.tsv$")
  
  targets <- setNames(files, labels)
  
} else {
  if (!file.exists(input_path)) stop("File not found: ", input_path)
  label <- if (is.null(label_in) || label_in == "") {
    tools::file_path_sans_ext(basename(input_path))
  } else label_in
  targets <- setNames(input_path, label)
}

for (lbl in names(targets)) {
  tryCatch({
    run_one(targets[[lbl]], lbl)
  }, error = function(e) {
    warning("Failed for ", lbl, ": ", e$message)
  })
}

message("\nAll done! Results saved in: ", normalizePath(out_dir))
