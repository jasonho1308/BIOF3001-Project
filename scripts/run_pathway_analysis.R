# Perform pathway analysis on DE results
# Runs GSEA using Hallmark gene sets from MSigDB on DESeq2 results

library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(tibble)
library(parallel)

set.seed(123)  # Ensures reproducible GSEA results

# Create output directory if it doesn't exist
dir.create("results/pathway_analysis", showWarnings = FALSE, recursive = TRUE)

# Define projects to analyze (matching DE.R)
projects <- c("BRCA", "THCA", "UCEC")

# Load Hallmark gene sets from MSigDB (H collection) for Homo sapiens
message("Loading Hallmark gene sets from MSigDB...")
hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene)
message(sprintf("Loaded %d gene set entries", nrow(hallmark)))

# Function to run GSEA for a single project
run_gsea_analysis <- function(project_name) {
  message(sprintf("\n========== Running GSEA for %s ==========", project_name))
  
  # Load DESeq2 results
  res_file <- sprintf("results/differential_analysis/%s_deseq2_results.rds", project_name)
  if (!file.exists(res_file)) {
    warning(sprintf("DESeq2 results not found for %s, skipping...", project_name))
    return(NULL)
  }
  
  res <- readRDS(res_file)
  message(sprintf("Loaded DESeq2 results: %d genes", nrow(res)))
  
  # 1) Build a ranked vector from DESeq2 results
  #    Use Wald statistic (signed, already scales by SE)
  res_tbl <- as_tibble(res, rownames = "gene_id") %>%
    filter(!is.na(stat))
  
  # Extract gene symbols from Ensembl IDs (format: ENSG00000000003.15)
  # The gene symbol needs to be mapped from the Ensembl ID
  res_tbl <- res_tbl %>%
    mutate(ensembl_base = gsub("\\..*", "", gene_id))  # Remove version number
  
  message(sprintf("Genes with valid stat: %d", nrow(res_tbl)))
  
  # 2) Map Ensembl IDs to Entrez IDs
  message("Mapping Ensembl IDs to Entrez IDs...")
  map_df <- clusterProfiler::bitr(res_tbl$ensembl_base,
                                  fromType = "ENSEMBL",
                                  toType   = c("ENTREZID", "SYMBOL"),
                                  OrgDb    = org.Hs.eg.db)
  
  message(sprintf("Successfully mapped %d genes", nrow(map_df)))
  
  # 3) Join and collapse duplicates (keep the entry with the largest |stat|)
  rank_df <- res_tbl %>%
    inner_join(map_df, by = c("ensembl_base" = "ENSEMBL")) %>%
    group_by(ENTREZID) %>%
    slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
    ungroup()
  
  message(sprintf("Unique Entrez IDs for GSEA: %d", nrow(rank_df)))
  
  # Create ranked gene list
  ranks <- rank_df$stat
  names(ranks) <- rank_df$ENTREZID
  ranks <- sort(ranks, decreasing = TRUE)
  
  message(sprintf("Rank range: %.2f to %.2f", min(ranks), max(ranks)))
  
  # 4) Run GSEA with Hallmark gene sets
  message("Running GSEA...")
  gsea_res <- clusterProfiler::GSEA(geneList     = ranks,
                                    TERM2GENE    = hallmark,
                                    minGSSize    = 10,
                                    maxGSSize    = 500,
                                    pvalueCutoff = 1,   # Keep all; filter later
                                    pAdjustMethod = "BH",
                                    seed = TRUE)
  
  # 5) Process results
  gsea_df <- as_tibble(gsea_res@result) %>%
    arrange(p.adjust)
  
  sig_gsea <- gsea_df %>% filter(p.adjust < 0.05)
  message(sprintf("Significant pathways (padj < 0.05): %d", nrow(sig_gsea)))
  
  # Save full and significant results
  output_prefix <- sprintf("results/pathway_analysis/%s", project_name)
  write.csv(gsea_df, paste0(output_prefix, "_gsea_all_results.csv"), row.names = FALSE)
  write.csv(sig_gsea, paste0(output_prefix, "_gsea_significant.csv"), row.names = FALSE)
  saveRDS(gsea_res, paste0(output_prefix, "_gsea_object.rds"))
  
  # 6) Generate visualizations if there are significant results
  if (nrow(sig_gsea) > 0) {
    # Bar Chart of significant pathways
    sig_gsea_plot <- sig_gsea %>%
      mutate(Description = gsub("^HALLMARK_", "", ID)) %>%
      arrange(p.adjust)
    
    p_bar <- ggplot(sig_gsea_plot,
                    aes(x = reorder(Description, NES),
                        y = NES,
                        fill = -log10(p.adjust))) +
      geom_col() +
      coord_flip() +
      scale_fill_gradient(name = "-log10(padj)",
                          low = "skyblue", high = "firebrick") +
      labs(title = paste0(project_name, ": Significant Hallmark Pathways (n = ", nrow(sig_gsea), ")"),
           subtitle = "Q4 (accelerated aging) vs Q1 (decelerated aging)",
           x = NULL,
           y = "Normalized Enrichment Score (NES)") +
      theme_minimal(base_size = 13)
    
    ggsave(paste0(output_prefix, "_gsea_barplot.png"), p_bar,
           width = 10, height = max(4, nrow(sig_gsea) * 0.3 + 2), dpi = 300)
    message(sprintf("Saved bar plot for %s", project_name))
    
    # Enrichment plot for top pathway
    top_term <- gsea_df$ID[1]
    ep <- enrichplot::gseaplot2(gsea_res, geneSetID = top_term, title = top_term)
    ggsave(paste0(output_prefix, "_gsea_top_enrichment.png"), ep,
           width = 8, height = 6, dpi = 300)
    message(sprintf("Saved enrichment plot for top pathway: %s", top_term))
    
    # Dot plot if enough pathways
    if (nrow(sig_gsea) >= 3) {
      p_dot <- dotplot(gsea_res, showCategory = min(20, nrow(sig_gsea)), 
                       title = paste0(project_name, ": GSEA Dotplot"))
      ggsave(paste0(output_prefix, "_gsea_dotplot.png"), p_dot,
             width = 10, height = 8, dpi = 300)
      message(sprintf("Saved dot plot for %s", project_name))
    }
  } else {
    message(sprintf("No significant pathways found for %s", project_name))
  }
  
  return(list(
    project = project_name,
    gsea_result = gsea_res,
    all_pathways = gsea_df,
    significant = sig_gsea
  ))
}

# Run GSEA for all projects
message("\n========== Starting Pathway Analysis ==========")
gsea_results <- list()

for (proj in projects) {
  result <- tryCatch(
    run_gsea_analysis(proj),
    error = function(e) {
      warning(sprintf("Error in %s: %s", proj, e$message))
      return(NULL)
    }
  )
  if (!is.null(result)) {
    gsea_results[[proj]] <- result
  }
}

# Save combined results
saveRDS(gsea_results, "results/pathway_analysis/all_gsea_results.rds")

# Print summary
message("\n========== Pathway Analysis Complete ==========")
for (proj in names(gsea_results)) {
  n_sig <- nrow(gsea_results[[proj]]$significant)
  message(sprintf("%s: %d significant Hallmark pathways", proj, n_sig))
}
message(sprintf("\nResults saved to: results/pathway_analysis/"))

