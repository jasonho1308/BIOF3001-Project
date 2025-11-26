# Perform differential expression analysis
# Analyzes BRCA, THCA (within-tissue Q1/Q4) and UCEC (pan-cancer Q1/Q4)

library(DESeq2)
library(ggplot2)
library(dplyr)
library(parallel)

# Detect number of cores (use half to avoid overwhelming the system)
n_cores <- max(1, floor(detectCores() / 2) + 2)
message(sprintf("Using %d cores for parallel processing", n_cores))

# Define projects to analyze
# BRCA and THCA use within-tissue quartiles
# UCEC uses pan-cancer quartiles
tissue_projects <- c("TCGA-BRCA", "TCGA-THCA")
pancan_projects <- c("TCGA-UCEC")
all_projects <- c(tissue_projects, pancan_projects)

# Load clock predictions (contains sample metadata but not counts)
gdc_pancan <- readRDS("results/rna/gdc_pancan_rna_predictions.rds")

if (nrow(gdc_pancan) == 0) {
  stop("No samples found for differential expression analysis!")
}

message(sprintf("Loaded %d samples with clock predictions", nrow(gdc_pancan)))

# Calculate pan-cancer quartiles from ALL samples BEFORE filtering to projects
# This ensures UCEC is compared against the true pan-cancer distribution
pancan_age_acc <- gdc_pancan$Horvath_residuals
pancan_q1 <- quantile(pancan_age_acc, 0.25, na.rm = TRUE)
pancan_q4 <- quantile(pancan_age_acc, 0.75, na.rm = TRUE)
message(sprintf("Pan-cancer Horvath_residuals Q1: %.2f, Q4: %.2f (from %d samples)", 
                pancan_q1, pancan_q4, sum(!is.na(pancan_age_acc))))

# Filter to only the projects of interest
gdc_pancan <- gdc_pancan[gdc_pancan$project %in% all_projects, ]
message(sprintf("Filtered to %d samples from %s", nrow(gdc_pancan), paste(all_projects, collapse = ", ")))

# Load RNA-seq queries to map sample IDs to file paths
rna_queries <- readRDS("data/rna/queries.rds")

# Function to read count data from a single file
read_count_file <- function(file_path) {
  counts <- read.delim(file_path, comment.char = "#", stringsAsFactors = FALSE)
  # Skip the first 4 rows which contain N_unmapped, N_multimapping, N_noFeature, N_ambiguous
  counts <- counts[!grepl("^N_", counts$gene_id), ]
  # Use unstranded counts (column 4)
  count_vec <- counts$unstranded
  names(count_vec) <- counts$gene_id
  return(count_vec)
}

# Build a mapping from sample_submitter_id to file paths (process projects in parallel)
message("Building sample to file path mapping...")

process_project_mapping <- function(proj, rna_queries) {
  query <- rna_queries[[proj]]
  if (is.null(query) || is.null(query$results[[1]])) return(NULL)
  
  results_df <- query$results[[1]]
  mappings <- list()
  
  for (i in seq_len(nrow(results_df))) {
    sample_id <- results_df$sample.submitter_id[i]
    file_id <- results_df$file_id[i]
    file_name <- results_df$file_name[i]
    
    file_path <- file.path("data/rna", proj, "Transcriptome_Profiling", 
                           "Gene_Expression_Quantification", file_id, file_name)
    
    if (file.exists(file_path)) {
      mappings[[sample_id]] <- file_path
    }
  }
  return(mappings)
}

# Process projects in parallel to build mapping
project_mappings <- mclapply(names(rna_queries), function(proj) {
  process_project_mapping(proj, rna_queries)
}, mc.cores = n_cores)

# Combine all project mappings into one list
sample_file_map <- do.call(c, project_mappings)

message(sprintf("Found %d samples with RNA-seq files", length(sample_file_map)))

# Filter predictions to samples with available count files
samples_with_counts <- gdc_pancan$sample_submitter_id %in% names(sample_file_map)
gdc_pancan <- gdc_pancan[samples_with_counts, ]
message(sprintf("After filtering: %d samples with both clock predictions and RNA-seq data", nrow(gdc_pancan)))

if (nrow(gdc_pancan) == 0) {
  stop("No samples found with both clock predictions and RNA-seq count files!")
}

# Read count data for all samples in parallel
message("Reading RNA-seq count files in parallel...")
sample_ids <- gdc_pancan$sample_submitter_id

count_list <- mclapply(sample_ids, function(sample_id) {
  file_path <- sample_file_map[[sample_id]]
  read_count_file(file_path)
}, mc.cores = n_cores)
names(count_list) <- sample_ids

message("Combining counts into matrix...")

# Combine into a count matrix (genes as rows, samples as columns)
all_genes <- unique(unlist(lapply(count_list, names)))
count_matrix <- matrix(0, nrow = length(all_genes), ncol = length(count_list),
                       dimnames = list(all_genes, names(count_list)))

for (sample_id in names(count_list)) {
  counts <- count_list[[sample_id]]
  count_matrix[names(counts), sample_id] <- counts
}

message(sprintf("Count matrix: %d genes x %d samples", nrow(count_matrix), ncol(count_matrix)))

# Function to run DE analysis for a specific project
run_de_analysis <- function(project_name, samples_df, count_mat, use_pancan_quartiles = FALSE) {
  message(sprintf("\n========== Processing %s ==========", project_name))
  
  # Filter to this project
  proj_samples <- samples_df[samples_df$project == project_name, ]
  proj_count_mat <- count_mat[, proj_samples$sample_submitter_id, drop = FALSE]
  
  message(sprintf("Project has %d samples", nrow(proj_samples)))
  
  # Calculate quartiles using age acceleration residuals
  # Use ageAcc2.Horvath for within-tissue, Horvath_residuals for pan-cancer
  if (use_pancan_quartiles) {
    age_acc <- proj_samples$Horvath_residuals
  } else {
    age_acc <- proj_samples$ageAcc2.Horvath
  }
  
  if (use_pancan_quartiles) {
    # Use pan-cancer quartiles
    q1 <- pancan_q1
    q4 <- pancan_q4
    message(sprintf("Using pan-cancer quartiles - Q1: %.2f, Q4: %.2f", q1, q4))
  } else {
    # Use within-tissue quartiles
    q1 <- quantile(age_acc, 0.25, na.rm = TRUE)
    q4 <- quantile(age_acc, 0.75, na.rm = TRUE)
    message(sprintf("Using within-tissue quartiles - Q1: %.2f, Q4: %.2f", q1, q4))
  }
  
  # Classify samples as Q1 (decelerated) or Q4 (accelerated)
  # Q1 = lower age acceleration (younger than expected)
  # Q4 = higher age acceleration (older than expected)
  age_acc_class <- ifelse(age_acc <= q1, "Q1",
                          ifelse(age_acc >= q4, "Q4", NA))
  
  # Filter to only Q1 and Q4 samples
  keep_samples <- !is.na(age_acc_class)
  proj_samples <- proj_samples[keep_samples, ]
  proj_count_mat <- proj_count_mat[, proj_samples$sample_submitter_id, drop = FALSE]
  age_acc_class <- age_acc_class[keep_samples]
  
  message(sprintf("After Q1/Q4 filtering: %d samples (Q1: %d, Q4: %d)", 
                  nrow(proj_samples), 
                  sum(age_acc_class == "Q1"), 
                  sum(age_acc_class == "Q4")))
  
  if (nrow(proj_samples) < 10) {
    warning(sprintf("Too few samples for %s, skipping...", project_name))
    return(NULL)
  }
  
  # Filter low count genes (keep genes with rowSums >= 10)
  proj_count_mat <- proj_count_mat[which(rowSums(proj_count_mat) >= 10), ]
  message(sprintf("After filtering low count genes: %d genes", nrow(proj_count_mat)))
  
  # Prepare conditions data frame (simple design, no age covariate)
  conditions <- data.frame(
    row.names = proj_samples$sample_submitter_id,
    Type = factor(age_acc_class, levels = c("Q1", "Q4"))
  )
  
  # Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(
    countData = proj_count_mat,
    colData = conditions,
    design = ~ Type
  )
  
  # Run DESeq2 analysis
  dds <- DESeq(dds)
  vsd <- vst(dds, blind = FALSE)
  
  # Create output file names
  proj_short <- gsub("TCGA-", "", project_name)
  output_prefix <- paste0("results/differential_analysis/", proj_short)
  
  # PCA plot
  pca_plot <- plotPCA(vsd, intgroup = c("Type"))
  pca_plot <- pca_plot + 
    scale_color_manual(values = c("Q1" = "grey60", "Q4" = "red")) +
    ggtitle(paste0("PCA: ", project_name)) +
    theme_minimal()
  ggsave(paste0(output_prefix, "_pca_plot.png"), plot = pca_plot, width = 8, height = 6)
  message(sprintf("PCA plot saved for %s", project_name))
  
  # Generate DE results
  res <- results(dds, contrast = c("Type", "Q4", "Q1"))
  res_df <- as.data.frame(res)
  res_sig <- subset(res_df, !is.na(padj) & padj < 0.05)
  res_sig <- res_sig[order(res_sig$padj), ]
  message(sprintf("Significant genes (padj < 0.05): %d", nrow(res_sig)))
  
  # Save results
  saveRDS(res, paste0(output_prefix, "_deseq2_results.rds"))
  write.csv(res_sig, paste0(output_prefix, "_significant_genes.csv"))
  message(sprintf("DESeq2 results saved for %s", project_name))
  
  # Volcano plot with color coding for significance
  volcano_df <- as.data.frame(res)
  volcano_df$gene <- rownames(volcano_df)
  
  log2FC_cutoff <- 1
  padj_cutoff <- 0.05
  
  volcano_df$color <- "Not significant"
  volcano_df$color[volcano_df$padj < padj_cutoff & abs(volcano_df$log2FoldChange) > log2FC_cutoff] <- "Significant"
  
  volcano_plot <- ggplot(volcano_df, aes(x = log2FoldChange, y = -log10(padj), color = color)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Not significant" = "grey", "Significant" = "red")) +
    geom_vline(xintercept = c(-log2FC_cutoff, log2FC_cutoff), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
    labs(title = paste0("Volcano Plot: ", project_name, " (Q4 vs Q1)"),
         x = expression(Log[2]~Fold~Change), 
         y = expression(-Log[10]~adjusted~p)) +
    theme_minimal()
  ggsave(paste0(output_prefix, "_volcano_plot.png"), plot = volcano_plot, width = 8, height = 6)
  message(sprintf("Volcano plot saved for %s", project_name))
  
  return(list(dds = dds, res = res, res_sig = res_sig, project = project_name))
}

# Run DE analysis for each project
results_list <- list()

# BRCA and THCA: within-tissue quartiles
for (proj in tissue_projects) {
  result <- run_de_analysis(proj, gdc_pancan, count_matrix, use_pancan_quartiles = FALSE)
  if (!is.null(result)) {
    results_list[[proj]] <- result
  }
}

# UCEC: pan-cancer quartiles
for (proj in pancan_projects) {
  result <- run_de_analysis(proj, gdc_pancan, count_matrix, use_pancan_quartiles = TRUE)
  if (!is.null(result)) {
    results_list[[proj]] <- result
  }
}

# Save combined results summary
saveRDS(results_list, "results/differential_analysis/all_de_results.rds")
message("\n========== All differential expression analyses completed ==========")
message(sprintf("Successfully analyzed: %s", paste(names(results_list), collapse = ", ")))
