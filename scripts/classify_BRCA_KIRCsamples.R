#!/usr/bin/env Rscript
# Simplified script for classifying and counting samples based on global epigenetic aging
# for specific tissues (TCGA-BRCA and TCGA-KIRC)

library(dplyr)

# Configuration
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
    clock_name <- args[1]
} else {
    clock_name <- "Horvath" # Default clock
}

message(sprintf("Classifying samples based on %s clock", clock_name))

# Load predictions
predictions <- readRDS("results/clock/gdc_pan/gdc_pancan_predictions.rds")
message(sprintf("Loaded %d samples with predictions", nrow(predictions)))

# Check if clock exists
if (!clock_name %in% colnames(predictions)) {
    stop(sprintf(
        "Clock '%s' not found in predictions. Available: %s",
        clock_name, paste(colnames(predictions), collapse = ", ")
    ))
}

# Apply global quartile-based classification using the entire dataset
residuals_vec <- predictions[[paste0(clock_name, "_residuals")]]
qs <- stats::quantile(residuals_vec, probs = c(0.25, 0.75), na.rm = TRUE)

# Classifying samples based on the entire dataset
predictions$acceleration_status <- NA_character_
predictions$acceleration_status[!is.na(residuals_vec) & residuals_vec < qs[1]] <- "Accelerated"
predictions$acceleration_status[!is.na(residuals_vec) & residuals_vec > qs[2]] <- "Decelerated"

# Filter for specific tissues AFTER classification
predictions_filtered <- predictions %>%
    filter(project %in% c("TCGA-BRCA", "TCGA-KIRC"))

# Count the samples for each tissue type
results <- predictions_filtered %>%
    group_by(project, acceleration_status) %>%
    summarise(count = n(), .groups = 'drop')

# Print results
for (tissue in unique(results$project)) {
    tissue_counts <- results %>% filter(project == tissue)
    
    # Handle counts safely
    acc_count <- ifelse(any(tissue_counts$acceleration_status == "Accelerated"), 
                        tissue_counts$count[tissue_counts$acceleration_status == "Accelerated"], 
                        0)
    
    dec_count <- ifelse(any(tissue_counts$acceleration_status == "Decelerated"), 
                        tissue_counts$count[tissue_counts$acceleration_status == "Decelerated"], 
                        0)
    
    cat(sprintf("\nTissue: %s\n", tissue))
    cat(sprintf("Accelerated samples: %d\n", acc_count))
    cat(sprintf("Decelerated samples: %d\n", dec_count))
}

# Set the relative output directory
output_dir <- "results/differential_analysis"

# Prepare sample lists for output
sample_list <- character()

# Ensure the directory exists
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

for (tissue in unique(results$project)) {
    tissue_counts <- results %>% filter(project == tissue)
    
    # Get accelerated and decelerated samples
    accelerated_samples <- predictions_filtered %>%
        filter(project == tissue, acceleration_status == "Accelerated") %>%
        pull(id)
    
    decelerated_samples <- predictions_filtered %>%
        filter(project == tissue, acceleration_status == "Decelerated") %>%
        pull(id)
    
    # Format and store the samples
    if (length(accelerated_samples) > 0) {
        sample_list <- c(sample_list, paste(tissue, "acc", sep = "_"))
        # Write to the specified relative output directory
        write.table(accelerated_samples, 
                    file = file.path(output_dir, paste0(tissue, "_acc.txt")), 
                    row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
    
    if (length(decelerated_samples) > 0) {
        sample_list <- c(sample_list, paste(tissue, "dec", sep = "_"))
        # Write to the specified relative output directory
        write.table(decelerated_samples, 
                    file = file.path(output_dir, paste0(tissue, "_dec.txt")), 
                    row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
}

# Print summary of stored samples
for (sample in sample_list) {
    cat(sprintf("Stored samples: %s\n", sample))
}

message("=== Classification Complete! ===")