#!/usr/bin/env Rscript
# Tissue type distribution analysis for accelerated vs decelerated epigenetic aging
# Shows distribution of each tissue type in Q1 (accelerated) vs Q4 (decelerated) quartiles

library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)

# Configuration
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  clock_name <- args[1]
} else {
  clock_name <- "Horvath" # Default clock
}

message(sprintf("Analyzing tissue type distribution for %s clock", clock_name))

# Create output directory
if (!dir.exists("results/tissue_analysis")) {
  dir.create("results/tissue_analysis", recursive = TRUE)
}

# ============================================================================
# Load and process data
# ============================================================================

# Load predictions
predictions <- readRDS("results/clock/gdc_pan/gdc_pancan_predictions.rds")
message(sprintf("Loaded %d samples with predictions", nrow(predictions)))

# Load phenotype data
pheno <- read.table("data/processed/gdc_pancan/normal_pheno.tsv",
  header = TRUE, sep = "\t", quote = "",
  stringsAsFactors = FALSE, na.strings = c("NA", "", "Not Reported", "not reported")
)

message(sprintf("Loaded phenotype data for %d samples", nrow(pheno)))

# Check if clock exists
if (!clock_name %in% colnames(predictions)) {
  stop(sprintf(
    "Clock '%s' not found in predictions. Available: %s",
    clock_name, paste(colnames(predictions), collapse = ", ")
  ))
}

# ============================================================================
# Apply quartile-based classification
# ============================================================================

# Classify acceleration by quartiles of residuals (same as clinical_correlation.R)
# Q1 (lowest 25%) -> Accelerated, Q4 (highest 25%) -> Decelerated
# Middle 50% are set to NA and excluded from analysis
residuals_vec <- predictions[[paste0(clock_name, "_residuals")]]
qs <- stats::quantile(residuals_vec, probs = c(0.25, 0.75), na.rm = TRUE)

predictions$acceleration_status <- NA_character_
predictions$acceleration_status[!is.na(residuals_vec) & residuals_vec <= qs[1]] <- "Accelerated"
predictions$acceleration_status[!is.na(residuals_vec) & residuals_vec >= qs[2]] <- "Decelerated"

message(sprintf(
  "Quartile classification: Q1=%.3f, Q4=%.3f | %d accelerated, %d decelerated, %d excluded",
  qs[1], qs[2],
  sum(predictions$acceleration_status == "Accelerated", na.rm = TRUE),
  sum(predictions$acceleration_status == "Decelerated", na.rm = TRUE),
  sum(is.na(predictions$acceleration_status))
))

# ============================================================================
# Merge with phenotype data
# ============================================================================

# Rename id to sample and match sample IDs
predictions <- predictions %>%
  rename(sample = id) %>%
  mutate(sample = substr(sample, 1, 16))

# Merge with phenotype data
merged <- predictions %>%
  select(sample, age, paste0(clock_name, "_residuals"), acceleration_status, all_of(clock_name)) %>%
  inner_join(pheno, by = "sample")

message(sprintf("Merged data: %d samples", nrow(merged)))

# Filter to only accelerated and decelerated samples (exclude middle 50%)
extreme_samples <- merged %>%
  filter(!is.na(acceleration_status))

message(sprintf("Analysis dataset: %d extreme samples (%d accelerated, %d decelerated)",
  nrow(extreme_samples),
  sum(extreme_samples$acceleration_status == "Accelerated"),
  sum(extreme_samples$acceleration_status == "Decelerated")
))

# ============================================================================
# Tissue type analysis
# ============================================================================

# Identify tissue type columns (project_id and sample_type are main indicators)
tissue_vars <- c("project.project_id", "samples.sample_type")

# Keep only existing columns
tissue_vars <- tissue_vars[tissue_vars %in% colnames(extreme_samples)]

if (length(tissue_vars) == 0) {
  stop("No tissue type variables found in data. Available columns: ", 
       paste(colnames(extreme_samples), collapse = ", "))
}

message(sprintf("Found %d tissue type variables: %s", 
                length(tissue_vars), paste(tissue_vars, collapse = ", ")))

# ============================================================================
# Analysis 1: Project-wise distribution
# ============================================================================

if ("project.project_id" %in% colnames(extreme_samples)) {
  project_summary <- extreme_samples %>%
    group_by(project.project_id, acceleration_status) %>%
    summarise(count = n(), .groups = "drop") %>%
    pivot_wider(names_from = acceleration_status, values_from = count, values_fill = 0) %>%
    mutate(
      total = Accelerated + Decelerated,
      pct_accelerated = round(100 * Accelerated / total, 1),
      pct_decelerated = round(100 * Decelerated / total, 1),
      ratio_acc_dec = round(Accelerated / pmax(Decelerated, 1), 2)
    ) %>%
    arrange(desc(total))
  
  # Save project summary
  write.csv(project_summary, 
    file = sprintf("results/tissue_analysis/%s_project_distribution.csv", clock_name),
    row.names = FALSE
  )
  
  message(sprintf("\nProject distribution summary (top 10):"))
  print(head(project_summary, 10))
  
  # Create project distribution plot
  plot_data <- project_summary %>%
    filter(total >= 5) %>%  # Only projects with >=5 extreme samples
    pivot_longer(cols = c(Accelerated, Decelerated), 
                 names_to = "status", values_to = "count") %>%
    mutate(
      project.project_id = reorder(project.project_id, count, sum),
      status = factor(status, levels = c("Accelerated", "Decelerated"))
    )
  
  p1 <- ggplot(plot_data, aes(x = project.project_id, y = count, fill = status)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Accelerated" = "#e74c3c", "Decelerated" = "#3498db")) +
    labs(
      title = sprintf("Tissue Distribution: %s Clock Acceleration", clock_name),
      subtitle = "Count of samples in Q1 (Accelerated) vs Q4 (Decelerated) by project",
      x = "TCGA Project",
      y = "Sample Count",
      fill = "Aging Status"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top"
    )
  
  ggsave(sprintf("results/tissue_analysis/%s_project_distribution.png", clock_name),
         plot = p1, width = 12, height = 8)
  
  # Create percentage stacked plot
  p2 <- plot_data %>%
    group_by(project.project_id) %>%
    mutate(percentage = 100 * count / sum(count)) %>%
    ggplot(aes(x = project.project_id, y = percentage, fill = status)) +
    geom_col() +
    scale_fill_manual(values = c("Accelerated" = "#e74c3c", "Decelerated" = "#3498db")) +
    labs(
      title = sprintf("Tissue Distribution: %s Clock Acceleration (Percentages)", clock_name),
      subtitle = "Proportion of Q1 (Accelerated) vs Q4 (Decelerated) within each project",
      x = "TCGA Project",
      y = "Percentage",
      fill = "Aging Status"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top"
    )
  
  ggsave(sprintf("results/tissue_analysis/%s_project_percentage.png", clock_name),
         plot = p2, width = 12, height = 8)
}

# ============================================================================
# Analysis 2: Sample type distribution (if available)
# ============================================================================

if ("samples.sample_type" %in% colnames(extreme_samples)) {
  sample_type_summary <- extreme_samples %>%
    filter(!is.na(samples.sample_type)) %>%
    group_by(samples.sample_type, acceleration_status) %>%
    summarise(count = n(), .groups = "drop") %>%
    pivot_wider(names_from = acceleration_status, values_from = count, values_fill = 0) %>%
    mutate(
      total = Accelerated + Decelerated,
      pct_accelerated = round(100 * Accelerated / total, 1),
      pct_decelerated = round(100 * Decelerated / total, 1),
      ratio_acc_dec = round(Accelerated / pmax(Decelerated, 1), 2)
    ) %>%
    arrange(desc(total))
  
  # Save sample type summary
  write.csv(sample_type_summary, 
    file = sprintf("results/tissue_analysis/%s_sample_type_distribution.csv", clock_name),
    row.names = FALSE
  )
  
  message(sprintf("\nSample type distribution:"))
  print(sample_type_summary)
  
  # Create sample type plot if there are multiple types
  if (nrow(sample_type_summary) > 1) {
    plot_data_st <- sample_type_summary %>%
      pivot_longer(cols = c(Accelerated, Decelerated), 
                   names_to = "status", values_to = "count") %>%
      mutate(
        samples.sample_type = reorder(samples.sample_type, count, sum),
        status = factor(status, levels = c("Accelerated", "Decelerated"))
      )
    
    p3 <- ggplot(plot_data_st, aes(x = samples.sample_type, y = count, fill = status)) +
      geom_col(position = "dodge") +
      scale_fill_manual(values = c("Accelerated" = "#e74c3c", "Decelerated" = "#3498db")) +
      labs(
        title = sprintf("Sample Type Distribution: %s Clock", clock_name),
        subtitle = "Count by sample type (Q1 vs Q4)",
        x = "Sample Type",
        y = "Sample Count",
        fill = "Aging Status"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top"
      )
    
    ggsave(sprintf("results/tissue_analysis/%s_sample_type_distribution.png", clock_name),
           plot = p3, width = 10, height = 6)
  }
}

# ============================================================================
# Statistical tests
# ============================================================================

# Test if tissue types have significantly different acceleration patterns
if ("project.project_id" %in% colnames(extreme_samples)) {
  # Chi-square test for project association
  project_table <- table(extreme_samples$project.project_id, extreme_samples$acceleration_status)
  
  # Only test if we have sufficient data
  if (all(dim(project_table) >= 2) && sum(project_table) >= 20) {
    chi_test <- chisq.test(project_table)
    
    message(sprintf("\nChi-square test for project association:"))
    message(sprintf("X-squared = %.3f, df = %d, p-value = %.3g", 
                    chi_test$statistic, chi_test$parameter, chi_test$p.value))
    
    if (chi_test$p.value < 0.05) {
      message("Significant association between tissue type and acceleration status (p < 0.05)")
    } else {
      message("No significant association between tissue type and acceleration status (p >= 0.05)")
    }
    
    # Save test results
    test_results <- data.frame(
      test = "chi_square_project_acceleration",
      statistic = chi_test$statistic,
      df = chi_test$parameter,
      p_value = chi_test$p.value,
      significant = chi_test$p.value < 0.05,
      stringsAsFactors = FALSE
    )
    
    write.csv(test_results,
      file = sprintf("results/tissue_analysis/%s_statistical_tests.csv", clock_name),
      row.names = FALSE
    )
  }
}

message(sprintf("\n=== Tissue distribution analysis complete! ==="))
message(sprintf("Results saved to results/tissue_analysis/%s_*", clock_name))