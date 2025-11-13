#!/usr/bin/env Rscript
# Clinical correlation analysis for accelerated/decelerated epigenetic aging
# Treats all variables as binary (yes/no) and computes associations

library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(gridExtra)

# Configuration
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  clock_name <- args[1]
} else {
  clock_name <- "Horvath" # Default clock
}

message(sprintf("Analyzing clinical correlations for %s clock", clock_name))

# Create output directory
if (!dir.exists("results/clinical")) {
  dir.create("results/clinical", recursive = TRUE)
}

# ============================================================================
# Load data
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

# ============================================================================
# Calculate age acceleration (binary: accelerated vs decelerated)
# ============================================================================

if (!clock_name %in% colnames(predictions)) {
  stop(sprintf(
    "Clock '%s' not found in predictions. Available: %s",
    clock_name, paste(colnames(predictions), collapse = ", ")
  ))
}

# Calculate residuals (age acceleration)
predictions$age_acceleration <- NA
ok_idx <- !is.na(predictions[[clock_name]]) & !is.na(predictions$age)

if (sum(ok_idx) >= 2) {
  model <- lm(predictions[[clock_name]][ok_idx] ~ predictions$age[ok_idx])
  predictions$age_acceleration[ok_idx] <- resid(model)
}

# Binarize: accelerated (positive residual) vs decelerated (negative residual)
predictions$accelerated <- ifelse(predictions$age_acceleration > 0, "Accelerated", "Decelerated")

message(sprintf(
  "Age acceleration: %d accelerated, %d decelerated, %d missing",
  sum(predictions$accelerated == "Accelerated", na.rm = TRUE),
  sum(predictions$accelerated == "Decelerated", na.rm = TRUE),
  sum(is.na(predictions$accelerated))
))

# ============================================================================
# Merge predictions with phenotype
# ============================================================================

# Predictions might have different sample ID format; try to match
# If predictions have sample_id column, use it; otherwise use rownames
if ("sample_id" %in% colnames(predictions)) {
  predictions$sample_barcode <- substr(predictions$sample_id, 1, 15)
} else {
  predictions$sample_barcode <- substr(rownames(predictions), 1, 15)
}

# Merge predictions with phenotype
merged <- predictions %>%
  select(sample, age, age_acceleration, accelerated, all_of(clock_name)) %>%
  inner_join(pheno, by = "sample")

message(sprintf("Merged data: %d samples", nrow(merged)))

# ============================================================================
# Binarize clinical variables
# ============================================================================

binarize_variable <- function(x) {
  # Convert to character
  x <- as.character(x)

  # Count unique non-NA values
  unique_vals <- unique(x[!is.na(x)])

  # If already binary-like (2 unique values), map to yes/no
  if (length(unique_vals) == 2) {
    # Common patterns
    if (all(tolower(unique_vals) %in% c("yes", "no"))) {
      return(ifelse(tolower(x) == "yes", "Yes", "No"))
    }
    if (all(tolower(unique_vals) %in% c("male", "female"))) {
      return(ifelse(tolower(x) == "male", "Male", "Female"))
    }
    if (all(tolower(unique_vals) %in% c("alive", "dead"))) {
      return(ifelse(tolower(x) == "alive", "Alive", "Dead"))
    }
    # Generic binary: just use the values as-is
    return(x)
  }

  # If numeric, split at median
  if (all(!is.na(as.numeric(x[!is.na(x)])))) {
    x_num <- as.numeric(x)
    if (sum(!is.na(x_num)) >= 2) {
      med <- median(x_num, na.rm = TRUE)
      return(ifelse(x_num > med, "High", "Low"))
    }
  }

  # If more than 2 categories, group rare ones as "Other"
  if (length(unique_vals) > 2) {
    # Find most common category
    tab <- table(x, useNA = "no")
    top_cat <- names(tab)[which.max(tab)]
    return(ifelse(x == top_cat, top_cat, "Other"))
  }

  # Otherwise return as-is
  return(x)
}

# Select clinical variables of interest
clinical_vars <- c(
  "demographic.gender",
  "demographic.race",
  "demographic.ethnicity",
  "demographic.vital_status",
  "diagnoses.prior_malignancy",
  "diagnoses.prior_treatment",
  "diagnoses.progression_or_recurrence",
  "diagnoses.tumor_grade",
  "diagnoses.tumor_stage",
  "exposures.alcohol_history",
  "exposures.bmi"
)

# Keep only variables that exist
clinical_vars <- clinical_vars[clinical_vars %in% colnames(merged)]

message(sprintf("Binarizing %d clinical variables", length(clinical_vars)))

# Binarize each variable
binarized_data <- merged %>%
  select(sample_barcode, accelerated, all_of(clinical_vars))

for (var in clinical_vars) {
  binarized_data[[var]] <- binarize_variable(merged[[var]])
}

# ============================================================================
# Statistical tests (Fisher's exact test for binary associations)
# ============================================================================

test_association <- function(outcome, predictor, outcome_name, predictor_name) {
  # Remove NAs
  complete_idx <- !is.na(outcome) & !is.na(predictor)

  if (sum(complete_idx) < 10) {
    return(NULL) # Not enough data
  }

  outcome <- outcome[complete_idx]
  predictor <- predictor[complete_idx]

  # Create contingency table
  tbl <- table(predictor, outcome)

  # Skip if any dimension is < 2
  if (nrow(tbl) < 2 || ncol(tbl) < 2) {
    return(NULL)
  }

  # Fisher's exact test
  test_result <- tryCatch(
    fisher.test(tbl),
    error = function(e) NULL
  )

  if (is.null(test_result)) {
    return(NULL)
  }

  # Calculate odds ratio (if 2x2 table)
  odds_ratio <- NA
  if (all(dim(tbl) == c(2, 2))) {
    odds_ratio <- (tbl[1, 1] * tbl[2, 2]) / (tbl[1, 2] * tbl[2, 1])
  }

  return(data.frame(
    clinical_variable = predictor_name,
    n_samples = sum(complete_idx),
    p_value = test_result$p.value,
    odds_ratio = odds_ratio,
    stringsAsFactors = FALSE
  ))
}

# Run tests for all clinical variables
results_list <- list()

for (var in clinical_vars) {
  result <- test_association(
    outcome = binarized_data$accelerated,
    predictor = binarized_data[[var]],
    outcome_name = "accelerated",
    predictor_name = var
  )

  if (!is.null(result)) {
    results_list[[var]] <- result
  }
}

# Combine results
results <- bind_rows(results_list)

# Adjust p-values for multiple testing (Bonferroni)
if (nrow(results) > 0) {
  results$p_adj <- p.adjust(results$p_value, method = "bonferroni")
  results$significant <- results$p_adj < 0.05

  # Sort by p-value
  results <- results %>% arrange(p_value)

  message("\n=== Clinical Correlation Results ===")
  print(results)

  # Save results
  write.csv(results,
    file = sprintf("results/clinical/%s_clinical_associations.csv", clock_name),
    row.names = FALSE
  )
  message(sprintf("\nResults saved to results/clinical/%s_clinical_associations.csv", clock_name))
}

# ============================================================================
# Visualizations
# ============================================================================

# 1. Odds ratio plot
if (nrow(results) > 0 && any(!is.na(results$odds_ratio))) {
  plot_data <- results %>%
    filter(!is.na(odds_ratio)) %>%
    mutate(
      clinical_variable = gsub("^(demographic|diagnoses|exposures)\\.", "", clinical_variable),
      log_odds_ratio = log(odds_ratio),
      sig_label = ifelse(significant, "*", "")
    )

  if (nrow(plot_data) > 0) {
    p1 <- ggplot(plot_data, aes(x = reorder(clinical_variable, odds_ratio), y = odds_ratio)) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
      geom_point(aes(color = significant), size = 4) +
      geom_errorbar(aes(ymin = odds_ratio * 0.8, ymax = odds_ratio * 1.2), width = 0.2) +
      geom_text(aes(label = sig_label), vjust = -1, size = 6) +
      scale_color_manual(values = c("FALSE" = "gray50", "TRUE" = "red")) +
      labs(
        title = sprintf("Clinical Associations with Age Acceleration (%s)", clock_name),
        subtitle = "Odds ratio for accelerated vs decelerated aging",
        x = "Clinical Variable",
        y = "Odds Ratio",
        color = "Significant\n(p < 0.05)"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    ggsave(
      filename = sprintf("results/clinical/%s_odds_ratio_plot.png", clock_name),
      plot = p1,
      width = 10,
      height = 6
    )
    message(sprintf("Odds ratio plot saved to results/clinical/%s_odds_ratio_plot.png", clock_name))
  }
}

# 2. Contingency table heatmaps
create_contingency_heatmap <- function(data, var_name) {
  tbl <- table(data[[var_name]], data$accelerated)

  if (nrow(tbl) < 2 || ncol(tbl) < 2) {
    return(NULL)
  }

  # Convert to proportions (row-wise)
  tbl_prop <- prop.table(tbl, margin = 1)

  # Convert to data frame for ggplot
  df <- as.data.frame(tbl_prop)
  colnames(df) <- c("Variable", "Acceleration", "Proportion")

  p <- ggplot(df, aes(x = Acceleration, y = Variable, fill = Proportion)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", Proportion)), color = "black") +
    scale_fill_gradient(low = "white", high = "steelblue") +
    labs(
      title = sprintf("%s vs Age Acceleration", gsub("^(demographic|diagnoses|exposures)\\.", "", var_name)),
      x = "Epigenetic Aging",
      y = ""
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 0))

  return(p)
}

# Create heatmaps for top 6 variables
top_vars <- head(results$clinical_variable, 6)
heatmap_plots <- list()

for (var in top_vars) {
  if (var %in% colnames(binarized_data)) {
    p <- create_contingency_heatmap(binarized_data, var)
    if (!is.null(p)) {
      heatmap_plots[[var]] <- p
    }
  }
}

if (length(heatmap_plots) > 0) {
  combined_heatmaps <- do.call(grid.arrange, c(heatmap_plots, ncol = 2))

  ggsave(
    filename = sprintf("results/clinical/%s_contingency_heatmaps.png", clock_name),
    plot = combined_heatmaps,
    width = 12,
    height = 8
  )
  message(sprintf("Contingency heatmaps saved to results/clinical/%s_contingency_heatmaps.png", clock_name))
}

# 3. Overall summary heatmap (all variables)
if (nrow(results) > 0) {
  # Create matrix of -log10(p-values)
  heatmap_data <- results %>%
    mutate(
      neg_log_p = -log10(p_value),
      clinical_variable = gsub("^(demographic|diagnoses|exposures)\\.", "", clinical_variable)
    ) %>%
    select(clinical_variable, neg_log_p)

  # Convert to matrix
  mat <- matrix(heatmap_data$neg_log_p, ncol = 1)
  rownames(mat) <- heatmap_data$clinical_variable
  colnames(mat) <- clock_name

  # Create heatmap
  png(sprintf("results/clinical/%s_pvalue_heatmap.png", clock_name),
    width = 600, height = 800
  )
  pheatmap(mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    display_numbers = TRUE,
    number_format = "%.2f",
    main = sprintf("Clinical Associations with %s Age Acceleration\n(-log10 p-value)", clock_name),
    color = colorRampPalette(c("white", "orange", "red"))(50),
    fontsize_row = 10,
    fontsize_col = 12
  )
  dev.off()
  message(sprintf("P-value heatmap saved to results/clinical/%s_pvalue_heatmap.png", clock_name))
}

message("\n=== Analysis complete! ===")
