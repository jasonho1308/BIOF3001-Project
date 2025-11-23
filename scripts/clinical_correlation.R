#!/usr/bin/env Rscript
# Clinical correlation analysis for accelerated/decelerated epigenetic aging
# Treats all variables as binary (yes/no) and computes associations

library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(gridExtra)
library(png)
library(grid)
library(stats)

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

# Classify acceleration by quartiles of residuals
# Q1 (lowest 25%) -> Accelerated, Q4 (highest 25%) -> Decelerated
# Middle 50% are set to NA and excluded from association tests/plots
residuals_vec <- predictions[[paste0(clock_name, "_residuals")]]
qs <- stats::quantile(residuals_vec, probs = c(0.25, 0.75), na.rm = TRUE)

predictions$accelerated <- NA_character_
predictions$accelerated[!is.na(residuals_vec) & residuals_vec <= qs[1]] <- "Accelerated"
predictions$accelerated[!is.na(residuals_vec) & residuals_vec >= qs[2]] <- "Decelerated"

predictions$accelerated <- factor(predictions$accelerated, levels = c("Accelerated", "Decelerated"))

message(sprintf(
  "Age acceleration (quartiles): %d accelerated (Q1), %d decelerated (Q4), %d mid/NA",
  sum(predictions$accelerated == "Accelerated", na.rm = TRUE),
  sum(predictions$accelerated == "Decelerated", na.rm = TRUE),
  sum(is.na(predictions$accelerated))
))

# ============================================================================
# Merge predictions with phenotype
# ============================================================================


# rename id to sample in predictions
predictions <- predictions %>%
  rename(sample = id)

# select first 15 characters of sample_barcode in pheno to match sample ids
predictions <- predictions %>%
  mutate(sample = substr(sample, 1, 16))

# Merge predictions with phenotype
merged <- predictions %>%
  select(sample, age, paste0(clock_name, "_residuals"), accelerated, all_of(clock_name)) %>%
  inner_join(pheno, by = "sample")

if (interactive()) View(merged)

message(sprintf("Merged data: %d samples", nrow(merged)))

# Derive a binary obesity flag from BMI > 30
if ("exposures.bmi" %in% colnames(merged)) {
  bmi_values <- suppressWarnings(as.numeric(merged[["exposures.bmi"]]))
  merged$exposures.obesity <- ifelse(is.na(bmi_values), NA_character_,
    ifelse(bmi_values > 30, "Obese", "Not Obese")
  )
} else {
  merged$exposures.obesity <- NA_character_
}

# ============================================================================
# Binarize clinical variables
# ============================================================================

# Helper function to check if values are numeric strings
# Normalize input by removing placeholder values
normalize_input <- function(x) {
  x[x %in% c("", "Not Reported", "not reported", "Not Available", "Unknown")] <- NA
  x
}

# Binarize numerical variables via median split
binarize_numerical <- function(x) {
  x <- normalize_input(x)
  x_chr <- as.character(x)
  x_num <- suppressWarnings(as.numeric(x_chr))

  # Check if we have enough valid numeric values
  if (sum(!is.na(x_num)) < 2) {
    rep(NA_character_, length(x_chr))
  }

  # Median split
  med <- median(x_num, na.rm = TRUE)
  result <- ifelse(is.na(x_num), NA, ifelse(x_num > med, "High", "Low"))
  result
}

# Binarize categorical variables (keep all categories as-is)
binarize_categorical <- function(x) {
  x <- normalize_input(x)
  x_chr <- as.character(x)
  x_chr
}

# Binarize multi-category variables (top category vs Other)
binarize_multicategory <- function(x) {
  x <- normalize_input(x)
  x_chr <- as.character(x)
  unique_vals <- unique(x_chr[!is.na(x_chr)])

  # If too few categories, treat as categorical
  if (length(unique_vals) <= 2) {
    binarize_categorical(x)
  }

  # Many categories - keep most frequent, others -> "Other"
  tab <- sort(table(x_chr[!is.na(x_chr)]), decreasing = TRUE)
  top_cat <- names(tab)[1]

  result <- ifelse(is.na(x_chr), NA,
    ifelse(x_chr == top_cat, top_cat, "Other")
  )
  result
}

# Categorize clinical variables by type
# Numerical/continuous variables
numerical_vars <- c(
  "demographic.age_at_index",
  "demographic.days_to_birth",
  "demographic.days_to_death",
  "demographic.year_of_birth",
  "demographic.year_of_death",
  "diagnoses.age_at_diagnosis",
  "diagnoses.days_to_diagnosis",
  "diagnoses.days_to_last_follow_up",
  "diagnoses.year_of_diagnosis",
  "exposures.cigarettes_per_day",
  "exposures.height",
  "exposures.pack_years_smoked",
  "exposures.weight",
  "exposures.years_smoked"
)

# Binary/categorical variables (yes/no, alive/dead, etc.)
categorical_vars <- c(
  "demographic.ethnicity",
  "demographic.gender",
  "demographic.race",
  "demographic.vital_status",
  "diagnoses.prior_malignancy",
  "diagnoses.prior_treatment",
  "diagnoses.progression_or_recurrence",
  "diagnoses.synchronous_malignancy",
  "diagnoses.tumor_grade",
  "diagnoses.tumor_stage",
  "exposures.alcohol_history",
  "exposures.obesity"
)

# Multi-category variables (diagnosis codes, tissue types, etc.)
multicategory_vars <- c(
  "diagnoses.classification_of_tumor",
  "diagnoses.icd_10_code",
  "diagnoses.last_known_disease_status",
  "diagnoses.morphology",
  "diagnoses.primary_diagnosis",
  "diagnoses.site_of_resection_or_biopsy",
  "diagnoses.tissue_or_organ_of_origin",
  "samples.sample_type"
)


# Combine useful clinical variables
clinical_vars <- c(numerical_vars, categorical_vars, multicategory_vars)

# Keep only variables that exist in merged data
clinical_vars <- clinical_vars[clinical_vars %in% colnames(merged)]

message(sprintf("Selected %d clinical variables for analysis:", length(clinical_vars)))
message(sprintf("  - Numerical: %d", sum(numerical_vars %in% clinical_vars)))
message(sprintf("  - Categorical: %d", sum(categorical_vars %in% clinical_vars)))
message(sprintf("  - Multi-category: %d", sum(multicategory_vars %in% clinical_vars)))

message(sprintf("Processing %d clinical variables", length(clinical_vars)))

# Process each variable type appropriately
binarized_data <- merged %>%
  select(sample, accelerated, all_of(clinical_vars))

for (var in clinical_vars) {
  original <- merged[[var]]
  before_na <- sum(is.na(original))

  # Determine variable type and apply appropriate binarization function
  if (var %in% numerical_vars) {
    transformed <- binarize_numerical(original)
    var_type <- "numerical"
  } else if (var %in% categorical_vars) {
    transformed <- binarize_categorical(original)
    var_type <- "categorical"
  } else if (var %in% multicategory_vars) {
    transformed <- binarize_multicategory(original)
    var_type <- "multi-category"
  } else {
    # Fallback to categorical for unknown types
    transformed <- binarize_categorical(original)
    var_type <- "unknown"
  }

  after_na <- sum(is.na(transformed))
  introduced <- max(0, after_na - before_na)

  # Log detailed info
  unique_after <- length(unique(transformed[!is.na(transformed)]))
  if (introduced > 0 || unique_after < 2) {
    message(sprintf(
      "[Info] %s (%s): %d unique values, %d NA added",
      var, var_type, unique_after, introduced
    ))
  }

  binarized_data[[var]] <- transformed
}

# ============================================================================
# Statistical tests (Fisher's exact test for binary associations)
# ============================================================================

test_association <- function(outcome, predictor, outcome_name, predictor_name, var_type = "unknown") {
  # Remove NAs
  complete_idx <- !is.na(outcome) & !is.na(predictor)

  if (sum(complete_idx) < 10) {
    message(sprintf("  Skipped %s: only %d complete cases", predictor_name, sum(complete_idx)))
    return(NULL) # Not enough data
  }

  outcome <- outcome[complete_idx]
  predictor <- predictor[complete_idx]

  # Create contingency table
  tbl <- table(predictor, outcome)
  print(tbl)

  # Skip if any dimension is < 2
  if (nrow(tbl) < 2 || ncol(tbl) < 2) {
    message(sprintf(
      "  Skipped %s: insufficient variation (dim: %dx%d)",
      predictor_name, nrow(tbl), ncol(tbl)
    ))
    return(NULL)
  }

  # Fisher's exact test
  test_result <- tryCatch(
    fisher.test(tbl),
    error = function(e) {
      message(sprintf("  Error in Fisher test for %s: %s", predictor_name, e$message))
      NULL
    }
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
    variable_type = var_type,
    n_samples = sum(complete_idx),
    n_levels = nrow(tbl),
    p_value = test_result$p.value,
    odds_ratio = odds_ratio,
    stringsAsFactors = FALSE
  ))
}

# Run tests for all clinical variables
results_list <- list()

message("\nTesting associations...")
for (var in clinical_vars) {
  # Determine variable type
  if (var %in% numerical_vars) {
    var_type <- "numerical"
  } else if (var %in% categorical_vars) {
    var_type <- "categorical"
  } else if (var %in% multicategory_vars) {
    var_type <- "multi-category"
  } else {
    var_type <- "unknown"
  }

  result <- test_association(
    outcome = binarized_data$accelerated,
    predictor = binarized_data[[var]],
    outcome_name = "accelerated",
    predictor_name = var,
    var_type = var_type
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
create_contingency_heatmap <- function(data, var_name, p_value = NULL) {
  tbl <- table(data[[var_name]], data$accelerated)

  if (nrow(tbl) < 2 || ncol(tbl) < 2) {
    return(NULL)
  }

  # Calculate total sample count (excluding NAs)
  n_samples <- sum(!is.na(data[[var_name]]) & !is.na(data$accelerated))

  # Convert to data frame so we can report counts and percentages
  df <- as.data.frame(tbl)
  colnames(df) <- c("BinaryValue", "Acceleration", "Count")
  df$BinaryValue <- as.character(df$BinaryValue)
  df$Acceleration <- as.character(df$Acceleration)
  row_totals <- rowSums(tbl)
  df$row_total <- row_totals[df$BinaryValue]
  df$Proportion <- ifelse(df$row_total == 0, NA, df$Count / df$row_total)
  percent_values <- ifelse(is.na(df$Proportion), 0, df$Proportion * 100)
  df$Label <- sprintf("%d (%.1f%%)", df$Count, percent_values)
  df$row_total <- NULL

  # Clean variable name for title
  clean_name <- gsub("^(demographic|diagnoses|exposures)\\.", "", var_name)

  # Format p-value and sample count for subtitle
  subtitle_text <- if (!is.null(p_value)) {
    sprintf("p = %.3g, N = %d", p_value, n_samples)
  } else {
    sprintf("N = %d", n_samples)
  }

  p <- ggplot(df, aes(x = Acceleration, y = BinaryValue, fill = Proportion)) +
    geom_tile(color = "white") +
    geom_text(aes(label = Label), color = "black") +
    scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 1)) +
    labs(
      title = sprintf("%s vs Age Acceleration", clean_name),
      subtitle = subtitle_text,
      x = "Epigenetic Aging",
      y = "Variable Value"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0),
      axis.text.y = element_text(size = 10)
    )

  return(p)
}

# Create heatmaps for all variables tested
if (nrow(results) > 0) {
  # Sort results by p-value to show most significant first
  results <- results %>% arrange(p_value)
  selected_vars <- results$clinical_variable
  message(sprintf("Creating heatmaps for all %d tested variables (sorted by p-value).", length(selected_vars)))

  heatmap_plots <- list()
  heatmap_files <- c()

  # Save each heatmap as a separate PNG
  for (i in seq_along(selected_vars)) {
    var <- selected_vars[i]
    if (var %in% colnames(binarized_data)) {
      # Get p-value for this variable
      p_val <- results$p_value[results$clinical_variable == var]
      p <- create_contingency_heatmap(binarized_data, var, p_value = p_val)
      if (!is.null(p)) {
        # Clean variable name for filename
        clean_var <- gsub("^(demographic|diagnoses|exposures)\\.", "", var)
        clean_var <- gsub("[^a-zA-Z0-9_]", "_", clean_var)

        filename <- sprintf("results/clinical/%s_%s_contingency.png", clock_name, clean_var)
        ggsave(filename, plot = p, width = 6, height = 4)
        heatmap_files <- c(heatmap_files, filename)
        heatmap_plots[[var]] <- p
      }
    }
  }

  if (length(heatmap_files) > 0) {
    message(sprintf("Saved %d individual contingency heatmaps", length(heatmap_files)))

    # Combine all PNGs into a single page PDF using gridExtra (horizontal layout)
    pdf(sprintf("results/clinical/%s_contingency_heatmaps.pdf", clock_name), width = 28, height = 20)

    # Calculate grid layout (prefer more columns for horizontal layout)
    n_plots <- length(heatmap_files)
    ncol <- ceiling(sqrt(n_plots * 1.4)) # Prefer horizontal layout
    nrow <- ceiling(n_plots / ncol)

    # Read all images
    img_list <- lapply(heatmap_files, function(f) {
      img <- png::readPNG(f)
      grid::rasterGrob(img, interpolate = TRUE)
    })

    # Arrange all plots on one page
    gridExtra::grid.arrange(grobs = img_list, ncol = ncol, nrow = nrow)

    dev.off()
    message(sprintf(
      "Combined %d heatmaps into single page: results/clinical/%s_contingency_heatmaps.pdf",
      n_plots, clock_name
    ))

    # Optionally clean up individual PNG files
    # Uncomment the next line if you want to delete individual PNGs after creating PDF
    # file.remove(heatmap_files)
  }
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
