#!/usr/bin/env Rscript
# Visualize clinical data availability and missingness patterns

library(ggplot2)
library(dplyr)
library(tidyr)

# Load predictions and phenotype data
predictions <- readRDS("results/clock/gdc_pan/gdc_pancan_predictions.rds")
pheno <- read.table("data/processed/gdc_pancan/methyl/normal_pheno.tsv",
  header = TRUE, sep = "\t", quote = "",
  stringsAsFactors = FALSE
)

message(sprintf("Loaded %d samples with predictions", nrow(predictions)))
message(sprintf("Loaded phenotype data for %d samples", nrow(pheno)))

# Prepare merged data (matching clinical_correlation.R logic)
predictions <- predictions %>%
  rename(sample = id) %>%
  mutate(sample = substr(sample, 1, 16))

merged <- predictions %>%
  select(sample, age, Horvath_residuals) %>%
  inner_join(pheno, by = "sample")

message(sprintf("Merged data: %d samples", nrow(merged)))

# Define clinical variables (from clinical_correlation.R)
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

all_vars <- c(numerical_vars, categorical_vars, multicategory_vars)
all_vars <- all_vars[all_vars %in% colnames(merged)]

# Calculate missingness and availability for each variable
availability_data <- data.frame()

for (var in all_vars) {
  values <- merged[[var]]
  
  # Normalize missing values
  values[values %in% c("", "Not Reported", "not reported", "Not Available", "Unknown")] <- NA
  
  n_missing <- sum(is.na(values))
  n_available <- sum(!is.na(values))
  pct_missing <- (n_missing / length(values)) * 100
  pct_available <- (n_available / length(values)) * 100
  
  # Count unique values (for categorical assessment)
  n_unique <- length(unique(values[!is.na(values)]))
  
  # Determine variable type
  var_type <- if (var %in% numerical_vars) {
    "Numerical"
  } else if (var %in% categorical_vars) {
    "Categorical"
  } else {
    "Multi-category"
  }
  
  availability_data <- rbind(availability_data, data.frame(
    variable = var,
    var_type = var_type,
    n_available = n_available,
    n_missing = n_missing,
    pct_available = pct_available,
    pct_missing = pct_missing,
    n_unique = n_unique,
    stringsAsFactors = FALSE
  ))
}

# Clean variable names for display
availability_data$variable_clean <- gsub("^(demographic|diagnoses|exposures)\\.", "", availability_data$variable)

# Sort by availability
availability_data <- availability_data %>%
  arrange(pct_available)

message("\n=== Clinical Variable Availability Summary ===")
print(availability_data %>% select(variable_clean, var_type, n_available, pct_available, n_unique))

# Create output directory
if (!dir.exists("results/clinical")) {
  dir.create("results/clinical", recursive = TRUE)
}

# 1. Stacked bar chart showing available vs missing data
availability_long <- availability_data %>%
  select(variable_clean, var_type, n_available, n_missing) %>%
  pivot_longer(cols = c(n_available, n_missing), 
               names_to = "status", 
               values_to = "count") %>%
  mutate(
    status = factor(status, levels = c("n_missing", "n_available"),
                   labels = c("Missing", "Available"))
  )

stacked_plot <- ggplot(availability_long, aes(x = reorder(variable_clean, count, sum), y = count, fill = status)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Missing" = "#d73027", "Available" = "#4575b4")) +
  labs(
    title = "Clinical Variable Data Availability",
    subtitle = sprintf("N = %d samples", nrow(merged)),
    x = "Clinical Variable",
    y = "Number of Samples",
    fill = "Status"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    plot.subtitle = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.position = "top"
  )

ggsave(
  filename = "results/clinical/clinical_data_availability.png",
  plot = stacked_plot,
  width = 16,
  height = 9,
  dpi = 300
)

message("\nStacked bar plot saved to: results/clinical/clinical_data_availability.png")

# 2. Percentage availability plot with threshold line
pct_plot <- ggplot(availability_data, aes(x = reorder(variable_clean, pct_available), y = pct_available)) +
  geom_bar(stat = "identity", aes(fill = var_type), alpha = 0.8) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red", linewidth = 1.5) +
  geom_text(aes(label = sprintf("%.0f%%", pct_available)), hjust = -0.1, size = 4) +
  coord_flip() +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Percentage of Available Clinical Data",
    subtitle = "Red line indicates 50% threshold",
    x = "Clinical Variable",
    y = "Data Availability (%)",
    fill = "Variable Type"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    plot.subtitle = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.position = "top"
  ) +
  ylim(0, 105)

ggsave(
  filename = "results/clinical/clinical_data_availability_pct.png",
  plot = pct_plot,
  width = 16,
  height = 9,
  dpi = 300
)

message("Percentage plot saved to: results/clinical/clinical_data_availability_pct.png")

# 3. Summary statistics
cat("\n=== Overall Data Quality Summary ===\n")
cat(sprintf("Total clinical variables analyzed: %d\n", nrow(availability_data)))
cat(sprintf("  - Numerical: %d\n", sum(availability_data$var_type == "Numerical")))
cat(sprintf("  - Categorical: %d\n", sum(availability_data$var_type == "Categorical")))
cat(sprintf("  - Multi-category: %d\n", sum(availability_data$var_type == "Multi-category")))

cat(sprintf("\nVariables with >50%% data: %d (%.1f%%)\n", 
            sum(availability_data$pct_available > 50),
            sum(availability_data$pct_available > 50) / nrow(availability_data) * 100))

cat(sprintf("Variables with >75%% data: %d (%.1f%%)\n", 
            sum(availability_data$pct_available > 75),
            sum(availability_data$pct_available > 75) / nrow(availability_data) * 100))

cat(sprintf("Variables with <25%% data: %d (%.1f%%)\n", 
            sum(availability_data$pct_available < 25),
            sum(availability_data$pct_available < 25) / nrow(availability_data) * 100))

cat("\nMost complete variables (top 5):\n")
print(head(availability_data %>% 
             arrange(desc(pct_available)) %>% 
             select(variable_clean, pct_available, n_unique), 5))

cat("\nLeast complete variables (bottom 5):\n")
print(head(availability_data %>% 
             arrange(pct_available) %>% 
             select(variable_clean, pct_available, n_unique), 5))

# Save summary table
write.csv(availability_data,
  file = "results/clinical/clinical_data_availability_summary.csv",
  row.names = FALSE
)

message("\nSummary table saved to: results/clinical/clinical_data_availability_summary.csv")
message("\n=== Analysis complete! ===")
