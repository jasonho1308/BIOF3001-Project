# Combine clock predictions for samples that have RNA-seq data

# Setup logging if running via Snakemake
if (exists("snakemake")) {
  log_file <- file(snakemake@log[[1]], open = "wt")
  sink(log_file, type = "output", split = FALSE)
  sink(log_file, type = "message", split = FALSE)
}

library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(rlang)

message("Combining pan-cancer predictions for samples with RNA-seq data...")

# Create output directory if needed
if (!dir.exists("results/rna")) {
  dir.create("results/rna", recursive = TRUE)
}

clock_cols <- c(
  "Horvath", "Hannum", "Levine", "BNN",
  "skinHorvath", "Wu", "TL", "BLUP", "EN", "PedBE"
)

build_clock_scatter <- function(prediction_df, clock) {
  if (!all(c("age", clock) %in% colnames(prediction_df))) {
    return(ggplot() +
      labs(title = sprintf("%s (R^2 = NA)", clock), x = "Chronological Age", y = "Predicted DNAmAge") +
      theme_void())
  }

  plot_df <- prediction_df[, c("age", clock), drop = FALSE]
  plot_df <- stats::na.omit(plot_df)

  r_squared <- NA_real_
  if (nrow(plot_df) >= 2) {
    model <- stats::lm(as.formula(sprintf("%s ~ age", clock)), data = plot_df)
    r_squared <- summary(model)$r.squared
  }

  title_text <- if (!is.na(r_squared)) {
    sprintf("%s (R^2 = %.2f)", clock, r_squared)
  } else {
    sprintf("%s (R^2 = NA)", clock)
  }

  ggplot(plot_df, aes(age, !!rlang::sym(clock))) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm", color = "red", formula = y ~ x) +
    labs(title = title_text, x = "Chronological Age", y = "Predicted DNAmAge") +
    theme_minimal()
}

# Load RNA-seq queries to get samples with RNA-seq data
rna_queries <- readRDS("data/rna/queries.rds")

if (length(rna_queries) == 0) {
  stop("No RNA-seq queries found!")
}

# Extract sample barcodes from RNA-seq queries
# The barcode in query results is sample.submitter_id (first 4 fields)
rna_sample_ids <- unlist(lapply(rna_queries, function(q) {
  if (!is.null(q) && !is.null(q$results[[1]])) {
    q$results[[1]]$sample.submitter_id
  } else {
    NULL
  }
}))

message(sprintf("Found %d samples with RNA-seq data", length(rna_sample_ids)))

# Load clock predictions
predictions <- readRDS("results/clock/gdc_pan/gdc_pancan_predictions.rds")
message(sprintf("Loaded %d samples with clock predictions", nrow(predictions)))

# Extract sample submitter ID from clock prediction id column (first 4 barcode fields)
predictions$sample_submitter_id <- sapply(predictions$id, function(bc) {
  parts <- strsplit(bc, "-", fixed = TRUE)[[1]]
  if (length(parts) >= 4) {
    paste(parts[1:4], collapse = "-")
  } else {
    substr(bc, 1, 16)
  }
})

# Filter predictions to only include samples with RNA-seq data
combined_store <- predictions %>%
  filter(sample_submitter_id %in% rna_sample_ids)

message(sprintf("Found %d samples with both methylation clock predictions and RNA-seq data", nrow(combined_store)))

if (nrow(combined_store) == 0) {
  stop("No samples found with both clock predictions and RNA-seq data!")
}

# Save combined predictions
saveRDS(combined_store, "results/rna/gdc_pancan_rna_predictions.rds")

# Create scatter plots for pancan data
scatter_plots <- lapply(clock_cols, function(clock) build_clock_scatter(combined_store, clock))
scatter_grid <- wrap_plots(scatter_plots, ncol = 3)
scatter_grid <- scatter_grid +
  plot_annotation(title = paste("DNAm Age Predictions for Samples with RNA-seq Data"))

ggsave(
  filename = "results/rna/gdc_pancan_rna_scatterplots.png",
  plot = scatter_grid,
  width = 15,
  height = 10
)

# Plot boxplots of residuals
residuals_long <- combined_store %>%
  mutate(sample_id = seq_len(n())) %>%
  pivot_longer(
    cols = all_of(clock_cols),
    names_to = "clock", values_to = "pred"
  ) %>%
  group_by(clock) %>%
  mutate(
    residual = {
      ok <- is.finite(pred) & is.finite(age)
      res_vec <- rep(NA_real_, n())
      if (sum(ok) >= 2) {
        res_vec[ok] <- stats::resid(stats::lm(pred[ok] ~ age[ok]))
      }
      res_vec
    }
  ) %>%
  ungroup() %>%
  mutate(clock = factor(clock, levels = clock_cols))

residual_plot <- ggplot(residuals_long, aes(x = clock, y = residual)) +
  geom_boxplot(na.rm = TRUE) +
  labs(
    title = paste("Residuals by clock for Samples with RNA-seq Data"),
    x = "Clock", y = "Residual (predicted - fitted)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = "results/rna/gdc_pancan_rna_residuals_boxplots.png",
  plot = residual_plot,
  width = 15,
  height = 10
)

message("Pan-cancer RNA-seq analysis complete!")

# Close log file if running via Snakemake
if (exists("snakemake")) {
  sink(type = "output")
  sink(type = "message")
  close(log_file)
}
