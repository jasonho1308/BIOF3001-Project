# Filter per-project clock predictions to samples with RNA-seq data and generate plots

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
library(parallel)

message("Processing per-project predictions for samples with RNA-seq data...")

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

# Load RNA-seq queries to get samples with RNA-seq data per project
rna_queries <- readRDS("data/rna/queries.rds")

if (length(rna_queries) == 0) {
  stop("No RNA-seq queries found!")
}

# Get list of projects with RNA-seq data
projects_with_rna <- names(rna_queries)
message(sprintf("Found %d projects with RNA-seq data", length(projects_with_rna)))

# Detect number of cores (use half to avoid overwhelming the system)
n_cores <- max(1, floor(detectCores() / 2))
message(sprintf("Using %d cores for parallel processing", n_cores))

# Function to process a single project
process_project <- function(proj, rna_queries) {
  result <- list(
    project = proj,
    n_samples = 0,
    success = FALSE
  )

  tryCatch({
    message(paste("Processing project:", proj))

    # Get RNA-seq sample IDs for this project
    rna_query <- rna_queries[[proj]]
    if (is.null(rna_query) || is.null(rna_query$results[[1]])) {
      message(sprintf("No RNA-seq query results for project: %s", proj))
      return(result)
    }

    rna_sample_ids <- rna_query$results[[1]]$sample.submitter_id
    message(sprintf("Found %d RNA-seq samples for project %s", length(rna_sample_ids), proj))

    # Load clock predictions for this project
    prediction_file <- paste0("results/clock/gdc_pan/", proj, "_predictions.rds")
    if (!file.exists(prediction_file)) {
      message(sprintf("No prediction file found for project: %s", proj))
      return(result)
    }

    predictions <- readRDS(prediction_file)
    if (nrow(predictions) == 0) {
      message(sprintf("Empty predictions for project: %s", proj))
      return(result)
    }

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
    filtered_predictions <- predictions %>%
      filter(sample_submitter_id %in% rna_sample_ids)

    if (nrow(filtered_predictions) == 0) {
      message(sprintf("No samples with both clock predictions and RNA-seq data for project: %s", proj))
      return(result)
    }

    message(sprintf("Found %d samples with both data for project %s", nrow(filtered_predictions), proj))

    # Save filtered predictions
    saveRDS(filtered_predictions, paste0("results/rna/", proj, "_rna_predictions.rds"))

    # Generate scatter plots
    scatter_plots <- lapply(clock_cols, function(clock) build_clock_scatter(filtered_predictions, clock))
    scatter_grid <- wrap_plots(scatter_plots, ncol = 3)
    scatter_grid <- scatter_grid +
      plot_annotation(title = paste("DNAm Age Predictions for", proj, "(RNA-seq samples)"))

    ggsave(
      filename = paste0("results/rna/", proj, "_rna_scatterplots.png"),
      plot = scatter_grid,
      width = 15,
      height = 10
    )

    # Generate residual boxplots
    residuals_long <- filtered_predictions %>%
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
        title = paste("Residuals by clock for", proj, "(RNA-seq samples)"),
        x = "Clock", y = "Residual (predicted - fitted)"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    ggsave(
      filename = paste0("results/rna/", proj, "_rna_residuals_boxplots.png"),
      plot = residual_plot,
      width = 12,
      height = 6
    )

    message(paste("Successfully processed project:", proj))

    result$n_samples <- nrow(filtered_predictions)
    result$success <- TRUE

  }, error = function(e) {
    message(paste("Error processing project:", proj))
    message(e)
  })

  result
}

# Process projects in parallel
results <- mclapply(
  projects_with_rna,
  process_project,
  rna_queries = rna_queries,
  mc.cores = n_cores,
  mc.preschedule = FALSE
)

# Summary
successful <- sum(sapply(results, function(r) if (is.list(r)) r$success else FALSE))
total_samples <- sum(sapply(results, function(r) if (is.list(r)) r$n_samples else 0))

message(sprintf("Processed %d/%d projects successfully", successful, length(projects_with_rna)))
message(sprintf("Total samples with both methylation and RNA-seq data: %d", total_samples))

message("Per-project RNA-seq analysis complete!")

# Close log file if running via Snakemake
if (exists("snakemake")) {
  sink(type = "output")
  sink(type = "message")
  close(log_file)
}
