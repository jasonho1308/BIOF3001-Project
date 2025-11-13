# Process individual TCGA project for epigenetic clock analysis

# Get project from Snakemake params (when run via Snakemake)
# or from command line args (when run directly)
if (exists("snakemake")) {
  proj <- snakemake@params$project
  # Redirect all output to log file only (no console output)
  log_file <- file(snakemake@log[[1]], open = "wt")
  sink(log_file, type = "output", split = FALSE)
  sink(log_file, type = "message", split = FALSE)
} else {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) {
    stop("Usage: Rscript process_project.R <PROJECT_NAME>")
  }
  proj <- args[1]
}

library(TCGAbiolinks)
library(methylclock)
library(methylclockData)
library(sesameData)
library(ggplot2)
library(SummarizedExperiment)
library(patchwork)
library(dplyr)
library(tidyr)
library(rlang)

message(paste("Processing project:", proj))

# Load queries
gdc_queries <- readRDS("data/processed/gdc_pancan/queries.rds")

# Create output directory
if (!dir.exists("results/clock/gdc_pan")) {
  dir.create("results/clock/gdc_pan", recursive = TRUE)
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

prepare_project_methylation <- function(proj, query, data_dir) {
  if (is.null(query$results) || length(query$results) == 0) {
    message(sprintf("No results stored for project: %s", proj))
    return(NULL)
  }

  results_tbl <- query$results[[1]]

  if (!"platform" %in% colnames(results_tbl)) {
    se_obj <- tryCatch(
      GDCprepare(query, directory = data_dir),
      error = function(e) {
        message(sprintf("Error preparing data for project: %s", proj))
        message(e)
        NULL
      }
    )
    if (is.null(se_obj)) {
      return(NULL)
    }
    return(list(`unknown platform` = se_obj))
  }

  platforms <- unique(results_tbl$platform)
  platform_ses <- vector("list", length(platforms))
  names(platform_ses) <- platforms

  for (platform_name in platforms) {
    message(sprintf("Preparing platform %s for project %s", platform_name, proj))
    platform_rows <- results_tbl$platform == platform_name
    if (!any(platform_rows)) {
      next
    }

    query_platform <- query
    query_platform$results[[1]] <- results_tbl[platform_rows, , drop = FALSE]

    platform_se <- tryCatch(
      GDCprepare(query_platform, directory = data_dir),
      error = function(e) {
        message(sprintf(
          "Error preparing platform %s for project %s",
          platform_name, proj
        ))
        message(e)
        NULL
      }
    )

    platform_ses[[platform_name]] <- platform_se
  }

  platform_ses <- Filter(Negate(is.null), platform_ses)

  if (!length(platform_ses)) {
    message(sprintf("No data available after platform preparation for project: %s", proj))
    return(NULL)
  }

  platform_ses
}

# Process the project
query_TCGA_meth <- gdc_queries[[proj]]

if (is.null(query_TCGA_meth)) {
  message(sprintf("No query found for project %s in queries.rds; skipping.", proj))
  saveRDS(data.frame(), paste0("results/clock/gdc_pan/", proj, "_predictions.rds"))
  quit(status = 0)
}

platform_ses <- prepare_project_methylation(proj, query_TCGA_meth, "data/raw/GDCdata/")
if (is.null(platform_ses) || !length(platform_ses)) {
  message(paste("No data available for project:", proj))
  # Save empty prediction to mark as processed
  saveRDS(data.frame(), paste0("results/clock/gdc_pan/", proj, "_predictions.rds"))
  quit(status = 0)
}

project_predictions <- list()

for (platform_name in names(platform_ses)) {
  se_obj <- platform_ses[[platform_name]]
  data <- SummarizedExperiment::assay(se_obj)
  if (!is.matrix(data)) {
    data <- as.matrix(data)
  }
  complete_rows <- stats::complete.cases(data)
  data <- data[complete_rows, , drop = FALSE]
  if (!nrow(data)) {
    next
  }

  age <- SummarizedExperiment::colData(se_obj)$age_at_index


  platform_prediction <- tryCatch(
    DNAmAge(data, age = age, cell.count = FALSE),
    error = function(e) {
      message(sprintf(
        "DNAmAge failed for platform %s in project %s", platform_name, proj
      ))
      message(e)
      NULL
    }
  )

  if (is.null(platform_prediction)) {
    next
  }

  platform_prediction$platform_source <- platform_name
  platform_prediction$project <- proj
  project_predictions[[platform_name]] <- platform_prediction
}

if (!length(project_predictions)) {
  message(sprintf("No predictions generated for project: %s", proj))
  # Save empty prediction to mark as processed
  saveRDS(data.frame(), paste0("results/clock/gdc_pan/", proj, "_predictions.rds"))
  quit(status = 0)
}

prediction <- dplyr::bind_rows(project_predictions)

# Save predictions for later combination
saveRDS(prediction, paste0("results/clock/gdc_pan/", proj, "_predictions.rds"))

# Generate scatter plots
scatter_plots <- lapply(clock_cols, function(clock) build_clock_scatter(prediction, clock))
scatter_grid <- wrap_plots(scatter_plots, ncol = 3)
scatter_grid <- scatter_grid +
  plot_annotation(title = paste("DNAm Age Predictions for", proj))

ggsave(
  filename = paste0("results/clock/gdc_pan/", proj, "_scatterplots.png"),
  plot = scatter_grid,
  width = 15,
  height = 10
)

# Generate residual boxplots
residuals_long <- prediction %>%
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
    title = paste("Residuals by clock for", proj),
    x = "Clock", y = "Residual (predicted - fitted)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = paste0("results/clock/gdc_pan/", proj, "_residuals_boxplots.png"),
  plot = residual_plot,
  width = 12,
  height = 6
)

message(paste("Successfully processed project:", proj))

# Close log file if running via Snakemake
if (exists("snakemake")) {
  sink(type = "output")
  sink(type = "message")
  close(log_file)
}
