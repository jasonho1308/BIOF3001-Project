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

gdc_queries <- readRDS("data/processed/gdc_pancan/queries.rds")

if (!dir.exists("results/clock/gdc_pan")) {
  dir.create("results/clock/gdc_pan", recursive = TRUE)
}

store_items <- c(
  "Horvath", "Hannum", "Levine", "BNN",
  "skinHorvath", "Wu", "TL", "BLUP", "EN", "PedBE", "age"
)

clock_cols <- c(
  "Horvath", "Hannum", "Levine", "BNN",
  "skinHorvath", "Wu", "TL", "BLUP", "EN", "PedBE"
)

prediction_store <- data.frame()

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

  ggplot(plot_df, aes_string(x = "age", y = clock)) +
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

# Per-tissue clock calculation and plotting
for (proj in names(gdc_queries)) {
  print(paste("Processing project:", proj))
  query_TCGA_meth <- gdc_queries[[proj]]

  platform_ses <- prepare_project_methylation(proj, query_TCGA_meth, "data/raw/GDCdata/")
  if (is.null(platform_ses) || !length(platform_ses)) {
    message(paste("No data available for project:", proj))
    next
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

    cpgs.missing <- checkClocks(data)

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

    platform_prediction$sample_id <- rownames(platform_prediction)
    if (is.null(platform_prediction$sample_id)) {
      platform_prediction$sample_id <- colnames(data)
    }
    platform_prediction$platform_source <- platform_name
    platform_prediction$project <- proj
    project_predictions[[platform_name]] <- platform_prediction
  }

  if (!length(project_predictions)) {
    message(sprintf("No predictions generated for project: %s", proj))
    next
  }

  prediction <- dplyr::bind_rows(project_predictions)

  # Merge elements for whole pancan analysis later
  prediction_store <- dplyr::bind_rows(prediction_store, prediction[, store_items])

  # Scatter plots with regression lines and R^2 annotation
  scatter_plots <- lapply(clock_cols, function(clock) build_clock_scatter(prediction, clock))

  scatter_grid <- wrap_plots(scatter_plots, ncol = 3)

  # add title to the combined plot
  scatter_grid <- scatter_grid +
    plot_annotation(title = paste("DNAm Age Predictions for", proj))

  ggsave(
    filename = paste0("results/clock/gdc_pan/", proj, "_scatterplots.png"),
    plot = scatter_grid,
    width = 15,
    height = 10
  )

  # Boxplots of residuals
  residual_cols <- c(
    "ageAcc2.Horvath", "ageAcc2.Hannum", "ageAcc2.Levine",
    "ageAcc2.BNN", "ageAcc2.skinHorvath", "ageAcc2.Wu",
    "ageAcc2.TL", "ageAcc2.BLUP", "ageAcc2.EN",
    "ageAcc2.PedBE"
  )
  residuals_long <- prediction %>%
    mutate(sample_id = seq_len(n())) %>%
    pivot_longer(
      cols = all_of(clock_cols),
      names_to = "clock", values_to = "pred"
    ) %>%
    group_by(clock) %>%
    mutate(
      residual = {
        ok <- !is.na(pred) & !is.na(age)
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
}

# Pancan-wide clock calculation and plotting
# Residuals calculation
# initialize empty data frame for residuals with em
residuals_store <- data.frame(id = seq_len(nrow(prediction_store)))
for (clock in colnames(prediction_store)) {
  if (clock == "age" || all(is.na(prediction_store[[clock]]))) {
    next
  }
  model <- lm(prediction_store[[clock]] ~ prediction_store$age)
  # insert NA for skipped samples
  residuals <- rep(NA, nrow(prediction_store))
  ok <- !is.na(prediction_store[[clock]]) & !is.na(prediction_store$age)
  residuals[ok] <- stats::resid(model)[ok]
  residuals_store <- cbind(
    residuals_store,
    setNames(
      data.frame(residuals),
      paste0(clock, "_residuals")
    )
  )
}
# Combine predictions and residuals
combined_store <- cbind(prediction_store, residuals_store)

# create scatter plots for pancan data in a grid layout
scatter_plots <- lapply(clock_cols, function(clock) build_clock_scatter(combined_store, clock))

scatter_grid <- wrap_plots(scatter_plots, ncol = 3)

# add title to the combined plot
scatter_grid <- scatter_grid +
  plot_annotation(title = paste("DNAm Age Predictions for GDC Pan-Cancer"))

ggsave(
  filename = paste0("results/clock/gdc_pan/gdc_pancan_scatterplots.png"),
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
      ok <- !is.na(pred) & !is.na(age)
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
    title = paste("Residuals by clock for GDC Pan-Cancer"),
    x = "Clock", y = "Residual (predicted - fitted)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = "results/clock/gdc_pan/gdc_pancan_residuals_boxplots.png",
  plot = residual_plot,
  width = 15,
  height = 10
)

# read from saved scatter plots dir and save all to a new page generated results to a pdf
png_files <- list.files("results/clock/gdc_pan/", pattern = ".png", full.names = TRUE)
pdf("results/clock/gdc_pan/gdc_pancan_methylclock.pdf")
# loop through each png file and add to pdf to separate pages
for (png_file in png_files) {
  img <- png::readPNG(png_file)
  grid::grid.raster(img)
  # don't create new page after the last image
  if (png_file != tail(png_files, n = 1)) {
    grid::grid.newpage()
  }
}
dev.off()
