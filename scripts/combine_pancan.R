# Combine all project predictions and generate pan-cancer wide analysis

library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(rlang)

message("Combining pan-cancer predictions...")

# Create output directory if needed
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

# Load all project predictions
prediction_files <- list.files("results/clock/gdc_pan/",
  pattern = "^TCGA-.*_predictions\\.rds$",
  full.names = TRUE
)

if (length(prediction_files) == 0) {
  stop("No prediction files found!")
}

prediction_list <- lapply(prediction_files, function(f) {
  pred <- readRDS(f)
  if (nrow(pred) == 0) {
    return(NULL)
  }
  pred
})

# Filter out empty predictions
prediction_list <- Filter(Negate(is.null), prediction_list)

if (length(prediction_list) == 0) {
  stop("All prediction files are empty!")
}

# Combine all predictions
prediction_store <- dplyr::bind_rows(prediction_list)
message(paste("Combined", nrow(prediction_store), "samples from", length(prediction_list), "projects"))

# Save combined predictions
saveRDS(prediction_store, "results/clock/gdc_pan/gdc_pancan_predictions.rds")

# Calculate residuals for pancan data
residuals_store <- data.frame(id = seq_len(nrow(prediction_store)))
for (clock in colnames(prediction_store)) {
  if (clock == "age" || all(is.na(prediction_store[[clock]]))) {
    next
  }
  ok <- !is.na(prediction_store[[clock]]) & !is.na(prediction_store$age)
  if (sum(ok) < 2) {
    next
  }
  model <- lm(prediction_store[[clock]][ok] ~ prediction_store$age[ok])
  residuals <- rep(NA, nrow(prediction_store))
  residuals[ok] <- stats::resid(model)
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

# Create scatter plots for pancan data
scatter_plots <- lapply(clock_cols, function(clock) build_clock_scatter(combined_store, clock))
scatter_grid <- wrap_plots(scatter_plots, ncol = 3)
scatter_grid <- scatter_grid +
  plot_annotation(title = paste("DNAm Age Predictions for GDC Pan-Cancer"))

ggsave(
  filename = "results/clock/gdc_pan/gdc_pancan_scatterplots.png",
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

message("Pan-cancer analysis complete!")
