library(TCGAbiolinks)
library(methylclock)
library(methylclockData)
library(sesameData)
library(ggplot2)
library(SummarizedExperiment)
library(patchwork)
library(dplyr)
library(tidyr)

gdc_queries <- readRDS("data/processed/gdc_pancan/queries.rds")

if (!dir.exists("results/clock/gdc_pan")) {
  dir.create("results/clock/gdc_pan", recursive = TRUE)
}

store_items <- c(
  "Horvath", "Hannum", "Levine", "BNN",
  "skinHorvath", "Wu", "TL", "BLUP", "EN", "PedBE", "age"
)
prediction_store <- data.frame()

# Per-tissue clock calculation and plotting
for (proj in names(gdc_queries)) {
  print(paste("Processing project:", proj))
  query_TCGA_meth <- gdc_queries[[proj]]

  me_TCGA <- NULL
  tryCatch(
    me_TCGA <- GDCprepare(query_TCGA_meth, directory = "data/raw/GDCdata/"),
    error = function(e) {
      message(paste("Error preparing data for project:", proj))
      message(e)
    }
  )
  if (is.null(me_TCGA)) {
    message(paste("No data available for project:", proj))
    next
  }
  dim(me_TCGA)

  # prepare data for methylclock

  data <- assay(me_TCGA) # converts to usable matrix
  data <- na.omit(data)

  age <- me_TCGA$age_at_index

  # checks for missing CpGs for the various clocks
  cpgs.missing <- checkClocks(data)

  # makes prediction based on methylation data of normal samples
  prediction <- DNAmAge(data, age = age, cell.count = FALSE)

  # Merge elements for whole pancan analysis later
  prediction_store <- rbind(prediction_store, prediction[, store_items])

  # Scatter plots with regression lines using tidy evaluation
  chron_age <- prediction$age
  clock_cols <- c(
    "Horvath", "Hannum", "Levine", "BNN",
    "skinHorvath", "Wu", "TL", "BLUP", "EN", "PedBE"
  )

  scatter_plots <- lapply(clock_cols, function(clock) {
    ggplot(prediction, aes(x = age, y = !!sym(clock))) +
      geom_point(alpha = 0.5) +
      geom_smooth(method = "lm", color = "red", formula = y ~ x) +
      labs(title = clock, x = "Chronological Age", y = "Predicted DNAmAge") +
      theme_minimal()
  })

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
scatter_plots <- lapply(clock_cols, function(clock) {
  ggplot(combined_store, aes(x = age, y = !!sym(clock))) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm", color = "red", formula = y ~ x) +
    labs(title = clock, x = "Chronological Age", y = "Predicted DNAmAge") +
    theme_minimal()
})

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
