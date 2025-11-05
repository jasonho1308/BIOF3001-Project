library(TCGAbiolinks)
library(methylclock)
library(methylclockData)
library(sesameData)
library(ggplot2)
library(SummarizedExperiment)
library(patchwork)

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
  # View query results
  res_TCGA_meth <- getResults(query_TCGA_meth)

  me_TCGA <- NULL
  tryCatch(
    me_TCGA <- GDCprepare(query_TCGA_meth, directory = "data/raw/GDCdata/"),
    error = function(e) {
      message(paste("Error preparing data for project:", proj))
      message(e)
    }
  )
  if (is.null(me_TCGA)) {
    next
  }
  dim(me_TCGA)

  # prepare data for methylclock

  data <- assay(me_TCGA) # converts to usable matrix
  data <- na.omit(data)
  samples <- colnames(data)

  age <- me_TCGA$age_at_index

  # checks for missing CpGs for the various clocks
  cpgs.missing <- checkClocks(data)

  # makes prediction based on methylation data of normal samples
  prediction <- DNAmAge(data, age = age, cell.count = FALSE)

  # Store elements for whole pancan analysis later
  prediction_store <- prediction[, store_items]

  # Scatter plots with regression lines using tidy evaluation
  chron_age <- prediction$age
  clock_cols <- c("Horvath", "Hannum", "Levine", "BNN",
                  "skinHorvath", "Wu", "TL", "BLUP", "EN", "PedBE")

  scatter_plots <- lapply(clock_cols, function(clock) {
    ggplot(prediction, aes(x = age, y = !!sym(clock))) +
      geom_point(alpha = 0.5) +
      geom_smooth(method = "lm", color = "red") +
      labs(title = clock, x = "Chronological Age", y = "Predicted DNAmAge") +
      theme_minimal()
  })

  scatter_grid <- wrap_plots(scatter_plots, ncol = 3)



  # add title to the combined plot
  scatter_grid <- scatter_grid +
    plot_annotation(title = paste("DNAm Age Predictions for", proj))

  ggsave(
    filename = paste0("results/clock/gdc_pan/", proj, "_methylclock_scatterplots.png"),
    plot = scatter_grid,
    width = 15,
    height = 10
  )

  # Boxplots of residuals
  residual_cols <- c("ageAcc2.Horvath", "ageAcc2.Hannum", "ageAcc2.Levine",
                     "ageAcc2.BNN", "ageAcc2.skinHorvath", "ageAcc2.Wu",
                     "ageAcc2.TL", "ageAcc2.BLUP", "ageAcc2.EN",
                     "ageAcc2.PedBE")
  residual_plots <- lapply(residual_cols, function(res_col) {
    # Skip when column not found
    if (!res_col %in% colnames(prediction)) {
      return(NULL)
    }
    ggplot(prediction, aes(x = "", y = !!sym(res_col))) +
      geom_boxplot() +
      labs(title = gsub("ageAcc2.", "", res_col),
           y = "Residuals") +
      theme_minimal()
  })
  # Remove NULL plots
  residual_plots <- Filter(Negate(is.null), residual_plots)
  residual_grid <- wrap_plots(residual_plots, ncol = 3) +
    plot_annotation(title = paste("DNAm Age Residuals for", proj))
  ggsave(
      filename = paste0("results/clock/gdc_pan/", proj, "_methylclock_residuals_boxplots.png"),
      plot = residual_grid,
      width = 15,
      height = 10
    )
}

# Pancan-wide clock calculation and plotting
# Residuals calculation
residuals_store <- data.frame()
for (clock in colnames(prediction_store)) {
  if (clock == "age") {
    next
  }
  model <- lm(prediction_store[[clock]] ~ prediction_store$age)
  residuals_store[[paste0(clock, "_residuals")]] <- resid(model)
}
# Combine predictions and residuals
combined_store <- cbind(prediction_store, residuals_store)

# Plot boxplots of residuals
residual_plots <- lapply(colnames(residuals_store), function(res_clock) {
  # Skip when column not found
  if (!res_clock %in% colnames(combined_store)) {
    return(NULL)
  }
  ggplot(combined_store, aes(x = "", y = !!sym(res_clock))) +
    geom_boxplot() +
    labs(title = paste("Residuals of", gsub("_residuals", "", res_clock)),
         y = "Residuals") +
    theme_minimal()
})
residual_grid <- wrap_plots(residual_plots, ncol = 3) +
  plot_annotation(title = "DNAm Age Residuals for GDC Pancan")
ggsave(
    filename = "results/clock/gdc_pan/gdc_pancan_methylclock_residuals_boxplots.png",
    plot = residual_grid,
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
