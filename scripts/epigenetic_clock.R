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

for (proj in names(gdc_queries)) {
  print(paste("Processing project:", proj))
  query_TCGA_meth <- gdc_queries[[proj]]
  # View query results
  res_TCGA_meth <- getResults(query_TCGA_meth)
  # View(res_TCGA_meth) # metadata

  # Download actual data
  GDCdownload(query_TCGA_meth, directory = "data/raw/GDCdata/")
  me_TCGA <- NULL
  tryCatch(
    me_TCGA <- GDCprepare(query_TCGA_meth, directory = "data/raw/GDCdata/"),
    error = function(e) {
      message(paste("Error preparing data for project:", proj))
      message(e)

      # skip to next iteration, outside of function scope
      return(NULL)
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

  # Scatter plots with regression lines using tidy evaluation
  chron_age <- prediction$age
  clock_cols <- colnames(prediction[, c(2, 5, 8, 12, 15, 18, 21, 24, 27)])

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
}

# read from saved scatter plots dir and save all to a new page generated results to a pdf
png_files <- list.files("results/clock/gdc_pan/", pattern = "_methylclock_scatterplots.png", full.names = TRUE)
pdf("results/clock/gdc_pan/gdc_pancan_methylclock_scatterplots.pdf")
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
