# Generate final PDF combining all plots for epigenetic clock analysis of GDC-PANCAN data

# Setup logging if running via Snakemake
if (exists("snakemake")) {
  log_file <- file(snakemake@log[[1]], open = "wt")
  sink(log_file, type = "output", split = FALSE)
  sink(log_file, type = "message", split = FALSE)
}

library(png)
library(grid)

message("Generating PDF report...")

# Get all PNG files
png_files <- list.files("results/clock/gdc_pan/",
  pattern = "\\.png$",
  full.names = TRUE
)

if (length(png_files) == 0) {
  stop("No PNG files found!")
}

# Sort files to have pancan plots at the end
png_files <- sort(png_files)
pancan_files <- grep("pancan", png_files, value = TRUE)
project_files <- setdiff(png_files, pancan_files)
png_files <- c(pancan_files, project_files)

message(paste("Creating PDF with", length(png_files), "plots"))

# Create PDF
pdf("results/clock/gdc_pan/gdc_pancan_methylclock.pdf", width = 15, height = 10)

# Loop through each PNG file and add to PDF on separate pages
for (i in seq_along(png_files)) {
  png_file <- png_files[i]
  message(paste("Adding:", basename(png_file)))

  img <- png::readPNG(png_file)
  grid::grid.raster(img)

  # Don't create new page after the last image
  if (i < length(png_files)) {
    grid::grid.newpage()
  }
}

dev.off()

message("PDF generation complete!")

# Close log file if running via Snakemake
if (exists("snakemake")) {
  sink(type = "output")
  sink(type = "message")
  close(log_file)
}
