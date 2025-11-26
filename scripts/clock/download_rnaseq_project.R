#!/usr/bin/env Rscript
# Download RNA-seq data for a single TCGA project
# Uses saved query from check_rnaseq_project.R

# Setup logging if running via Snakemake
if (exists("snakemake")) {
  log_file <- file(snakemake@log[[1]], open = "wt")
  sink(log_file, type = "output", split = FALSE)
  sink(log_file, type = "message", split = FALSE)

  project <- snakemake@params$project
  query_path <- snakemake@input$query
  output_dir <- snakemake@output$data_dir
} else {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 3) {
    stop("Usage: Rscript download_rnaseq_project.R <project> <query_path> <output_dir>")
  }
  project <- args[1]
  query_path <- args[2]
  output_dir <- args[3]
}

suppressPackageStartupMessages({
  library(TCGAbiolinks)
})

message(sprintf("[RNA-seq] Downloading data for %s", project))

if (!file.exists(query_path)) {
  stop(sprintf("Query file not found: %s", query_path))
}

# Load saved query
message(sprintf("Loading query from %s", query_path))
query <- readRDS(query_path)

if (is.null(query)) {
  stop(sprintf("Query is NULL for project %s", project))
}

# Create output directory
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Download RNA-seq data
message(sprintf("Downloading RNA-seq data to %s", output_dir))
tryCatch(
  {
    GDCdownload(query, directory = output_dir, files.per.chunk = 20)
    message(sprintf("Successfully downloaded RNA-seq data for %s", project))
  },
  error = function(e) {
    message(sprintf("Error downloading RNA-seq data for %s: %s", project, e$message))
    stop(e)
  }
)

message(sprintf("RNA-seq download complete for %s", project))

# Close log file if running via Snakemake
if (exists("snakemake")) {
  sink(type = "output")
  sink(type = "message")
  close(log_file)
}
