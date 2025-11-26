#!/usr/bin/env Rscript
# Combine per-project RNA-seq availability lists into a single TSV
# Called by Snakemake after all check_rnaseq_project rules complete

# Setup logging if running via Snakemake
if (exists("snakemake")) {
  log_file <- file(snakemake@log[[1]], open = "wt")
  sink(log_file, type = "output", split = FALSE)
  sink(log_file, type = "message", split = FALSE)

  input_files <- snakemake@input$project_lists
  output_file <- snakemake@output$combined
} else {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2) {
    stop("Usage: Rscript combine_rnaseq_availability.R <output_file> <input_file1> [<input_file2> ...]")
  }
  output_file <- args[1]
  input_files <- args[-1]
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

message(sprintf("Combining %d project RNA-seq availability lists", length(input_files)))

# Read and combine all project lists
all_data <- list()
for (f in input_files) {
  if (!file.exists(f)) {
    warning(sprintf("File not found: %s", f))
    next
  }

  tryCatch(
    {
      df <- read_tsv(f, show_col_types = FALSE)
      if (nrow(df) > 0) {
        all_data[[f]] <- df
      }
    },
    error = function(e) {
      warning(sprintf("Failed to read %s: %s", f, e$message))
    }
  )
}

if (length(all_data) == 0) {
  warning("No data found in any input files")
  combined <- data.frame(sample = character(0), project = character(0))
} else {
  combined <- bind_rows(all_data) %>%
    distinct(sample, project)
}

message(sprintf("Combined total: %d unique samples with RNA-seq", nrow(combined)))

# Write combined output
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write_tsv(combined, output_file)
message(sprintf("Wrote combined availability list to %s", output_file))

# Close log file if running via Snakemake
if (exists("snakemake")) {
  sink(type = "output")
  sink(type = "message")
  close(log_file)
}
