#!/usr/bin/env Rscript
# Check RNA-seq availability for a single TCGA project
# Called by Snakemake per project for parallel processing

# Setup logging if running via Snakemake
if (exists("snakemake")) {
  log_file <- file(snakemake@log[[1]], open = "wt")
  sink(log_file, type = "output", split = FALSE)
  sink(log_file, type = "message", split = FALSE)

  project <- snakemake@params$project
  pheno_path <- snakemake@input$pheno
  output_file <- snakemake@output$available
  query_output <- snakemake@output$query
} else {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 4) {
    stop("Usage: Rscript check_rnaseq_project.R <project> <pheno_path> <output_file> <query_output>")
  }
  project <- args[1]
  pheno_path <- args[2]
  output_file <- args[3]
  query_output <- args[4]
}

suppressPackageStartupMessages({
  library(dplyr)
  library(TCGAbiolinks)
  library(readr)
})

message(sprintf("[RNA-seq] Checking availability for %s", project))

workflow_type <- "STAR - Counts"

if (!file.exists(pheno_path)) {
  stop(sprintf("Phenotype file not found: %s", pheno_path))
}

# Read phenotype and filter for this project
pheno <- read.table(pheno_path,
  header = TRUE, sep = "\t", quote = "",
  stringsAsFactors = FALSE
)

project_samples <- pheno %>%
  filter(project.project_id == project) %>%
  pull(sample) %>%
  unique()

if (length(project_samples) == 0) {
  message(sprintf("No samples found for project %s", project))
  # Write empty output
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  write_tsv(data.frame(sample = character(0), project = character(0)), output_file)
  saveRDS(NULL, file = query_output)
  quit(save = "no", status = 0)
}

message(sprintf("Checking %d samples for %s", length(project_samples), project))

# Build query
q <- NULL
tryCatch(
  {
    q <- GDCquery(
      project       = project,
      data.category = "Transcriptome Profiling",
      data.type     = "Gene Expression Quantification",
      workflow.type = workflow_type,
      barcode       = project_samples
    )
  },
  error = function(e) {
    message(sprintf("Query failed: %s", e$message))
  }
)

if (is.null(q)) {
  # Write empty output
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  write_tsv(data.frame(sample = character(0), project = character(0)), output_file)
  saveRDS(NULL, file = query_output)
  quit(save = "no", status = 0)
}

# Get available samples
available <- character(0)
rr <- try(getResults(q), silent = TRUE)
if (!inherits(rr, "try-error") && is.data.frame(rr) && "sample.submitter_id" %in% names(rr)) {
  available <- unique(rr$sample.submitter_id)
}

message(sprintf("Found %d/%d samples with RNA-seq available", length(available), length(project_samples)))

# Save query for later download
dir.create(dirname(query_output), recursive = TRUE, showWarnings = FALSE)
saveRDS(q, file = query_output)
message(sprintf("Saved query to %s", query_output))

# Write output
result <- data.frame(
  sample = available,
  project = project,
  stringsAsFactors = FALSE
)

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write_tsv(result, output_file)

message(sprintf("Wrote availability list to %s", output_file))

# Close log file if running via Snakemake
if (exists("snakemake")) {
  sink(type = "output")
  sink(type = "message")
  close(log_file)
}
