#!/usr/bin/env Rscript
# Title: RNA-seq downloader for extreme acceleration cohorts
# Description: Reads the Accelerated/Decelerated sample sheets produced by
#   classify_extreme_samples.R, extracts the TCGA sample IDs (first four barcode
#   fields), and invokes get_brca_rnaseq_by_barcode.R to download STAR count
#   matrices plus metadata for each cohort. Defaults cover LIHC (global), BRCA,
#   and KIRC tissue-specific extremes, but a custom TSV (label,project,sample_file)
#   can be provided to override or extend the cases.

suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

get_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_flag <- "--file="
  file_arg <- cmd_args[startsWith(cmd_args, file_flag)]
  if (length(file_arg)) {
    return(normalizePath(sub(file_flag, "", file_arg[length(file_arg)])))
  }
  normalizePath(getwd())
}

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1 && grepl("\\.tsv$|\\.csv$", args[[1]], ignore.case = TRUE)) {
  args[[1]]
} else {
  NA_character_
}
output_dir <- if (length(args) >= 2) args[[2]] else file.path("results", "differential_analysis")
id_dir <- if (length(args) >= 3) args[[3]] else output_dir
workflow_type <- if (length(args) >= 4) args[[4]] else "STAR - Counts"
chunk_size <- if (length(args) >= 5) as.integer(args[[5]]) else 50L

script_dir <- dirname(get_script_path())
helper_script <- file.path(script_dir, "get_brca_rnaseq_by_barcode.R")
if (!file.exists(helper_script)) {
  stop(sprintf("Helper script not found at %s", helper_script))
}

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(id_dir)) dir.create(id_dir, recursive = TRUE, showWarnings = FALSE)

default_cases <- data.table(
  label = c(
    "TCGA-LIHC_balanced_global",
    "TCGA-BRCA_tissue_quartiles",
    "TCGA-KIRC_tissue_quartiles"
  ),
  project = c("TCGA-LIHC", "TCGA-BRCA", "TCGA-KIRC"),
  sample_file = file.path(
    "results", "clock", "derived",
    sprintf("%s_extreme_samples.txt", c(
      "TCGA-LIHC_balanced_global",
      "TCGA-BRCA_tissue_quartiles",
      "TCGA-KIRC_tissue_quartiles"
    ))
  )
)

cases <- if (!is.na(config_path)) {
  if (!file.exists(config_path)) {
    stop(sprintf("Config file not found: %s", config_path))
  }
  dt <- fread(config_path)
  required <- c("label", "project", "sample_file")
  missing <- setdiff(required, names(dt))
  if (length(missing)) {
    stop(sprintf("Config file missing columns: %s", paste(missing, collapse = ", ")))
  }
  dt
} else {
  default_cases
}

read_sample_ids <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Sample sheet not found: %s", path))
  }
  dt <- fread(path)
  if ("sample_id" %in% names(dt)) {
    ids <- dt[["sample_id"]]
  } else {
    ids <- dt[[1]]
  }
  ids <- trimws(ids)
  ids <- ids[nchar(ids) > 0]
  unique(ids)
}

failed_cases <- list()

for (i in seq_len(nrow(cases))) {
  label <- cases$label[i]
  project <- cases$project[i]
  sample_file <- cases$sample_file[i]
  message(sprintf("[%s] Extracting barcodes from %s", label, sample_file))
  sample_ids <- read_sample_ids(sample_file)
  if (!length(sample_ids)) {
    stop(sprintf("No sample IDs found in %s", sample_file))
  }
  id_path <- file.path(id_dir, sprintf("%s_extreme_sample_ids.txt", label))
  writeLines(sample_ids, con = id_path)
  counts_path <- file.path(output_dir, sprintf("%s_extreme_rnaseq_counts.tsv", label))
  metadata_path <- file.path(output_dir, sprintf("%s_extreme_rnaseq_metadata.tsv", label))
  message(sprintf("[%s] Downloading counts for %d samples", label, length(sample_ids)))
  status <- system2(
    "Rscript",
    args = c(
      helper_script,
      id_path,
      project,
      counts_path,
      metadata_path,
      as.character(chunk_size),
      shQuote(workflow_type)
    )
  )
  if (!identical(status, 0L)) {
    warning(sprintf("[%s] Download failed (exit code %s)", label, as.character(status)))
    failed_cases[[label]] <- status
    next
  }
  message(sprintf("[%s] Counts saved to %s", label, counts_path))
}

if (length(failed_cases)) {
  warning(sprintf(
    "Completed with %d failure(s): %s",
    length(failed_cases),
    paste(names(failed_cases), collapse = ", ")
  ))
} else {
  message("All downloads complete.")
}
