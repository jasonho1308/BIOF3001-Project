#!/usr/bin/env Rscript
# Title: Extreme acceleration status classifier
# Description: Reads TCGA pan-cancer clock predictions, assigns acceleration
#   statuses using Q1/Q4 Horvath residual cutoffs, and exports per-project
#   sample lists (LIHC global, BRCA/KIRC tissue-specific) limited to the first
#   four barcode fields. Accelerated samples represent the lower tail, while
#   Decelerated samples come from the upper tail. Outputs one TSV per cohort
#   listing the sample ID and its classification.

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
predictions_path <- if (length(args) >= 1) args[[1]] else file.path(
  "results", "clock", "gdc_pan", "gdc_pancan_predictions.rds"
)
output_dir <- if (length(args) >= 2) args[[2]] else file.path(
  "results", "clock", "derived"
)
lower_quantile <- if (length(args) >= 3) as.numeric(args[[3]]) else 0.25
upper_quantile <- if (length(args) >= 4) as.numeric(args[[4]]) else 0.75

if (!file.exists(predictions_path)) {
  stop(sprintf("Predictions file not found: %s", predictions_path))
}

if (is.na(lower_quantile) || is.na(upper_quantile) || lower_quantile >= upper_quantile) {
  stop("Quantile cutoffs must be numeric and lower < upper.")
}

message(sprintf("Reading predictions from %s", predictions_path))
raw_obj <- readRDS(predictions_path)
if (!is.data.frame(raw_obj)) {
  stop("Predictions RDS must deserialize to a data.frame/data.table object.")
}
dt <- data.table::as.data.table(raw_obj)

required_cols <- c("project", "id", "Horvath_residuals")
missing <- setdiff(required_cols, names(dt))
if (length(missing)) {
  stop(sprintf(
    "Missing required columns in predictions file: %s",
    paste(missing, collapse = ", ")
  ))
}
dt <- dt[, ..required_cols]

sample_id_from_barcode <- function(barcodes) {
  vapply(barcodes, function(bc) {
    parts <- strsplit(bc, "-", fixed = TRUE)[[1]]
    if (length(parts) >= 4) {
      paste(parts[1:4], collapse = "-")
    } else {
      substr(bc, 1, 16)
    }
  }, character(1))
}

assign_status_vector <- function(values, lower, upper) {
  status <- rep(NA_character_, length(values))
  keep <- !is.na(values)
  if (!any(keep)) return(status)
  qs <- stats::quantile(values[keep], probs = c(lower, upper), na.rm = TRUE)
  status[keep & values < qs[1]] <- "Accelerated"
  status[keep & values > qs[2]] <- "Decelerated"
  status
}

write_status_file <- function(data, label, out_dir) {
  if (!nrow(data)) {
    stop(sprintf("No extreme samples found for %s", label))
  }
  out_path <- file.path(out_dir, sprintf("%s_extreme_samples.txt", label))
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  data <- data[order(data$classification, data$sample_id), ]
  data.table::fwrite(data, out_path, sep = "\t")
  message(sprintf("Wrote %d samples to %s", nrow(data), out_path))
}

# Global LIHC classification (pan-cancer residual quantiles)
dt$global_status <- assign_status_vector(dt$Horvath_residuals, lower_quantile, upper_quantile)
lihc_mask <- dt$project == "TCGA-LIHC" & !is.na(dt$global_status)
lihc_subset <- data.table(
  sample_id = sample_id_from_barcode(dt$id[lihc_mask]),
  classification = dt$global_status[lihc_mask]
)
lihc_subset <- unique(lihc_subset, by = "sample_id")
write_status_file(lihc_subset, "TCGA-LIHC_balanced_global", output_dir)

classify_tissue <- function(project_code) {
  mask <- dt$project == project_code & !is.na(dt$Horvath_residuals)
  if (!any(mask)) {
    stop(sprintf("No Horvath residuals available for %s", project_code))
  }
  residuals <- dt$Horvath_residuals[mask]
  statuses <- assign_status_vector(residuals, lower_quantile, upper_quantile)
  keep <- !is.na(statuses)
  if (!any(keep)) {
    stop(sprintf("No extreme samples found for %s", project_code))
  }
  subset <- data.table(
    sample_id = sample_id_from_barcode(dt$id[mask][keep]),
    classification = statuses[keep]
  )
  unique(subset, by = "sample_id")
}

brca_subset <- classify_tissue("TCGA-BRCA")
write_status_file(brca_subset, "TCGA-BRCA_tissue_quartiles", output_dir)

kirc_subset <- classify_tissue("TCGA-KIRC")
write_status_file(kirc_subset, "TCGA-KIRC_tissue_quartiles", output_dir)

message("Classification complete.")
