#!/usr/bin/env Rscript
# Orchestrates two differential-expression analyses:
# 1) TCGA-LIHC samples classified globally and balanced between acceleration groups.
# 2) TCGA-BRCA and TCGA-KIRC samples classified within each tissue separately.

suppressPackageStartupMessages({
  library(data.table)
  library(TCGAbiolinks)
})

acceleration_status <- project <- residual_value <- global_status <- id <- NULL

get_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_flag <- "--file="
  file_arg <- cmd_args[startsWith(cmd_args, file_flag)]
  if (length(file_arg)) {
    return(normalizePath(sub(file_flag, "", file_arg[length(file_arg)])))
  }
  normalizePath(getwd())
}

script_dir <- dirname(get_script_path())
repo_root <- normalizePath(file.path(script_dir, ".."))
get_counts_script <- file.path(repo_root, "scripts", "get_brca_rnaseq_by_barcode.R")
de_script <- file.path(repo_root, "scripts", "run_accel_de_analysis.R")

params <- list(
  predictions_path = file.path(repo_root, "results", "clock", "gdc_pan", "gdc_pancan_predictions.rds"),
  clock_residual_col = "Horvath_residuals",
  lower_quantile = 0.25,
  upper_quantile = 0.75,
  lihc_project = "TCGA-LIHC",
  brca_project = "TCGA-BRCA",
  kirc_project = "TCGA-KIRC",
  balanced_seed = 42L,
  output_dir = file.path(repo_root, "results", "differential_analysis"),
  classification_dir = file.path(repo_root, "results", "clock", "derived"),
  chunk_size = 50L,
  default_workflow = "STAR - Counts",
  workflow_overrides = list(),
  default_sample_codes = c("01"),
  sample_code_overrides = list(
    "TCGA-LIHC" = c("11", "01"),
    "TCGA-BRCA" = c("11", "01"),
    "TCGA-KIRC" = c("11", "01")
  )
)

override_params <- function(arg_vec, defaults) {
  if (!length(arg_vec)) return(defaults)
  for (arg in arg_vec) {
    if (!startsWith(arg, "--") || !grepl("=", arg)) next
    key_val <- sub("^--", "", arg)
    parts <- strsplit(key_val, "=", fixed = TRUE)[[1]]
    if (length(parts) != 2) next
    key <- parts[1]
    val <- parts[2]
    if (!nzchar(key) || !key %in% names(defaults)) next
    if (key %in% c("lower_quantile", "upper_quantile")) {
      defaults[[key]] <- as.numeric(val)
    } else if (key == "balanced_seed") {
      defaults[[key]] <- as.integer(val)
    } else {
      defaults[[key]] <- val
    }
  }
  defaults
}

assign_status_vector <- function(values, lower, upper) {
  status <- rep(NA_character_, length(values))
  ok <- !is.na(values)
  if (!any(ok)) return(status)
  qs <- stats::quantile(values[ok], probs = c(lower, upper), na.rm = TRUE)
  status[ok & values < qs[1]] <- "Accelerated"
  status[ok & values > qs[2]] <- "Decelerated"
  status
}

balance_status_groups <- function(dt, seed) {
  status_counts <- table(dt$acceleration_status)
  if (length(status_counts) < 2) {
    stop(sprintf(
      "Cannot balance groups: need both statuses, found: %s",
      paste(names(status_counts), collapse = ", ")
    ))
  }
  target <- min(status_counts)
  set.seed(seed)
  keep_idx <- unlist(lapply(names(status_counts), function(state) {
    idx <- which(dt$acceleration_status == state)
    if (length(idx) <= target) {
      idx
    } else {
      sample(idx, target)
    }
  }))
  dt[sort(keep_idx)]
}

write_vector_file <- function(values, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(values, con = path)
}

create_classification_rds <- function(ids, statuses, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  table <- data.table(
    id = ids,
    classification_score = ifelse(statuses == "Accelerated", 1, -1)
  )
  saveRDS(table, path)
}

subset_sample_types <- function(dt, allowed_codes) {
  if (!nrow(dt)) return(dt)
  dt[substr(dt$id, 14, 15) %in% allowed_codes]
}

resolve_sample_codes <- function(project_code, cfg) {
  if (!is.null(cfg$sample_code_overrides[[project_code]])) {
    return(cfg$sample_code_overrides[[project_code]])
  }
  cfg$default_sample_codes
}

make_paths <- function(label, params) {
  list(
    barcode = file.path(params$output_dir, sprintf("%s_sample_ids.txt", label)),
    counts = file.path(params$output_dir, sprintf("%s_rnaseq_counts.tsv", label)),
    metadata = file.path(params$output_dir, sprintf("%s_rnaseq_metadata.tsv", label)),
    classification = file.path(params$classification_dir, sprintf("%s_classification.rds", label))
  )
}

run_script <- function(script_path, args) {
  status <- system2("Rscript", args = c(script_path, args))
  if (!identical(status, 0L)) {
    stop(sprintf("Command failed for %s", script_path))
  }
}

resolve_workflow <- function(project_code, cfg) {
  if (!is.null(cfg$workflow_overrides[[project_code]])) {
    return(cfg$workflow_overrides[[project_code]])
  }
  cfg$default_workflow
}

available_cache <- new.env(parent = emptyenv())

fetch_available_samples <- function(project_code, workflow_type) {
  key <- sprintf("%s::%s", project_code, workflow_type)
  if (exists(key, envir = available_cache, inherits = FALSE)) {
    return(get(key, envir = available_cache, inherits = FALSE))
  }
  message(sprintf("[%s] Fetching available samples for '%s'", project_code, workflow_type))
  query <- TCGAbiolinks::GDCquery(
    project = project_code,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = workflow_type
  )
  results <- TCGAbiolinks::getResults(query)
  available_ids <- unique(results[["sample.submitter_id"]])
  available_ids <- available_ids[!is.na(available_ids)]
  if (!length(available_ids)) {
    stop(sprintf("No RNA-seq samples found for %s with workflow %s", project_code, workflow_type))
  }
  assign(key, available_ids, envir = available_cache)
  available_ids
}

intersect_with_available_samples <- function(case_dt, project_code, workflow_type) {
  if (!nrow(case_dt)) return(case_dt)
  available_ids <- fetch_available_samples(project_code, workflow_type)
  filtered_dt <- data.table::copy(case_dt)
  filtered_dt[, sample_short := substr(id, 1, 16)]
  keep <- filtered_dt$sample_short %in% available_ids
  dropped <- sum(!keep)
  if (dropped > 0) {
    message(sprintf(
      "[%s] Dropping %d sample(s) without '%s' counts",
      project_code,
      dropped,
      workflow_type
    ))
  }
  filtered_dt <- filtered_dt[keep]
  if (!nrow(filtered_dt)) {
    stop(sprintf(
      "No samples remain for %s after enforcing '%s' availability",
      project_code,
      workflow_type
    ))
  }
  filtered_dt[, sample_short := NULL]
  filtered_dt
}

run_case <- function(case_dt, label, project_code, params, workflow_type) {
  if (!nrow(case_dt)) {
    stop(sprintf("No samples available for %s", label))
  }
  case_dt <- unique(case_dt, by = "id")
  unique_status <- unique(case_dt$acceleration_status)
  if (length(unique_status) < 2) {
    stop(sprintf("Need both statuses for %s, found: %s", label, paste(unique_status, collapse = ", ")))
  }
  paths <- make_paths(label, params)
  write_vector_file(case_dt$id, paths$barcode)
  create_classification_rds(case_dt$id, case_dt$acceleration_status, paths$classification)
  message(sprintf("[%s] Downloading counts for %d samples", label, length(case_dt$id)))
  run_script(
    get_counts_script,
    c(
      paths$barcode,
      project_code,
      paths$counts,
      paths$metadata,
      as.character(params$chunk_size),
      shQuote(workflow_type)
    )
  )
  message(sprintf("[%s] Running differential analysis", label))
  run_script(de_script, c(paths$counts, paths$classification, label, "classification_score", params$output_dir))
}

execute_workflow <- function(arg_vec = commandArgs(trailingOnly = TRUE)) {
  cfg <- override_params(arg_vec, params)
  if (!file.exists(cfg$predictions_path)) {
    stop(sprintf("Predictions file not found: %s", cfg$predictions_path))
  }
  dt <- as.data.table(readRDS(cfg$predictions_path))
  if (!cfg$clock_residual_col %in% names(dt)) {
    stop(sprintf("Residual column '%s' missing in predictions", cfg$clock_residual_col))
  }
  dt <- dt[!is.na(dt$project)]
  dt$residual_value <- dt[[cfg$clock_residual_col]]
  dt$global_status <- assign_status_vector(dt$residual_value, cfg$lower_quantile, cfg$upper_quantile)

  # LIHC global balanced case
  lihc_idx <- dt$project == cfg$lihc_project & !is.na(dt$global_status)
  lihc_dt <- dt[lihc_idx, c("id", "global_status"), with = FALSE]
  if (!nrow(lihc_dt)) {
    stop("No LIHC samples available after global classification.")
  }
  data.table::setnames(lihc_dt, "global_status", "acceleration_status")
  lihc_dt <- subset_sample_types(lihc_dt, resolve_sample_codes(cfg$lihc_project, cfg))
  if (!nrow(lihc_dt)) {
    stop("No LIHC samples available after filtering to primary tumor codes.")
  }
  lihc_workflow <- resolve_workflow(cfg$lihc_project, cfg)
  lihc_dt <- intersect_with_available_samples(lihc_dt, cfg$lihc_project, lihc_workflow)
  balanced_lihc <- balance_status_groups(lihc_dt, cfg$balanced_seed)
  run_case(
    balanced_lihc,
    sprintf("%s_balanced_global", cfg$lihc_project),
    cfg$lihc_project,
    cfg,
    lihc_workflow
  )

  # BRCA per tissue
  brca_idx <- dt$project == cfg$brca_project & !is.na(dt$residual_value)
  brca_dt <- dt[brca_idx, c("id", "residual_value"), with = FALSE]
  brca_dt$acceleration_status <- assign_status_vector(brca_dt$residual_value, cfg$lower_quantile, cfg$upper_quantile)
  brca_case <- brca_dt[!is.na(brca_dt$acceleration_status), c("id", "acceleration_status"), with = FALSE]
  brca_case <- subset_sample_types(brca_case, resolve_sample_codes(cfg$brca_project, cfg))
  if (!nrow(brca_case)) {
    stop("No BRCA samples remain after filtering to primary tumor codes.")
  }
  brca_workflow <- resolve_workflow(cfg$brca_project, cfg)
  brca_case <- intersect_with_available_samples(brca_case, cfg$brca_project, brca_workflow)
  run_case(
    brca_case,
    sprintf("%s_tissue_quartiles", cfg$brca_project),
    cfg$brca_project,
    cfg,
    brca_workflow
  )

  # KIRC per tissue
  kirc_idx <- dt$project == cfg$kirc_project & !is.na(dt$residual_value)
  kirc_dt <- dt[kirc_idx, c("id", "residual_value"), with = FALSE]
  kirc_dt$acceleration_status <- assign_status_vector(kirc_dt$residual_value, cfg$lower_quantile, cfg$upper_quantile)
  kirc_case <- kirc_dt[!is.na(kirc_dt$acceleration_status), c("id", "acceleration_status"), with = FALSE]
  kirc_case <- subset_sample_types(kirc_case, resolve_sample_codes(cfg$kirc_project, cfg))
  if (!nrow(kirc_case)) {
    stop("No KIRC samples remain after filtering to primary tumor codes.")
  }
  kirc_workflow <- resolve_workflow(cfg$kirc_project, cfg)
  kirc_case <- intersect_with_available_samples(kirc_case, cfg$kirc_project, kirc_workflow)
  run_case(
    kirc_case,
    sprintf("%s_tissue_quartiles", cfg$kirc_project),
    cfg$kirc_project,
    cfg,
    kirc_workflow
  )

  message("Workflows complete.")
}

execute_workflow()
