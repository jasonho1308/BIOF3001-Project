#!/usr/bin/env Rscript
# Title: RNA-seq downloader for extreme acceleration cohorts
# Description: Reads the Accelerated/Decelerated sample sheets produced by
#   classify_extreme_samples.R, extracts the TCGA sample IDs (first four barcode
#   fields), and retrieves STAR (or chosen workflow) count matrices plus metadata
#   directly via TCGAbiolinks. All raw downloads are cached under data/GDCdata to
#   keep GEO/GDC assets within the data directory tree.

suppressPackageStartupMessages({
  library(data.table)
  library(TCGAbiolinks)
  library(SummarizedExperiment)
})

utils::globalVariables(c("sample_id", "ord"))
preferred_sample_types <- c("01", "03", "06", "11", "12", "02", "10", "13")
available_barcode_cache <- new.env(parent = emptyenv())

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
gdc_storage_dir <- if (length(args) >= 6) args[[6]] else file.path("data", "GDCrnaseq")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(id_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gdc_storage_dir, recursive = TRUE, showWarnings = FALSE)

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
  ids <- if ("sample_id" %in% names(dt)) dt[["sample_id"]] else dt[[1]]
  ids <- trimws(ids)
  ids <- ids[nchar(ids) > 0]
  unique(ids)
}

flatten_lists <- function(dt) {
  list_cols <- names(dt)[vapply(dt, is.list, logical(1L))]
  if (!length(list_cols)) return(dt)
  for (col in list_cols) {
    dt[[col]] <- vapply(dt[[col]], function(x) {
      if (length(x) == 0 || all(is.na(x))) return(NA_character_)
      paste(unique(as.character(unlist(x))), collapse = ";")
    }, character(1))
  }
  dt
}

get_available_barcodes <- function(project, workflow) {
  key <- sprintf("%s::%s", project, workflow)
  if (exists(key, envir = available_barcode_cache, inherits = FALSE)) {
    return(get(key, envir = available_barcode_cache, inherits = FALSE))
  }
  message(sprintf("Fetching available barcodes for %s (%s)", project, workflow))
  query <- TCGAbiolinks::GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = workflow
  )
  res <- TCGAbiolinks::getResults(query)
  ids <- unique(res$cases)
  ids <- ids[nchar(ids) > 0]
  assign(key, ids, envir = available_barcode_cache)
  ids
}

rank_barcode <- function(barcode) {
  parts <- strsplit(barcode, "-", fixed = TRUE)[[1]]
  if (length(parts) < 4) {
    return(length(preferred_sample_types) + 1L)
  }
  sample_type <- substr(parts[4], 1, 2)
  rank <- match(sample_type, preferred_sample_types)
  if (is.na(rank)) length(preferred_sample_types) + 1L else rank
}

choose_barcode <- function(matches) {
  if (!length(matches)) return(NA_character_)
  ranks <- vapply(matches, rank_barcode, integer(1L))
  matches[order(ranks, matches)][1]
}

expand_to_available_barcodes <- function(sample_ids, project, workflow) {
  prefixes <- unique(substr(sample_ids, 1, 16))
  available <- get_available_barcodes(project, workflow)
  resolved <- vapply(prefixes, function(prefix) {
    matches <- available[substr(available, 1, 16) == prefix]
    choose_barcode(matches)
  }, character(1))
  list(
    expanded = resolved[!is.na(resolved)],
    missing = prefixes[is.na(resolved)],
    mapping = data.table(
      sample_id = prefixes,
      aliquot_barcode = resolved
    )
  )
}

prepare_with_fallback <- function(query) {
  ns <- asNamespace("TCGAbiolinks")
  original_coldata <- get("colDataPrepare", envir = ns)
  safe_wrapper <- function(barcode) {
    tryCatch(
      original_coldata(barcode),
      error = function(e) {
        message("colDataPrepare fallback applied: ", e$message)
        S4Vectors::DataFrame(
          barcode = barcode,
          patient = substr(barcode, 1, 12),
          sample = substr(barcode, 1, 16)
        )
      }
    )
  }
  unlockBinding("colDataPrepare", ns)
  on.exit({
    unlockBinding("colDataPrepare", ns)
    assign("colDataPrepare", original_coldata, envir = ns)
    lockBinding("colDataPrepare", ns)
  }, add = TRUE)
  assign("colDataPrepare", safe_wrapper, envir = ns)
  lockBinding("colDataPrepare", ns)
  TCGAbiolinks::GDCprepare(query, directory = gdc_storage_dir)
}

download_rnaseq_for_samples <- function(sample_ids, project, workflow, chunk_size) {
  cleaned_ids <- trimws(sample_ids)
  cleaned_ids <- cleaned_ids[nchar(cleaned_ids) > 0]
  if (!length(cleaned_ids)) {
    stop("No barcodes detected in input file.")
  }
  sample_ids_unique <- unique(substr(cleaned_ids, 1, 16))
  message(sprintf("Detected %d unique sample IDs", length(sample_ids_unique)))

  expanded <- expand_to_available_barcodes(sample_ids_unique, project, workflow)
  matched_barcodes <- unique(expanded$expanded)
  if (!length(matched_barcodes)) {
    stop(sprintf("No matching aliquot barcodes were found for %s", project))
  }
  if (length(expanded$missing)) {
    warning(sprintf(
      "%s: %d sample(s) lacked matching aliquots (first few: %s)",
      project,
      length(expanded$missing),
      paste(head(expanded$missing, 5), collapse = ", ")
    ))
  }
  message(sprintf(
    "Matched %d/%d sample IDs to aliquot barcodes",
    length(matched_barcodes),
    length(sample_ids_unique)
  ))

  chunks <- split(matched_barcodes, ceiling(seq_along(matched_barcodes) / chunk_size))
  message(sprintf("Processing %d chunk(s) for project %s", length(chunks), project))
  se_list <- vector("list", length(chunks))

  for (idx in seq_along(chunks)) {
    ids <- chunks[[idx]]
    message(sprintf("Querying chunk %d/%d (%d samples)", idx, length(chunks), length(ids)))
    query <- tryCatch(
      TCGAbiolinks::GDCquery(
        project = project,
        data.category = "Transcriptome Profiling",
        data.type = "Gene Expression Quantification",
        workflow.type = workflow,
        barcode = ids
      ),
      error = function(e) {
        stop(sprintf("GDCquery failed for chunk %d: %s", idx, e$message))
      }
    )
    TCGAbiolinks::GDCdownload(query, directory = gdc_storage_dir)
    se_list[[idx]] <- tryCatch(
      TCGAbiolinks::GDCprepare(query, directory = gdc_storage_dir),
      error = function(err) {
        message(sprintf(
          "Standard GDCprepare failed for chunk %d: %s. Retrying with minimal clinical data.",
          idx,
          err$message
        ))
        prepare_with_fallback(query)
      }
    )
  }

  if (!length(se_list) || all(vapply(se_list, is.null, logical(1L)))) {
    stop("No SummarizedExperiment objects retrieved; aborting.")
  }

  assay_list <- lapply(se_list, SummarizedExperiment::assay)
  gene_order <- rownames(assay_list[[1]])
  assay_list <- lapply(assay_list, function(mat) {
    if (!identical(rownames(mat), gene_order)) {
      mat <- mat[gene_order, , drop = FALSE]
    }
    mat
  })
  counts_mat <- do.call(cbind, assay_list)
  if (is.null(rownames(counts_mat))) {
    stop("Counts matrix is missing gene identifiers.")
  }

  coldata_list <- lapply(se_list, function(se) {
    df <- as.data.frame(SummarizedExperiment::colData(se))
    df$sample_id <- colnames(se)
    df
  })
  coldata_dt <- rbindlist(coldata_list, fill = TRUE, use.names = TRUE)
  order_idx <- match(coldata_dt$sample_id, colnames(counts_mat))
  coldata_dt <- coldata_dt[order(order_idx)]
  setcolorder(coldata_dt, c("sample_id", setdiff(names(coldata_dt), "sample_id")))
  coldata_dt <- flatten_lists(coldata_dt)

  counts_dt <- as.data.table(counts_mat, keep.rownames = "gene_id")
  list(
    counts = counts_dt,
    metadata = coldata_dt,
    barcodes = matched_barcodes,
    mapping = expanded$mapping
  )
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
  barcode_path <- file.path(id_dir, sprintf("%s_extreme_sample_barcodes.txt", label))
  mapping_path <- file.path(id_dir, sprintf("%s_extreme_sample_mapping.tsv", label))
  counts_path <- file.path(output_dir, sprintf("%s_extreme_rnaseq_counts.tsv", label))
  metadata_path <- file.path(output_dir, sprintf("%s_extreme_rnaseq_metadata.tsv", label))
  message(sprintf("[%s] Downloading counts for %d samples", label, length(sample_ids)))

  result <- tryCatch(
    download_rnaseq_for_samples(sample_ids, project, workflow_type, chunk_size),
    error = function(e) {
      warning(sprintf("[%s] Download failed: %s", label, e$message))
      NULL
    }
  )

  if (is.null(result)) {
    failed_cases[[label]] <- NA_character_
    next
  }

  fwrite(result$counts, counts_path, sep = "\t")
  fwrite(result$metadata, metadata_path, sep = "\t")
  writeLines(result$barcodes, con = barcode_path)
  result$mapping[, status := ifelse(is.na(aliquot_barcode), "missing", "matched")]
  fwrite(result$mapping, mapping_path, sep = "\t")
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
