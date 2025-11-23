#!/usr/bin/env Rscript
# Utility script for downloading TCGA RNA-seq counts for a set of BRCA barcodes

ensure_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE)
  )
}

packages <- c("data.table", "TCGAbiolinks", "SummarizedExperiment")
invisible(lapply(packages, ensure_pkg))

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
  TCGAbiolinks::GDCprepare(query)
}

args <- commandArgs(trailingOnly = TRUE)
barcode_file <- if (length(args) >= 1) args[[1]] else "results/differential_analysis/TCGA-BRCA_dec.txt"
project_code <- if (length(args) >= 2) args[[2]] else "TCGA-BRCA"
counts_out <- if (length(args) >= 3) args[[3]] else "results/differential_analysis/TCGA-BRCA_dec_rnaseq_counts.tsv"
metadata_out <- if (length(args) >= 4) args[[4]] else "results/differential_analysis/TCGA-BRCA_dec_rnaseq_metadata.tsv"
chunk_size <- if (length(args) >= 5) as.integer(args[[5]]) else 50
workflow_type <- if (length(args) >= 6) args[[6]] else "HTSeq - Counts"

if (!file.exists(barcode_file)) {
  stop(sprintf("Barcode file not found: %s", barcode_file))
}

barcodes <- trimws(readLines(barcode_file, warn = FALSE))
barcodes <- barcodes[nchar(barcodes) > 0]
if (!length(barcodes)) {
  stop("No barcodes detected in input file.")
}

sample_ids <- unique(substr(barcodes, 1, 16))
message(sprintf("Detected %d unique sample IDs", length(sample_ids)))

chunks <- split(sample_ids, ceiling(seq_along(sample_ids) / chunk_size))
message(sprintf("Processing %d chunk(s)", length(chunks)))
message(sprintf("Using workflow type: %s", workflow_type))

fetch_chunk <- function(ids, idx, total) {
  message(sprintf("Querying chunk %d/%d (%d samples)", idx, total, length(ids)))
  query <- tryCatch({
    TCGAbiolinks::GDCquery(
      project = project_code,
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      workflow.type = workflow_type,
      barcode = ids
    )
  }, error = function(e) {
    stop(sprintf("GDCquery failed for chunk %d: %s", idx, e$message))
  })
  TCGAbiolinks::GDCdownload(query)
  tryCatch(
    TCGAbiolinks::GDCprepare(query),
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

se_list <- vector("list", length(chunks))
for (i in seq_along(chunks)) {
  se_list[[i]] <- fetch_chunk(chunks[[i]], i, length(chunks))
}

if (length(se_list) == 0) {
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
coldata_dt <- data.table::rbindlist(coldata_list, fill = TRUE, use.names = TRUE)
coldata_dt[, ord := match(sample_id, colnames(counts_mat))]
setorder(coldata_dt, ord)
coldata_dt[, ord := NULL]
setcolorder(coldata_dt, c("sample_id", setdiff(names(coldata_dt), "sample_id")))

flatten_lists <- function(dt) {
  list_cols <- names(dt)[vapply(dt, is.list, logical(1L))]
  if (!length(list_cols)) return(dt)
  for (col in list_cols) {
    dt[[col]] <- vapply(dt[[col]], function(x) {
      if (length(x) == 0 || all(is.na(x))) {
        return(NA_character_)
      }
      paste(unique(as.character(unlist(x))), collapse = ";")
    }, character(1))
  }
  dt
}

coldata_dt <- flatten_lists(coldata_dt)

counts_dt <- data.table::as.data.table(counts_mat, keep.rownames = "gene_id")

dir.create(dirname(counts_out), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(metadata_out), recursive = TRUE, showWarnings = FALSE)

fwrite(counts_dt, counts_out, sep = "\t")
fwrite(coldata_dt, metadata_out, sep = "\t")

message(sprintf("Counts written to %s", counts_out))
message(sprintf("Metadata written to %s", metadata_out))
message("Download complete.")
