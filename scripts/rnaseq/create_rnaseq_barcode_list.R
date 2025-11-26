# To download RNA-seq Gene Expression data for normal tissue samples from GDC.

# Setup logging if running via Snakemake
if (exists("snakemake")) {
  log_file <- file(snakemake@log[[1]], open = "wt")
  sink(log_file, type = "output", split = FALSE)
  sink(log_file, type = "message", split = FALSE)
}

library(dplyr)
library(TCGAbiolinks)
library(parallel)

# Load predictions (samples that have methylation data with clock predictions)
predictions <- readRDS("results/clock/gdc_pan/gdc_pancan_predictions.rds")
message(sprintf("Loaded %d samples with methylation clock predictions", nrow(predictions)))

# Extract sample submitter id from id
predictions$sample <- sapply(predictions$id, function(bc) {
  parts <- strsplit(bc, "-", fixed = TRUE)[[1]]
  if (length(parts) >= 4) {
    paste(parts[1:4], collapse = "-")
  } else {
    substr(bc, 1, 16)
  }
})

# Split samples by project
normal_samples_by_project <- split(predictions, predictions$project)

# count samples by project
sample_counts <- sapply(normal_samples_by_project, nrow)
print(sample_counts)

# Detect number of cores (use half to avoid overwhelming the system)
n_cores <- max(1, floor(detectCores() / 2))
message(sprintf("Using %d cores for parallel processing", n_cores))

# make sure output directory exists
if (!dir.exists("data/rna")) {
  dir.create("data/rna", recursive = TRUE)
}

# Function to process a single project
process_project <- function(proj, samples_by_project) {
  ids <- unlist(samples_by_project[[proj]]["sample"])

  result <- list(
    project = proj,
    query = NULL,
    missing_ids = NULL,
    success = FALSE
  )

  tryCatch(
    {
      message(paste("Processing project:", proj))

      # Construct GDC query for RNA-seq Gene Expression data
      query_TCGA_rnaseq <- GDCquery(
        project       = proj,
        data.category = "Transcriptome Profiling",
        data.type     = "Gene Expression Quantification",
        workflow.type = "STAR - Counts",
        barcode       = ids
      )

      result$query <- query_TCGA_rnaseq

      # Download data files from GDC
      # files.per.chunk = 20 to avoid timeout issues
      GDCdownload(query_TCGA_rnaseq, directory = "data/rna", files.per.chunk = 20)
      message(paste("Downloaded data for project:", proj))

      result$success <- TRUE
    },
    error = function(e) {
      message(paste("Error processing project:", proj))
      message(e)
      
      # used <<- to assign to outer scope
      result$missing_ids <<- ids
    }
  )

  result
}

# Process projects in parallel
# mc.preschedule = FALSE ensures each job completes independently
results <- mclapply(
  names(normal_samples_by_project),
  process_project,
  samples_by_project = normal_samples_by_project,
  mc.cores = n_cores,
  mc.preschedule = FALSE
)

# Collect queries and missing IDs from results
gdc_queries <- list()
missing_ids <- list()

for (result in results) {
  # Skip if result is an error or NULL
  if (inherits(result, "try-error") || is.null(result)) {
    next
  }

  # Check if result is a list with expected structure
  if (is.list(result) && !is.null(result$project)) {
    if (!is.null(result$query)) {
      gdc_queries[[result$project]] <- result$query
    }
    if (!is.null(result$missing_ids)) {
      missing_ids[[result$project]] <- result$missing_ids
    }
  }
}

# Total number of normal samples
total_samples <- sum(sapply(gdc_queries, function(q) {
  if (!is.null(q) && !is.null(q$results[[1]])) nrow(q$results[[1]]) else 0
}))
cat("Total normal samples with RNA-seq downloaded:", total_samples, "\n")

# Samples per project
samples_per_project <- sapply(gdc_queries, function(q) {
  if (!is.null(q) && !is.null(q$results[[1]])) nrow(q$results[[1]]) else 0
})
print(samples_per_project)

# save missing ids, gdc queries to processed data
saveRDS(missing_ids, file = "data/rna/missing_ids.rds")
saveRDS(gdc_queries, file = "data/rna/queries.rds")

# Close log file if running via Snakemake
if (exists("snakemake")) {
  sink(type = "output")
  sink(type = "message")
  close(log_file)
}


