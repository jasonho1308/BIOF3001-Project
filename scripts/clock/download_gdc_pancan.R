# To download GDC-PANCAN phenotype data from UCSC Xena and
# then download DNA Methylation Beta Value data for normal tissue samples from GDC.

# Setup logging if running via Snakemake
if (exists("snakemake")) {
  log_file <- file(snakemake@log[[1]], open = "wt")
  sink(log_file, type = "output", split = FALSE)
  sink(log_file, type = "message", split = FALSE)
}

# Get output paths from snakemake or use defaults
if (exists("snakemake")) {
  pheno_output <- snakemake@output$pheno
  counts_output <- snakemake@output$counts
  queries_output <- snakemake@output$queries
  missing_output <- snakemake@output$missing
} else {
  pheno_output <- "data/processed/gdc_pancan/normal_pheno.tsv"
  counts_output <- "data/processed/gdc_pancan/normal_counts_by_project.csv"
  queries_output <- "data/processed/gdc_pancan/queries.rds"
  missing_output <- "data/processed/gdc_pancan/missing_ids.rds"
}

library(UCSCXenaTools)
library(dplyr)
library(TCGAbiolinks)
library(parallel)

# Generate phenotype query
xe_pheno <- XenaGenerate(subset = XenaHostNames == "gdcHub") %>%
  XenaFilter(filterDatasets = "GDC-PANCAN") %>%
  XenaFilter(filterDatasets = "TCGA_phenotype") %>%
  XenaFilter(filterDatasets = "tsv")

# Create query object
pheno_query <- XenaQuery(xe_pheno)

# Download and prepare phenotype data
pheno_download <- XenaDownload(pheno_query, destdir = "data/raw/phenotype")
pheno_data <- XenaPrepare(pheno_download)

print(dim(pheno_data))

# Filter for normal tissue, finding if "normal" in sample type
normal_samples <- pheno_data %>% filter(grepl("normal", samples.sample_type, ignore.case = TRUE))


# make sure output directory exists
output_dir <- dirname(pheno_output)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# save normal samples into tsv
write.table(normal_samples,
  file = pheno_output, sep = "\t",
  row.names = FALSE, quote = FALSE
)

# function to split tibble by column
split_tibble <- function(tibble, col = "col") tibble %>% split(., .[, col])

normal_samples_by_project <- split_tibble(normal_samples, "project.project_id")

# count samples by project
sample_counts <- sapply(normal_samples_by_project, nrow)
print(sample_counts)

write.csv(as.data.frame(sample_counts),
  file = counts_output,
  row.names = TRUE,
  quote = FALSE
)

# Detect number of cores (use half to avoid overwhelming the system)
n_cores <- max(1, floor(detectCores() / 2))
message(sprintf("Using %d cores for parallel processing", n_cores))

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

      # Construct GDC query for DNA Methylation Beta Value data
      query_TCGA_meth <- GDCquery(
        project       = proj,
        data.category = "DNA Methylation",
        data.type     = "Methylation Beta Value",
        barcode       = ids
      )

      result$query <- query_TCGA_meth

      # Download data files from GDC
      # files.per.chunk = 20 to avoid timeout issues
      GDCdownload(query_TCGA_meth, directory = "data/raw/GDCdata", files.per.chunk = 20)
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

# save missing ids, gdc queries to processed data
saveRDS(missing_ids, file = missing_output)
saveRDS(gdc_queries, file = queries_output)

# Close log file if running via Snakemake
if (exists("snakemake")) {
  sink(type = "output")
  sink(type = "message")
  close(log_file)
}
