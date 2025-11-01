library(UCSCXenaTools)
library(dplyr)
library(TCGAbiolinks)


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
if (!dir.exists("data/processed/gdc_pancan")) {
  dir.create("data/processed/gdc_pancan", recursive = TRUE)
}

# save normal samples into tsv
write.table(normal_samples,
  file = "data/processed/gdc_pancan/normal_pheno.tsv", sep = "\t",
  row.names = FALSE, quote = FALSE
)

# function to split tibble by column
split_tibble <- function(tibble, col = "col") tibble %>% split(., .[, col])

normal_samples_by_project <- split_tibble(normal_samples, "project.project_id")

# count each samples by project
sample_counts <- sapply(normal_samples_by_project, nrow)
print(sample_counts)

write.csv(as.data.frame(sample_counts),
  file = "data/processed/gdc_pancan/normal_counts_by_project.csv",
  row.names = TRUE,
  quote = FALSE
)

gdc_queries <- list()
missing_ids <- list()
for (proj in names(normal_samples_by_project)) {
  ids <- unlist(normal_samples_by_project[[proj]]["sample"])
  tryCatch(
    {
      print(paste("Processing project:", proj))

      # Construct GDC query for DNA Methylation Beta Value data for the given project and sample barcodes
      query_TCGA_meth <- GDCquery(
        project       = proj,
        data.category = "DNA Methylation",
        data.type     = "Methylation Beta Value",
        barcode       = ids
      )

      gdc_queries[[proj]] <- query_TCGA_meth

      # To download the actual data files from the GDC based on the previously constructed query:
      # files.per.chunk = 20 to avoid timeout issues
      GDCdownload(query_TCGA_meth, directory = "data/raw/GDCdata", files.per.chunk = 20)
      print(paste("Downloaded data for project:", proj))
    },
    error = function(e) {
      message(paste("Error processing project:", proj))
      message(e)

      # save all missing ids to a list, for logging by project
      # <<- to modify the variable outside the function scope
      missing_ids[[proj]] <<- ids
    }
  )
}

# save missing ids, gdc queries to processed data
saveRDS(missing_ids, file = "data/processed/gdc_pancan/missing_ids.rds")
saveRDS(gdc_queries, file = "data/processed/gdc_pancan/queries.rds")