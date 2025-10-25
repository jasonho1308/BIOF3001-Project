library(UCSCXenaTools)
library(dplyr)
library(TCGAbiolinks)


# Generate phenotype query
XenaGenerate(subset = XenaHostNames == "gdcHub") %>%
  XenaFilter(filterDatasets = "PAN") %>%
  XenaFilter(filterDatasets = "TCGA_phenotype") %>%
  XenaFilter(filterDatasets = "tsv") -> xe_pheno

pheno_query <- XenaQuery(xe_pheno)

pheno_download <- XenaDownload(pheno_query, destdir = "data/raw/phenotype") 
pheno_data <- XenaPrepare(pheno_download)

print(dim(pheno_data))

# Filter for normal tissue, finding if "normal" in sample type
normal_samples <- pheno_data %>% filter(grepl("normal", samples.sample_type, ignore.case = TRUE))

# save normal samples into tsv
write.table(normal_samples,
  file = "data/processed/gdc_pancan_normal_pheno.tsv", sep = "\t",
  row.names = FALSE, quote = FALSE
)

split_tibble <- function(tibble, col = 'col') tibble %>% split(., .[, col])

normal_samples_by_project <- split_tibble(normal_samples, 'project.project_id')

# initialize missing ids list
missing_ids <- list()

for (proj in names(normal_samples_by_project)) {
  ids <- unlist(normal_samples_by_project[[proj]]["sample"])
  tryCatch({
    print(paste("Processing project:", proj))
    query_TCGA_meth <- GDCquery(
      project       = proj, 
      data.category = "DNA Methylation", 
      data.type     = "Methylation Beta Value",
      barcode       = ids
    )

    # To download the actual data files from the GDC based on the previously constructed query: 
    GDCdownload(query_TCGA_meth, directory = "data/raw/GDCdata", files.per.chunk = 20) 
    print(paste("Downloaded data for project:", proj))               
  }, error = function(e) {
    message(paste("Error processing project:", proj))
    message(e)

    # save all missing ids to a list, for logging by project
    missing_ids[[proj]] <- ids
  })
}

# save missing ids log
saveRDS(missing_ids, file = "data/processed/gdc_pancan_missing_ids.rds")
