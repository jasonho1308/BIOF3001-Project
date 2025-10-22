barcodes_small <- c("TCGA-E9-A1N8-11A", "TCGA-E2-A1BC-11A", "TCGA-E9-A1ND-11A", "TCGA-AC-A2FM-11B")

# This line searches in the database and builds the query, selecting the beta values data from the TCGA project, pulling only the files for the 4 samples. This will take 2-3 minutes
query_TCGA_meth <- GDCquery(
  project       = "TCGA-BRCA",
  data.category = "DNA Methylation",
  data.type     = "Methylation Beta Value",
  barcode       = barcodes_small # <— restrict to just these 4 samples
)

# To extract and view the metadata of the data we'll be downloading:
res_TCGA_meth <- getResults(query_TCGA_meth)
View(res_TCGA_meth) # metadata

# To download the actual data files from the GDC based on the previously constructed query:
GDCdownload(query_TCGA_meth, directory = "data/raw/GDCdata/")
me_TCGA <- GDCprepare(query_TCGA_meth, directory = "data/raw/GDCdata/")
dim(me_TCGA) # methylation x 4 samples

# prepare data for methylclock

data <- assay(me_TCGA)
data <- na.omit(data)
samples <- colnames(data)

# ordered.metadata <- metadata[match(samples, metadata$Sample), ]
age <- me_TCGA$age_at_index

# checks for missing CpGs for the various clocks
cpgs.missing <- checkClocks(data) # approx 1 min to run

# makes prediction based on methylation data of normal samples
prediction <- DNAmAge(data, age = age, cell.count = FALSE) # approx <2 mins to run

# Scatter plots with regression lines using tidy evaluation
chron_age <- prediction$age
clock_cols <- colnames(prediction[, c(2, 5, 8, 12, 15, 18, 21, 24, 27)])

scatter_plots <- lapply(clock_cols, function(clock) {
  ggplot(prediction, aes(x = age, y = !!sym(clock))) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm", color = "red") +
    labs(title = clock, x = "Chronological Age", y = "Predicted DNAmAge") +
    theme_minimal()
})

scatter_grid <- wrap_plots(scatter_plots, ncol = 3)

print(scatter_grid)
