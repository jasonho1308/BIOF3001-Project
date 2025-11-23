# load required libraries
packages <- c("BiocManager", "data.table", "tidyverse", "dplyr", "ggplot2", "patchwork", "DT", "tidyr")
not_installed <- setdiff(packages, rownames(installed.packages()))
if (length(not_installed)) install.packages(not_installed)
lapply(packages, library, character.only = TRUE)

if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install(version = "3.21")

if (!requireNamespace("DESeq2", quietly = TRUE)) {
    BiocManager::install("DESeq2")
}

if (!requireNamespace("methylclock", quietly = TRUE)) {
    BiocManager::install("methylclock")
    BiocManager::install("methylclockData")
}

if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    BiocManager::install("clusterProfiler")
}
if (!requireNamespace("msigdbr", quietly = TRUE)) {
    BiocManager::install("msigdbr")
}
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    BiocManager::install("org.Hs.eg.db")
}
if (!requireNamespace("enrichplot", quietly = TRUE)) {
    BiocManager::install("enrichplot")
}
if (!requireNamespace("sesame", quietly = TRUE)) {
    BiocManager::install("sesame")
}
if (!requireNamespace("sesameData", quietly = TRUE)) {
    BiocManager::install("sesameData")
}
if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) BiocManager::install("TCGAbiolinks")
if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) BiocManager::install("SummarizedExperiment")
if (!requireNamespace("UCSCXenaTools", quietly = TRUE)) {
    install.packages("UCSCXenaTools")
}
if (!requireNamespace("pheatmap", quietly = TRUE)) {
    install.packages("pheatmap")
}
