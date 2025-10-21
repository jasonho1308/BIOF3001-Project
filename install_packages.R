#!/usr/bin/env Rscript
# Install Required Packages for BIOF3001 Project

cat("Installing required R packages for epigenetic aging analysis...\n\n")

# Install CRAN packages
cran_packages <- c(
  "dplyr",
  "tidyr",
  "readr",
  "readxl",
  "ggplot2",
  "corrplot",
  "pheatmap",
  "ggpubr",
  "car",
  "broom",
  "lme4",
  "multcomp",
  "here",
  "rmarkdown",
  "knitr",
  "rstudioapi"
)

cat("Installing CRAN packages...\n")
for (pkg in cran_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    cat(paste("Installed:", pkg, "\n"))
  } else {
    cat(paste("Already installed:", pkg, "\n"))
  }
}

# Install Bioconductor packages
cat("\nInstalling Bioconductor packages...\n")
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c(
  "minfi",
  "ChAMP",
  "wateRmelon",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"
)

for (pkg in bioc_packages) {
  if (!require(pkg, character.only = TRUE)) {
    BiocManager::install(pkg, update = FALSE)
    cat(paste("Installed:", pkg, "\n"))
  } else {
    cat(paste("Already installed:", pkg, "\n"))
  }
}

# Optional: Install methylclock package for epigenetic clocks
cat("\nInstalling methylclock package...\n")
if (!require("methylclock", quietly = TRUE)) {
  BiocManager::install("methylclock", update = FALSE)
}

cat("\nPackage installation complete!\n")
cat("You can now run the analysis scripts.\n")
