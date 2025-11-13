# Snakemake workflow for epigenetic clock analysis
# This workflow processes TCGA projects independently to avoid reprocessing all clocks

import os

# Get all TCGA projects from the queries file
PROJECTS = [
    "TCGA-BLCA", "TCGA-BRCA", "TCGA-CESC", "TCGA-CHOL", "TCGA-COAD",
    "TCGA-ESCA", "TCGA-GBM", "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP",
    "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-OV", "TCGA-PAAD",
    "TCGA-PCPG", "TCGA-PRAD", "TCGA-READ", "TCGA-SARC", "TCGA-SKCM",
    "TCGA-STAD", "TCGA-THCA", "TCGA-THYM", "TCGA-UCEC"
]

os.makedirs("logs/clock/gdc_pan", exist_ok=True)

# Final output
rule all:
    input:
        "results/clock/gdc_pan/gdc_pancan_methylclock.pdf"

# Process individual project and generate predictions + plots
rule process_project:
    input:
        queries = "data/processed/gdc_pancan/queries.rds",
        data = "data/raw/GDCdata/"
    output:
        scatter = "results/clock/gdc_pan/{project}_scatterplots.png",
        residuals = "results/clock/gdc_pan/{project}_residuals_boxplots.png",
        predictions = "results/clock/gdc_pan/{project}_predictions.rds"
    log:
        "logs/clock/gdc_pan/{project}.log"
    shell:
        """
        Rscript scripts/process_project.R {wildcards.project} > {log} 2>&1
        """

# Combine all project predictions and generate pancan-wide analysis
rule combine_pancan:
    input:
        predictions = expand("results/clock/gdc_pan/{project}_predictions.rds", project=PROJECTS)
    output:
        scatter = "results/clock/gdc_pan/gdc_pancan_scatterplots.png",
        residuals = "results/clock/gdc_pan/gdc_pancan_residuals_boxplots.png",
        combined = "results/clock/gdc_pan/gdc_pancan_predictions.rds"
    log:
        "logs/clock/gdc_pan/pancan.log"
    shell:
        """
        Rscript scripts/combine_pancan.R > {log} 2>&1
        """

# Generate final PDF combining all plots
rule generate_pdf:
    input:
        scatter_plots = expand("results/clock/gdc_pan/{project}_scatterplots.png", project=PROJECTS),
        residual_plots = expand("results/clock/gdc_pan/{project}_residuals_boxplots.png", project=PROJECTS),
        pancan_scatter = "results/clock/gdc_pan/gdc_pancan_scatterplots.png",
        pancan_residuals = "results/clock/gdc_pan/gdc_pancan_residuals_boxplots.png"
    output:
        "results/clock/gdc_pan/gdc_pancan_methylclock.pdf"
    log:
        "logs/clock/gdc_pan/generate_pdf.log"
    shell:
        """
        Rscript scripts/generate_pdf.R > {log} 2>&1
        """
