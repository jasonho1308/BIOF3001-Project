# BIOF3001-Project

Investigating the molecular basis of accelerated epigenetic aging

## Tasks distribution

| Name | Tasks |
| --- | --- |
| Jason | - the download of GDC-PANCAN data <br> - initial clock scripts, fix plot <br> - clinical correlation analysis <br> - snakemake, parallel processing and dvc implementation | 
| Jacob | |
| William | |

## Running the Analysis

This project uses Snakemake to efficiently process epigenetic clock data for multiple TCGA projects. Snakemake tracks which outputs have been generated and only reprocesses missing or out-of-date files.

### Prerequisites

Run `library.R` to download all required libraries for R

Run `uv sync` to install all required libraries for python

For data management, we have used [dvc](https://dvc.org/) to store and share our data, please contact Jason Ho for the server address.

For workflow, we have used [snakemake](https://snakemake.readthedocs.io/en/stable/) to ensure reproducible results.

### Usage

Run the entire workflow:

```bash
uv run snakemake --cores 1
```
N: for number of cores

Process specific projects only:

```bash
uv run snakemake results/clock/gdc_pan/TCGA-BRCA_scatterplots.png --cores 1
```

Dry run to see what would be executed:

```bash
uv run snakemake -n
```

### Workflow Structure

The workflow consists of the following steps:

1. **Process individual projects** (`scripts/clock/process_project.R`): Calculates epigenetic clocks and generates plots for each TCGA project independently
2. **Combine pan-cancer analysis** (`scripts/clock/combine_pancan.R`): Aggregates all project predictions and generates pan-cancer wide plots
3. **Generate PDF report** (`scripts/clock/generate_pdf.R`): Combines all plots into a single PDF document
4. **Clinical correlation analysis** (`scripts/clinical_correlation.R`): Performs clinical correlation analysis

### Outputs

- `results/clock/gdc_pan/{PROJECT}_scatterplots.png`: Scatter plots for each project
- `results/clock/gdc_pan/{PROJECT}_residuals_boxplots.png`: Residual boxplots for each project
- `results/clock/gdc_pan/{PROJECT}_predictions.rds`: Prediction data for each project
- `results/clock/gdc_pan/gdc_pancan_scatterplots.png`: Combined pan-cancer scatter plots
- `results/clock/gdc_pan/gdc_pancan_residuals_boxplots.png`: Combined pan-cancer residual boxplots
- `results/clock/gdc_pan/gdc_pancan_methylclock.pdf`: Final PDF report with all plots
- `results/clinical/*`: Heat map for clinical correlation

> some `*.rds` was saved to prevent excess api usage and avoid recomputing of variables
