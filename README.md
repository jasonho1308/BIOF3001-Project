# BIOF3001-Project

Investigating the molecular basis of accelerated epigenetic aging

## Running the Analysis

This project uses Snakemake to efficiently process epigenetic clock data for multiple TCGA projects. Snakemake tracks which outputs have been generated and only reprocesses missing or out-of-date files.

### Prerequisites

- R with required packages: TCGAbiolinks, methylclock, methylclockData, sesameData, ggplot2, SummarizedExperiment, patchwork, dplyr, tidyr, rlang, png, grid
- Snakemake (install with: `pip install snakemake` or `uv pip install snakemake`)

### Usage

Run the entire workflow:

```bash
snakemake --cores 1
```

Process specific projects only:

```bash
snakemake results/clock/gdc_pan/TCGA-BRCA_scatterplots.png --cores 1
```

Run with multiple cores (parallel processing):

```bash
snakemake --cores 4
```

Dry run to see what would be executed:

```bash
snakemake --dry-run
```

### Workflow Structure

The workflow consists of three main steps:

1. **Process individual projects** (`scripts/process_project.R`): Calculates epigenetic clocks and generates plots for each TCGA project independently
2. **Combine pan-cancer analysis** (`scripts/combine_pancan.R`): Aggregates all project predictions and generates pan-cancer wide plots
3. **Generate PDF report** (`scripts/generate_pdf.R`): Combines all plots into a single PDF document

### Outputs

- `results/clock/gdc_pan/{PROJECT}_scatterplots.png`: Scatter plots for each project
- `results/clock/gdc_pan/{PROJECT}_residuals_boxplots.png`: Residual boxplots for each project
- `results/clock/gdc_pan/{PROJECT}_predictions.rds`: Prediction data for each project
- `results/clock/gdc_pan/gdc_pancan_scatterplots.png`: Combined pan-cancer scatter plots
- `results/clock/gdc_pan/gdc_pancan_residuals_boxplots.png`: Combined pan-cancer residual boxplots
- `results/clock/gdc_pan/gdc_pancan_methylclock.pdf`: Final PDF report with all plots

### Benefits of Snakemake

- **Incremental processing**: Only reprocesses projects with missing outputs
- **Parallel execution**: Process multiple projects simultaneously with `--cores`
- **Automatic dependency tracking**: Ensures outputs are regenerated when inputs change
- **Resume capability**: Can resume interrupted analyses without restarting from scratch
