# BIOF3001-Project

Investigating the molecular basis of accelerated epigenetic aging

## Tasks distribution

| Name | Tasks |
| --- | --- |
| Jason | - the download of GDC-PANCAN data <br> - initial clock scripts, fix plot <br> - clinical correlation analysis <br> - snakemake, parallel processing and dvc implementation | 
| Jacob | - plot boxplots and scatterplots for different epigenetic clocks <br> - fix some projects having multiple platforms causing error when running clock <br> - other error fixes|
| William | - classify samples into accelerated and decelerated aging <br> - download corresponding rna sequence data <br> - differential expression analysis <br> - pathway analysis <br> |

## Running the Analysis

This project uses Snakemake to efficiently process epigenetic clock data for multiple TCGA projects. Snakemake tracks which outputs have been generated and only reprocesses missing or out-of-date files.

### Prerequisites

Run `library.R` to download all required libraries for R

Run `uv sync` to install all required libraries for python

For data management, we have used [dvc](https://dvc.org/) to store and share our raw and processed data, please contact Jason Ho for the server address. (just the raw data account for 1XGB, it takes at least 3 hours to download)

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

## Big Data Approach

This project implements several big data strategies to efficiently handle large-scale genomic data from TCGA (The Cancer Genome Atlas):

### Data Management

- **DVC (Data Version Control)**: Large raw data files (`.dvc` files) are tracked separately from the git repository, enabling efficient versioning and sharing of multi-gigabyte methylation datasets without bloating the repository
- **Selective Processing**: Only normal tissue samples are extracted from the full TCGA-PANCAN dataset, reducing computational overhead
- **Cached Intermediate Results**: Key data structures (`.rds` files) are saved to prevent redundant API calls and expensive recomputations

### Parallel Processing

- **Snakemake Workflow**: Orchestrates parallel execution across 24 TCGA projects, automatically managing dependencies and utilizing multiple CPU cores
- **Multi-core Downloads**: The `download_gdc_pancan.R` script uses R's `parallel::mclapply()` to download data from multiple projects concurrently, significantly reducing wall-clock time
- **Project-level Parallelism**: Each TCGA project is processed independently, allowing Snakemake to distribute work across available cores without waiting for sequential completion
- **Adaptive Core Usage**: Scripts automatically detect available CPU cores and use half to balance performance with system stability

### Scalability Features

- **Incremental Execution**: Snakemake tracks which outputs exist and only reprocesses missing or outdated files, enabling efficient iterative development
- **Chunk-based Downloads**: GDC downloads use `files.per.chunk = 20` to prevent timeout issues when fetching large batches of methylation files
- **Memory-efficient Processing**: Data is processed project-by-project rather than loading the entire pan-cancer dataset into memory at once
- **Error Recovery**: Parallel processing includes robust error handling to continue processing remaining projects even if individual projects fail

### Performance Optimization

- **Binary Serialization**: R's `.rds` format provides fast, compressed storage for intermediate results
- **Lazy Evaluation**: Phenotype and methylation data are only loaded when needed for specific analyses
- **Distributed I/O**: Multiple projects write to separate output files simultaneously, avoiding I/O bottlenecks
