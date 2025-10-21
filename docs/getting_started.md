# Getting Started Guide

## Quick Start

This guide will help you get started with the BIOF3001 epigenetic aging analysis project.

## Prerequisites

1. **Install R and RStudio**
   - Download R from: https://cran.r-project.org/
   - Download RStudio from: https://www.rstudio.com/products/rstudio/download/
   - Recommended: R version 4.0 or higher

2. **Clone the Repository**
   ```bash
   git clone https://github.com/jasonho1308/BIOF3001-Project.git
   cd BIOF3001-Project
   ```

3. **Open the Project in RStudio**
   - Double-click `BIOF3001-Project.Rproj` to open in RStudio
   - Or in RStudio: File > Open Project > Select BIOF3001-Project.Rproj

## Step-by-Step Setup

### 1. Install Required Packages

Run the installation script to install all required R packages:

```r
source("install_packages.R")
```

This will install:
- Data manipulation packages (dplyr, tidyr, readr)
- Visualization packages (ggplot2, corrplot, pheatmap)
- Statistical packages (car, broom, lme4)
- Bioconductor packages (minfi, ChAMP, wateRmelon)
- Epigenetic clock packages (methylclock)

**Note:** Installation may take 15-30 minutes depending on your system.

### 2. Prepare Your Data

Place your raw methylation data in the appropriate directory:

```
data/
├── raw/              # Put your raw data files here
│   ├── methylation_raw.csv
│   ├── sample_metadata.csv
│   └── ...
└── processed/        # Processed data will be saved here
```

**Supported data formats:**
- CSV files with beta values or M-values
- IDAT files from Illumina arrays (450K or EPIC)
- Pre-normalized data matrices

### 3. Run the Analysis Pipeline

#### Option A: Run the Complete Pipeline

Execute the main analysis script:

```r
source("scripts/main_analysis.R")
```

#### Option B: Run Step-by-Step

For more control, run each step individually:

1. **Preprocess Data**
   ```r
   source("scripts/01_preprocess_data.R")
   ```

2. **Exploratory Analysis**
   ```r
   source("scripts/02_exploratory_analysis.R")
   ```

3. **Statistical Analysis**
   ```r
   source("scripts/03_statistical_analysis.R")
   ```

### 4. Generate Reports

Create an analysis report using the R Markdown template:

```r
rmarkdown::render("docs/analysis_report.Rmd")
```

The rendered HTML report will be saved in the `docs/` directory.

## Data Requirements

Your data should include:

1. **Methylation Data**
   - Beta values (0-1) or M-values
   - Rows = CpG sites, Columns = Samples
   - Minimum: Contains clock CpG sites (Horvath, Hannum, etc.)

2. **Sample Metadata**
   - Sample IDs (matching methylation data)
   - Chronological age (years)
   - Sex (Male/Female)
   - Optional: Additional clinical/phenotypic variables

### Example Data Format

**methylation_data.csv:**
```
CpG,Sample1,Sample2,Sample3
cg00000029,0.83,0.81,0.85
cg00000165,0.45,0.47,0.44
...
```

**sample_metadata.csv:**
```
SampleID,Age,Sex,Group
Sample1,45,M,Control
Sample2,52,F,Case
Sample3,38,M,Control
...
```

## Customizing the Analysis

### Modifying Scripts

All analysis scripts are in the `scripts/` directory and can be customized:

- `helper_functions.R`: Add your own functions
- `01_preprocess_data.R`: Adjust QC thresholds and normalization methods
- `02_exploratory_analysis.R`: Add custom visualizations
- `03_statistical_analysis.R`: Modify statistical tests

### Using Different Epigenetic Clocks

To use specific clocks, modify the relevant sections in `main_analysis.R`:

```r
# For Horvath clock
horvath_age <- calculate_horvath_age(methylation_data)

# For Hannum clock
hannum_age <- calculate_hannum_age(methylation_data)

# For PhenoAge
phenoage <- calculate_phenoage(methylation_data, clinical_data)
```

## Troubleshooting

### Common Issues

1. **Package Installation Fails**
   - Ensure you have the latest version of R
   - Try installing BiocManager first: `install.packages("BiocManager")`
   - Check your internet connection

2. **Memory Issues**
   - Increase R memory limit: `memory.limit(size = 16000)` (Windows)
   - Process data in batches
   - Use more efficient data structures (e.g., data.table)

3. **Missing CpG Sites**
   - Some clocks require specific CpG sites
   - Check if your array platform includes the required sites
   - Use imputation methods for missing sites (with caution)

4. **Batch Effects**
   - Check for batch effects using PCA
   - Apply batch correction using ComBat or similar methods
   - Include batch as a covariate in your models

## Getting Help

### Documentation

- **Analysis Guidelines**: See `docs/analysis_guidelines.md`
- **R Markdown Template**: See `docs/analysis_report.Rmd`
- **Package Documentation**: Use `?function_name` in R console

### Resources

- **Bioconductor Support**: https://support.bioconductor.org/
- **minfi User Guide**: https://bioconductor.org/packages/release/bioc/vignettes/minfi/inst/doc/minfi.html
- **ChAMP Tutorial**: https://bioconductor.org/packages/release/bioc/vignettes/ChAMP/inst/doc/ChAMP.html

### Key Papers

1. Horvath, S. (2013). DNA methylation age of human tissues and cell types. *Genome Biology*.
2. Hannum, G. et al. (2013). Genome-wide methylation profiles reveal quantitative views of human aging rates. *Molecular Cell*.
3. Levine, M. E. et al. (2018). An epigenetic biomarker of aging for lifespan and healthspan. *Aging*.

## Next Steps

1. Read the analysis guidelines in `docs/analysis_guidelines.md`
2. Familiarize yourself with the helper functions in `scripts/helper_functions.R`
3. Prepare your data according to the format requirements
4. Start with exploratory data analysis
5. Consult the literature for appropriate analysis methods

## Tips for Success

- **Start Small**: Test your pipeline with a subset of data first
- **Document Everything**: Keep notes on analysis decisions and parameters
- **Version Control**: Commit changes frequently with descriptive messages
- **Validate Results**: Cross-check findings with published literature
- **Reproducibility**: Set random seeds and document software versions

Happy analyzing!
