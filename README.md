# BIOF3001-Project
Investigating the molecular basis of accelerated epigenetic aging

## Overview

This repository contains R scripts and analysis workflows for investigating the molecular mechanisms underlying accelerated epigenetic aging. The project focuses on DNA methylation-based age estimators (epigenetic clocks) and factors that contribute to age acceleration.

## Project Structure

```
BIOF3001-Project/
├── data/
│   ├── raw/              # Raw data files (not tracked in git)
│   ├── processed/        # Processed data files
│   └── README.md
├── scripts/
│   ├── main_analysis.R              # Main analysis pipeline
│   ├── helper_functions.R           # Helper functions
│   ├── 01_preprocess_data.R         # Data preprocessing
│   ├── 02_exploratory_analysis.R    # Exploratory data analysis
│   └── 03_statistical_analysis.R    # Statistical analysis
├── results/              # Analysis results (not tracked in git)
├── figures/              # Generated figures (not tracked in git)
├── docs/                 # Documentation
├── install_packages.R    # Package installation script
├── requirements.R        # Required packages list
└── README.md
```

## Getting Started

### Prerequisites

- R (version 4.0 or higher recommended)
- RStudio (optional but recommended)

### Installation

1. Clone this repository:
```bash
git clone https://github.com/jasonho1308/BIOF3001-Project.git
cd BIOF3001-Project
```

2. Install required R packages:
```r
source("install_packages.R")
```

### Usage

1. **Data Preparation**: Place your raw methylation data in `data/raw/`

2. **Preprocessing**: Run the preprocessing script
```r
source("scripts/01_preprocess_data.R")
```

3. **Exploratory Analysis**: Explore data patterns
```r
source("scripts/02_exploratory_analysis.R")
```

4. **Statistical Analysis**: Perform statistical tests
```r
source("scripts/03_statistical_analysis.R")
```

5. **Main Analysis**: Run the complete pipeline
```r
source("scripts/main_analysis.R")
```

## Research Questions

1. What are the molecular markers associated with accelerated epigenetic aging?
2. How does epigenetic age differ from chronological age in various populations?
3. What factors contribute to age acceleration (genetics, lifestyle, environment)?
4. Can we identify biomarkers that predict accelerated aging?

## Epigenetic Clocks

This project utilizes several DNA methylation-based age estimators:
- **Horvath Clock**: Multi-tissue age predictor (353 CpG sites)
- **Hannum Clock**: Blood-based age predictor (71 CpG sites)
- **PhenoAge**: Predicts mortality and healthspan
- **GrimAge**: Predicts time to death

## Data Requirements

The analysis expects DNA methylation data from platforms such as:
- Illumina HumanMethylation450 BeadChip
- Illumina MethylationEPIC BeadChip

Data should include:
- Beta values or M-values for CpG sites
- Sample metadata (age, sex, clinical variables)
- Quality control metrics

## Key Analyses

1. **Quality Control**: Sample and probe filtering
2. **Normalization**: Quantile or functional normalization
3. **Age Prediction**: Calculate epigenetic age using various clocks
4. **Age Acceleration**: Compute difference from chronological age
5. **Association Studies**: Identify factors related to age acceleration
6. **Visualization**: Generate plots and heatmaps

## Contributing

This is a research project for BIOF3001. Contributions and suggestions are welcome.

## References

- Horvath, S. (2013). DNA methylation age of human tissues and cell types. *Genome Biology*, 14(10), R115.
- Hannum, G. et al. (2013). Genome-wide methylation profiles reveal quantitative views of human aging rates. *Molecular Cell*, 49(2), 359-367.
- Levine, M. E. et al. (2018). An epigenetic biomarker of aging for lifespan and healthspan. *Aging*, 10(4), 573-591.

## License

This project is for educational and research purposes.

## Contact

For questions or collaboration, please open an issue on GitHub.
