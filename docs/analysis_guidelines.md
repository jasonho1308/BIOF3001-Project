# Analysis Guidelines

## Data Preprocessing Best Practices

### Quality Control
1. **Sample QC**
   - Remove samples with >5% failed probes
   - Check for sample mix-ups using SNP probes
   - Verify sex prediction matches metadata
   - Identify and handle outliers

2. **Probe QC**
   - Filter probes with detection p-value > 0.01
   - Remove cross-reactive probes
   - Remove probes with SNPs at CpG sites
   - Filter sex chromosome probes (optional)

### Normalization
- **Recommended**: Functional normalization (for EPIC arrays) or BMIQ
- **Alternative**: Quantile normalization, SWAN
- Normalize within batches if batch effects present
- Consider using ComBat for batch effect correction

## Epigenetic Clock Calculations

### Horvath Clock (2013)
- Uses 353 CpG sites
- Applicable to multiple tissues
- Normalized age transformation required
- R package: `methylclock` or custom implementation

### Hannum Clock (2013)
- Uses 71 CpG sites
- Specific to blood samples
- Linear model

### PhenoAge (2018)
- Incorporates clinical biomarkers
- Better predictor of healthspan
- Requires additional phenotypic data

### GrimAge (2019)
- Best predictor of mortality
- Uses DNAm surrogates of plasma proteins
- Requires specific CpG sites

## Statistical Analysis Considerations

### Age Acceleration Calculation
```r
# Simple age acceleration
age_accel = epigenetic_age - chronological_age

# Residual age acceleration (preferred)
model <- lm(epigenetic_age ~ chronological_age)
age_accel_residual <- residuals(model)
```

### Multiple Testing Correction
- Use FDR (Benjamini-Hochberg) for exploratory analyses
- Use Bonferroni for confirmatory analyses
- Report both raw and adjusted p-values

### Covariates to Consider
- Sex
- Cell type composition (blood samples)
- Batch/plate effects
- Technical covariates (array row/column)
- Smoking status
- BMI

## Visualization Guidelines

### Required Plots
1. QC plots (pre and post normalization)
2. Age correlation plot (epigenetic vs chronological)
3. Age acceleration distribution
4. Methylation heatmaps for clock CpGs
5. PCA plots colored by key variables

### Figure Quality
- Use consistent color schemes
- Include error bars/confidence intervals
- Label axes clearly with units
- Use appropriate scales (linear/log)
- Export at high resolution (300 DPI minimum)

## Reproducibility

### Documentation
- Document all analysis steps
- Record software versions
- Note any deviations from standard protocols
- Keep analysis logs

### Code Organization
- Use relative paths (with `here` package)
- Write modular, reusable functions
- Comment complex operations
- Use meaningful variable names

### Version Control
- Commit after each major analysis step
- Write descriptive commit messages
- Don't commit large data files
- Use .gitignore appropriately

## Common Pitfalls to Avoid

1. **Overfitting**: Don't use training data for validation
2. **Batch Effects**: Always check for and correct batch effects
3. **Cell Type Heterogeneity**: Adjust for cell composition in blood
4. **Technical Artifacts**: Remove probes prone to technical issues
5. **Multiple Testing**: Always correct for multiple comparisons

## Resources

### Key Papers
- Horvath (2013) - Original multi-tissue clock
- Hannum et al. (2013) - Blood-based clock
- Levine et al. (2018) - PhenoAge
- Lu et al. (2019) - GrimAge

### Software Packages
- `minfi`: Methylation array analysis
- `ChAMP`: Comprehensive analysis pipeline
- `wateRmelon`: QC and normalization
- `methylclock`: Epigenetic age calculation
- `ENmix`: Preprocessing and QC

### Online Tools
- EWAS Atlas: Database of EWAS results
- MethBase: Methylation database
- GEO: Gene Expression Omnibus (for public datasets)
