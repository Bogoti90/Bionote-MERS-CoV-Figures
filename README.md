# MERS-CoV Bionote Rapid Antigen Test Validation Analysis

[![R](https://img.shields.io/badge/R-%3E%3D4.0-blue.svg)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

## Overview

This repository contains the complete statistical analysis and figure generation code for the manuscript:

**"On-site Detection of MERS-CoV Infections in a Camel Slaughterhouse in Kenya Using a Commercial Rapid Antigen Test"**

The code generates two main publication figures:
- **Figure 1**: Clade comparison and validation analysis
- **Figure 2**: Four-panel validation analysis (ROC, LOD, sensitivity)

## Quick Start

### Prerequisites
```r
# Required R packages
install.packages(c("tidyverse", "ggplot2", "cowplot", "scales", "pROC"))
```

### Usage
1. Clone this repository
2. Place your data files in the working directory:
   - `clade_bionote_analysis.csv`
   - `Bionote_results.csv`
3. Run the analysis:
```r
source("Final_Bionotefigures_Code.R")
```


## Generated Figures

### Figure 1: Clade Comparison and Validation Analysis
- **Panel A**: MERS-CoV clade performance comparison (RNA copies/mL and TCID₅₀/mL)
- **Panel B**: Field validation with viral isolation correlation

### Figure 2: Four-Panel Validation Analysis  
- **Panel A**: Sensitivity vs viral load threshold
- **Panel B**: ROC curve analysis (AUC calculation)
- **Panel C**: Limit of detection determination
- **Panel D**: Detection probability modeling


## Statistical Methods

- **Non-parametric testing**: Mann-Whitney U tests for group comparisons
- **ROC analysis**: pROC package for performance evaluation  
- **Confidence intervals**: Wilson score method for proportions
- **LOD calculation**: Binomial exact tests with Youden's J statistic

## Data Requirements

### clade_bionote_analysis.csv
Required columns:
- `Sample.ID`: Unique sample identifier
- `Virus.dose..expected.TCID50.mL.`: Virus dose information
- `Rapid.Ag.Test.result..15.min.`: Rapid test results
- `TCID50.mL..VeroE6.T.titration.`: TCID₅₀ measurements
- `genome.copies.mL`: RNA copy measurements

### Bionote_results.csv  
Required columns:
- `Animal_ID`: Animal identifier
- `Final_result`: PCR results (Positive/Negative)
- `MERS_RNA_copiesPerml`: Viral load measurements
- `Bionote_result`: Rapid test results (Positive/Negative)
- `Viral.Isolation`: Virus isolation results

## Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| tidyverse | ≥1.3.0 | Data manipulation and visualization |
| ggplot2 | ≥3.3.0 | Statistical graphics |
| cowplot | ≥1.1.0 | Figure composition |
| scales | ≥1.1.0 | Scale formatting |
| pROC | ≥1.17.0 | ROC analysis |


## Author

**Brian Ogoti**  
Department of Medical Microbiology and Immunology  
University of Nairobi, Kenya  
📧 brian.ogoti@cema.africa

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---
