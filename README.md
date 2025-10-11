# CSB195
Code and data for CSB195 Computational Biology Foundations, University of Toronto
Test commit Mon Oct  6 12:49:30 EDT 2025
GCM OK Mon Oct  6 12:57:31 EDT 2025
Another test Mon Oct  6 13:15:12 EDT 2025

# CSB195 – Report 1

This repository contains the source files for **Report 1** in CSB195 Computational Biology Foundations, University of Toronto.

## Contents
- `report1.qmd` – Main Quarto document for the report.  
- `dat/aaSim.4.1.Rds` – Data file used for similarity scoring.  
- `src/R/check_benchmark.R` – Script to verify the benchmark score.  
- `src/R/random_codes.R` – Script to generate random genetic codes under different null models.  
- `out/` – Folder with simulation results (CSV tables and figures).  

## How to Reproduce
1. Clone this repository.  
2. Open `report1.qmd` in RStudio or run:  
   ```bash
   quarto render report1.qmd
