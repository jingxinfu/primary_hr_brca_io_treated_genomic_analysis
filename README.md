# Genomic analysis scripts for "Efficacy, safety, and comprehensive biomarker investigation in patients with early-stage hormone receptor-positive/HER2-negative breast cancer treated with neoadjuvant nab-paclitaxel and pembrolizumab"

## Analysis environment set up
0. Download `mamba` following instruction: [mamba](https://github.com/mamba-org/mamba) 
1. clone this repo:`git clone `
2. switch to the repo directory: `cd primary_hr_brca_io_treated_genomic_analysis`
3. Enter command: `make create_environment`, this will create a conda environment called `hr_brca_16_466`
   
## Preprocess raw data from workflows
Raw data will be shared upon request. The raw data should be stored under `data/raw` folder.
1. preprocess data from bulk WES profiles: `make wes_data`
   - the process detail can be found at `src/dataset/wesDataset.py`
2.  preprocess data from bulk RNA-Seq profiles: `make rna_data`
    - the process detail can be found at `src/dataset/rnaDataset.py`

## Analyze steps
After having preprocessed the bulk RNA-Seq and WES data, we could use following commands to reproduce the analysis been done in our study:
1. perform somatic mutation analysis: `make analyze_mutation`
2. perform differential gene expression and pathway enrichment analysis: `make analyze_rna`
3. perform cytokine and immune infiltration analysis: `make analyze_feature`

The related scripts can be found under `src/models` folder.

## Visualization steps
All figures presented on the manuscript can be reproduced by running the notebooks under `notebooks/` folder.