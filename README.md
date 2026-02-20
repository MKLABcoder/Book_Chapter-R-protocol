# Book_Chapter-R-protocol
An R-based Machine Learning Protocol to Predict Diagnostic Biomarkers in Breast Cancer using RNA-sequencing Data

On this page, we provide the processed RNA-seq raw counts datasets as well as the R-based code that can be used to download the TCGA-BRCA data 

# Requisites for data downloading
*TCGAbiolinks*, *SummarizedExperiment*, *tidyverse*, *dplyr*

# Processed datasets

1. **raw_counts_brca.csv** - Raw TCGA processed counts for breast cancer samples
2. **raw_counts_normal.csv** - Raw TCGA processed counts for normal samples
3. **GSE71651.csv** - Independent processed dataset taken from BARRA:CuRDa database

# Scripts

TCGA_data_download.R - Contains a step by step script on how TCGA-BRCA data can be downloaded from R 
