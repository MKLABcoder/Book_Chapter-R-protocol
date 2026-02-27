# Book_Chapter-R-protocol
An R-based Machine Learning Protocol to Predict Diagnostic Biomarkers in Breast Cancer using RNA-sequencing Data

On this page, we provide the processed RNA-seq raw counts datasets as well as the R-based code that can be used to download the TCGA-BRCA data.
We also provide a DEG extraction protocol as well as the RMarkdown pdf for machine learning section. 

# Requisites for data downloading
*TCGAbiolinks*, *SummarizedExperiment*, *tidyverse*, *dplyr*

# Requisites for DEG extraction
*DESeq2*, *dplyr*, *caret*

# Processed datasets
raw_counts_degs_all_samples.csv - contains the top 200 DEGs raw count data for all 1224 TCGA-BRCA samples

# Scripts

TCGA_data_download.R - Contains a step-by-step script on how TCGA-BRCA data can be downloaded from R 

# PDFs

Final_rmark.pdf - RMarkdown pdf protocol
