# Book_Chapter-R-protocol
An R-based Machine Learning Protocol to Predict Diagnostic Biomarkers in Breast Cancer using RNA-sequencing Data

On this page, we provide the processed RNA-seq raw counts datasets as well as the R-based code that can be used to download the TCGA-BRCA data.
We also provide a DEG extraction protocol as well as the RMarkdown pdf for machine learning section. 

# Requisites for data downloading
*TCGAbiolinks*, *SummarizedExperiment*, *tidyverse*, *dplyr*

# Requisites for DEG extraction
*DESeq2*, *dplyr*, *caret*, *stringr*

# Processed datasets
This file needs to be downloaded to follow the protocol
lncRNA_DEGs_ml_input.csv - Contains the significant differentially expressed lncRNA raw counts 

# Scripts

TCGA_data_download.R - Contains a step-by-step script on how TCGA-BRCA data can be downloaded from R 
DEG_extraction_protocol.R - Contains the step-by-step script on extracting DEGs and making the processed file provided

# PDFs

Final_rmark_cp_1.pdf - RMarkdown pdf protocol
