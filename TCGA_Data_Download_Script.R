# Downloading raw data from TCGA
# Using GDC Query
library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)
library(dplyr)

# Setting up the query variable
query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

print(query$results[[1]]$cases)

# Downloading in chunks - large file hnadling 
download_tcga_safely <- function(query, max_retries = 5, files_per_chunk = 5) {
  # Set longer timeout
  options(timeout = 600)  # 10 minutes
  
  for (attempt in 1:max_retries) {
    message(paste("\n=== Attempt", attempt, "of", max_retries, "==="))
    
    tryCatch({
      # Try with very small chunks
      GDCdownload(
        query = query,
        method = "api",
        directory = "GDCdata",
        files.per.chunk = files_per_chunk
      )
      message("✓ Download successful!")
      return(TRUE)
      
    }, error = function(e) {
      message(paste("✗ Attempt failed:", e$message))
      
      # Increase chunk size slightly each retry
      files_per_chunk <- min(files_per_chunk * 2, 20)
      message(paste("Next attempt with", files_per_chunk, "files per chunk"))
      
      # Clean partial downloads
      if (attempt < max_retries) {
        message("Cleaning failed downloads...")
        tar_files <- list.files(pattern = "\\.tar\\.gz$")
        if (length(tar_files) > 0) file.remove(tar_files)
        Sys.sleep(30)  # Wait before retry
      }
    })
  }
  return(FALSE)
}

download_query_safely(query)

# Prepaing counts and metadata files 

data <- GDCprepare(query, directory = "GDCdata")
counts <- assay(data, "unstranded")
counts_df <- as.data.frame(counts)
counts_df <- cbind(gene = rownames(counts_df), counts_df)

metadata <- colData(data)
metadata_df <- as.data.frame(metadata)

# Saving the files 
fwrite(metadata_df, "metadata.csv")
fwrite(counts_df, "counts.csv")

#counts_df <- read.delim("counts.csv", header = T, sep = ",")

# Creating separate Tumor and Normal Counts files
# Matching barcodes between counts_df and metadata
metadata$barcode <- gsub("-",".", metadata$barcode)

counts_df <- as.data.frame(t(counts_df))
counts_df$barcode <- rownames(counts_df)
colnames(counts_df) <- counts_df[1, ]
counts_df <- counts_df [-1, ]
# Joining with metadata to get tumor and normal sample separately

counts_df <- counts_df %>% 
  left_join(metadata[, c("barcode", "shortLetterCode")], # shortLetterCode contains tumor and normal info
            by = "barcode") 
rownames(counts_df) <- counts_df$barcode
counts_df <- counts_df %>% select(shortLetterCode, everything())

## Visualizing the number of Tumor and Normal Samples 
counts_df %>% 
  group_by(shortLetterCode) %>% 
  summarise(n = n())

## Splitting the file into raw_tumor and raw_normal
tumor_raw <- counts_df %>% filter(shortLetterCode == "TP")
head(tumor_raw[,1:5])
write.table(tumor_raw, "raw_counts_brca.csv", row.names = T, col.names = T, sep = ",", quote = F)

normal_raw <- counts_df %>% filter(shortLetterCode == "NT")
write.table(normal_raw, "raw_counts_normal.csv", row.names = T, col.names = T, sep = ",", quote = F)
###---- Proceed to book chapter protocol ---###