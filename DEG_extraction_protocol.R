# Loading Libraries 
library(dplyr)
library(DESeq2)
library(caret)
library(ggplot2)


# Loading Raw Counts for Tumor and Normal
tumor_raw  <- read.csv("raw_counts_brca.csv", check.names = FALSE)
normal_raw <- read.csv("raw_counts_normal.csv", check.names = FALSE)

head(tumor_raw[, 1:5])

# Loaded Data should be a genes x samples dataframe
tumor_raw <- as.data.frame(t(tumor_raw))
normal_raw <- as.data.frame(t(normal_raw))

# Binding into a single dataframe 
counts <- cbind(tumor_raw, normal_raw)
dim(counts)

# Storing the labels separately
counts_labels <- counts[1, ]
counts <- counts[-1, ]

# Filtering against Low read counts (Rows with less than 10 counts total are removed)
keep <- rowSums(counts >= 10) >= 10
counts_filtered <- counts[keep, ]
dim(counts_filtered)

# Removing Low Variance Genes
counts_filtered <- as.data.frame(
  apply(counts_filtered, 2, function(x) as.numeric(trimws(x))),
  row.names = rownames(counts_filtered)
)
gene_variance <- apply(counts_filtered, 1, var)
# Threshold out the bottom 25% 
threshold <- quantile(counts_filtered, 0.5, na.rm = T)
# Only keep the genes that are above the threshold
counts_filtered <- counts_filtered[gene_variance > threshold, ]
# Remove NA rownames 
counts_filtered <- counts_filtered[!is.na(rownames(counts_filtered)), ]
dim(counts_filtered)


# Stratified Train-Test Split(70:30)
ml_data <- as.data.frame(t(counts_filtered))
ml_data$type <- as.factor(as.vector(unlist(counts_labels)))
ml_data <- ml_data %>% dplyr::select(type, everything())

set.seed(345)
# Creating the index to split
trainIndex <- createDataPartition(
  ml_data$type,
  p = 0.7,
  list = FALSE
)

trainData <- ml_data[trainIndex, ]
testData  <- ml_data[-trainIndex, ]

table(trainData$type)
table(testData$type)

# Performing DESeq2 on the training data 
# Preparing Training data for DESeq2
trainData <- as.data.frame(t(trainData)) # Getting a genes x samples 
train_group <- as.vector(unlist(trainData[1, ]))# Storing the counts 
train_counts <- trainData[-1, ] # Storing the labels separately (TP or NT)
train_counts <- as.data.frame(
  apply(train_counts, 2, function(x) as.numeric(trimws(x))),
  row.names = rownames(train_counts)
)
train_counts[is.na(train_counts)] <- 0 # Converting NA counts to 0
train_counts <- as.matrix(train_counts) # DESeq2 requires a matrix as input for counts

col_data <- data.frame(
  row.names = colnames(train_counts),
  condition = factor(train_group)
)
# Creating a DESeq2 object
dds_train <- DESeqDataSetFromMatrix(
  countData = train_counts,
  colData = col_data,
  design = ~ condition
)
# Running DESeq2
dds_train <- DESeq(dds_train)

# Obtaining DEseq2 results
res_train <- results(dds_train)
res_df_train <- as.data.frame(res_train)

# Significance Thresholding 
res_sig <- res_df_train %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  arrange(padj)
cat("Number of significant DEGs:", nrow(res_sig), "\n")

write.csv(res_sig, "Significant_DEGs_all.csv", row.names = T)

# Selecting the top 200 DEGs
top_genes <- rownames(res_sig)[1:200]

# Getting the DEGs raw counts from the main counts_filtered file
# And the sample information

counts <- cbind(tumor_raw, normal_raw)
dim(counts)
counts <- as.data.frame(t(counts))


counts_labels <- counts$shortLetterCode
counts_degs <- counts[, top_genes]
counts_degs$type <- as.factor(as.vector(unlist(counts_labels)))
counts_degs <- counts_degs %>% dplyr::select(type, everything())

write.csv(counts_degs, "raw_counts_degs_all_samples.csv", row.names = T, col.names = T)
