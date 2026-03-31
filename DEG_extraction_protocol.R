# Loading Libraries 
library(dplyr)
library(DESeq2)
library(caret)
library(ggplot2)
library(stringr)

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
# Threshold out the bottom 50% 
threshold <- quantile(counts_filtered, 0.5, na.rm = T)
# Only keep the genes that are above the threshold
counts_filtered <- counts_filtered[gene_variance > threshold, ]
# Remove NA rownames 
counts_filtered <- counts_filtered[!is.na(rownames(counts_filtered)), ]
dim(counts_filtered)

# Stratified Train-Test Split(70:30)
ml_data <- as.data.frame(t(counts_filtered))
counts_labels <- as.character(unlist(counts_labels))
ml_data$type <- factor(counts_labels)

ml_data <- ml_data %>%
  mutate(type = counts_labels)
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

top_genes <- rownames(res_sig)

# Download the lncRNA file from gencode v32
# Subsetting the DEGs for differentially expressed lncRNA only 
gencode.v32.long_noncoding_RNAs <- read.delim("gencode.v32.long_noncoding_RNAs.gtf", header=FALSE, comment.char="#")

lncrna_genes <- gencode.v32.long_noncoding_RNAs %>% filter(V3 == "gene")
lncrna_genes <- gencode.v32.long_noncoding_RNAs %>%
  filter(V3 == "gene") %>%
  mutate(gene_id = str_remove(str_extract(V9, 'ENSG[^"]+'), '\\..*'))
lncrna_gene_ids <- as.data.frame(unique(lncrna_genes$gene_id))

write.csv(lncrna_gene_ids, "lncRNA_Gene_IDs.csv")


lncrna_degs <- intersect(lncrna_gene_ids$`unique(lncrna_genes$gene_id)`, top_genes)

counts_lncrna_degs <- counts[lncrna_degs, ]
counts_lncrna_degs <- as.data.frame(t(counts_lncrna_degs))
counts_lncrna_degs$type <- factor(counts_labels)

counts_lncrna_degs <- counts_lncrna_degs %>%
  mutate(type = counts_labels)
counts_lncrna_degs <- counts_lncrna_degs %>% dplyr::select(type, everything())

write.csv(counts_lncrna_degs, "lncRNA_DEGs_ml_input.csv", col.names = T, row.names = T)
### --- Proceed to Book chapter protocol --- ###
