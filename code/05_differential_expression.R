
library(DESeq2)

# Import count data
counts <- read.csv(
  "C:/Users/shaat/OneDrive/Desktop/DR_Bioinformatics_ICA3/data/GSE160306_human_retina_DR_totalRNA_counts.txt",
  check.names = FALSE
)

count_matrix <- as.matrix(counts[, -1])
rownames(count_matrix) <- counts$ensemblID
storage.mode(count_matrix) <- "integer"

# Import verified sample metadata
sample_data <- read.csv(
  "C:/Users/shaat/OneDrive/Desktop/DR_Bioinformatics_ICA3/data/sample_metadata.csv",
  stringsAsFactors = TRUE
)

# Match metadata order to count matrix columns
rownames(sample_data) <- sample_data$sample
sample_data <- sample_data[colnames(count_matrix), , drop = FALSE]

stopifnot(identical(rownames(sample_data), colnames(count_matrix)))

sample_data$group <- relevel(
  factor(sample_data$group),
  ref = "Control"
)

# Filter low-count genes
keep <- rowSums(count_matrix >= 10) >= 2
count_matrix <- count_matrix[keep, ]

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_data,
  design = ~ group
)

# Run differential expression analysis
dds <- DESeq(dds)

# Example: NPDR compared with Control
# Change NPDR to an actual verified group name
res <- results(
  dds,
  contrast = c("group", "NPDR", "Control"),
  alpha = 0.05
)

res_df <- as.data.frame(res)
res_df$ensemblID <- rownames(res_df)

# Save results
write.csv(
  res_df,
  "C:/Users/shaat/OneDrive/Desktop/DR_Bioinformatics_ICA3/results/tables/NPDR_vs_Control_DEGs.csv",
  row.names = FALSE
)

# Significant genes
significant_genes <- subset(
  res_df,
  !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1
)

write.csv(
  significant_genes,
  "C:/Users/shaat/OneDrive/Desktop/DR_Bioinformatics_ICA3/results/tables/NPDR_vs_Control_significant_genes.csv",
  row.names = FALSE
)

