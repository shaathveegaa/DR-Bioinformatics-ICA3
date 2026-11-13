
# ============================================================
# BIOINFORMATICS ICA 3
# ANNOTATED GENE EXPRESSION HEATMAP
# Dataset: GSE160306
# ============================================================

# 1. Load required libraries
library(DESeq2)
library(pheatmap)

# 2. Set project working directory
setwd("C:/Users/shaat/OneDrive/Desktop/DR_Bioinformatics_ICA3")

# 3. Create output directory
dir.create(
  "results/plots",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# STEP 1: IMPORT RAW COUNT DATA
# ============================================================
# Import GEO metadata using a new variable name
geo_metadata <- read.csv(
  "results/tables/GSE160306_metadata.csv",
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Check that it is a data frame
class(geo_metadata)

# Display column names
colnames(geo_metadata)

# Display first 6 rows
head(geo_metadata)

counts <- read.csv(
  "data/GSE160306_human_retina_DR_totalRNA_counts.txt",
  check.names = FALSE
)

# Prepare gene expression count matrix
count_matrix <- as.matrix(counts[, -1])

rownames(count_matrix) <- counts$ensemblID

storage.mode(count_matrix) <- "integer"

# Check dimensions
cat("Count matrix dimensions:\n")
print(dim(count_matrix))

# ============================================================
# STEP 2: IMPORT GEO SAMPLE METADATA
# ============================================================

metadata <- read.csv(
  "results/tables/GSE160306_metadata.csv",
  row.names = 1,
  check.names = FALSE
)

# Check required metadata columns
required_columns <- c(
  "sampleID:ch1",
  "sample_group:ch1",
  "sample_site:ch1"
)

if (!all(required_columns %in% colnames(metadata))) {
  stop("Required metadata columns are missing.")
}

# Extract sample information
metadata$Sample_ID <- metadata[["sampleID:ch1"]]

metadata$Group <- metadata[["sample_group:ch1"]]

metadata$Site <- metadata[["sample_site:ch1"]]

# Check missing values
if (any(is.na(metadata$Sample_ID)) ||
    any(is.na(metadata$Group)) ||
    any(is.na(metadata$Site))) {
  stop("Missing sample information in metadata.")
}

# Check duplicate sample IDs
if (anyDuplicated(metadata$Sample_ID)) {
  stop("Duplicate sample IDs found in metadata.")
}

# Set sample IDs as row names
rownames(metadata) <- metadata$Sample_ID

# ============================================================
# STEP 3: MATCH SAMPLES
# ============================================================

missing_samples <- setdiff(
  colnames(count_matrix),
  rownames(metadata)
)

if (length(missing_samples) > 0) {
  cat("Samples missing from metadata:\n")
  print(missing_samples)
  
  stop("Sample matching failed.")
}

# Reorder metadata according to count matrix
sample_data <- metadata[
  colnames(count_matrix),
  ,
  drop = FALSE
]

# Confirm correct matching
stopifnot(
  identical(
    rownames(sample_data),
    colnames(count_matrix)
  )
)

cat("\nSample matching successful!\n")

cat("Total samples:", ncol(count_matrix), "\n")

# Display group information
cat("\nSamples in each group:\n")
print(table(sample_data$Group))

cat("\nSamples at each retinal site:\n")
print(table(sample_data$Site))

# ============================================================
# STEP 4: CREATE DESEQ2 DATASET
# ============================================================

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_data,
  design = ~ 1
)

# ============================================================
# STEP 5: VARIANCE STABILIZING TRANSFORMATION
# ============================================================

vsd <- vst(
  dds,
  blind = TRUE
)

# Extract transformed expression values
vst_matrix <- assay(vsd)

cat("\nVST transformation completed!\n")

# ============================================================
# STEP 6: IDENTIFY TOP 30 VARIABLE GENES
# ============================================================

# Calculate variance for each gene
gene_variance <- apply(
  vst_matrix,
  1,
  var
)

# Sort genes by variance
sorted_variance <- sort(
  gene_variance,
  decreasing = TRUE
)

# Select top 30 genes
top_30_genes <- names(
  sorted_variance
)[1:30]

# Extract expression values
heatmap_data <- vst_matrix[
  top_30_genes,
  ,
  drop = FALSE
]

cat("\nTop 30 variable genes selected!\n")

# Save selected gene information
gene_table <- data.frame(
  Ensembl_ID = top_30_genes,
  Variance = gene_variance[top_30_genes]
)

write.csv(
  gene_table,
  "results/tables/top_30_variable_genes.csv",
  row.names = FALSE
)

# ============================================================
# STEP 7: CREATE SAMPLE ANNOTATIONS
# ============================================================

annotation_col <- data.frame(
  Disease_Group = factor(sample_data$Group),
  Retinal_Site = factor(sample_data$Site)
)

rownames(annotation_col) <- rownames(sample_data)

# ============================================================
# STEP 8: GENERATE ANNOTATED HEATMAP
# ============================================================

pheatmap(
  heatmap_data,
  
  # Standardize expression values by gene
  scale = "row",
  
  # Hierarchical clustering
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  
  # Sample annotations
  annotation_col = annotation_col,
  
  # Display labels
  show_rownames = TRUE,
  show_colnames = TRUE,
  
  # Heatmap title
  main = "Top 30 Most Variable Genes in Human Retina",
  
  # Font sizes
  fontsize = 10,
  fontsize_row = 8,
  fontsize_col = 7,
  
  # Sample label angle
  angle_col = 45,
  
  # Legends
  legend = TRUE,
  annotation_legend = TRUE,
  annotation_names_col = TRUE,
  
  # Appearance
  border_color = NA,
  cellwidth = 12,
  cellheight = 12,
  
  # Save heatmap
  filename = "results/plots/01_annotated_heatmap.png",
  
  width = 18,
  height = 12
)

# ============================================================
# STEP 9: SAVE TRANSFORMED EXPRESSION MATRIX
# ============================================================

write.csv(
  heatmap_data,
  "results/tables/top_30_gene_expression.csv"
)

# ============================================================
# STEP 10: COMPLETION MESSAGE
# ============================================================

cat("\n====================================\n")
cat("HEATMAP GENERATED SUCCESSFULLY!\n")
cat("====================================\n")

cat(
  "\nHeatmap: results/plots/01_annotated_heatmap.png\n"
)

cat(
  "\nGene table: results/tables/top_30_variable_genes.csv\n"
)

cat(
  "\nExpression table: results/tables/top_30_gene_expression.csv\n"
)

# Verify output files
cat("\nOutput file verification:\n")

print(
  file.exists(
    "results/plots/01_annotated_heatmap.png"
  )
)

print(
  file.exists(
    "results/tables/top_30_variable_genes.csv"
  )
)

print(
  file.exists(
    "results/tables/top_30_gene_expression.csv"
  )
)


