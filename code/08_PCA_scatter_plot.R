
# ============================================================
# BIOINFORMATICS ICA 3
# 08. PCA SCATTER PLOT
# Dataset: GSE160306
# ============================================================

# 1. Load libraries
library(DESeq2)
library(ggplot2)
library(ggrepel)

# 2. Set working directory
setwd("C:/Users/shaat/OneDrive/Desktop/DR_Bioinformatics_ICA3")

# 3. Create output directory
dir.create(
  "results/plots",
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# STEP 1: IMPORT COUNT DATA
# ============================================================

counts <- read.csv(
  "data/GSE160306_human_retina_DR_totalRNA_counts.txt",
  check.names = FALSE
)

count_matrix <- as.matrix(counts[, -1])

rownames(count_matrix) <- counts$ensemblID

storage.mode(count_matrix) <- "integer"

# ============================================================
# STEP 2: IMPORT GEO METADATA
# ============================================================

geo_metadata <- read.csv(
  "results/tables/GSE160306_metadata.csv",
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Extract metadata columns
geo_metadata$Sample_ID <-
  geo_metadata[["sampleID:ch1"]]

geo_metadata$Group <-
  geo_metadata[["sample_group:ch1"]]

geo_metadata$Site <-
  geo_metadata[["sample_site:ch1"]]

# Check missing values
if (any(is.na(geo_metadata$Sample_ID)) ||
    any(is.na(geo_metadata$Group))) {
  stop("Missing sample information in metadata.")
}

# Check duplicate sample IDs
if (anyDuplicated(geo_metadata$Sample_ID)) {
  stop("Duplicate sample IDs found.")
}

rownames(geo_metadata) <- geo_metadata$Sample_ID

# ============================================================
# STEP 3: MATCH SAMPLES
# ============================================================

missing_samples <- setdiff(
  colnames(count_matrix),
  rownames(geo_metadata)
)

if (length(missing_samples) > 0) {
  print(missing_samples)
  stop("Sample matching failed.")
}

sample_data <- geo_metadata[
  colnames(count_matrix),
  ,
  drop = FALSE
]

stopifnot(
  identical(rownames(sample_data), colnames(count_matrix))
)

cat("Sample matching successful!\n")

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

# ============================================================
# STEP 6: PERFORM PCA
# ============================================================

pca_data <- plotPCA(
  vsd,
  intgroup = c("Group", "Site"),
  returnData = TRUE
)

# Calculate percentage of variance
percent_variance <- round(
  100 * attr(pca_data, "percentVar"),
  2
)

# Add sample labels
pca_data$Sample <- rownames(pca_data)

# ============================================================
# STEP 7: CREATE PCA SCATTER PLOT
# ============================================================

pca_plot <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    color = Group,
    shape = Site
  )
) +
  
  geom_point(
    size = 3.5,
    alpha = 0.85
  ) +
  
  geom_text_repel(
    aes(label = Sample),
    size = 2.5,
    max.overlaps = 25
  ) +
  
  labs(
    title = "PCA of Retinal Gene Expression",
    subtitle = "Principal Component Analysis of GSE160306 Samples",
    x = paste0(
      "PC1 (",
      percent_variance[1],
      "% variance)"
    ),
    y = paste0(
      "PC2 (",
      percent_variance[2],
      "% variance)"
    ),
    color = "Disease Group",
    shape = "Retinal Site"
  ) +
  
  theme_minimal(base_size = 13) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15
    ),
    axis.title = element_text(face = "bold"),
    legend.position = "right"
  )

# ============================================================
# STEP 8: SAVE PCA PLOT
# ============================================================

ggsave(
  "results/plots/02_PCA_scatter_plot.png",
  plot = pca_plot,
  width = 12,
  height = 8,
  dpi = 300
)

# Display plot
print(pca_plot)

cat("\nPCA scatter plot generated successfully!\n")
