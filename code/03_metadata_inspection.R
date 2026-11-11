
# Install GEOquery if necessary
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("GEOquery", quietly = TRUE)) {
  BiocManager::install("GEOquery")
}

library(GEOquery)

# Download GEO metadata
gse <- getGEO("GSE160306", GSEMatrix = TRUE)

# Extract metadata
metadata <- pData(gse[[1]])

# Inspect metadata
dim(metadata)
colnames(metadata)

# Display sample titles
metadata$title

# Display sample characteristics
characteristic_columns <- grep(
  "characteristics|title",
  colnames(metadata),
  ignore.case = TRUE,
  value = TRUE
)

metadata[, characteristic_columns, drop = FALSE]

write.csv(
  metadata,
  "C:/Users/shaat/OneDrive/Desktop/DR_Bioinformatics_ICA3/results/tables/GSE160306_metadata.csv",
  row.names = FALSE
)

file.exists("results/tables/GSE160306_metadata.csv")

View(metadata)