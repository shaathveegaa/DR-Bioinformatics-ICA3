
library(DESeq2)
library(ggplot2)

# Import count data
counts <- read.csv(
  "data/GSE160306_human_retina_DR_totalRNA_counts.txt",
  check.names = FALSE
)

# Extract gene IDs
gene_ids <- counts$ensemblID

# Create count matrix
count_matrix <- as.matrix(counts[, -1])
rownames(count_matrix) <- gene_ids

# Ensure integer counts
storage.mode(count_matrix) <- "integer"

# Remove low-expression genes
keep <- rowSums(count_matrix >= 10) >= 2
filtered_counts <- count_matrix[keep, ]

# Inspect library sizes
library_sizes <- colSums(filtered_counts)

print(library_sizes)

# Save filtered counts
write.csv(
  filtered_counts,
  "C:/Users/shaat/OneDrive/Desktop/DR_Bioinformatics_ICA3/results/tables/filtered_counts.csv"
)

# Plot library sizes
library_data <- data.frame(
  Sample = names(library_sizes),
  LibrarySize = as.numeric(library_sizes)
)

p <- ggplot(library_data, aes(Sample, LibrarySize)) +
  geom_col() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

ggsave(
  "C:/Users/shaat/OneDrive/Desktop/DR_Bioinformatics_ICA3/results/plots/library_sizes.png",
  p,
  width = 12,
  height = 6
)
