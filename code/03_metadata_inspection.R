
# ============================================================
# 03_metadata_inspection.R
# Prepare verified sample metadata
# Dataset: GSE160306
# ============================================================

# 1. Set working directory
setwd("C:/Users/shaat/OneDrive/Desktop/DR_Bioinformatics_ICA3")


# 2. Load count data
counts <- read.csv(
  "data/GSE160306_human_retina_DR_totalRNA_counts.txt",
  check.names = FALSE
)

count_matrix <- as.matrix(counts[, -1])
colnames(count_matrix) <- colnames(counts)[-1]


# 3. Load original GEO metadata
geo_metadata <- read.csv(
  "C:/Users/shaat/OneDrive/Desktop/DR_Bioinformatics_ICA3/data/sample_metadata_original.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# 4. Check required columns
if (!all(c("title", "disease_group:ch1") %in%
         colnames(geo_metadata))) {
  stop("Required metadata columns are missing.")
}


# 5. Match sample names
sample_names <- colnames(count_matrix)

metadata_index <- match(
  sample_names,
  geo_metadata$title
)

if (any(is.na(metadata_index))) {
  cat("Unmatched samples:\n")
  print(sample_names[is.na(metadata_index)])
  stop("Sample matching failed.")
}


# 6. Reorder metadata
matched_metadata <- geo_metadata[
  metadata_index,
  ,
  drop = FALSE
]


# 7. Create the required two-column metadata
sample_metadata <- data.frame(
  sample = sample_names,
  group = matched_metadata$`disease_group:ch1`
)


# 8. Verify the metadata
if (any(is.na(sample_metadata$group)) ||
    any(trimws(sample_metadata$group) == "")) {
  stop("Missing disease group information.")
}

if (anyDuplicated(sample_metadata$sample)) {
  stop("Duplicate sample names found.")
}

stopifnot(
  identical(sample_metadata$sample, colnames(count_matrix))
)


# 9. Display group distribution
cat("\nDisease group distribution:\n")
print(table(sample_metadata$group))


# 10. Save the final metadata file
write.csv(
  sample_metadata,
  "data/sample_metadata.csv",
  row.names = FALSE
)


# 11. Final verification
cat("\nMetadata verification successful!\n")
cat("Total samples:", nrow(sample_metadata), "\n")
cat("Total columns:", ncol(sample_metadata), "\n")

print(head(sample_metadata))
