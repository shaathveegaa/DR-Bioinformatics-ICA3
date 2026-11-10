
# 1. Load dataset
file_path <- "data/GSE160306_human_retina_DR_totalRNA_counts.txt"

counts <- read.csv(
  file_path,
  header = TRUE,
  check.names = FALSE
)

# 2. Dataset dimensions
cat("Dataset dimensions:\n")
print(dim(counts))

# 3. Column names
cat("Sample column names:\n")
print(colnames(counts))

# 4. First few rows
head(counts)

# 5. Data structure
str(counts)

# 6. Missing values
cat("Missing values:\n")
print(sum(is.na(counts)))

# 7. Duplicate gene IDs
cat("Duplicate gene IDs:\n")
print(sum(duplicated(counts$ensemblID)))

# 8. Summary statistics
summary(counts)

# 9. Check count values
count_matrix <- counts[, -1]

cat("Minimum count:\n")
print(min(as.matrix(count_matrix), na.rm = TRUE))

cat("Maximum count:\n")
print(max(as.matrix(count_matrix), na.rm = TRUE))


# Create results/tables folder
dir.create(
  "results/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

# Create dataset exploration summary
exploration_summary <- data.frame(
  Total_Genes = nrow(counts),
  Total_Columns = ncol(counts),
  Total_Samples = ncol(counts) - 1,
  Missing_Values = sum(is.na(counts)),
  Duplicate_Gene_IDs = sum(duplicated(counts$ensemblID)),
  Minimum_Count = min(as.matrix(count_matrix)),
  Maximum_Count = max(as.matrix(count_matrix))
)

# Display summary
print(exploration_summary)

# Save summary as CSV
write.csv(
  exploration_summary,
  "results/tables/dataset_exploration_summary.csv",
  row.names = FALSE
)

# Verify file creation
file.exists("results/tables/dataset_exploration_summary.csv")

