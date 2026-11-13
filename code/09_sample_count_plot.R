
# ============================================================
# BIOINFORMATICS ICA 3
# 09. SAMPLE COUNT PLOT
# Dataset: GSE160306
# ============================================================

# 1. Load library
library(ggplot2)

# 2. Set working directory
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
# STEP 1: IMPORT GEO METADATA
# ============================================================

geo_metadata <- read.csv(
  "results/tables/GSE160306_metadata.csv",
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Extract sample information
geo_metadata$Sample_ID <-
  geo_metadata[["sampleID:ch1"]]

geo_metadata$Group <-
  geo_metadata[["sample_group:ch1"]]

geo_metadata$Site <-
  geo_metadata[["sample_site:ch1"]]

# Check missing values
if (any(is.na(geo_metadata$Sample_ID)) ||
    any(is.na(geo_metadata$Group)) ||
    any(is.na(geo_metadata$Site))) {
  stop("Missing sample information in metadata.")
}

# Check duplicate samples
if (anyDuplicated(geo_metadata$Sample_ID)) {
  stop("Duplicate sample IDs found.")
}

# ============================================================
# STEP 2: COUNT SAMPLES IN EACH DISEASE GROUP
# ============================================================

group_counts <- as.data.frame(
  table(geo_metadata$Group)
)

colnames(group_counts) <- c(
  "Disease_Group",
  "Number_of_Samples"
)

# Remove empty groups
group_counts <- group_counts[
  group_counts$Number_of_Samples > 0,
]

# Display count table
print(group_counts)

# Save count table
write.csv(
  group_counts,
  "results/tables/sample_group_counts.csv",
  row.names = FALSE
)

# ============================================================
# STEP 3: CREATE SAMPLE COUNT BAR PLOT
# ============================================================

count_plot <- ggplot(
  group_counts,
  aes(
    x = Disease_Group,
    y = Number_of_Samples,
    fill = Disease_Group
  )
) +
  
  geom_col(
    width = 0.65
  ) +
  
  geom_text(
    aes(label = Number_of_Samples),
    vjust = -0.5,
    size = 5,
    fontface = "bold"
  ) +
  
  labs(
    title = "Distribution of Retinal Samples Across Groups",
    subtitle = "Sample counts based on GEO metadata",
    x = "Disease Group",
    y = "Number of Samples",
    fill = "Disease Group"
  ) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.15))
  ) +
  
  theme_minimal(base_size = 13) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15
    ),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    ),
    legend.position = "none"
  )

# ============================================================
# STEP 4: SAVE PLOT
# ============================================================

ggsave(
  "results/plots/03_sample_count_plot.png",
  plot = count_plot,
  width = 11,
  height = 7,
  dpi = 300
)

# Display plot
print(count_plot)

cat("\nSample count plot generated successfully!\n")
