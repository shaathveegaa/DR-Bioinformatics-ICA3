
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)

# Load differential expression results
res <- read.csv(
  "C:/Users/shaat/OneDrive/Desktop/DR_Bioinformatics_ICA3/results/tables/NPDR_vs_Control_DEGs.csv"
)

# Volcano plot
res$Significance <- "Not significant"

res$Significance[
  !is.na(res$padj) &
    res$padj < 0.05 &
    res$log2FoldChange > 1
] <- "Upregulated"

res$Significance[
  !is.na(res$padj) &
    res$padj < 0.05 &
    res$log2FoldChange < -1
] <- "Downregulated"

volcano <- ggplot(
  res,
  aes(x = log2FoldChange, y = -log10(padj),
      color = Significance)
) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: NPDR vs Control",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  )

ggsave(
  "C:/Users/shaat/OneDrive/Desktop/DR_Bioinformatics_ICA3/results/plots/volcano_NPDR_Control.png",
  volcano,
  width = 8,
  height = 6
)

# Convert significant Ensembl IDs to Entrez IDs
sig_ids <- res$ensemblID[
  !is.na(res$padj) &
    res$padj < 0.05 &
    abs(res$log2FoldChange) > 1
]

converted <- bitr(
  sig_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# GO Biological Process enrichment
go_results <- enrichGO(
  gene = unique(converted$ENTREZID),
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

write.csv(
  as.data.frame(go_results),
  "C:/Users/shaat/OneDrive/Desktop/DR_Bioinformatics_ICA3/results/tables/GO_enrichment.csv",
  row.names = FALSE
)

# KEGG enrichment
kegg_results <- enrichKEGG(
  gene = unique(converted$ENTREZID),
  organism = "hsa",
  pvalueCutoff = 0.05
)

write.csv(
  as.data.frame(kegg_results),
  "C:/Users/shaat/OneDrive/Desktop/DR_Bioinformatics_ICA3/results/tables/KEGG_enrichment.csv",
  row.names = FALSE
)
