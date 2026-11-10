
# Comparative Transcriptomic Analysis of Diabetic Retinopathy

## 1. Project Overview

This project investigates gene expression changes associated with different stages of Diabetic Retinopathy (DR) using publicly available transcriptomic data.

The study aims to identify differentially expressed genes between healthy control samples and samples representing different stages of diabetic retinopathy.

## 2. Research Aim

To perform a comparative transcriptomic analysis of gene expression changes across different stages of diabetic retinopathy using the GSE160306 dataset.

## 3. Research Questions

1. Which genes show significant expression differences between healthy controls and diabetic retinopathy samples?
2. How do gene expression patterns differ across the available stages of diabetic retinopathy?
3. Which biological processes and pathways are associated with the identified differentially expressed genes?

## 4. Dataset Information

- **Dataset:** GSE160306
- **Source:** NCBI Gene Expression Omnibus (GEO)
- **Organism:** Homo sapiens
- **Data Type:** RNA sequencing gene expression data
- **Data Format:** Gene count matrix

Dataset link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160306

## 5. Tools and Technologies

- R
- RStudio
- GEOquery
- DESeq2
- ggplot2
- pheatmap
- clusterProfiler
- org.Hs.eg.db
- GitHub

## 6. Methodology

The analysis will follow these steps:

1. Dataset acquisition and import.
2. Sample metadata inspection.
3. Data preprocessing and quality control.
4. Differential gene expression analysis.
5. Gene expression visualization.
6. Functional enrichment and pathway analysis.
7. Biological interpretation and reporting.

## 7. Project Structure

```text
DR_Bioinformatics_ICA3/
│
├── data/
│   └── README.txt
│
├── code/
│   └── R scripts
│
├── results/
│   └── Analysis outputs
│
├── report/
│   └── Final report
│
├── README.md
└── .gitignore
```

## 8. Current Progress

- Project topic selected.
- GSE160306 dataset identified.
- Gene count matrix downloaded.
- Count matrix imported into RStudio.
- Project repository structure created.

## 9. Expected Outcomes

The study is expected to identify differentially expressed genes and biological pathways associated with diabetic retinopathy progression.

## 10. References

1. NCBI Gene Expression Omnibus (GEO). GSE160306.
   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160306
