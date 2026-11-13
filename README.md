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
- **Organism:** Homo sapiens (Human)
- **Data Type:** RNA sequencing gene expression data
- **Data Format:** Raw gene count matrix
- **Number of genes:** 58,051
- **Number of samples:** 79
- **Dataset file:** `GSE160306_human_retina_DR_totalRNA_counts.txt`

Dataset link: [NCBI GEO – GSE160306](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160306)

The raw count file is stored in the `data/` folder locally. It is excluded from GitHub because of its file size. Refer to `data/README.txt` for download instructions.

## 5. Tools and Technologies

- R
- RStudio
- GEOquery
- DESeq2
- ggplot2
- pheatmap
- clusterProfiler
- org.Hs.eg.db
- Git and GitHub

## 6. Methodology

The analysis will follow these steps:

1. Dataset acquisition and import.
2. Dataset exploration and summary statistics.
3. Sample metadata inspection using GEO information.
4. Data preprocessing and quality control.
5. Differential gene expression analysis.
6. Gene expression visualization using appropriate plots, including PCA and sample count plots.
7. Functional enrichment and pathway analysis.
8. Biological interpretation and final reporting.

**Note:** Expression boxplots are not included in the current analysis plan.

## 7. Project Structure

```text
DR_Bioinformatics_ICA3/
│
├── data/
│   └── README.txt
│
├── code/
│   ├── 01_import_data.R
│   ├── 02_dataset_exploration.R
│   ├── 03_metadata_inspection.R
│   ├── 04_preprocessing_QC.R
│   ├── 05_differential_expression.R
│   └── 06_visualization_enrichment.R
│
├── results/
│   ├── plots/
│   └── tables/
│
├── report/
│
├── README.md
└── .gitignore
```

## 8. Current Progress

- Project topic selected.
- GSE160306 dataset identified from NCBI GEO.
- Raw gene count matrix downloaded.
- Count matrix successfully imported into RStudio.
- Dataset exploration completed.
- Identified 58,051 genes and 79 samples.
- Confirmed that the dataset contains no missing values.
- Confirmed that there are no duplicate Ensembl gene IDs.
- Dataset exploration summary prepared.
- GitHub repository created and initial project files pushed.
- Dataset documentation added to the `data/` folder.

## 9. Expected Outcomes

The study is expected to identify differentially expressed genes and biological pathways associated with diabetic retinopathy progression.

The results may provide insights into molecular changes associated with different stages of diabetic retinopathy.

## 10. References

1. NCBI Gene Expression Omnibus (GEO). **GSE160306**. Available at: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160306
