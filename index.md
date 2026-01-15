**Workshop website:**
<https://phipsonlab.github.io/single_cell_workshop/>

## Overview

Single-cell RNA sequencing (scRNA-seq) has revolutionised our ability to
study gene expression at the resolution of individual cells, enabling
the discovery of novel cell types and providing insights into the
cellular composition of complex tissues. This workshop provides a
comprehensive introduction to the computational analysis of scRNA-seq
data using R and Bioconductor.

We analyse single-nucleus RNA-sequencing (snRNA-seq) data from human
heart tissue across three developmental stages: foetal, young, and
adult. The dataset originates from [Sim et
al. (2021)](https://doi.org/10.1161/CIRCULATIONAHA.120.051921) examining
sex-specific control of human heart maturation (*Circulation*).

## Pre-requisites

This workshop is designed for researchers and students who:

- Have basic familiarity with R programming (data manipulation,
  plotting)
- Are interested in single-cell transcriptomics analysis
- Want to understand best practices for scRNA-seq data processing

No prior experience with single-cell analysis or Bioconductor is
required. All concepts are introduced from first principles with
detailed explanations.

## System Requirements

| Resource   | Minimum   | Recommended |
|------------|-----------|-------------|
| RAM        | 8 GB      | 16 GB       |
| Disk space | 5 GB free | 10 GB free  |
| R version  | 4.3+      | 4.5.2       |
| RStudio    | 2023.06+  | Latest      |

## Workshop Outline

### Session 1: Core Single Cell Analysis (Morning, ~3 hours)

| Module       | Topic                       | Duration |
|--------------|-----------------------------|----------|
| **Module 1** | Quality Control             | 45 min   |
|              | *Break*                     | 10 min   |
| **Module 2** | Normalisation & Integration | 50 min   |
|              | *Break*                     | 10 min   |
| **Module 3** | Cell Type Annotation        | 20 min   |
| **Module 4** | Differential Expression     | 55 min   |
|              | Wrap-up & Q&A               | 10 min   |

### Session 2: Downstream Analysis (Afternoon, ~3 hours)

*Coming soon* - Additional downstream analyses using the same heart
development dataset.

## Learning Objectives

By the end of this workshop, participants will be able to:

- Load and explore 10X Genomics scRNA-seq data in R using Seurat
- Calculate and interpret per-cell quality control metrics
- Apply appropriate filtering thresholds to remove low-quality cells
- Normalise data using SCTransform and correct batch effects with
  Harmony
- Perform graph-based clustering and visualise results with UMAP
- Annotate cell types using canonical marker genes
- Understand the pseudoreplication problem in single-cell differential
  expression
- Perform statistically rigorous differential expression analysis using
  pseudobulk methods
- Analyse cell type composition changes using propeller

## Dataset

The workshop uses snRNA-seq data from human heart tissue (Sim et al.,
2021):

| Group  | Samples | Age Range   | Description          |
|--------|---------|-------------|----------------------|
| Foetal | 3       | 19-20 weeks | Developing heart     |
| Young  | 3       | 4-14 years  | Postnatal maturation |
| Adult  | 3       | 35-42 years | Mature heart         |

**Total**: 9 samples, ~47,000 nuclei after quality control

## Methods Covered

| Analysis Step            | Method                      | Package           |
|--------------------------|-----------------------------|-------------------|
| Quality control          | Per-cell metrics, filtering | Seurat            |
| Normalisation            | SCTransform v2              | Seurat, glmGamPoi |
| Batch correction         | Harmony                     | harmony           |
| Dimensionality reduction | PCA, UMAP                   | Seurat            |
| Clustering               | Louvain algorithm           | Seurat            |
| Cell type annotation     | Marker-based (manual)       | Seurat            |
| Differential expression  | Pseudobulk + limma-voom     | edgeR, limma      |
| Composition analysis     | propeller                   | speckle           |

## Quick Start

**Please complete setup at least one day before the workshop.**

1.  **Clone or download** this repository
2.  **Open** `single_cell_workshop.Rproj` in RStudio
3.  **Follow** [Module 0:
    Setup](https://phipsonlab.github.io/single_cell_workshop/articles/00_setup.html)
    for detailed instructions

The setup involves: - Installing packages with
[`renv::restore()`](https://rstudio.github.io/renv/reference/restore.html)
(~10-15 minutes) - Downloading data from Zenodo (~420 MB, ~5 minutes)

## Key Package Versions

This workshop uses pinned package versions for reproducibility:

| Package      | Version | Package      | Version |
|--------------|---------|--------------|---------|
| R            | 4.5.2   | Bioconductor | 3.22    |
| Seurat       | 5.4.0   | edgeR        | 4.8.2   |
| SeuratObject | 5.3.0   | limma        | 3.66.0  |
| harmony      | 1.2.4   | speckle      | 1.10.0  |
| glmGamPoi    | 1.22.0  |              |         |

## Workshop Materials

### Session 1: Core Single Cell Analysis

| Module                                                                                                 | Topic           | Description                                 |
|--------------------------------------------------------------------------------------------------------|-----------------|---------------------------------------------|
| [Module 0](https://phipsonlab.github.io/single_cell_workshop/articles/00_setup.html)                   | Setup           | Environment setup and data download         |
| [Module 1](https://phipsonlab.github.io/single_cell_workshop/articles/01_quality_control.html)         | Quality Control | QC metrics, cell filtering                  |
| [Module 2](https://phipsonlab.github.io/single_cell_workshop/articles/02_integration_clustering.html)  | Integration     | Normalisation, batch correction, clustering |
| [Module 3](https://phipsonlab.github.io/single_cell_workshop/articles/03_cell_type_annotation.html)    | Annotation      | Marker genes and cell type assignment       |
| [Module 4](https://phipsonlab.github.io/single_cell_workshop/articles/04_differential_expression.html) | DE Analysis     | Pseudobulk DE and composition analysis      |

### Session 2: Downstream Analysis

*Coming soon* - Additional downstream analyses using the same heart
development dataset.

## Citation

If you use materials from this workshop, please cite:

**Original dataset:**

> Sim CB, Phipson B, Ziemann M, et al. Sex-Specific Control of Human
> Heart Maturation by the Progesterone Receptor. *Circulation*.
> 2021;143(10):1614-1628. <doi:10.1161/CIRCULATIONAHA.120.051921>

## Acknowledgements

This workshop was developed by the [Phipson
Lab](https://www.phipsonlab.org/) using data from the Porrello and
Hewitt laboratories. We thank the original authors for making their data
publicly available.

## License

This project is licensed under the MIT License - see the
[LICENSE](https://phipsonlab.github.io/single_cell_workshop/LICENSE)
file for details.
