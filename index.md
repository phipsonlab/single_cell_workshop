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

## Workshop Outline (~3 hours)

| Module       | Topic                       | Duration |
|--------------|-----------------------------|----------|
| **Module 1** | Quality Control             | 45 min   |
|              | *Break*                     | 10 min   |
| **Module 2** | Normalisation & Integration | 50 min   |
|              | *Break*                     | 10 min   |
| **Module 3** | Cell Type Annotation        | 20 min   |
| **Module 4** | Differential Expression     | 55 min   |
|              | Wrap-up & Q&A               | 10 min   |

## Learning Objectives

By the end of this workshop, participants will be able to:

- Load and explore 10X Genomics scRNA-seq data in R using Seurat
- Calculate and interpret per-cell quality control metrics
- Identify and remove doublets using computational methods
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

**Total**: 9 samples, ~43,000 nuclei after quality control

## Methods Covered

| Analysis Step            | Method                      | Package      |
|--------------------------|-----------------------------|--------------|
| Quality control          | Per-cell metrics, filtering | Seurat       |
| Doublet detection        | Simulation-based            | scDblFinder  |
| Normalisation            | SCTransform v2              | Seurat       |
| Batch correction         | Harmony                     | harmony      |
| Dimensionality reduction | PCA, UMAP                   | Seurat       |
| Clustering               | Louvain algorithm           | Seurat       |
| Cell type annotation     | Marker-based (manual)       | Seurat       |
| Differential expression  | Pseudobulk + limma-voom     | edgeR, limma |
| Composition analysis     | propeller                   | speckle      |

## Installation

### Package Versions

This workshop uses pinned package versions for reproducibility. The
tutorial outputs were generated with these exact versions:

| Package      | Version | Package      | Version |
|--------------|---------|--------------|---------|
| R            | 4.5.2   | Bioconductor | 3.22    |
| Seurat       | 5.4.0   | scDblFinder  | 1.24.0  |
| SeuratObject | 5.3.0   | edgeR        | 4.8.2   |
| harmony      | 1.2.4   | limma        | 3.66.0  |
| ggplot2      | 4.0.1   | speckle      | 1.10.0  |

### Step 1: Install Required Packages

``` r
# Install remotes and BiocManager
install.packages(c("remotes", "BiocManager"))

# Set Bioconductor version
BiocManager::install(version = "3.22", ask = FALSE)

# Install Bioconductor packages
BiocManager::install(c(
    "scDblFinder",
    "SingleCellExperiment",
    "edgeR",
    "limma",
    "org.Hs.eg.db",
    "AnnotationDbi",
    "speckle"
))

# Install CRAN packages with specific versions
remotes::install_version("Seurat", version = "5.4.0")
remotes::install_version("SeuratObject", version = "5.3.0")
remotes::install_version("harmony", version = "1.2.4")
remotes::install_version("ggplot2", version = "4.0.1")
remotes::install_version("patchwork", version = "1.3.2")
remotes::install_version("dplyr", version = "1.1.4")
remotes::install_version("tidyr", version = "1.3.2")
remotes::install_version("RColorBrewer", version = "1.1.3")
remotes::install_version("clustree", version = "0.5.1")
remotes::install_version("pheatmap", version = "1.0.13")
```

### Step 2: Download Workshop Data

The workshop data (~420 MB) is hosted on Zenodo:

``` r
# Download data from Zenodo
zenodo_record <- "18237749"
base_url <- paste0("https://zenodo.org/records/", zenodo_record, "/files/")

dir.create("data", showWarnings = FALSE)
for (f in c("heart-counts.Rds", "cellinfo_updated.Rds")) {
    download.file(paste0(base_url, f, "?download=1"), file.path("data", f), mode = "wb")
}
```

### Step 3: Verify Installation

``` r
# Check that key packages load correctly
packages <- c("Seurat", "harmony", "scDblFinder", "edgeR", "limma", "speckle")
sapply(packages, requireNamespace, quietly = TRUE)
```

## Workshop Materials

| Module                                                                                                 | Topic           | Description                                   |
|--------------------------------------------------------------------------------------------------------|-----------------|-----------------------------------------------|
| [Module 0](https://phipsonlab.github.io/single_cell_workshop/articles/00_setup.html)                   | Setup           | Environment setup and package installation    |
| [Module 1](https://phipsonlab.github.io/single_cell_workshop/articles/01_quality_control.html)         | Quality Control | QC metrics, doublet detection, cell filtering |
| [Module 2](https://phipsonlab.github.io/single_cell_workshop/articles/02_integration_clustering.html)  | Integration     | Normalisation, batch correction, clustering   |
| [Module 3](https://phipsonlab.github.io/single_cell_workshop/articles/03_cell_type_annotation.html)    | Annotation      | Marker genes and cell type assignment         |
| [Module 4](https://phipsonlab.github.io/single_cell_workshop/articles/04_differential_expression.html) | DE Analysis     | Pseudobulk DE and composition analysis        |

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
