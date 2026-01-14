# Module 0: Environment Setup

## Introduction

Welcome to the Single Cell RNA-Sequencing Workshop. This document guides
you through setting up your R environment before the workshop begins.
Single cell RNA-sequencing (scRNA-seq) has revolutionised our ability to
study gene expression at the resolution of individual cells, enabling
the discovery of novel cell types and providing insights into the
cellular composition of complex tissues.

In this workshop, we will analyse single nuclei RNA-sequencing
(snRNA-seq) data from human heart tissue across three developmental
stages: foetal, young, and adult. The dataset originates from a study
examining sex-specific control of human heart maturation (Sim et al.,
2021, *Circulation*).

### Why is Setup Important?

Single cell analysis requires several specialised R packages that work
together. Installing these packages correctly before the workshop
ensures that:

1.  You can follow along with the live coding without delays
2.  Package dependencies are resolved correctly
3.  Any installation issues can be addressed before the session

**Please complete this setup at least one day before the workshop.** The
installation process typically takes 15-30 minutes, depending on your
internet connection and system configuration.

## System Requirements

Before proceeding, it is important to verify that your system meets the
minimum requirements for running single cell analyses. Single cell
datasets can be large, and insufficient resources may cause R to crash
or analyses to run very slowly.

| Resource   | Minimum   | Recommended | Notes                                           |
|------------|-----------|-------------|-------------------------------------------------|
| RAM        | 8 GB      | 16 GB       | More RAM allows analysis of larger datasets     |
| Disk space | 5 GB free | 10 GB free  | Needed for packages and intermediate files      |
| R version  | 4.3+      | 4.4+        | Earlier versions may lack required features     |
| RStudio    | 2023.06+  | Latest      | Provides the integrated development environment |

### Checking Your R Version

To check which version of R you have installed, open RStudio and run the
following command in the console:

``` r
R.version.string
```

If your R version is older than 4.3, we recommend updating R before
installing the workshop packages. You can download the latest version
from [CRAN](https://cran.r-project.org/).

## Understanding R Package Sources

R packages come from three main sources, each serving different
purposes:

1.  **CRAN (Comprehensive R Archive Network)**: The official repository
    for general-purpose R packages. Packages here undergo basic quality
    checks and are easy to install using
    [`install.packages()`](https://rdrr.io/r/utils/install.packages.html).

2.  **Bioconductor**: A specialised repository for bioinformatics
    packages. Bioconductor packages follow stricter development
    guidelines and are designed to work together for genomic data
    analysis. We install these using `BiocManager::install()`.

3.  **GitHub**: A code hosting platform where developers share packages
    that may be in active development or not yet submitted to
    CRAN/Bioconductor. We install these using
    `remotes::install_github()`.

We install packages in a specific order (CRAN → Bioconductor → GitHub)
to ensure that dependencies are resolved correctly.

## Automatic Package Installation (Recommended)

The following code chunk automatically detects which packages are
missing from your system and installs them. This is the easiest way to
set up your environment - simply run this chunk and wait for it to
complete.

``` r
# =============================================================================
# Automatic Package Installation
# This chunk detects missing packages and installs them automatically
# =============================================================================

# Define all required packages with their sources
cran_packages <- c(
    "Seurat",
    "harmony",
    "ggplot2",
    "patchwork",
    "dplyr",
    "tidyr",
    "RColorBrewer",
    "clustree",
    "pheatmap",
    "remotes"
)

bioc_packages <- c(
    "SingleCellExperiment",
    "scDblFinder",
    "edgeR",
    "limma",
    "org.Hs.eg.db",
    "AnnotationDbi",
    "speckle"
)

# Function to check if a package is installed
is_installed <- function(pkg) {
    requireNamespace(pkg, quietly = TRUE)
}

# --------------------------------------------------------------------------
# Step 1: Ensure BiocManager is available
# --------------------------------------------------------------------------
if (!is_installed("BiocManager")) {
    message("Installing BiocManager...")
    install.packages("BiocManager", quiet = TRUE)
}

# --------------------------------------------------------------------------
# Step 2: Check and install CRAN packages
# --------------------------------------------------------------------------
missing_cran <- cran_packages[!sapply(cran_packages, is_installed)]

if (length(missing_cran) > 0) {
    message("\n=== Installing ", length(missing_cran), " missing CRAN package(s) ===")
    message("Packages: ", paste(missing_cran, collapse = ", "))
    install.packages(missing_cran, quiet = TRUE)
} else {
    message("\n=== All CRAN packages are already installed ===")
}

# --------------------------------------------------------------------------
# Step 3: Check and install Bioconductor packages
# --------------------------------------------------------------------------
missing_bioc <- bioc_packages[!sapply(bioc_packages, is_installed)]

if (length(missing_bioc) > 0) {
    message("\n=== Installing ", length(missing_bioc), " missing Bioconductor package(s) ===")
    message("Packages: ", paste(missing_bioc, collapse = ", "))
    BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
} else {
    message("\n=== All Bioconductor packages are already installed ===")
}

# --------------------------------------------------------------------------
# Final verification
# --------------------------------------------------------------------------
all_packages <- c(cran_packages, bioc_packages)
final_check <- sapply(all_packages, is_installed)
missing_final <- names(final_check)[!final_check]

if (length(missing_final) == 0) {
    message("\n", paste(rep("=", 50), collapse = ""))
    message("SUCCESS! All ", length(all_packages), " packages are installed.")
    message(paste(rep("=", 50), collapse = ""))
    message("\nYour environment is ready for the workshop!")
} else {
    message("\n", paste(rep("=", 50), collapse = ""))
    message("WARNING: The following packages could not be installed:")
    for (pkg in missing_final) {
        message("  - ", pkg)
    }
    message(paste(rep("=", 50), collapse = ""))
    message("\nPlease see the Troubleshooting section below or try")
    message("installing these packages manually using the steps that follow.")
}
```

If you prefer to install packages manually, or if the automatic
installation fails for some packages, follow the step-by-step
instructions below.

------------------------------------------------------------------------

## Manual Package Installation

### Step 1: Install BiocManager

*BiocManager* is a package that manages installation of Bioconductor
packages. It ensures that you get compatible versions of all
Bioconductor packages. If you do not already have BiocManager installed,
run the following code:

``` r
# Check if BiocManager is installed; if not, install it
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
```

The [`require()`](https://rdrr.io/r/base/library.html) function attempts
to load a package and returns `TRUE` if successful, `FALSE` otherwise.
The `quietly = TRUE` argument suppresses messages if the package is not
found.

### Step 2: Install Bioconductor Packages

The following Bioconductor packages provide essential functionality for
our single cell analysis:

| Package                | Purpose                                                           |
|------------------------|-------------------------------------------------------------------|
| `SingleCellExperiment` | Data structure for single cell data, required by many tools       |
| `scDblFinder`          | Computational detection of doublets (two cells captured together) |
| `edgeR`                | Differential expression analysis using negative binomial models   |
| `limma`                | Linear models for differential expression, works with edgeR       |
| `org.Hs.eg.db`         | Human gene annotation database                                    |
| `AnnotationDbi`        | Interface to annotation databases                                 |
| `speckle`              | Statistical methods for cell type composition analysis            |

``` r
BiocManager::install(c(
    "SingleCellExperiment",
    "scDblFinder",
    "edgeR",
    "limma",
    "org.Hs.eg.db",
    "AnnotationDbi",
    "speckle"
))
```

When prompted “Update all/some/none? \[a/s/n\]:”, you can type `a` to
update all packages, or `n` if you prefer not to update existing
packages.

### Step 3: Install CRAN Packages

The following CRAN packages form the core of our analysis toolkit:

| Package        | Purpose                                                                          |
|----------------|----------------------------------------------------------------------------------|
| `Seurat`       | Comprehensive toolkit for single cell analysis (clustering, visualisation, etc.) |
| `harmony`      | Batch effect correction and data integration                                     |
| `ggplot2`      | Grammar of graphics for creating publication-quality plots                       |
| `patchwork`    | Combining multiple ggplot2 plots into one figure                                 |
| `dplyr`        | Data manipulation (filtering, summarising, etc.)                                 |
| `tidyr`        | Data tidying (reshaping data frames)                                             |
| `RColorBrewer` | Colour palettes for visualisation                                                |
| `clustree`     | Visualising how clusters change across different resolutions                     |
| `pheatmap`     | Creating heatmaps with clustering                                                |
| `remotes`      | Installing packages from GitHub (useful for development versions)                |

``` r
install.packages(c(
    "Seurat",
    "harmony",
    "ggplot2",
    "patchwork",
    "dplyr",
    "tidyr",
    "RColorBrewer",
    "clustree",
    "pheatmap",
    "remotes"
))
```

## Verify Installation

Once all packages are installed, it is essential to verify that they
load correctly. A package may fail to install silently, or there may be
version conflicts that prevent loading. The following code attempts to
load each required package and reports any failures:

``` r
# List of required packages
packages <- c(
    "Seurat",
    "harmony",
    "scDblFinder",
    "SingleCellExperiment",
    "edgeR",
    "limma",
    "speckle",
    "org.Hs.eg.db",
    "ggplot2",
    "patchwork",
    "dplyr",
    "tidyr",
    "pheatmap",
    "clustree"
)

# Check which packages can be loaded
check_results <- sapply(packages, function(pkg) {
    requireNamespace(pkg, quietly = TRUE)
})
```

    ## 

``` r
# Report results
missing <- names(check_results)[!check_results]

if (length(missing) > 0) {
    message("=== WARNING ===")
    message("The following packages failed to install or load:")
    for (pkg in missing) {
        message("  - ", pkg)
    }
    message("\nPlease try reinstalling these packages before the workshop.")
} else {
    message("=== SUCCESS ===")
    message("All required packages are installed and can be loaded.")
    message("Your environment is ready for the workshop!")
}
```

    ## === SUCCESS ===

    ## All required packages are installed and can be loaded.

    ## Your environment is ready for the workshop!

### Check Package Versions

Different versions of packages may behave differently. Here, we print
the versions of key packages to help with troubleshooting if issues
arise during the workshop:

``` r
# Print versions of key packages
key_packages <- c("Seurat", "harmony", "edgeR", "limma", "scDblFinder")

cat("Key package versions:\n")
```

    ## Key package versions:

``` r
cat(paste(rep("-", 40), collapse = ""), "\n")
```

    ## ----------------------------------------

``` r
for (pkg in key_packages) {
    if (requireNamespace(pkg, quietly = TRUE)) {
        cat(sprintf("  %-15s : %s\n", pkg, as.character(packageVersion(pkg))))
    } else {
        cat(sprintf("  %-15s : NOT INSTALLED\n", pkg))
    }
}
```

    ##   Seurat          : 5.4.0
    ##   harmony         : 1.2.4
    ##   edgeR           : 4.8.2
    ##   limma           : 3.66.0
    ##   scDblFinder     : 1.24.0

## Download Workshop Data

The workshop uses single-nucleus RNA-seq data from human heart tissue.
The data files are hosted on Zenodo and must be downloaded before
running the workshop modules.

### About the Dataset

The dataset contains ~54,000 nuclei from 9 human heart samples across
three developmental stages:

| Group  | Samples | Age Range             |
|--------|---------|-----------------------|
| Foetal | 3       | 19-20 weeks gestation |
| Young  | 3       | 4-14 years            |
| Adult  | 3       | 35-42 years           |

Data source: Sim et al. (2021) “Sex-Specific Control of Human Heart
Maturation by the Progesterone Receptor”, *Circulation*. DOI:
[10.1161/CIRCULATIONAHA.120.051921](https://doi.org/10.1161/CIRCULATIONAHA.120.051921)

### Download the Data

Run the following code to download the workshop data files from Zenodo.
This will create a `data/` folder in your working directory containing
the required files. **Note:** The total download size is approximately
200-300 MB. Download time depends on your internet connection.

``` r
# Download workshop data from Zenodo
# This will create a data/ folder with the required files

# Zenodo record ID
zenodo_record <- "18237749"
base_url <- paste0("https://zenodo.org/records/", zenodo_record, "/files/")

# Files to download
files <- c("heart-counts.Rds", "cellinfo_updated.Rds")

# Create data directory
if (!dir.exists("data")) {
    dir.create("data")
    message("Created data/ directory")
}

# Download each file
for (f in files) {
    dest_file <- file.path("data", f)
    if (file.exists(dest_file)) {
        message(f, " already exists, skipping...")
        next
    }
    message("Downloading ", f, "...")
    download.file(
        url = paste0(base_url, f, "?download=1"),
        destfile = dest_file,
        mode = "wb"
    )
    message("  Done!")
}

message("\nDownload complete! Data saved to data/ folder.")
```

## Test Data Loading

As a final verification step, we check that the workshop data files are
accessible and can be loaded correctly. The workshop uses two data
files:

1.  **heart-counts.Rds**: A sparse matrix containing gene expression
    counts (genes × cells)
2.  **cellinfo_updated.Rds**: A data frame containing metadata for each
    cell (sample, group, etc.)

``` r
# Set path to data directory (relative to tutorials folder)
data_dir <- "../data"

# Check that data files exist
counts_file <- file.path(data_dir, "heart-counts.Rds")
cellinfo_file <- file.path(data_dir, "cellinfo_updated.Rds")

cat("Checking for workshop data files...\n")
```

    ## Checking for workshop data files...

``` r
cat(paste(rep("-", 40), collapse = ""), "\n")
```

    ## ----------------------------------------

``` r
if (file.exists(counts_file)) {
    cat("  heart-counts.Rds    : FOUND\n")

    # Load and report dimensions
    counts <- readRDS(counts_file)
    cat(sprintf("    - Dimensions: %d genes x %d cells\n", nrow(counts), ncol(counts)))
    rm(counts)  # Free memory immediately

} else {
    cat("  heart-counts.Rds    : NOT FOUND\n")
}
```

    ##   heart-counts.Rds    : FOUND
    ##     - Dimensions: 33939 genes x 54140 cells

``` r
if (file.exists(cellinfo_file)) {
    cat("  cellinfo_updated.Rds: FOUND\n")
} else {
    cat("  cellinfo_updated.Rds: NOT FOUND\n")
}
```

    ##   cellinfo_updated.Rds: FOUND

``` r
if (file.exists(counts_file) && file.exists(cellinfo_file)) {
    cat("\n=== Data files are ready! ===\n")
} else {
    cat("\n=== WARNING: Some data files are missing ===\n")
    cat("Please check the data directory path.\n")
}
```

    ## 
    ## === Data files are ready! ===

## Troubleshooting Common Issues

If you encounter problems during installation, this section provides
solutions to the most common issues.

### Issue 1: “Package ‘X’ is not available for R version Y.Z”

**Cause**: Your R version is too old for the package.

**Solution**: Update R to version 4.3 or later from
[CRAN](https://cran.r-project.org/), then restart RStudio and retry the
installation.

### Issue 2: “There is no package called ‘X’”

**Cause**: The package failed to install, possibly due to missing
dependencies or network issues.

**Solution**: Try installing the specific package again:

``` r
# For CRAN packages:
install.packages("package_name")

# For Bioconductor packages:
BiocManager::install("package_name")

# For GitHub packages:
remotes::install_github("username/package_name")
```

### Issue 3: Compilation Errors on macOS

**Cause**: Missing Xcode Command Line Tools, which provide compilers
needed to build some packages.

**Solution**: Open Terminal (not R) and run:

``` bash
xcode-select --install
```

Follow the prompts to install, then restart RStudio and retry package
installation.

### Issue 4: Memory Errors (“cannot allocate vector of size…”)

**Cause**: Insufficient RAM available.

**Solution**: 1. Close other applications to free memory 2. Restart R
(Session → Restart R in RStudio) 3. Try installing packages one at a
time rather than all at once

### Issue 5: Permission Errors on Windows

**Cause**: R does not have permission to write to the package library.

**Solution**: Run RStudio as Administrator (right-click RStudio icon →
“Run as administrator”).

### Issue 6: Package Loading Conflicts

**Cause**: Multiple versions of a package are installed, or packages
have conflicting dependencies.

**Solution**: Remove and reinstall the problematic package:

``` r
# Remove the package
remove.packages("package_name")

# Restart R, then reinstall
install.packages("package_name")
```

## Getting Help

If you encounter issues that you cannot resolve using this
troubleshooting guide, please:

1.  **Note the exact error message** - Copy the complete error text
2.  **Record your R version** - Run `R.version.string`
3.  **List installed package versions** - Run
    [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html)
4.  **Contact the workshop organisers** before the session

Having setup issues resolved before the workshop begins ensures that
everyone can participate fully in the hands-on exercises.

## What to Expect in the Workshop

Once your environment is set up, you will be ready for the workshop,
which covers:

1.  **Module 1: Quality Control** - Loading data, calculating QC
    metrics, filtering low-quality cells
2.  **Module 2: Integration & Clustering** - Normalisation, batch
    correction, dimensionality reduction, and clustering
3.  **Module 3: Cell Type Annotation** - Identifying cell types using
    marker genes
4.  **Module 4: Differential Expression** - Finding genes and cell types
    that differ between conditions

We look forward to seeing you at the workshop!

## Session Information

The output below records your R environment configuration, which is
useful for reproducibility and troubleshooting:

``` r
sessionInfo()
```

    ## R version 4.5.2 (2025-10-31)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.3 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    ##  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    ##  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    ## [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] fs_1.6.6                    matrixStats_1.5.0          
    ##   [3] spatstat.sparse_3.1-0       bitops_1.0-9               
    ##   [5] httr_1.4.7                  RColorBrewer_1.1-3         
    ##   [7] tools_4.5.2                 sctransform_0.4.3          
    ##   [9] R6_2.6.1                    lazyeval_0.2.2             
    ##  [11] uwot_0.2.4                  withr_3.0.2                
    ##  [13] sp_2.2-0                    gridExtra_2.3              
    ##  [15] progressr_0.18.0            cli_3.6.5                  
    ##  [17] Biobase_2.70.0              textshaping_1.0.4          
    ##  [19] spatstat.explore_3.6-0      fastDummies_1.7.5          
    ##  [21] sass_0.4.10                 Seurat_5.4.0               
    ##  [23] S7_0.2.1                    spatstat.data_3.1-9        
    ##  [25] ggridges_0.5.7              pbapply_1.7-4              
    ##  [27] pkgdown_2.2.0               Rsamtools_2.26.0           
    ##  [29] systemfonts_1.3.1           harmony_1.2.4              
    ##  [31] scater_1.38.0               parallelly_1.46.1          
    ##  [33] limma_3.66.0                RSQLite_2.4.5              
    ##  [35] generics_0.1.4              BiocIO_1.20.0              
    ##  [37] ica_1.0-3                   spatstat.random_3.4-3      
    ##  [39] dplyr_1.1.4                 Matrix_1.7-4               
    ##  [41] ggbeeswarm_0.7.3            S4Vectors_0.48.0           
    ##  [43] abind_1.4-8                 lifecycle_1.0.5            
    ##  [45] yaml_2.3.12                 edgeR_4.8.2                
    ##  [47] SummarizedExperiment_1.40.0 SparseArray_1.10.8         
    ##  [49] Rtsne_0.17                  grid_4.5.2                 
    ##  [51] blob_1.2.4                  promises_1.5.0             
    ##  [53] dqrng_0.4.1                 crayon_1.5.3               
    ##  [55] miniUI_0.1.2                speckle_1.10.0             
    ##  [57] lattice_0.22-7              beachmat_2.26.0            
    ##  [59] cowplot_1.2.0               cigarillo_1.0.0            
    ##  [61] KEGGREST_1.50.0             pillar_1.11.1              
    ##  [63] knitr_1.51                  metapod_1.18.0             
    ##  [65] GenomicRanges_1.62.1        rjson_0.2.23               
    ##  [67] xgboost_3.1.3.1             future.apply_1.20.1        
    ##  [69] codetools_0.2-20            glue_1.8.0                 
    ##  [71] spatstat.univar_3.1-5       data.table_1.18.0          
    ##  [73] vctrs_0.6.5                 png_0.1-8                  
    ##  [75] spam_2.11-3                 gtable_0.3.6               
    ##  [77] cachem_1.1.0                xfun_0.55                  
    ##  [79] S4Arrays_1.10.1             mime_0.13                  
    ##  [81] tidygraph_1.3.1             Seqinfo_1.0.0              
    ##  [83] survival_3.8-3              SingleCellExperiment_1.32.0
    ##  [85] pheatmap_1.0.13             statmod_1.5.1              
    ##  [87] bluster_1.20.0              fitdistrplus_1.2-4         
    ##  [89] ROCR_1.0-11                 nlme_3.1-168               
    ##  [91] bit64_4.6.0-1               RcppAnnoy_0.0.23           
    ##  [93] GenomeInfoDb_1.46.2         bslib_0.9.0                
    ##  [95] irlba_2.3.5.1               vipor_0.4.7                
    ##  [97] KernSmooth_2.23-26          otel_0.2.0                 
    ##  [99] BiocGenerics_0.56.0         DBI_1.2.3                  
    ## [101] tidyselect_1.2.1            bit_4.6.0                  
    ## [103] compiler_4.5.2              curl_7.0.0                 
    ## [105] BiocNeighbors_2.4.0         desc_1.4.3                 
    ## [107] DelayedArray_0.36.0         plotly_4.11.0              
    ## [109] rtracklayer_1.70.1          scales_1.4.0               
    ## [111] lmtest_0.9-40               stringr_1.6.0              
    ## [113] digest_0.6.39               goftest_1.2-3              
    ## [115] spatstat.utils_3.2-1        rmarkdown_2.30             
    ## [117] XVector_0.50.0              htmltools_0.5.9            
    ## [119] pkgconfig_2.0.3             MatrixGenerics_1.22.0      
    ## [121] fastmap_1.2.0               rlang_1.1.7                
    ## [123] htmlwidgets_1.6.4           UCSC.utils_1.6.1           
    ## [125] shiny_1.12.1                farver_2.1.2               
    ## [127] jquerylib_0.1.4             zoo_1.8-15                 
    ## [129] jsonlite_2.0.0              BiocParallel_1.44.0        
    ## [131] BiocSingular_1.26.1         RCurl_1.98-1.17            
    ## [133] magrittr_2.0.4              scuttle_1.20.0             
    ## [135] dotCall64_1.2               patchwork_1.3.2            
    ## [137] Rcpp_1.1.1                  viridis_0.6.5              
    ## [139] reticulate_1.44.1           stringi_1.8.7              
    ## [141] ggraph_2.2.2                MASS_7.3-65                
    ## [143] plyr_1.8.9                  org.Hs.eg.db_3.22.0        
    ## [145] parallel_4.5.2              listenv_0.10.0             
    ## [147] ggrepel_0.9.6               deldir_2.0-4               
    ## [149] scDblFinder_1.24.0          graphlayouts_1.2.2         
    ## [151] Biostrings_2.78.0           splines_4.5.2              
    ## [153] tensor_1.5.1                locfit_1.5-9.12            
    ## [155] clustree_0.5.1              igraph_2.2.1               
    ## [157] spatstat.geom_3.6-1         RcppHNSW_0.6.0             
    ## [159] reshape2_1.4.5              stats4_4.5.2               
    ## [161] ScaledMatrix_1.18.0         XML_3.99-0.20              
    ## [163] evaluate_1.0.5              SeuratObject_5.3.0         
    ## [165] scran_1.38.0                tweenr_2.0.3               
    ## [167] httpuv_1.6.16               RANN_2.6.2                 
    ## [169] tidyr_1.3.2                 purrr_1.2.1                
    ## [171] polyclip_1.10-7             future_1.68.0              
    ## [173] scattermore_1.2             ggplot2_4.0.1              
    ## [175] ggforce_0.5.0               rsvd_1.0.5                 
    ## [177] xtable_1.8-4                restfulr_0.0.16            
    ## [179] RSpectra_0.16-2             later_1.4.5                
    ## [181] viridisLite_0.4.2           ragg_1.5.0                 
    ## [183] tibble_3.3.1                memoise_2.0.1              
    ## [185] beeswarm_0.4.0              AnnotationDbi_1.72.0       
    ## [187] GenomicAlignments_1.46.0    IRanges_2.44.0             
    ## [189] cluster_2.1.8.1             globals_0.18.0
