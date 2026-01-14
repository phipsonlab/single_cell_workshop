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

| Resource | Minimum | Recommended | Notes |
|----|----|----|----|
| RAM | 8 GB | 16 GB | More RAM allows analysis of larger datasets |
| Disk space | 5 GB free | 10 GB free | Needed for packages and intermediate files |
| R version | 4.3+ | 4.4+ | Earlier versions may lack required features |
| RStudio | 2023.06+ | Latest | Provides the integrated development environment |

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
    analysis. We install these using
    [`BiocManager::install()`](https://bioconductor.github.io/BiocManager/reference/install.html).

3.  **GitHub**: A code hosting platform where developers share packages
    that may be in active development or not yet submitted to
    CRAN/Bioconductor. We install these using
    [`remotes::install_github()`](https://remotes.r-lib.org/reference/install_github.html).

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

| Package | Purpose |
|----|----|
| `SingleCellExperiment` | Data structure for single cell data, required by many tools |
| `scDblFinder` | Computational detection of doublets (two cells captured together) |
| `edgeR` | Differential expression analysis using negative binomial models |
| `limma` | Linear models for differential expression, works with edgeR |
| `org.Hs.eg.db` | Human gene annotation database |
| `AnnotationDbi` | Interface to annotation databases |
| `speckle` | Statistical methods for cell type composition analysis |

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

| Package | Purpose |
|----|----|
| `Seurat` | Comprehensive toolkit for single cell analysis (clustering, visualisation, etc.) |
| `harmony` | Batch effect correction and data integration |
| `ggplot2` | Grammar of graphics for creating publication-quality plots |
| `patchwork` | Combining multiple ggplot2 plots into one figure |
| `dplyr` | Data manipulation (filtering, summarising, etc.) |
| `tidyr` | Data tidying (reshaping data frames) |
| `RColorBrewer` | Colour palettes for visualisation |
| `clustree` | Visualising how clusters change across different resolutions |
| `pheatmap` | Creating heatmaps with clustering |
| `remotes` | Installing packages from GitHub (useful for development versions) |

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

    ##   Seurat          : 5.3.0
    ##   harmony         : 1.2.4
    ##   edgeR           : 4.6.2
    ##   limma           : 3.64.0
    ##   scDblFinder     : 1.22.0

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
    ## Platform: aarch64-apple-darwin20
    ## Running under: macOS Tahoe 26.2
    ## 
    ## Matrix products: default
    ## BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: Australia/Melbourne
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] fs_1.6.6                    matrixStats_1.5.0          
    ##   [3] spatstat.sparse_3.1-0       bitops_1.0-9               
    ##   [5] httr_1.4.7                  RColorBrewer_1.1-3         
    ##   [7] tools_4.5.2                 sctransform_0.4.2          
    ##   [9] R6_2.6.1                    lazyeval_0.2.2             
    ##  [11] uwot_0.2.3                  withr_3.0.2                
    ##  [13] sp_2.2-0                    gridExtra_2.3              
    ##  [15] progressr_0.15.1            cli_3.6.5                  
    ##  [17] Biobase_2.68.0              textshaping_1.0.1          
    ##  [19] spatstat.explore_3.4-3      fastDummies_1.7.5          
    ##  [21] sass_0.4.10                 Seurat_5.3.0               
    ##  [23] spatstat.data_3.1-6         ggridges_0.5.6             
    ##  [25] pbapply_1.7-2               pkgdown_2.2.0              
    ##  [27] Rsamtools_2.24.1            systemfonts_1.2.3          
    ##  [29] harmony_1.2.4               scater_1.35.0              
    ##  [31] parallelly_1.44.0           limma_3.64.0               
    ##  [33] RSQLite_2.3.11              generics_0.1.4             
    ##  [35] BiocIO_1.18.0               ica_1.0-3                  
    ##  [37] spatstat.random_3.4-1       dplyr_1.1.4                
    ##  [39] Matrix_1.7-4                ggbeeswarm_0.7.2           
    ##  [41] S4Vectors_0.46.0            abind_1.4-8                
    ##  [43] lifecycle_1.0.4             yaml_2.3.10                
    ##  [45] edgeR_4.6.2                 SummarizedExperiment_1.38.1
    ##  [47] SparseArray_1.8.0           Rtsne_0.17                 
    ##  [49] grid_4.5.2                  blob_1.2.4                 
    ##  [51] promises_1.3.2              dqrng_0.4.1                
    ##  [53] crayon_1.5.3                miniUI_0.1.2               
    ##  [55] speckle_1.8.0               lattice_0.22-7             
    ##  [57] beachmat_2.24.0             cowplot_1.1.3              
    ##  [59] KEGGREST_1.48.1             pillar_1.10.2              
    ##  [61] knitr_1.50                  metapod_1.16.0             
    ##  [63] GenomicRanges_1.60.0        rjson_0.2.23               
    ##  [65] xgboost_3.1.3.1             future.apply_1.11.3        
    ##  [67] codetools_0.2-20            glue_1.8.0                 
    ##  [69] spatstat.univar_3.1-3       data.table_1.17.2          
    ##  [71] vctrs_0.6.5                 png_0.1-8                  
    ##  [73] spam_2.11-1                 gtable_0.3.6               
    ##  [75] cachem_1.1.0                xfun_0.52                  
    ##  [77] S4Arrays_1.8.0              mime_0.13                  
    ##  [79] tidygraph_1.3.1             survival_3.8-3             
    ##  [81] SingleCellExperiment_1.30.1 pheatmap_1.0.13            
    ##  [83] statmod_1.5.0               bluster_1.18.0             
    ##  [85] fitdistrplus_1.2-2          ROCR_1.0-11                
    ##  [87] nlme_3.1-168                bit64_4.6.0-1              
    ##  [89] RcppAnnoy_0.0.22            GenomeInfoDb_1.44.0        
    ##  [91] bslib_0.9.0                 irlba_2.3.5.1              
    ##  [93] vipor_0.4.7                 KernSmooth_2.23-26         
    ##  [95] colorspace_2.1-1            BiocGenerics_0.54.0        
    ##  [97] DBI_1.2.3                   tidyselect_1.2.1           
    ##  [99] bit_4.6.0                   compiler_4.5.2             
    ## [101] curl_7.0.0                  BiocNeighbors_2.2.0        
    ## [103] desc_1.4.3                  DelayedArray_0.34.1        
    ## [105] plotly_4.10.4               rtracklayer_1.68.0         
    ## [107] scales_1.4.0                lmtest_0.9-40              
    ## [109] stringr_1.5.1               digest_0.6.37              
    ## [111] goftest_1.2-3               spatstat.utils_3.1-4       
    ## [113] rmarkdown_2.29              XVector_0.48.0             
    ## [115] htmltools_0.5.8.1           pkgconfig_2.0.3            
    ## [117] MatrixGenerics_1.20.0       fastmap_1.2.0              
    ## [119] rlang_1.1.6                 htmlwidgets_1.6.4          
    ## [121] UCSC.utils_1.4.0            shiny_1.10.0               
    ## [123] farver_2.1.2                jquerylib_0.1.4            
    ## [125] zoo_1.8-14                  jsonlite_2.0.0             
    ## [127] BiocParallel_1.42.0         BiocSingular_1.24.0        
    ## [129] RCurl_1.98-1.17             magrittr_2.0.3             
    ## [131] scuttle_1.18.0              GenomeInfoDbData_1.2.14    
    ## [133] dotCall64_1.2               patchwork_1.3.0            
    ## [135] Rcpp_1.0.14                 viridis_0.6.5              
    ## [137] reticulate_1.42.0           stringi_1.8.7              
    ## [139] ggraph_2.2.1                MASS_7.3-65                
    ## [141] plyr_1.8.9                  org.Hs.eg.db_3.21.0        
    ## [143] parallel_4.5.2              listenv_0.9.1              
    ## [145] ggrepel_0.9.6               deldir_2.0-4               
    ## [147] scDblFinder_1.22.0          graphlayouts_1.2.2         
    ## [149] Biostrings_2.76.0           splines_4.5.2              
    ## [151] tensor_1.5                  locfit_1.5-9.12            
    ## [153] clustree_0.5.1              igraph_2.1.4               
    ## [155] spatstat.geom_3.4-1         RcppHNSW_0.6.0             
    ## [157] reshape2_1.4.4              stats4_4.5.2               
    ## [159] ScaledMatrix_1.16.0         XML_3.99-0.20              
    ## [161] evaluate_1.0.3              SeuratObject_5.1.0         
    ## [163] scran_1.36.0                tweenr_2.0.3               
    ## [165] httpuv_1.6.16               RANN_2.6.2                 
    ## [167] tidyr_1.3.1                 purrr_1.0.4                
    ## [169] polyclip_1.10-7             future_1.49.0              
    ## [171] scattermore_1.2             ggplot2_3.5.2              
    ## [173] ggforce_0.4.2               rsvd_1.0.5                 
    ## [175] xtable_1.8-4                restfulr_0.0.16            
    ## [177] RSpectra_0.16-2             later_1.4.2                
    ## [179] viridisLite_0.4.2           ragg_1.4.0                 
    ## [181] tibble_3.2.1                memoise_2.0.1              
    ## [183] beeswarm_0.4.0              AnnotationDbi_1.70.0       
    ## [185] GenomicAlignments_1.44.0    IRanges_2.42.0             
    ## [187] cluster_2.1.8.1             globals_0.18.0
