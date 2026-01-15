# Module 0: Environment Setup

## Introduction

Welcome to the Single Cell RNA-Sequencing Workshop. This document guides
you through setting up your R environment **before the workshop
begins**.

**Please complete this setup at least one day before the workshop.**
Package installation takes approximately 10-15 minutes, and data
download takes approximately 5 minutes depending on your internet
connection.

### What We Will Set Up

1.  **R packages** using `renv` for reproducibility
2.  **Workshop data** downloaded from Zenodo (~420 MB)

### System Requirements

| Resource   | Minimum   | Recommended |
|------------|-----------|-------------|
| RAM        | 8 GB      | 16 GB       |
| Disk space | 5 GB free | 10 GB free  |
| R version  | 4.3+      | 4.5+        |
| RStudio    | 2023.06+  | Latest      |

Check your R version:

``` r
R.version.string
```

    ## [1] "R version 4.5.2 (2025-10-31)"

If your R version is older than 4.3, please update from
[CRAN](https://cran.r-project.org/) before proceeding.

## Step 1: Clone or Download the Repository

First, obtain the workshop materials:

**Option A: Clone with Git (Recommended)**

``` bash
git clone https://github.com/phipsonlab/single_cell_workshop.git
cd single_cell_workshop
```

**Option B: Download ZIP**

1.  Visit <https://github.com/phipsonlab/single_cell_workshop>
2.  Click the green “Code” button
3.  Select “Download ZIP”
4.  Extract the ZIP file

Then open `single_cell_workshop.Rproj` in RStudio.

## Step 2: Install Packages with renv

This workshop uses **renv** for reproducible package management. All
package versions are locked in `renv.lock`, ensuring everyone has
identical environments.

When you open the project in RStudio, renv should automatically
bootstrap itself. If prompted to install renv, select “Yes”.

Then run:

``` r
# Install all packages with exact versions from renv.lock
# This may take 10-15 minutes on first run
renv::restore()
```

When prompted “Do you want to proceed?”, type `y` and press Enter.

### What renv Does

- Installs packages to a project-specific library (not your global R
  library)
- Uses exact package versions from `renv.lock`
- Ensures reproducibility across different machines

### Key Packages Installed

| Package   | Version | Purpose                     |
|-----------|---------|-----------------------------|
| Seurat    | 5.4.0   | Core single cell analysis   |
| harmony   | 1.2.4   | Batch correction            |
| glmGamPoi | 1.22.0  | Fast SCTransform (critical) |
| edgeR     | 4.8.2   | Differential expression     |
| limma     | 3.66.0  | Linear models               |
| speckle   | 1.10.0  | Composition analysis        |

## Step 3: Verify Package Installation

After
[`renv::restore()`](https://rstudio.github.io/renv/reference/restore.html)
completes, verify that all critical packages are installed:

``` r
packages <- c("Seurat", "harmony", "glmGamPoi", "edgeR", "limma", "speckle")

check <- sapply(packages, requireNamespace, quietly = TRUE)
```

    ## Warning: replacing previous import 'S4Arrays::makeNindexFromArrayViewport' by
    ## 'DelayedArray::makeNindexFromArrayViewport' when loading 'SummarizedExperiment'

``` r
missing <- names(check)[!check]

if (length(missing) == 0) {
    message("All critical packages installed successfully!")
} else {
    message("WARNING: Missing packages: ", paste(missing, collapse = ", "))
    message("Try running renv::restore() again.")
}
```

    ## All critical packages installed successfully!

## Step 4: Download Workshop Data

The workshop data files are hosted on Zenodo. Run the following code to
download them:

``` r
# Zenodo record for workshop data
zenodo_url <- "https://zenodo.org/records/18237749/files/"
files <- c("heart-counts.Rds", "cellinfo_updated.Rds")

# Create data directory if needed
if (!dir.exists("data")) dir.create("data")

# Download each file
for (f in files) {
    dest <- file.path("data", f)
    if (file.exists(dest)) {
        message(f, " already exists, skipping")
        next
    }
    message("Downloading ", f, "...")
    download.file(
        url = paste0(zenodo_url, f, "?download=1"),
        destfile = dest,
        mode = "wb"
    )
}
message("Download complete!")
```

**Note:** Total download size is approximately 420 MB.

## Step 5: Final Verification

Run this final check to ensure everything is ready:

``` r
cat("=== Workshop Setup Verification ===\n\n")
```

    ## === Workshop Setup Verification ===

``` r
# Check R version
cat("R Version:", R.version.string, "\n\n")
```

    ## R Version: R version 4.5.2 (2025-10-31)

``` r
# Check critical packages
cat("Package Status:\n")
```

    ## Package Status:

``` r
packages <- c("Seurat", "harmony", "glmGamPoi", "edgeR", "limma", "speckle")
for (pkg in packages) {
    status <- if (requireNamespace(pkg, quietly = TRUE)) "OK" else "MISSING"
    cat(sprintf("  %-12s %s\n", pkg, status))
}
```

    ##   Seurat       OK
    ##   harmony      OK
    ##   glmGamPoi    OK
    ##   edgeR        OK
    ##   limma        OK
    ##   speckle      OK

``` r
# Check data files
cat("\nData Files:\n")
```

    ## 
    ## Data Files:

``` r
data_dir <- if (dir.exists("data")) "data" else "../data"
for (f in c("heart-counts.Rds", "cellinfo_updated.Rds")) {
    path <- file.path(data_dir, f)
    status <- if (file.exists(path)) "FOUND" else "NOT FOUND"
    cat(sprintf("  %-25s %s\n", f, status))
}
```

    ##   heart-counts.Rds          FOUND
    ##   cellinfo_updated.Rds      FOUND

``` r
cat("\n================================\n")
```

    ## 
    ## ================================

## Troubleshooting

### renv::restore() fails

If
[`renv::restore()`](https://rstudio.github.io/renv/reference/restore.html)
encounters errors:

1.  **Restart R** (Session \> Restart R)
2.  Run
    [`renv::restore()`](https://rstudio.github.io/renv/reference/restore.html)
    again
3.  If specific packages fail, try installing them manually then run
    [`renv::restore()`](https://rstudio.github.io/renv/reference/restore.html)

### Memory errors

If you encounter “cannot allocate vector” errors:

1.  Close other applications
2.  Restart R
3.  Try again

### macOS compilation errors

Install Xcode Command Line Tools:

``` bash
xcode-select --install
```

### Package conflicts

If you have existing packages causing conflicts:

``` r
# Use a clean renv library
renv::rebuild()
```

## Getting Help

If you cannot resolve setup issues:

1.  Note the exact error message
2.  Run [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html) and
    save the output
3.  Contact the workshop organisers before the session

## Session Information

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
    ## [1] stats     graphics  grDevices datasets  utils     methods   base     
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] RColorBrewer_1.1-3          jsonlite_2.0.0             
    ##   [3] magrittr_2.0.4              spatstat.utils_3.2-1       
    ##   [5] farver_2.1.2                rmarkdown_2.30             
    ##   [7] fs_1.6.6                    ragg_1.5.0                 
    ##   [9] vctrs_0.6.5                 ROCR_1.0-11                
    ##  [11] spatstat.explore_3.6-0      S4Arrays_1.10.1            
    ##  [13] htmltools_0.5.9             SparseArray_1.10.8         
    ##  [15] sass_0.4.10                 sctransform_0.4.3          
    ##  [17] parallelly_1.46.1           KernSmooth_2.23-26         
    ##  [19] bslib_0.9.0                 htmlwidgets_1.6.4          
    ##  [21] desc_1.4.3                  ica_1.0-3                  
    ##  [23] plyr_1.8.9                  speckle_1.10.0             
    ##  [25] plotly_4.11.0               zoo_1.8-15                 
    ##  [27] cachem_1.1.0                igraph_2.2.1               
    ##  [29] mime_0.13                   lifecycle_1.0.5            
    ##  [31] pkgconfig_2.0.3             Matrix_1.7-4               
    ##  [33] R6_2.6.1                    fastmap_1.2.0              
    ##  [35] MatrixGenerics_1.22.0       fitdistrplus_1.2-4         
    ##  [37] future_1.68.0               shiny_1.12.1               
    ##  [39] digest_0.6.39               patchwork_1.3.2            
    ##  [41] S4Vectors_0.48.0            Seurat_5.4.0               
    ##  [43] tensor_1.5.1                RSpectra_0.16-2            
    ##  [45] irlba_2.3.5.1               GenomicRanges_1.62.1       
    ##  [47] textshaping_1.0.4           beachmat_2.26.0            
    ##  [49] progressr_0.18.0            spatstat.sparse_3.1-0      
    ##  [51] httr_1.4.7                  polyclip_1.10-7            
    ##  [53] abind_1.4-8                 compiler_4.5.2             
    ##  [55] S7_0.2.1                    fastDummies_1.7.5          
    ##  [57] MASS_7.3-65                 DelayedArray_0.36.0        
    ##  [59] tools_4.5.2                 lmtest_0.9-40              
    ##  [61] otel_0.2.0                  httpuv_1.6.16              
    ##  [63] future.apply_1.20.1         goftest_1.2-3              
    ##  [65] glmGamPoi_1.22.0            glue_1.8.0                 
    ##  [67] nlme_3.1-168                promises_1.5.0             
    ##  [69] grid_4.5.2                  Rtsne_0.17                 
    ##  [71] cluster_2.1.8.1             reshape2_1.4.5             
    ##  [73] generics_0.1.4              gtable_0.3.6               
    ##  [75] spatstat.data_3.1-9         tidyr_1.3.2                
    ##  [77] data.table_1.18.0           XVector_0.50.0             
    ##  [79] sp_2.2-0                    BiocGenerics_0.56.0        
    ##  [81] spatstat.geom_3.6-1         RcppAnnoy_0.0.23           
    ##  [83] ggrepel_0.9.6               RANN_2.6.2                 
    ##  [85] pillar_1.11.1               stringr_1.6.0              
    ##  [87] limma_3.66.0                spam_2.11-3                
    ##  [89] RcppHNSW_0.6.0              later_1.4.5                
    ##  [91] splines_4.5.2               dplyr_1.1.4                
    ##  [93] lattice_0.22-7              renv_1.1.5                 
    ##  [95] survival_3.8-3              deldir_2.0-4               
    ##  [97] tidyselect_1.2.1            SingleCellExperiment_1.32.0
    ##  [99] locfit_1.5-9.12             miniUI_0.1.2               
    ## [101] pbapply_1.7-4               knitr_1.51                 
    ## [103] gridExtra_2.3               Seqinfo_1.0.0              
    ## [105] IRanges_2.44.0              edgeR_4.8.2                
    ## [107] SummarizedExperiment_1.40.0 scattermore_1.2            
    ## [109] stats4_4.5.2                xfun_0.55                  
    ## [111] Biobase_2.70.0              statmod_1.5.1              
    ## [113] matrixStats_1.5.0           stringi_1.8.7              
    ## [115] lazyeval_0.2.2              yaml_2.3.12                
    ## [117] evaluate_1.0.5              codetools_0.2-20           
    ## [119] tibble_3.3.1                BiocManager_1.30.27        
    ## [121] cli_3.6.5                   uwot_0.2.4                 
    ## [123] xtable_1.8-4                reticulate_1.44.1          
    ## [125] systemfonts_1.3.1           jquerylib_0.1.4            
    ## [127] harmony_1.2.4               Rcpp_1.1.1                 
    ## [129] globals_0.18.0              spatstat.random_3.4-3      
    ## [131] png_0.1-8                   spatstat.univar_3.1-5      
    ## [133] parallel_4.5.2              pkgdown_2.2.0              
    ## [135] ggplot2_4.0.1               dotCall64_1.2              
    ## [137] listenv_0.10.0              viridisLite_0.4.2          
    ## [139] scales_1.4.0                ggridges_0.5.7             
    ## [141] SeuratObject_5.3.0          purrr_1.2.1                
    ## [143] rlang_1.1.7                 cowplot_1.2.0
