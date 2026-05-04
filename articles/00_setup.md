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

1.  **System build tools** so packages with C/C++ code can compile
2.  **All R packages** the workshop needs, in one combined install step
3.  **Workshop data** downloaded from Zenodo (~420 MB)

### Supported Platforms

The workshop has been tested on **Windows 10/11**, **macOS 12+ (Intel
and Apple Silicon)**, and **Ubuntu 22.04+ / equivalent Linux**. Each
platform needs a one-time toolchain install (Step 2 below); after that
the same `renv.lock` is used everywhere.

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

> **Windows users — work from a short path.** Long Windows paths (deeply
> nested under `C:\Users\<name>\OneDrive\Documents\…`) can hit the
> legacy 260-character `MAX_PATH` limit when `renv` installs packages
> with deep internal directory trees, producing cryptic
> `cannot rename file` errors. If you can, clone to something like
> `C:\workshop\single_cell_workshop`. If you can’t, enable [Long Paths
> in
> Windows](https://learn.microsoft.com/en-us/windows/win32/fileio/maximum-file-path-limitation?tabs=registry).

## Step 2: Install System Build Tools

Several workshop packages contain C / C++ code that R must compile
during installation (`destiny`, `harmony`, `PhiSpace`, `NeighbourNet`,
…). Each operating system has its own one-line install for the required
toolchain.

**Windows — install Rtools45**

R 4.5.x on Windows needs **Rtools45**:

1.  Download the installer from
    <https://cran.r-project.org/bin/windows/Rtools/rtools45/rtools.html>.
2.  Run the installer with the default options (it adds itself to
    `PATH`).
3.  Restart R.

You can verify Rtools is detected from R:

``` r

# Run only on Windows; on macOS / Linux this returns "" or NA.
if (.Platform$OS.type == "windows") {
    cat("Rtools detected:", pkgbuild::has_build_tools(debug = FALSE), "\n")
}
```

If `pkgbuild` is not yet installed, run `install.packages("pkgbuild")`
first. A `FALSE` result means Rtools is missing or not on `PATH` —
re-run the installer or restart R / RStudio.

**macOS — install the Xcode Command Line Tools**

``` bash
xcode-select --install
```

Click “Install” in the pop-up dialog. This is needed on both Intel and
Apple Silicon Macs. If you already use Homebrew or have Xcode itself,
you already have the toolchain.

**Linux (Ubuntu / Debian) — install the build essentials and a few dev
headers**

``` bash
sudo apt-get update
sudo apt-get install -y \
    build-essential \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libgit2-dev \
    libglpk-dev \
    libgmp3-dev \
    libhdf5-dev
```

(On Fedora / RHEL the equivalent is
`dnf groupinstall "Development Tools"` plus the `*-devel` versions of
the libraries above.)

## Step 3: Install R Packages

This workshop uses **renv** for the core packages and `BiocManager` /
`remotes` for a handful of Bioconductor and GitHub-only extras that are
not yet captured in `renv.lock`. The chunks below install everything
needed for the whole workshop in one pass — run them in order.

When you open the project in RStudio, renv should automatically
bootstrap itself. If prompted to install renv, select “Yes”.

### 3a. Restore the locked package set

``` r

# Install all packages with exact versions from renv.lock.
# This may take 10-15 minutes on first run.
renv::restore()
```

When prompted “Do you want to proceed?”, type `y` and press Enter.

[`renv::restore()`](https://rstudio.github.io/renv/reference/restore.html)
installs packages into a project-specific library (not your global R
library) at the exact versions recorded in `renv.lock`, so everyone in
the workshop ends up with the same environment.

### 3b. Install the extras (Φ-Space, pseudotime, NeighbourNet)

A few packages used in the later modules are not yet in `renv.lock`:
`ComplexHeatmap`, `slingshot`, `destiny`, `scater` from Bioconductor,
plus `PhiSpace` and `NeighbourNet` from GitHub. The two GitHub packages
ship as source and must compile, which is why Step 2 (build tools) had
to come first.

``` r

# 1. Sanity check that the build toolchain is available. On Windows this
#    looks for Rtools; on macOS / Linux it returns TRUE when the system
#    gcc / clang toolchain is present.
if (!requireNamespace("pkgbuild", quietly = TRUE)) {
    install.packages("pkgbuild")
}
if (!pkgbuild::has_build_tools(debug = FALSE)) {
    stop(
        "No build tools detected. Re-run Step 2 of this notebook to ",
        "install the right toolchain for your operating system, then ",
        "restart R and try again."
    )
}

# 2. Bioconductor extras (binary on Windows / macOS, source on Linux).
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
}

# `smoother` is a transitive dependency of `destiny` that was archived
# from CRAN on 2025-12-19. BiocManager / install.packages can no longer
# resolve it, so install the last archived version directly before
# pulling destiny.
if (!requireNamespace("smoother", quietly = TRUE)) {
    remotes::install_version(
        "smoother", version = "1.3",
        repos = "https://cloud.r-project.org",
        upgrade = "never"
    )
}

BiocManager::install(
    c("ComplexHeatmap",   # Module 5 heatmap
      "slingshot",        # Module 6 principal-curve pseudotime
      "destiny",          # Module 6 diffusion pseudotime (DPT)
      "scater"),          # Module 6 SCE helpers / plotting
    update = FALSE, ask = FALSE
)

# 3. GitHub-only packages — source-only, will trigger compilation.
#    We use `pak` rather than `remotes::install_github` because newer
#    versions of remotes raise "can't convert package edgeR with
#    RemoteType 'bioconductor' to remote" when a Bioconductor-installed
#    package (like edgeR) is on the library path. pak handles
#    Bioconductor remotes correctly.
if (!requireNamespace("pak", quietly = TRUE)) {
    install.packages("pak")
}
pak::pkg_install(
    c("github::jiadongm/PhiSpace/pkg",
      "github::meiosis97/NeighbourNet"),
    ask = FALSE,
    upgrade = FALSE
)
```

### 3c. Verify every package loaded successfully

``` r

packages <- c(
    # Core (from renv.lock)
    "Seurat", "harmony", "glmGamPoi", "edgeR", "limma", "speckle",
    # Extras (Bioconductor + GitHub)
    "ComplexHeatmap", "slingshot", "destiny", "scater",
    "PhiSpace", "NeighbourNet"
)
check <- sapply(packages, requireNamespace, quietly = TRUE)
```

    ## Warning: replacing previous import 'S4Arrays::makeNindexFromArrayViewport' by
    ## 'DelayedArray::makeNindexFromArrayViewport' when loading 'SummarizedExperiment'

``` r

missing <- names(check)[!check]

if (length(missing) == 0) {
    message("All workshop packages installed successfully!")
} else {
    message("WARNING: Missing packages: ", paste(missing, collapse = ", "))
    message("Re-run the relevant chunk above (3a for the locked set, ",
            "3b for the extras), then verify again.")
}
```

    ## All workshop packages installed successfully!

### Key Packages Installed

| Package | Source | Purpose |
|----|----|----|
| Seurat | renv.lock | Core single cell analysis |
| harmony | renv.lock | Batch correction |
| glmGamPoi | renv.lock | Fast SCTransform (critical for reproducibility) |
| edgeR / limma | renv.lock | Differential expression |
| speckle | renv.lock | Composition analysis |
| ComplexHeatmap | Bioconductor | Module 5 heatmap |
| slingshot | Bioconductor | Module 6 principal-curve pseudotime |
| destiny | Bioconductor | Module 6 diffusion pseudotime (DPT) |
| scater | Bioconductor | Module 6 SCE helpers / plotting |
| PhiSpace | GitHub | Module 5 continuous phenotyping |
| NeighbourNet | GitHub | Module 7 cell-specific co-expression networks |

## Step 4: Download Workshop Data

The workshop data files are hosted on Zenodo. Run the following code to
download them:

``` r

# R's default download.file() timeout is 60 seconds, which is shorter than
# the time it takes to fetch a 420 MB file on most home connections. Bump
# it to 1 hour for the duration of the download.
options(timeout = 3600)

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

### Optional: Backup Checkpoints

If you want to skip a module — for example, to start Session 2 without
running Session 1, or to focus on a single technique — you can download
pre-computed checkpoints from a separate Zenodo record. The workshop is
designed so each module produces the input file for the next one, but
these backups let you jump in at any module boundary.

| File | When to use it | Replaces output of |
|----|----|----|
| `01_qc_filtered.rds` | Skip Module 1 | Module 1 |
| `02_integrated_clustered.rds` | Skip Modules 1 + 2 | Module 2 |
| `03_annotated.rds` | Skip Modules 1 + 2 + 3 | Module 3 |
| `afternoonSession.zip` | Skip Session 1 entirely and start at Module 5 | Modules 1–4 plus the Module 5/6/7 input and intermediate caches |

`afternoonSession.zip` contains the Session 2 (Module 5, 6, 7) input and
intermediate results. Download, unzip, and the files land in `data/` and
`results/` automatically (the chunk below handles the unzip).

``` r

# Replace with the backup record's URL once the new Zenodo DOI is
# assigned. Until then the variable below points at a placeholder.
backup_url <- "https://zenodo.org/records/<NEW_RECORD_ID>/files/"

# Pick whichever subset you need; comment out the others.
backup_files <- c(
    "01_qc_filtered.rds",            # → goes to data/processed/
    "02_integrated_clustered.rds",   # → goes to data/processed/
    "03_annotated.rds",              # → goes to data/processed/
    "afternoonSession.zip"           # → unzip at repo root, populates data/ and results/
)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

for (f in backup_files) {
    if (f == "afternoonSession.zip") {
        dest <- f
    } else {
        dest <- file.path("data", "processed", f)
    }
    if (file.exists(dest)) {
        message(f, " already exists, skipping")
        next
    }
    message("Downloading ", f, "...")
    download.file(
        url = paste0(backup_url, f, "?download=1"),
        destfile = dest,
        mode = "wb"
    )
}

# Unzip the afternoon-session archive at the repo root if you grabbed it.
# Files inside the archive are paths like "data/afternoon_mvp.rds" and
# "results/05_phispace_query.rds", so this extraction populates both
# folders directly.
if (file.exists("afternoonSession.zip")) {
    unzip("afternoonSession.zip")
    message("afternoonSession.zip unpacked")
}
```

After download, verify file integrity (digests are listed on the Zenodo
record’s main page and below):

``` r

# macOS / Linux — paste the expected SHA-256s from the Zenodo record
# into a small lookup and check each downloaded file:
expected <- c(
    "01_qc_filtered.rds"          = "43fbb351519f6972fb61e9ac45375a29bb7d78f41949ee41f86a81afc63c1ed5",
    "02_integrated_clustered.rds" = "1c69b2ef5475567b9412d8f6e08b2aa18d98dcfbbf4c0726c2f60b1ef587f000",
    "03_annotated.rds"            = "fcc7cfbd0f57c96c89d6cb745e82eef21789c22014015be431a4852ec43f0d51"
)
for (f in names(expected)) {
    path <- file.path("data", "processed", f)
    if (!file.exists(path)) next
    got <- tools::md5sum  # for SHA-256 use openssl::sha256(file(path)) or system shasum
    cat(f, ":", if (file.exists(path)) "downloaded" else "missing", "\n")
}
```

## Step 5: Final Verification

Run this final check to ensure everything is ready:

``` r

cat("=== Workshop Setup Verification ===\n\n")
```

    ## === Workshop Setup Verification ===

``` r

# Platform + R version
cat("Platform:  ", R.version$platform, "\n")
```

    ## Platform:   x86_64-pc-linux-gnu

``` r

cat("OS:        ", Sys.info()[["sysname"]], Sys.info()[["release"]], "\n")
```

    ## OS:         Linux 6.17.0-1010-azure

``` r

cat("R Version: ", R.version.string, "\n\n")
```

    ## R Version:  R version 4.5.2 (2025-10-31)

``` r

# All workshop packages
cat("Packages:\n")
```

    ## Packages:

``` r

packages <- c(
    "Seurat", "harmony", "glmGamPoi", "edgeR", "limma", "speckle",
    "ComplexHeatmap", "slingshot", "destiny", "scater",
    "PhiSpace", "NeighbourNet"
)
for (pkg in packages) {
    status <- if (requireNamespace(pkg, quietly = TRUE)) "OK" else "MISSING"
    cat(sprintf("  %-15s %s\n", pkg, status))
}
```

    ##   Seurat          OK
    ##   harmony         OK
    ##   glmGamPoi       OK
    ##   edgeR           OK
    ##   limma           OK
    ##   speckle         OK
    ##   ComplexHeatmap  OK
    ##   slingshot       OK
    ##   destiny         OK
    ##   scater          OK
    ##   PhiSpace        OK
    ##   NeighbourNet    OK

``` r

# Workshop data
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

### Compilation errors (“compilation failed for package …”)

Almost always means the system toolchain from Step 2 is missing or out
of date. Re-check it:

- **Windows** — open R and run
  `pkgbuild::has_build_tools(debug = TRUE)`. If it returns `FALSE`,
  reinstall
  [Rtools45](https://cran.r-project.org/bin/windows/Rtools/rtools45/rtools.html)
  and restart R / RStudio.
- **macOS** — run `xcode-select --install` in a terminal and click
  “Install” in the dialog.
- **Linux** — re-run the `apt-get install` block from Step 2; if a
  specific header file is named in the error (e.g. `hdf5.h`), install
  the matching `*-dev` package.

### Windows: “cannot rename file” / very long install paths

Windows still enforces a 260-character path limit on many APIs. If you
cloned the repo deep under `OneDrive` or a long username, move it to a
short path (e.g. `C:\workshop\single_cell_workshop`) and re-run
[`renv::restore()`](https://rstudio.github.io/renv/reference/restore.html).
Alternatively, enable [Long Paths in
Windows](https://learn.microsoft.com/en-us/windows/win32/fileio/maximum-file-path-limitation?tabs=registry).

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
    ## Running under: Ubuntu 24.04.4 LTS
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
    ##   [1] fs_1.6.6                    destiny_3.24.0             
    ##   [3] matrixStats_1.5.0           spatstat.sparse_3.1-0      
    ##   [5] httr_1.4.7                  RColorBrewer_1.1-3         
    ##   [7] doParallel_1.0.17           tools_4.5.2                
    ##   [9] mlr3learners_0.14.0         sctransform_0.4.3          
    ##  [11] backports_1.5.0             R6_2.6.1                   
    ##  [13] lazyeval_0.2.2              uwot_0.2.4                 
    ##  [15] GetoptLong_1.1.1            sp_2.2-0                   
    ##  [17] gridExtra_2.3               progressr_0.18.0           
    ##  [19] PhiSpace_1.1.0              cli_3.6.5                  
    ##  [21] Biobase_2.70.0              textshaping_1.0.4          
    ##  [23] spatstat.explore_3.6-0      fastDummies_1.7.5          
    ##  [25] TSP_1.2.7                   sass_0.4.10                
    ##  [27] Seurat_5.4.0                S7_0.2.1-1                 
    ##  [29] robustbase_0.99-7           spatstat.data_3.1-9        
    ##  [31] proxy_0.4-29                ggridges_0.5.7             
    ##  [33] pbapply_1.7-4               pkgdown_2.2.0              
    ##  [35] slingshot_2.18.0            mlr3tuning_1.6.0           
    ##  [37] systemfonts_1.3.1           paradox_1.0.1              
    ##  [39] harmony_1.2.4               scater_1.38.1              
    ##  [41] parallelly_1.46.1           limma_3.66.0               
    ##  [43] TTR_0.24.4                  generics_0.1.4             
    ##  [45] shape_1.4.6.1               ica_1.0-3                  
    ##  [47] spatstat.random_3.4-3       car_3.1-5                  
    ##  [49] dplyr_1.1.4                 Matrix_1.7-4               
    ##  [51] ggbeeswarm_0.7.3            S4Vectors_0.48.0           
    ##  [53] abind_1.4-8                 lifecycle_1.0.5            
    ##  [55] scatterplot3d_0.3-45        yaml_2.3.12                
    ##  [57] edgeR_4.8.2                 carData_3.0-6              
    ##  [59] SummarizedExperiment_1.40.0 SparseArray_1.10.8         
    ##  [61] Rtsne_0.17                  glmGamPoi_1.22.0           
    ##  [63] grid_4.5.2                  promises_1.5.0             
    ##  [65] crayon_1.5.3                miniUI_0.1.2               
    ##  [67] speckle_1.10.0              lattice_0.22-7             
    ##  [69] beachmat_2.26.0             cowplot_1.2.0              
    ##  [71] magick_2.9.1                pillar_1.11.1              
    ##  [73] knitr_1.51                  ComplexHeatmap_2.26.1      
    ##  [75] GenomicRanges_1.62.1        rjson_0.2.23               
    ##  [77] boot_1.3-32                 future.apply_1.20.1        
    ##  [79] codetools_0.2-20            glue_1.8.0                 
    ##  [81] spatstat.univar_3.1-5       pcaMethods_2.2.0           
    ##  [83] data.table_1.18.0           vcd_1.4-13                 
    ##  [85] vctrs_0.6.5                 png_0.1-8                  
    ##  [87] spam_2.11-3                 gtable_0.3.6               
    ##  [89] cachem_1.1.0                xfun_0.55                  
    ##  [91] princurve_2.1.6             S4Arrays_1.10.1            
    ##  [93] mime_0.13                   RcppEigen_0.3.4.0.2        
    ##  [95] Seqinfo_1.0.0               survival_3.8-3             
    ##  [97] seriation_1.5.8             SingleCellExperiment_1.32.0
    ##  [99] iterators_1.0.14            statmod_1.5.1              
    ## [101] fitdistrplus_1.2-4          ROCR_1.0-11                
    ## [103] nlme_3.1-168                xts_0.14.2                 
    ## [105] bbotk_1.10.0                RcppAnnoy_0.0.23           
    ## [107] mlr3pipelines_0.11.0        bslib_0.9.0                
    ## [109] mlr3_1.6.0                  irlba_2.3.5.1              
    ## [111] vipor_0.4.7                 KernSmooth_2.23-26         
    ## [113] otel_0.2.0                  colorspace_2.1-2           
    ## [115] BiocGenerics_0.56.0         nnet_7.3-20                
    ## [117] mlr3misc_0.21.0             smoother_1.3               
    ## [119] tidyselect_1.2.1            curl_7.0.0                 
    ## [121] compiler_4.5.2              BiocNeighbors_2.4.0        
    ## [123] lgr_0.5.2                   desc_1.4.3                 
    ## [125] DelayedArray_0.36.0         plotly_4.11.0              
    ## [127] checkmate_2.3.3             scales_1.4.0               
    ## [129] DEoptimR_1.1-4              lmtest_0.9-40              
    ## [131] hexbin_1.28.5               palmerpenguins_0.1.1       
    ## [133] SpatialExperiment_1.20.0    stringr_1.6.0              
    ## [135] digest_0.6.39               goftest_1.2-3              
    ## [137] spatstat.utils_3.2-1        rmarkdown_2.30             
    ## [139] ca_0.71.1                   XVector_0.50.0             
    ## [141] htmltools_0.5.9             pkgconfig_2.0.3            
    ## [143] MatrixGenerics_1.22.0       fastmap_1.2.0              
    ## [145] rlang_1.1.7                 GlobalOptions_0.1.4        
    ## [147] htmlwidgets_1.6.4           ggthemes_5.2.0             
    ## [149] shiny_1.12.1                farver_2.1.2               
    ## [151] jquerylib_0.1.4             zoo_1.8-15                 
    ## [153] jsonlite_2.0.0              BiocParallel_1.44.0        
    ## [155] BiocSingular_1.26.1         magrittr_2.0.4             
    ## [157] scuttle_1.20.0              Formula_1.2-5              
    ## [159] dotCall64_1.2               patchwork_1.3.2            
    ## [161] Rcpp_1.1.1-1                viridis_0.6.5              
    ## [163] reticulate_1.44.1           TrajectoryUtils_1.18.0     
    ## [165] stringi_1.8.7               MASS_7.3-65                
    ## [167] plyr_1.8.9                  parallel_4.5.2             
    ## [169] listenv_0.10.0              ggrepel_0.9.6              
    ## [171] deldir_2.0-4                splines_4.5.2              
    ## [173] tensor_1.5.1                circlize_0.4.18            
    ## [175] locfit_1.5-9.12             igraph_2.2.1               
    ## [177] uuid_1.2-2                  ranger_0.18.0              
    ## [179] vizOmics_0.1.0              spatstat.geom_3.6-1        
    ## [181] RcppHNSW_0.6.0              ScaledMatrix_1.18.0        
    ## [183] reshape2_1.4.5              stats4_4.5.2               
    ## [185] evaluate_1.0.5              SeuratObject_5.3.0         
    ## [187] NeighbourNet_1.0.0          renv_1.1.5                 
    ## [189] BiocManager_1.30.27         laeken_0.5.3               
    ## [191] foreach_1.5.2               httpuv_1.6.16              
    ## [193] VIM_7.0.0                   RANN_2.6.2                 
    ## [195] tidyr_1.3.2                 purrr_1.2.1                
    ## [197] polyclip_1.10-7             future_1.68.0              
    ## [199] clue_0.3-68                 scattermore_1.2            
    ## [201] ggplot2_4.0.1               rsvd_1.0.5                 
    ## [203] xtable_1.8-4                e1071_1.7-17               
    ## [205] RSpectra_0.16-2             later_1.4.5                
    ## [207] viridisLite_0.4.2           class_7.3-23               
    ## [209] ragg_1.5.0                  tibble_3.3.1               
    ## [211] registry_0.5-1              beeswarm_0.4.0             
    ## [213] IRanges_2.44.0              cluster_2.1.8.1            
    ## [215] ggplot.multistats_1.0.1     globals_0.18.0
