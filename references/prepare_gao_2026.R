# ---------------------------------------------------------------------------
# prepare_gao_2026.R
#
# One-off preprocessing script: downloads the Gao et al. 2026 pan-developmental
# human left ventricle snRNA-seq atlas from GEO, filters to non-diseased
# donors, and saves a subsampled SingleCellExperiment for use as the PhiSpace
# reference in Module 5.
#
# Gao W. et al. (2026) "An integrative single-nucleus multiomic atlas of the
# human left ventricle identifies gene regulatory network dynamics across
# cardiac development, aging, and disease", Genome Biology.
# DOI: 10.1186/s13059-026-04061-7
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE290367
#
# Data: GSE290367_integrated_RNA_adata.h5ad.gz (13.9 GB compressed)
#   - adata.X          : log-normalised counts
#   - adata.layers['counts'] : raw counts
#   - 90 samples spanning fetal (8-18 wk) through adult (22-91 y)
#   - 13 annotated cell types + developmental stage, age, sex, disease status
#
# This script is run once, offline, by an instructor. Participants load the
# processed .rds directly from Zenodo during the workshop.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(zellkonverter)   # readH5AD()
    library(PhiSpace)        # subsample()
    library(scran)           # computeSumFactors(), for scranTransf
})

# --- Paths -----------------------------------------------------------------

dest_dir <- file.path("data", "processed")
raw_h5ad <- file.path(dest_dir, "gao_2026_raw.h5ad.gz")
out_rds  <- file.path(dest_dir, "reference_gao2026.rds")

dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

# --- 1. Download -----------------------------------------------------------
# 13.9 GB compressed; uncompressed will be considerably larger.
# Requires ~40 GB free disk space and ~32 GB RAM to load in full.

geo_url <- paste0(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE290nnn/GSE290367/suppl/",
    "GSE290367_integrated_RNA_adata.h5ad.gz"
)

if (!file.exists(raw_h5ad)) {
    message("Downloading GSE290367 RNA h5ad (~13.9 GB compressed) ...")
    options(timeout = 7200)
    download.file(geo_url, raw_h5ad, mode = "wb")
}

# --- 2. Load ---------------------------------------------------------------
# Use HDF5Array backend to avoid loading everything into RAM at once.
# TODO: test with use_hdf5 = TRUE; if downstream steps require in-memory,
# switch to use_hdf5 = FALSE on a machine with sufficient RAM (>=64 GB).

message("Loading h5ad (this may take several minutes) ...")
sce <- readH5AD(raw_h5ad, use_hdf5 = TRUE)

cat("Loaded SCE:\n"); print(sce)
cat("\ncolData columns:\n"); print(colnames(colData(sce)))
cat("\nAssay names:\n"); print(assayNames(sce))

# --- 3. Inspect metadata ---------------------------------------------------
# TODO: confirm exact column names after loading; expected:
#   - cell_type / celltype / leiden  : cell type labels
#   - development_stage / stage      : fetal / adult
#   - age / age_group                : numeric or binned age
#   - sex / gender
#   - disease / condition            : non-diseased, DCM, HCM, AMI, etc.
#   - donor_id / sample

# Quick look
cat("\nCell type distribution:\n")
# print(table(colData(sce)$cell_type))   # TODO: use correct column name
cat("\nDisease distribution:\n")
# print(table(colData(sce)$disease))

# --- 4. Filter: non-diseased donors only -----------------------------------
# TODO: confirm 'disease' column name and the value for non-diseased
# (likely "normal", "healthy", "non-diseased", or "control").

# is_healthy <- colData(sce)$disease %in% c("normal", "non-diseased", "control")
# cat("Non-diseased nuclei:", sum(is_healthy), "of", ncol(sce), "\n")
# sce <- sce[, is_healthy]

# --- 5. Filter: LV only (if multi-region) ---------------------------------
# Sim 2021 query is apical LV snRNA-seq. Restricting to LV keeps
# reference phenotypes directly comparable; include all regions if LV-only
# gives too few cells for rare types.
# TODO: confirm tissue/region column and values.

# is_lv <- colData(sce)$tissue %in% c("left ventricle", "LV")
# sce <- sce[, is_lv]

# --- 6. Subsample for balance and tractable PhiSpace training --------------
# Full non-diseased LV could be >500k nuclei. We subsample to a balanced
# set: enough to train PLS robustly (~50-200 cells per cell type), but
# manageable for a laptop-based workshop.
# TODO: inspect cell type counts after filtering, then tune proportion/minCellNum.

# sce <- subsample(
#     sce,
#     key        = "cell_type",   # TODO: use confirmed column name
#     proportion = 0.05,
#     minCellNum = 100,
#     seed       = 5202056
# )
# cat("After subsampling:", ncol(sce), "nuclei\n")

# --- 7. Assay setup --------------------------------------------------------
# PhiSpace needs a normalised assay. The h5ad X slot is already log-normalised
# so we rename it "logcounts" and use that as refAssay in PhiSpace().
# Raw counts are in layers[['counts']] if scran normalisation is preferred.

# assayNames(sce)[assayNames(sce) == "X"] <- "logcounts"
# TODO: verify this rename is needed or if zellkonverter already names it correctly.

# --- 8. Save ---------------------------------------------------------------

# saveRDS(sce, out_rds)
# message("Saved processed Gao 2026 reference to: ", out_rds)
# message("Nuclei: ", ncol(sce), " | Cell types: ",
#         length(unique(colData(sce)$cell_type)))
