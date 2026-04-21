# ---------------------------------------------------------------------------
# prepare_kanemaru_2023.R
#
# One-off preprocessing script: downloads the Kanemaru et al. 2023 (Nature)
# adult human heart snRNA-seq reference from CELLxGENE, filters to the
# snRNA-seq subset, and saves a SingleCellExperiment object for use as a
# PhiSpace reference in Module 5.
#
# Kanemaru K. et al. (2023) "Spatially resolved multiomics of human cardiac
# niches", Nature. DOI: 10.1038/s41586-023-06311-1
#
# This script is run once, offline, by an instructor. Participants load the
# processed .rds directly from Zenodo during the workshop.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(zellkonverter)   # readH5AD()
    library(PhiSpace)        # subsample()
})

# --- Paths -----------------------------------------------------------------

dest_dir  <- file.path("data", "processed")
raw_h5ad  <- file.path(dest_dir, "kanemaru_2023_raw.h5ad")
out_rds   <- file.path(dest_dir, "reference_kanemaru.rds")

dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

# --- 1. Download -----------------------------------------------------------
# CELLxGENE-hosted h5ad (~5 GB). This URL points to the global multiomic
# object; the snRNA-seq subset is filtered out below.

cellxgene_url <- "https://datasets.cellxgene.cziscience.com/1a7a9fb0-aee1-437f-8a7c-9d132253a4db.h5ad"

if (!file.exists(raw_h5ad)) {
    message("Downloading Kanemaru 2023 h5ad (~5 GB) from CELLxGENE ...")
    options(timeout = 3600)
    download.file(cellxgene_url, raw_h5ad, mode = "wb")
}

# --- 2. Load as SCE --------------------------------------------------------

sce <- readH5AD(raw_h5ad, use_hdf5 = TRUE)
# TODO: verify assay name (likely "counts" + "X" as logcounts); rename if needed
# TODO: verify reducedDims transferred; drop unused ones

cat("Loaded Kanemaru SCE:\n")
print(sce)
cat("colData columns:\n"); print(colnames(colData(sce)))

# --- 3. Filter to snRNA-seq only ------------------------------------------
# The global object mixes scRNA-seq and snRNA-seq. Filter to nuclei only.
# TODO: confirm the correct column after loading — candidates:
#   - colData(sce)$cell_source   (likely "nuclei" vs "cells")
#   - colData(sce)$assay         (e.g. "10x 3' v3" — same for sc and sn)
# After the first run, replace the placeholder below with the confirmed column.

# is_nuclei <- colData(sce)$cell_source == "nuclei"
# sce <- sce[, is_nuclei]
# cat("After snRNA-seq filter: ", ncol(sce), " nuclei\n")

# --- 4. Harmonise metadata -------------------------------------------------
# TODO: standardise column names the query will join on, e.g.:
#   - cell_type     : main cell type annotation (keep Kanemaru's curated labels)
#   - region        : anatomical region (LV/RV/LA/RA/SAN/AVN/...)
#   - donor         : donor id
#   - sex, age      : donor metadata
#   - stage         : add "adult" to harmonise with Nicin labels

# colData(sce)$stage <- "adult"

# --- 5. Subsample for balance and speed ------------------------------------
# The full nuclei subset is hundreds of thousands of cells. PhiSpace trains
# faster on a balanced subsample; we keep enough per cell type to retain
# the phenotypic richness.
# TODO: pick proportion / minCellNum after inspecting cell-type counts.

# sce <- subsample(
#     sce,
#     key        = "cell_type",
#     proportion = 0.2,
#     minCellNum = 100,
#     seed       = 5202056
# )

# --- 6. Save ---------------------------------------------------------------

saveRDS(sce, out_rds)
message("Saved processed Kanemaru reference to: ", out_rds)
