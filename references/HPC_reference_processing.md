# HPC Instructions: Processing the Gao 2026 Reference Dataset

Reference data processing for Module 5 (Φ-Space annotation). Run once,
offline, on HPC. The output `.rds` is then uploaded to Zenodo for participants
to download.

**Script:** `references/prepare_gao_2026.R`
**Dataset:** Gao et al. 2026, *Genome Biology*, GSE290367
**File to download:** `GSE290367_integrated_RNA_adata.h5ad.gz` (13.9 GB compressed)

---

## 1. System Requirements

| Resource | Minimum | Recommended |
|---|---|---|
| RAM | 64 GB | 128 GB |
| Disk (scratch) | 60 GB free | 100 GB free |
| R version | 4.3+ | 4.5+ |
| Time (walltime) | 4 h | 6 h |

---

## 2. Getting the Repository on HPC

```bash
git clone https://github.com/phipsonlab/single_cell_workshop.git
cd single_cell_workshop
```

---

## 3. Install Required R Packages

These packages are **not** in the workshop `renv.lock` (which covers only
what participants need). Install them separately into your HPC R library.

```r
# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "SingleCellExperiment",
    "zellkonverter",       # readH5AD()
    "HDF5Array",           # delayed arrays for large h5ad
    "scran",               # scranTransf
    "scuttle"              # logNormCounts
))

# PhiSpace (subsample function)
BiocManager::install("jiadongm/PhiSpace/pkg")

# Utilities
install.packages(c("rhdf5", "Matrix"))
```

Verify:

```r
pkgs <- c("SingleCellExperiment", "zellkonverter", "HDF5Array",
          "scran", "PhiSpace")
ok <- sapply(pkgs, requireNamespace, quietly = TRUE)
stopifnot(all(ok))
```

---

## 4. SLURM Job Script

Save as `references/run_gao_processing.slurm`, then submit:

```bash
sbatch references/run_gao_processing.slurm
```

```bash
#!/bin/bash
#SBATCH --job-name=gao2026_ref
#SBATCH --output=references/gao2026_%j.log
#SBATCH --error=references/gao2026_%j.err
#SBATCH --time=06:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --partition=regular        # TODO: change to your cluster's partition name

# Load R module (adjust version to what is available)
module load R/4.5.2                # TODO: check with `module avail R`

# Run from repo root
cd $SLURM_SUBMIT_DIR
Rscript references/prepare_gao_2026.R
```

Check available R modules with:

```bash
module avail R
```

---

## 5. Two-Phase Processing

Because the h5ad contains ~2.27M nuclei and is 13.9 GB compressed, it is
best to run the script in two phases:

### Phase 1 — Inspect metadata only (interactive or short job)

Run this interactively (or as a short 1 h, 64 GB job) to confirm column names
**before** enabling the filtering and subsampling TODOs.

```r
library(zellkonverter)
library(HDF5Array)

# Load with HDF5 backend — delays loading the full expression matrix
sce <- readH5AD("data/processed/gao_2026_raw.h5ad.gz", use_hdf5 = TRUE)

# Inspect
dim(sce)                              # genes × nuclei
colnames(colData(sce))                # metadata column names
assayNames(sce)                       # assay names

# Key distributions to check
table(colData(sce)$cell_type)         # TODO: replace with actual column name
table(colData(sce)$disease)           # TODO: values for non-diseased
table(colData(sce)$development_stage) # TODO: fetal / adult breakdown
table(colData(sce)$tissue)            # TODO: LV / RV / LA / RA etc.

# Check assay names — we expect something like "X" (logcounts) or "counts"
print(assayNames(sce))
```

**Fill in the TODOs in `prepare_gao_2026.R`** based on what you see:

| TODO in script | What to fill in |
|---|---|
| `cell_type` column | Exact column name for cell-type labels |
| `disease` column + "healthy" value | Column name + exact string(s) for non-diseased |
| `tissue` column + "LV" value | Column name + exact string for left ventricle |
| Assay rename | Whether `assayNames(sce)` returns `"X"` (needs rename) or `"logcounts"` already |
| `proportion` and `minCellNum` in `subsample()` | Based on cell-type count table |

### Phase 2 — Full processing (submit as SLURM job)

Once TODOs are filled, uncomment the filter + subsample + save blocks in
`prepare_gao_2026.R` and submit the full job:

```bash
sbatch references/run_gao_processing.slurm
```

Expected output file: `data/processed/reference_gao2026.rds`

---

## 6. Subsampling Guidance

PhiSpace trains PLS on the reference. A well-balanced subsample is better
than the full dataset. Target:

- ≥ 200 cells per cell type (for rare types, set `minCellNum = 100`)
- ≤ ~10,000 total nuclei for fast in-workshop fitting
- Stratify by `cell_type`; optionally by `development_stage` within `cell_type`
  if you want equal representation of fetal and adult cardiomyocytes

```r
# Suggested starting point — tune based on Phase 1 cell-type counts
sce_ref <- PhiSpace::subsample(
    sce,
    key        = "cell_type",   # use confirmed column name
    proportion = 0.005,         # ~0.5% of 2.27M = ~11k; adjust
    minCellNum = 100,
    seed       = 5202056
)
cat("Reference nuclei after subsample:", ncol(sce_ref), "\n")
print(table(sce_ref$cell_type))
```

After subsampling, inspect the stage breakdown:

```r
table(sce_ref$development_stage, sce_ref$cell_type)
```

If the fetal:adult ratio for cardiomyocytes is very uneven, run a
two-level subsample (stage within cell type). Adjust as needed.

---

## 7. Verifying the Output

Before saving, run a quick sanity check:

```r
# Dimensions
cat("Reference dimensions:", nrow(sce_ref), "genes x", ncol(sce_ref), "nuclei\n")

# Cell types
cat("\nCell type counts:\n")
print(sort(table(sce_ref$cell_type), decreasing = TRUE))

# Stage distribution
cat("\nDevelopmental stage × cell type:\n")
print(table(sce_ref$development_stage, sce_ref$cell_type))

# Assay
cat("\nAssay names:\n"); print(assayNames(sce_ref))
stopifnot("logcounts" %in% assayNames(sce_ref))
```

Then save:

```r
saveRDS(sce_ref, "data/processed/reference_gao2026.rds")
```

---

## 8. Transfer and Zenodo Upload

Once the `.rds` is generated, transfer it back to your local machine:

```bash
# From local machine
scp user@hpc.cluster.edu:/path/to/single_cell_workshop/data/processed/reference_gao2026.rds \
    data/processed/reference_gao2026.rds
```

Then upload to the existing Zenodo record (record `18237749`) or a new
record alongside the query data files (`heart-counts.Rds`,
`cellinfo_updated.Rds`). Update `R/download_data.R` to include
`reference_gao2026.rds` in the `files` vector.

---

## 9. Troubleshooting

**Out of memory when loading h5ad**
: Always use `readH5AD(..., use_hdf5 = TRUE)` to keep arrays on disk. If
you still hit limits, try loading only the `obs` (metadata) layer first via
`rhdf5::h5read("file.h5ad", "obs")` to inspect columns without touching the
expression matrix.

**`gzip` / decompression error**
: GEO serves `.h5ad.gz`; `readH5AD` handles gzipped h5ad automatically in
recent `zellkonverter` versions. If it fails, decompress manually first:
`R.utils::gunzip("gao_2026_raw.h5ad.gz")` then re-read.

**Slow download on HPC**
: GEO FTP is sometimes throttled. Alternative: download to local machine
first, then `scp` to HPC scratch. Or use `wget` with resumption:
`wget -c <url> -O data/processed/gao_2026_raw.h5ad.gz`.

**Module not found: `R/4.5.2`**
: Check available versions with `module avail R` and update the `module load`
line in the SLURM script accordingly.

**`PhiSpace::subsample` not found**
: Reinstall PhiSpace: `BiocManager::install("jiadongm/PhiSpace/pkg")`.
