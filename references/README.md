# `references/` — reference data processing and internal docs

This folder holds materials that support development of the afternoon session
but are not themselves part of what participants run live.

## Contents

| File | Purpose |
|---|---|
| `PhiSpace.pdf` | Mao, Deng & Lê Cao 2025, *Genome Biology* — the Φ-Space paper |
| `NNet.pdf` | NeighbourNet paper |
| `PhiSpace_Guide_for_VibeCoding.md` | Internal API reference for PhiSpace |
| `prepare_kanemaru_2023.R` | Download + process the Kanemaru 2023 adult-heart snRNA-seq reference |
| `prepare_nicin_2022.R` | Download + process the Nicin 2022 foetal/paediatric heart snRNA-seq reference |

## PhiSpace reference strategy

Module 5 uses a **single pan-developmental reference**:

**Gao et al. 2026** (*Genome Biology*, GSE290367) — integrative snRNA-seq atlas
of the human left ventricle spanning fetal (8-18 wk) through old adult (≥60 y),
299 donors, 13 annotated cell types, with developmental stage, age, sex, and
disease metadata. We filter to non-diseased donors only.

This single reference covers the full maturation axis of the Sim 2021 query
(foetal → young → adult) without needing a multi-reference approach.

A Kanemaru 2023 adult reference is retained as a fallback
([prepare_kanemaru_2023.R](prepare_kanemaru_2023.R)).

## How to use the processing scripts

These scripts are **run once, offline**, by an instructor, not by workshop
participants. They produce a processed `SingleCellExperiment` object that is
uploaded to Zenodo alongside the query data, and loaded directly by the
Module 5 vignette.

```r
# From the repo root — requires ~64 GB RAM and ~40 GB free disk
source("references/prepare_gao_2026.R")
```

Output (not tracked in git):

```
data/processed/reference_gao2026.rds
```

## Open issues

- **Metadata column names**: exact colData column names for cell type, disease,
  stage, and tissue need to be confirmed after first loading the h5ad.
  TODOs are marked in `prepare_gao_2026.R`.
- **Subsampling**: need to inspect cell-type counts after filtering to
  non-diseased LV, then tune `proportion` and `minCellNum` in `subsample()`.
- **Memory**: the h5ad is 13.9 GB compressed; loading uncompressed will require
  a machine with ≥64 GB RAM. The HDF5Array backend (`use_hdf5 = TRUE`) may
  help defer materialisation.
