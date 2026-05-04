# `references/` — reference data processing and internal docs

This folder holds materials that support development of the afternoon session
but are not themselves part of what participants run live.

## Contents

| File | Purpose |
|---|---|
| `PhiSpace.pdf` | Mao, Deng & Lê Cao 2025, *Genome Biology* — the Φ-Space paper |
| `NNet.pdf` | Deng, Mao, Choi & Lê Cao 2026, *Genome Research* — the NeighbourNet (NNet) paper, "Scalable cell-specific coexpression networks for granular regulatory pattern discovery with NeighbourNet" |
| `PhiSpace_Guide_for_VibeCoding.md` | Internal API reference for PhiSpace |
| `NeighbourNet_Guide_for_VibeCoding.md` | Internal API reference for NeighbourNet (NNet) |
| `prepare_kanemaru_2023.R` | Download + process the Kanemaru 2023 adult-heart snRNA-seq reference (fallback only) |
| `HPC_reference_processing.md` | HPC recipe for processing the Gao 2026 reference (the script itself is to be written — see "Open issues" below) |

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

The Gao 2026 processing script (`prepare_gao_2026.R`) is not yet committed —
follow the recipe in [`HPC_reference_processing.md`](HPC_reference_processing.md)
to assemble it on an HPC node with ~64 GB RAM and ~40 GB free disk. The
expected output (not tracked in git) is:

```
data/processed/reference_gao2026.rds
```

## Open issues

- **Script not yet committed**: `prepare_gao_2026.R` itself is yet to be
  written; the HPC instructions describe the intended workflow but are not a
  runnable file.
- **Metadata column names**: exact colData column names for cell type, disease,
  stage, and tissue need to be confirmed after first loading the h5ad. TODOs
  to add to `prepare_gao_2026.R` are listed in `HPC_reference_processing.md`.
- **Subsampling**: need to inspect cell-type counts after filtering to
  non-diseased LV, then tune `proportion` and `minCellNum` in `subsample()`.
- **Memory**: the h5ad is 13.9 GB compressed; loading uncompressed will require
  a machine with ≥64 GB RAM. The HDF5Array backend (`use_hdf5 = TRUE`) may
  help defer materialisation.
