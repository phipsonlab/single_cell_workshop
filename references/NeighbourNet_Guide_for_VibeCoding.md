# NeighbourNet (NNet) — Guide for Vibe-Coding

Internal reference for working with the `NeighbourNet` R package
([meiosis97/NeighbourNet](https://github.com/meiosis97/NeighbourNet),
docs: <https://meiosis97.github.io/NeighbourNet/>).

This is **not** the phylogenetic NeighbourNet algorithm by Bryant &
Moulton (2004). It is a single-cell method for inferring
**cell-specific co-expression networks (CSNs)** and aggregated
**meta-networks**, from the Lê Cao lab (the same group that produced
PhiSpace).

**Citation**: Deng, Mao, Choi & Lê Cao (2026), "Scalable cell-specific
coexpression networks for granular regulatory pattern discovery with
NeighbourNet", *Genome Research* 36:785–801,
[10.1101/gr.281171.125](https://doi.org/10.1101/gr.281171.125). Yidi
Deng and Jiadong Mao are co-first authors.

The method, in one diagram (Figure 1 of the paper):

1. PCA on the Seurat expression matrix → low-rank approximation
   (LRA).
2. kNN graph in PC space defines each cell's neighbourhood.
3. Per-cell PC regression: a *predictor* gene `q` is embedded in the
   PCs of cell `n`'s kNN; the regression model predicts the
   *response* gene `p`. Co-expression `CSN[n,p,q]` = how strongly `q`
   contributes through the PCs to predicting `p`.
4. Repeat for every (response, predictor, cell) triple → a 3-D
   tensor (the "cube").
5. **nPCA** on the cube along the cell axis → soft cell clusters and
   their **meta-networks** (the typical co-expression patterns each
   cell-cluster shows). The same trick along the gene axis yields
   meta-genes. (The paper text says "NMF"; the package implementation
   is nPCA via `npca()` — non-negative PCA with deflation, which gives
   an ordered basis like PCA but with non-negativity constraints.)
6. Optional: integrate prior knowledge from OmniPath / NicheNet to
   contextualise edges as gene-regulation evidence and to infer
   upstream signaling pathway activity per cell.

---

## 1. What NNet does (and doesn't)

**Does:**
- Builds a per-cell network of **TF → target gene co-expression**, using
  PCs of nearest neighbours as predictors and permutation feature
  importance as the co-expression signal.
- Aggregates the per-cell networks via **non-negative PCA (nPCA)**
  with deflation, exported as `npca()`, into a small set of
  **meta-networks** that expose recurring co-expression patterns.
  (The paper text talks about "NMF" loosely; the implementation is
  nPCA. nPCA differs from classical NMF in that components are
  extracted sequentially with a deflation step, giving an ordered
  basis where meta-network 1 captures the strongest pattern, 2 the
  next, and so on — analogous to PCA component ordering, but with
  non-negativity.)
- Optionally maps cell-cell-network variation to receptor activities
  using a built-in receptor-TF signaling prior (PageRank-derived).

**Doesn't:**
- Infer a *gene regulatory network* (GRN) in the directed/causal sense.
  Edges are co-expression associations. Treat the output as a
  **scaffold** that may *contain* GRN structure rather than as a GRN
  itself. (Cf. Saint-Antoine & Singh, *Nat Rev Genet* 2026, on the
  GRN-vs-co-expression distinction.)
- Provide directionality or kinetics; no edge "from TF to target" is
  evidence of regulation in the mechanistic sense.

The framing the package authors use: *individual cell-specific networks
are noisy, but aggregating across many of them via meta-networks
denoises the signal.*

---

## 2. Installation

CRAN-hosted dependencies plus a GitHub install:

```r
remotes::install_github("meiosis97/NeighbourNet")
```

Hard dependencies (loaded by the package): `Seurat`, `ggplot2`,
`ggraph`, `igraph`, `Matrix`. Optional/extension: `patchwork` for
combining plots.

---

## 3. Input format

NNet operates on a **Seurat object** with normalised data
(`NormalizeData()` already run; the `data` layer is what NNet uses
during scaling).

If your data is in a `SingleCellExperiment` you must convert before
using NNet:

```r
library(Seurat)
seu <- as.Seurat(sce, counts = "counts", data = "logcounts")
```

Note: `as.Seurat()` requires `counts` and a normalised assay; the
`logcounts` from scran is acceptable as the `data` layer (NNet expects
log-normalised, not raw).

---

## 4. Built-in prior knowledge: `gene.list`, `gr.graph`, `sig.graph`, `receptor.ppr`

The package ships four data objects compiled from integrated prior
knowledge:

| Object | Purpose |
|---|---|
| `gene.list` | Named list with `tfs`, `targets`, `receptors`, `ligands` — character vectors of human gene symbols |
| `gr.graph` | Weighted directed igraph: TF → target gene regulation prior |
| `sig.graph` | Weighted directed igraph: receptor → TF signaling prior |
| `receptor.ppr` | Personalised PageRank scores propagated through `sig.graph` |

**Species**: human gene symbols throughout. For mouse data, you must
supply your own TF/target lists via the `tfs`/`targets` arguments to
`select.gene()` (case-mapping to mouse symbols is up to you).

---

## 5. Canonical workflow

The full pipeline runs on a Seurat object and stores everything back
in `obj@misc`. Order matters; each step depends on the previous.

```r
library(Seurat)
library(NeighbourNet)

# 1) Filter genes and pick TF / target lists
genes <- select.gene(obj, min.cells = 10)            # uses built-in gene.list

# 2) Scale + PCA on the filtered gene set
obj <- prepare.seurat(obj, genes = genes$genes, npcs = 100)

# 3) kNN graph in PC space
obj <- prepare.graph(obj, knn = 30)

# 4) (Optional) balanced cell sub-sample for compute reasons
obj <- select.cell(obj, p = 0.1)

# 5) Set up regression scaffolding (responses + predictors)
obj <- prepare.reg(
    obj,
    predictors = genes$tfs,
    responses  = genes$targets
)

# 6) Run the per-cell PC regression
obj <- run.nn.reg(obj, responses = genes$targets, return.p.val = TRUE)

# 7) Aggregate cell-specific networks into meta-networks
obj <- build.meta.network(obj, n.net = 20)

# 8) Pick hub genes for visualisation
ctr.genes <- select.central.genes(obj, n.per.component = 4, k = 2)

# 9) Prepare visualisation context
obj <- prepare.visualise(obj, central.genes = ctr.genes, n.clu = 4)

# 10) Plot per-cell or meta-network
visualise.network(obj, i = 1,
                  radius     = c(0.4, 0.7, 0.85, 1),
                  pie.radius = 0.04,
                  text.size  = 5)
```

---

## 6. Function reference

### 6.1 `select.gene()` — gene QC + TF/target list resolution

```r
select.gene(seurat.obj, tfs = NULL, targets = NULL, bgs = NULL,
            min.cells = 20)
```

- `min.cells` — drop genes detected in fewer than this many cells.
- `tfs`, `targets`, `bgs` — character vectors. If `NULL`, defaults to
  `gene.list$tfs`, `gene.list$targets`, and `gene.list$ligands`
  (background) respectively.

Returns a list with `tfs`, `targets`, `bgs`, and a unified `genes`
vector. Use `genes$genes` as the input to `prepare.seurat()`.

### 6.2 `prepare.seurat()` — scale + PCA

```r
prepare.seurat(seurat.obj, genes,
               npcs = 100, truncated = TRUE,
               ScaleData.ctrl = list(), RunPCA.ctrl = list())
```

- `genes` — character vector to restrict scaling/PCA to.
- `npcs` — max PCs to compute.
- `truncated = TRUE` — keeps only PCs above a noise threshold (uses
  `find.significant.pcs()` internally).
- `ScaleData.ctrl`, `RunPCA.ctrl` — pass-through args to Seurat.

Stores `pca` reduction back on the Seurat object plus an
`NNet.setting` list in `@misc`.

### 6.3 `prepare.graph()` — kNN graph in PC space

```r
prepare.graph(seurat.obj, knn = 30)
```

Builds an undirected kNN graph (default k = 30) on the PCs from step
6.2. Stores three pieces under `obj@misc$NNet.setting`:

- `p` — sparse cell × cell affinity matrix.
- `nn.idx` — cell × k matrix of neighbour indices.
- `nn.w` — cell × k matrix of edge weights.

### 6.4 `select.cell()` — balanced sub-sampling (optional)

```r
select.cell(seurat.obj, p = 0.1, n = NULL, all = FALSE, ...)
```

- `p` — fraction of cells to keep (default 10%).
- `n` — explicit count (overrides `p`).
- `all = TRUE` — keep all cells.

Uses k-means in PC space + nearest-cell-to-centroid selection rather
than uniform random sampling. Skip this step on small datasets
(< ~5k cells) — the per-cell regression runs fast enough.

### 6.5 `prepare.reg()` — regression scaffolding

```r
prepare.reg(seurat.obj,
            responses = NULL, predictors = NULL,
            cells = NULL, check.expressed = FALSE)
```

Pre-computes local variances and low-rank approximations needed by
`run.nn.reg`. Just records gene lists; no per-cell regression yet.

### 6.6 `run.nn.reg()` — the regression itself

```r
run.nn.reg(seurat.obj,
           responses = NULL, Y = NULL, predictors = NULL,
           t = 3, k = NULL,
           remove.self.loops = TRUE,
           f = function(x) 2 * x^2,
           assay = c("effect", "p.val"),
           prune = TRUE, cutoff = 0.95,
           return.p.val = FALSE,
           return.smooth = TRUE, return.prune = FALSE)
```

The core computation. For every cell, regresses each response gene's
expression on PCs of its neighbourhood, then computes
**permutation feature importance** for each predictor (TF) as the
co-expression signal.

Stores everything under `obj@misc$NNet.mod`:

- `effect` — `responses × predictors × cells` tensor of co-expression
  signals.
- `p.val` — same shape, only if `return.p.val = TRUE`.
- Plus smoothing/pruning state, noise distributions, and the
  Laplacian operator pieces.

Compute scales roughly with `cells × responses × predictors`. On the
default Nestorowa example (~1.6k cells, ~10 targets, ~hundreds of
TFs) the run takes seconds.

**Key knobs:**
- `t` — Laplacian power; controls how aggressively neighbourhood
  smoothing is applied. Default 3.
- `return.p.val = TRUE` — enables the meta-network step's significance
  filtering; almost always wanted.

### 6.7 `build.meta.network()` — aggregate to meta-networks

```r
build.meta.network(seurat.obj = NULL, network = NULL,
                   k = 100, cutoff = NULL, big.memory = FALSE,
                   scale = TRUE, truncated = TRUE,
                   n.net = 20,
                   non.neg = TRUE,
                   max.iter = 1000, tol = 1e-10,
                   return.p.val = TRUE)
```

Vectorises each cell's `responses × predictors` adjacency, runs PCA on
the `cells × edges` matrix, then non-negative PCA (`nPCA`) for an
interpretable basis. Adds under `NNet.mod$meta.network`:

- `meta.network` — `responses × predictors × n.net` tensor.
- `p.val` — significance per meta-network.
- `pcs` — `cells × n.net` cell embedding (each row is a cell, each
  column its weight on a meta-network).
- `pca.loadings`, `npca.loadings`, `pca.sd` — embedding internals.

Sister function `build.meta.response()` does the same trick on
`cell × predictor` profiles per response gene rather than full
networks; useful when you want to find groups of *responses* that
share predictors rather than groups of *cell states* that share
network structure.

### 6.8 `get.network()` — extract sub-tensors

```r
get.network(seurat.obj, i = NULL,
            assay = NULL, remove.self.loops = NULL,
            responses = NULL, predictors = NULL,
            f = NULL, drop = TRUE, cutoff = NULL)
```

Convenience getter. Returns:

- single cell `i` → `responses × predictors` matrix.
- single response → `predictors × cells` matrix.
- single predictor → `responses × cells` matrix.

`assay = "effect"` returns the regression-importance tensor (default
applies `f(x) = 2x^2`); `assay = "p.val"` returns p-values.

### 6.9 `select.central.genes()` — pick hubs

```r
select.central.genes(seurat.obj = NULL, network = NULL,
                     n.net = NULL,
                     k = 1, n.per.component = 4,
                     keep.responses = FALSE)
```

Per meta-network, runs SVD on the `responses × predictors` adjacency
and picks the top genes by absolute loading on the leading singular
vectors — i.e. eigenvector centrality on the bipartite graph. Returns
a character vector. Defaults to 4 genes × 1 vector per meta-network;
bumping `k` and `n.per.component` widens the set.

### 6.10 `prepare.visualise()` — set up the plot context

```r
prepare.visualise(seurat.obj,
                  n.clu = 4, central.genes = NULL,
                  check.gr.evidence = TRUE, t = 2, p = NULL,
                  as.g2 = c("predictors", "responses"),
                  g1 = NULL, g2 = NULL, receptors = NULL)
```

Picks the gene clusters that will form the concentric layers in
`visualise.network()`. With `central.genes` from step 6.9 and
`n.clu = 4`, you get 4 clusters of related genes around 4 hub genes.
Stores `NNet.visual.setting` in `@misc`.

### 6.11 `visualise.network()` — the figure

```r
visualise.network(seurat.obj, i,
                  meta.network = FALSE,
                  fix.cluster = TRUE, hubs = NULL, n.clu = 4,
                  cutoff = NULL, show.pathways = TRUE,
                  change.receptors = TRUE,
                  receptor.activity = c("cprod", "dist"),
                  check.receptor.expression = TRUE,
                  scale.ppr = FALSE, scale.network = FALSE,
                  scale.signifiance = FALSE,
                  swap.layers = FALSE,
                  k = 2, n.per.component = 10,
                  radius = NULL, pie.radius = 0.05, text.size = 4)
```

Returns a `ggplot` object (built with `ggraph` and `ggplot2`).
Concentric layers: outer = target gene clusters, middle = TFs, inner
= receptors. Edge thickness = significance, pie charts = node
properties.

- `i` — cell index (when `meta.network = FALSE`) or meta-network
  component (when `TRUE`).
- `radius = c(0.4, 0.7, 0.85, 1)` is a sensible default for a
  3-layer plot with a tight inner ring.
- `pie.radius = 0.04` keeps node pies legible without overlap.
- `cutoff` — p-value threshold for edge visibility. The
  vignette example computes a per-meta-network cutoff via
  `mean(apply(meta.p, 1, max))`.

### 6.12 `receptor.activity()` — receptors → TF → target propagation

```r
receptor.activity(seurat.obj, i = NULL, meta.network = FALSE,
                  cutoff = NULL,
                  check.receptor.expression = TRUE,
                  check.gr.evidence = TRUE,
                  t = 2, p = NULL,
                  scale.ppr = TRUE, scale.network = TRUE,
                  as.tfs = c("predictors", "responses"),
                  receptors = NULL, tfs = NULL, targets = NULL,
                  receptor.activity = c("cprod", "dist"))
```

Combines `receptor.ppr` (receptor → TF prior) with the cell-specific
TF → target networks to estimate **receptor activity per cell**.
Returns a list with `receptor.act`, `tf.act`, `target.act` matrices
when called over multiple cells; a smaller list with the activity
matrix and the propagation pieces for a single cell.

This is the route from co-expression networks back to upstream
**signaling pathway hypotheses**: high `NOTCH2` activity in a cell
means the propagation signal flowing through the prior signaling
graph and the cell's specific TF→target network coincides on the
NOTCH2 receptor node.

---

## 7. Output structures: where everything lives

After the canonical workflow, `obj@misc` carries:

```
obj@misc$NNet.setting        # gene lists, kNN graph, local variances
obj@misc$NNet.mod
   ├── effect                 # responses × predictors × cells
   ├── p.val                  # same shape, optional
   └── meta.network
          ├── meta.network    # responses × predictors × n.net
          ├── p.val           # same shape
          ├── pcs             # cells × n.net
          ├── pca.loadings    # for the SVD step
          ├── npca.loadings   # for the nPCA step
          └── pca.sd
obj@misc$NNet.visual.setting  # g1, g2, hubs, evidence, ppr
```

Use `Seurat::Misc(obj, "NNet.mod") %>% names()` to inspect.

---

## 8. Quickstart in 25 lines

(From the package vignette, lightly compressed.)

```r
library(Seurat); library(NeighbourNet); library(ggplot2)

obj <- readRDS("my_seurat.rds")            # already NormalizeData()'d
genes <- select.gene(obj, min.cells = 10)

obj <- obj |>
    prepare.seurat(genes = genes$genes) |>
    prepare.graph() |>
    prepare.reg(predictors = genes$tfs, responses = genes$targets)

obj <- run.nn.reg(obj, responses = genes$targets,
                  return.p.val = TRUE) |>
       build.meta.network()

ctr <- select.central.genes(obj, n.per.component = 4, k = 2)
obj <- prepare.visualise(obj, central.genes = ctr, n.clu = 4)

# A meta-network plot:
i  <- 1
mp <- Seurat::Misc(obj, "NNet.mod")$meta.network$p.val[,, i]
visualise.network(obj, i, meta.network = TRUE,
                  cutoff     = mean(apply(mp, 1, max)),
                  radius     = c(0.4, 0.7, 0.85, 1),
                  pie.radius = 0.04,
                  text.size  = 5)
```

---

## 9. Integrating with PhiSpace / pseudotime workflows

The natural play with PhiSpace + Module 6 pseudotime:

1. Run PhiSpace and a pseudotime method (e.g. Module 6's
   `pt_phi_score = adult − fetal`). This gives you a continuous
   per-cell scalar.
2. Use `PhiSpace::rankFeatures(method = "PLS")` with `pt_phi_score`
   as the response and `source = "assay"` to rank genes by their
   association with the maturation axis.
3. Pass the **top-ranked targets** to `select.gene(targets = …)` so
   that NNet builds co-expression networks among the genes that
   actually move along the trajectory of interest, rather than over
   the full default target list.
4. Use the meta-network `pcs` embedding to colour cells by pseudotime
   — cells in the same region of meta-network space share a
   co-expression program; whether that program shifts smoothly along
   pseudotime tells you whether the network rewires during
   maturation.

The TF list (`gene.list$tfs`) is generic prior; this stays unchanged.

---

## 10. Gotchas

- **Seurat object only.** Convert from SCE before starting.
- **Human gene symbols** in the built-in priors. Mouse needs
  user-supplied lists (`select.gene(tfs = …, targets = …)`).
- **`run.nn.reg` is the bottleneck.** Sub-sampling via `select.cell()`
  helps if the cell count is large. Otherwise compute scales linearly
  with cells × responses × predictors.
- **`return.p.val = TRUE`** is required if you want to filter edges
  in `visualise.network()` or build meaningful meta-networks. It
  doubles memory but is worth it.
- **`prepare.visualise` needs `central.genes`.** Without them you'll
  see an awkward default layout. Always run `select.central.genes()`
  first.
- **The `cutoff` argument to `visualise.network()` is a p-value**, so
  *lower* values keep more edges. Confusingly, p-values are inverted
  inside the package's significance scoring — sanity check by trying
  `cutoff = 0.05` vs `cutoff = 0.95` and seeing which gives the
  denser graph.
- **`gene.list` is not exhaustive.** Some perfectly good TF or
  target genes are missing from the prior. Pre-filter your data and
  inspect `genes$tfs` / `genes$targets` before assuming coverage.

---

## 11. References

- Package: <https://github.com/meiosis97/NeighbourNet>
- Docs: <https://meiosis97.github.io/NeighbourNet/>
- Cell-cycle vignette:
  <https://meiosis97.github.io/NeighbourNet/articles/cell-cycle.html>
- Manuscript: Deng, Mao, Choi & Lê Cao, *Genome Research* 36:785–801
  (2026), DOI 10.1101/gr.281171.125. PDF in
  [`references/NNet.pdf`](NNet.pdf). See LinkedIn post for
  plain-English motivation.
- Prior-knowledge sources used by the package:
  [OmniPath](https://omnipathdb.org/) (signaling) and
  [NicheNet](https://github.com/saeyslab/nichenetr) (regulatory
  potential).
- Co-expression vs GRN distinction: Saint-Antoine & Singh, *Nature
  Reviews Genetics* (2026) — <https://www.nature.com/articles/s41576-026-00939-1>
