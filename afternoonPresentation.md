# Afternoon Session — Mini-presentation Plan

Each afternoon module opens with a ~10-15 min mini-presentation before
participants touch the vignette. Goal: set up the statistical ideas and
the biological motivation, so the code that follows is "confirm what you
just heard" rather than "decode what it does".

---

## Module 5 — Continuous Phenotyping with Φ-Space (~15 min)

### Learning outcomes

By the end of the talk participants should be able to say:

- why a single hard label throws away biological information for this
  dataset;
- what the Φ-Space score matrix is (cells × reference phenotypes, per-cell
  per-phenotype real number);
- where the scores come from (PLS regression from reference expression
  onto a dummy phenotype matrix);
- why scores can be collapsed to labels but not vice versa.

### Slide structure

**1. Motivation — the Sim et al. heart data (~2 min)**

- Three developmental stages (foetal / young / adult), ~54k nuclei.
- Show a UMAP with hard CM labels, then the same UMAP coloured by
  continuous age.
- Punchline: the hard label hides the maturation gradient we want to
  trace.

**2. Hard vs soft annotation (~2 min)**

- Hard: clusters → marker lookup → one label per cluster. Module 3.
- Soft: for each cell and each reference phenotype, a real-valued
  "how much does this look like X?" score.
- Cells are allowed to be *partly* fetal-CM and *partly* adult-CM at the
  same time — this is the whole point.

**3. What Φ-Space does, in one picture (~3 min)**

- Reference: annotated SCE with phenotypes in `colData` (e.g.
  `cell_type`, `age_group`).
- Encode phenotype columns as a dummy matrix Y (cells × levels,
  ±1 coding).
- Fit PLS regression: expression X → Y. (Why PLS, not regression: it
  finds components that maximise covariance with Y, so noise orthogonal
  to the phenotypes is down-weighted for free — this is how Φ-Space
  avoids needing explicit batch correction between reference and query.)
- Apply the fitted model to the query → a scores matrix in
  `reducedDim(query, "PhiSpace")`.
- Visual: expression-space cartoon → PLS arrow → phenotype-space cartoon
  with the same cells.

**4. Why this dataset in particular (~2 min)**

- Developmental continuum → continuous annotation is natural.
- Multiple phenotype layers at once (cell type *and* age).
- Bulk references are allowed — relevant because there is no single
  cross-age snRNA-seq atlas of the human heart; bulk developmental
  time-courses are.
- Scores are a `reducedDim` → they drop straight into slingshot in
  Module 6.

**5. Design choices we made (~3 min)**

- Reference: Gao 2026 LV atlas, pre-balanced.
- **Splitting CM and Fibroblast into subtypes** — populous classes carry
  the biology we care about; splitting also helps class balance.
- **Rare cell type filter** (`cellTypeThreshold`) — PLS components for
  cell types with <20 reference cells are mostly noise.
- **scran normalisation** on both sides. Rank transform is an
  alternative (robust across platforms, but throws away count-level
  variation that PLS can exploit when both sides are snRNA-seq).

**6. What the vignette will do (~2 min)**

- Run `PhiSpace()` on the prepared reference + query.
- Inspect the score matrix via heatmap (cells × phenotypes).
- PCA biplot of the score matrix coloured by stage and by broad cell
  type.
- Collapse to hard labels with `getClass()` as a sanity check.
- Save scores for Module 6.

---

## Module 6 — Pseudotime (~15 min)

### Learning outcomes

Participants should be able to explain:

- what pseudotime is (a scalar per cell, monotone along the
  trajectory);
- the three statistical families and which one gives directionality
  "for free";
- why slingshot depends on both an embedding and a clustering;
- why we compare PCA vs. Φ-Space embeddings and slingshot vs. DPT
  method families;
- why velocity is discussed but not run on this dataset.

### Slide structure

**1. The problem statement (~2 min)**

- Clusters give discrete identities; often biology is a continuum.
- Pseudotime = a scalar per cell such that cells that are close in the
  process have close values.
- Input: an embedding. Output: one (or several, for branching
  trajectories) real numbers per cell.
- This is *not* real time, and it is *not* identified — it only orders
  cells. Direction must come from somewhere else.

**2. Three statistical families (~5 min)**

A single table, then one slide per family.

| Family              | Idea                                                        | Examples                 | Directionality | Topology it assumes |
|---------------------|-------------------------------------------------------------|--------------------------|----------------|---------------------|
| Principal curve/graph | Fit a skeleton through cluster centroids; pseudotime = arc length | slingshot, Monocle 2/3, TSCAN | External       | Line / tree         |
| Diffusion / random walk | Markov chain on cell-cell graph; pseudotime = accumulated / absorption probability | DPT, Palantir, Wishbone  | External (root cell) | Arbitrary graph    |
| Velocity-based      | Splicing dynamics give local time derivative; pseudotime from integrating or Markov-weighting velocities | scVelo, CellRank         | Intrinsic       | Arbitrary           |

Per-family slide covers: **core statistical object, what you have to
pick, main failure mode.**

**3. Slingshot — the main method for this module (~5 min)**

Principal-curve pseudotime is the default in most single-cell
workflows; we spend the majority of the talk here.

- Inputs: a reducedDim + a clustering. We'll compare PCA of scran
  logcounts against the Φ-Space score matrix, and compare two
  clusterings: unsupervised Louvain clusters and the fine PhiCellType
  labels from Module 5.
- Algorithm in three steps:
  1. Cluster centroids in the chosen embedding.
  2. Minimum spanning tree across the centroids → a tree of
     lineages.
  3. Fit a smooth principal curve through each lineage and project
     every cell onto it. Arc length along the projection = pseudotime
     for that lineage.
- Direction: external. `start.clus` (and optionally `end.clus`) pin the
  lineage orientation. We use the CM cluster with the highest mean
  `fetal` Φ-score as `start.clus`. This is the payoff of Module 5 —
  the age-group scores give us a principled, biological anchor rather
  than a hand-picked root cell.
- Outputs: one pseudotime per lineage, per cell, plus a cell-weight
  matrix assigning cells softly to lineages. Slots directly into
  `tradeSeq::fitGAM` for DE along pseudotime.
- **Comparison built into this module**: run slingshot four ways —
  two embeddings crossed with two clusterings. Same fetal-score anchor.
  The point is to show that clustering is not a preprocessing detail:
  it changes the MST that slingshot fits.

**4. DPT via destiny — the diffusion cross-check (~3 min)**

A different statistical family gives a useful sanity check. We run DPT
on both PCA and Φ-Space so the final comparison crosses two methods
with two embeddings.

- Why destiny: the only serious diffusion-pseudotime package that is
  pure R. Palantir / scFates / CellRank are Python; PAGA lives inside
  scanpy. Destiny keeps Module 6 self-contained.
- Statistical object, one slide:
  - Gaussian-kernel affinity matrix on the chosen embedding → row-
    normalise to a Markov transition matrix P.
  - Eigendecompose P. Top eigenvectors (diffusion components) are
    noise-smoothed coordinates of the trajectory.
  - DPT(i, j) = accumulated transition probability between cells. Fix
    a root → scalar pseudotime per cell.
- Three-line API: `DiffusionMap(data)` → (pick root) → `DPT(dm, ...)`.
- We run destiny on **PCA of scran logcounts** and on **Φ-Space
  scores**. That gives the same embedding comparison as slingshot, but
  with a graph/diffusion method rather than a principal curve.
- Direction: external. We pin the root cell to the cell with the
  highest fetal Φ-score, matching the slingshot anchor logic and keeping
  the pseudotimes comparable.
- Caveats to flag: O(n²) affinity matrix (fine at 10k cells, not at
  100k); destiny returns one scalar plus branch labels, not
  per-lineage arc lengths.

**5. Velocity — why we discuss but don't run it (~2 min)**

- Would give directionality from the data, not a hand-picked root.
- Requires unspliced counts — typically discarded by our upstream
  pipeline; for this workshop they are not available.
- On snRNA-seq specifically the unspliced fraction is unusually high
  and the standard velocity assumptions break down; even with counts
  this would not be a slam-dunk.

**6. Sensitivity and four-way comparison (~2 min)**

Six pseudotimes on the same CM subset:

- slingshot × PCA/Φ-Space × unsupervised/PhiCellType clusters;
- destiny-DPT × PCA/Φ-Space.

Two figures are worth articulating for participants:

- **Slingshot sensitivity grid** — shows how changing clustering moves
  the trajectory even when the embedding is fixed.
- **Canonical four-way comparison** — PhiCellType-clustered slingshot
  vs DPT, each on PCA and Φ-Space. Rows isolate embedding effects;
  columns isolate method effects.

Strong agreement across the canonical four is the case we hope for: the
trajectory is robust to both embedding and method. Disagreement is still
useful because it teaches participants how to sanity-check a pseudotime
result rather than treating one trajectory as truth.

---

## Module 7 — Cell-specific co-expression networks (NeighbourNet) (~15 min)

### Learning outcomes

Participants should be able to explain:

- why co-expression networks are useful but are not causal GRNs;
- how NeighbourNet turns noisy per-cell networks into denoised
  meta-networks;
- why we use `pt_phi_score` to select maturation-associated target genes;
- how to read the NeighbourNet meta-network plot without over-claiming
  directionality.

### Slide structure

**1. The framing: scaffold, not causal GRN (~3 min)**

- A true GRN is directional, mechanistic and causal.
- scRNA-seq expression alone mostly gives associations, not perturbation
  evidence.
- What we can build is a co-expression scaffold: a structured set of
  TF-target associations that may contain regulatory hypotheses.
- This distinction matters because the figure will look directional, but
  the direction comes partly from prior knowledge, not from the expression
  data alone.

**2. Why cell-specific networks? (~3 min)**

- A global network averages over fetal-like and adult-like cardiomyocytes.
- A per-cell network asks whether the local gene-gene relationships change
  across the maturation continuum.
- Individual per-cell networks are noisy, so NeighbourNet borrows strength
  from k-nearest neighbours in PC space.
- The key idea: many noisy local networks can be aggregated into stable
  recurring meta-networks.

**3. What NeighbourNet does (~4 min)**

- Input: a Seurat object with normalised expression.
- Choose targets and TFs. For this module, targets are the top
  maturation-associated genes from `PhiSpace::rankFeatures()` using
  `pt_phi_score`.
- PCA + kNN graph defines each cell's local neighbourhood.
- For each cell, regress target-gene expression on local PC structure and
  score TF-target co-expression.
- Stack all per-cell networks into a cell × edge matrix.
- Run non-negative PCA to get meta-networks and per-cell meta-network
  scores.

**4. Why `pt_phi_score`, not slingshot pseudotime? (~2 min)**

- `pt_phi_score = adult-like Φ-Space score - fetal Φ-Space score`.
- It is trajectory-free: no slingshot curve, no DPT graph, no clustering
  decision.
- That keeps Module 7 focused on network rewiring, while Module 6 remains
  the place where we compare trajectory methods.

**5. Reading the output (~3 min)**

- First plot: meta-network score space coloured by donor stage and
  `pt_phi_score`.
- Pick the maturation meta-network by correlation with `pt_phi_score`.
- Final plot: receptors / TFs / target clusters arranged by the prior
  knowledge graph; edge strength reflects co-expression evidence.
- Take-home: the plot suggests hypotheses about maturation-associated
  regulatory programs, but perturbation or external evidence is needed
  before calling any edge causal.
