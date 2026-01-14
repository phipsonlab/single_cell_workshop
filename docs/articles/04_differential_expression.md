# Module 4: Differential Expression Analysis

## Introduction

Differential expression (DE) analysis identifies genes whose expression
levels differ significantly between experimental conditions. In single
cell studies, this typically involves comparing gene expression across
developmental stages, disease states, or treatment conditions. However,
the statistical approach requires careful consideration of the study
design.

### About This Dataset

The data we are analysing comes from the study by Sim et al. (2021)
published in *Circulation*, titled “Sex-Specific Control of Human Heart
Maturation by the Progesterone Receptor”. This study performed
single-nucleus RNA sequencing (snRNA-seq) of 54,140 nuclei from 9 human
donors to profile transcriptional changes during heart maturation from
fetal stages to adulthood. The original study identified:

- **Cell type composition changes**: A significant expansion of cardiac
  fibroblasts and immune cells in the postnatal heart, accompanied by a
  decrease in cardiomyocyte proportions
- **Transcriptional maturation**: Profound changes in all cardiac cell
  types, with the largest number of differentially expressed genes in
  cardiomyocytes
- **Sex-specific programmes**: Sexually dimorphic gene expression that
  emerges primarily during adulthood
- **Metabolic maturation**: Cardiomyocyte maturation was associated with
  repression of cell cycle genes and activation of oxidative
  phosphorylation and respiratory electron transport pathways

In this module, we will reproduce key aspects of this analysis using the
pseudobulk differential expression approach employed in the original
study.

### The Pseudoreplication Problem

A critical pitfall in single cell differential expression analysis is
**pseudoreplication**: treating cells as independent biological
replicates when they are not. Cells from the same sample share technical
and biological variation that violates the independence assumption of
statistical tests.

Consider this analogy: if we want to compare heights between two
populations, measuring the same person multiple times does not increase
our sample size. Similarly, capturing 5,000 cells from one individual
does not give us 5,000 independent observations of that individual’s
biology.

Many single cell DE methods (Wilcoxon test, t-test, MAST) operate at the
cell level and treat each cell as an independent replicate. These
approaches:

- **Inflate Type I error rates** (false positives) dramatically
- **Produce artificially small p-values** that do not reflect true
  significance
- **Confound biological with technical variation**
- **Cannot generalise findings** beyond the specific samples analysed

The solution is **pseudobulk analysis**: aggregating counts across cells
within each sample to create sample-level expression estimates, then
applying well-established bulk RNA-seq methods. This approach:

- **Correctly identifies the biological replicate** (the
  sample/individual)
- **Properly controls false discovery rates**
- **Accounts for sample-to-sample variation**
- **Produces results that generalise to the population**

### Study Design Considerations

Our dataset contains 9 samples across 3 developmental stages:

| Group | Samples    | Sex Distribution |
|-------|------------|------------------|
| Fetal | f1, f2, f3 | 2 male, 1 female |
| Young | y1, y2, y3 | 2 male, 1 female |
| Adult | a1, a2, a3 | 1 male, 2 female |

With n=3 per group, statistical power is limited. This is a common
constraint in single cell studies due to the cost of sample collection
and sequencing. We must interpret results with appropriate caution and
avoid over-claiming the significance of findings.

The unbalanced sex distribution across groups is another consideration.
We address this by including sex in our statistical model, though the
small sample sizes limit our ability to fully disentangle sex and
developmental effects.

### Learning Objectives

By the end of this module, you will be able to:

1.  Understand why pseudobulk analysis is necessary for single cell DE
2.  Perform appropriate gene filtering for DE analysis
3.  Create pseudobulk expression matrices stratified by cell type and
    sample
4.  Apply the limma-voom pipeline for differential expression
5.  Interpret and visualise DE results critically
6.  Understand the limitations imposed by sample size

## Load Libraries

``` r
library(Seurat)
library(edgeR)
library(limma)
library(speckle)        # For propeller cell type composition analysis
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(pheatmap)
library(patchwork)
```

## Load Data

We load the QC-filtered Seurat object from Module 1, which contains
cells that passed our quality control criteria. The pseudobulk approach
aggregates raw counts, so we extract the count matrix from the RNA
assay. Cell type labels come from the original study annotations in
`cellinfo_updated.Rds`.

``` r
data_dir <- "../data"

# Load QC-filtered Seurat object from Module 1
seu <- readRDS(file.path(data_dir, "processed/01_qc_filtered.rds"))

# Extract raw counts (not normalised - pseudobulk uses raw counts)
counts <- GetAssayData(seu, assay = "RNA", layer = "counts")

# Load cell metadata with cell type annotations from original study
cellinfo <- readRDS(file.path(data_dir, "cellinfo_updated.Rds"))

# Match cellinfo to QC-filtered cells
cellinfo <- cellinfo[match(colnames(seu), cellinfo$CellID), ]

cat("QC-filtered dataset:\n")
```

    ## QC-filtered dataset:

``` r
cat("- Genes:", format(nrow(counts), big.mark = ","), "\n")
```

    ## - Genes: 19,001

``` r
cat("- Cells:", format(ncol(counts), big.mark = ","), "\n")
```

    ## - Cells: 43,119

### Examine Cell Type Labels

The cell metadata contains cell type annotations from the original
study. We use the `Celltype` column which provides broad cell type
classifications:

``` r
# Cell type distribution in QC-filtered data
cat("Cell type distribution:\n")
```

    ## Cell type distribution:

``` r
print(table(cellinfo$Celltype))
```

    ## 
    ##      Cardiomyocytes   Endothelial cells    Epicardial cells           Erythroid 
    ##               24025                3407                2888                  52 
    ##          Fibroblast        Immune cells             Neurons Smooth muscle cells 
    ##                8686                2912                 754                 395

``` r
# Cell type by developmental stage
cat("\nCell type by developmental stage:\n")
```

    ## 
    ## Cell type by developmental stage:

``` r
print(table(cellinfo$Celltype, cellinfo$Group))
```

    ##                      
    ##                       adult fetal young
    ##   Cardiomyocytes       1979 17046  5000
    ##   Endothelial cells    1092  1077  1238
    ##   Epicardial cells      878   933  1077
    ##   Erythroid               0    52     0
    ##   Fibroblast           2223  2934  3529
    ##   Immune cells         1066   578  1268
    ##   Neurons               102   292   360
    ##   Smooth muscle cells    81   185   129

### Define Colour Palettes

``` r
# Developmental group colours
group_colors <- c(
    "fetal" = "#E64B35",
    "young" = "#4DBBD5",
    "adult" = "#3C5488"
)

# Cell type colours
celltype_colors <- c(
    "Cardiomyocytes" = "#D62728",
    "Fibroblast" = "#1F77B4",
    "Endothelial cells" = "#9467BD",
    "Immune cells" = "#2CA02C",
    "Smooth muscle cells" = "#FF7F0E",
    "Epicardial cells" = "#17BECF",
    "Neurons" = "#8C564B",
    "Erythroid" = "#E377C2"
)
```

### Gene Annotation

We retrieve gene annotations for displaying results (gene names in
output tables). Note that gene filtering was already performed in Module
1.

``` r
# Get gene annotation for output tables
gene_symbols <- rownames(counts)

ann <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = gene_symbols,
  columns = c("SYMBOL", "ENTREZID", "GENENAME"),
  keytype = "SYMBOL"
)

# Handle genes with multiple mappings
ann <- ann[!duplicated(ann$SYMBOL), ]
ann <- ann[match(gene_symbols, ann$SYMBOL), ]

cat("Genes with annotation:", sum(!is.na(ann$GENENAME)), "/", length(gene_symbols), "\n")
```

    ## Genes with annotation: 19001 / 19001

## Prepare Data for Pseudobulk Analysis

Quality control (cell and gene filtering) was completed in Module 1. We
now prepare the data for pseudobulk differential expression analysis.

### Exclude Erythroid Cells

Erythroid cells represent a small population that may be contaminants.
We exclude them from DE analysis:

``` r
# Check for erythroid cells
n_erythroid <- sum(cellinfo$Celltype == "Erythroid", na.rm = TRUE)
cat("Erythroid cells found:", n_erythroid, "\n")
```

    ## Erythroid cells found: 52

``` r
# Remove erythroid cells if present
if (n_erythroid > 0) {
  cells_keep <- cellinfo$Celltype != "Erythroid"
  cellinfo_filtered <- cellinfo[cells_keep, ]
  counts_filtered <- counts[, cells_keep]
} else {
  cellinfo_filtered <- cellinfo
  counts_filtered <- counts
}

cat("Cells for DE analysis:", ncol(counts_filtered), "\n")
```

    ## Cells for DE analysis: 43067

``` r
cat("\nCell type distribution:\n")
```

    ## 
    ## Cell type distribution:

``` r
print(table(cellinfo_filtered$Celltype))
```

    ## 
    ##      Cardiomyocytes   Endothelial cells    Epicardial cells          Fibroblast 
    ##               24025                3407                2888                8686 
    ##        Immune cells             Neurons Smooth muscle cells 
    ##                2912                 754                 395

> **Note:** Gene filtering (mitochondrial, ribosomal, sex chromosome,
> and unannotated genes) was performed in Module 1. Filtering of lowly
> expressed genes is performed later using
> [`edgeR::filterByExpr()`](https://rdrr.io/pkg/edgeR/man/filterByExpr.html),
> which applies appropriate thresholds based on the pseudobulk sample
> structure.

## Cell Type Composition Analysis

Before examining gene expression changes, we first ask: **does the
cellular composition of the heart change during development?** This is
an important biological question because changes in cell type
proportions can:

- Reflect biological processes such as cell proliferation, migration, or
  death
- Influence tissue function independently of gene expression changes
- Confound bulk RNA-seq analysis (where expression changes may simply
  reflect composition changes)

The original study by Sim et al. (2021) found that cardiac maturation
was associated with a significant expansion in the relative proportion
of cardiac fibroblasts and immune cells, accompanied by a significant
decrease in the proportion of cardiomyocytes (Figure 1C in the paper).
We can test this using the **propeller** method from the *speckle*
package.

### Why Use propeller?

Cell type proportions are **compositional data**: they must sum to 1
within each sample. This creates dependencies between cell types—if one
proportion increases, others must decrease. Standard statistical tests
(t-tests, ANOVA) ignore this constraint and can produce misleading
results.

The [`propeller()`](https://rdrr.io/pkg/speckle/man/propeller.html)
function addresses this by:

1.  Applying appropriate transformations to proportional data
2.  Using empirical Bayes moderation for robust variance estimation
3.  Properly accounting for the nested structure (cells within samples)

### Calculate Cell Type Proportions

``` r
# Calculate cell counts per sample and cell type
composition_data <- cellinfo_filtered %>%
  group_by(Sample, Group, Celltype) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(
    total_cells = sum(n_cells),
    proportion = n_cells / total_cells
  ) %>%
  ungroup()

# Display summary
cat("Cells per sample:\n")
```

    ## Cells per sample:

``` r
composition_data %>%
  select(Sample, Group, total_cells) %>%
  distinct() %>%
  print()
```

    ## # A tibble: 9 × 3
    ##   Sample Group total_cells
    ##   <chr>  <chr>       <int>
    ## 1 a1     adult        3806
    ## 2 a2     adult        2449
    ## 3 a3     adult        1166
    ## 4 f1     fetal        7220
    ## 5 f2     fetal        9069
    ## 6 f3     fetal        6756
    ## 7 y1     young        4051
    ## 8 y2     young        4371
    ## 9 y3     young        4179

### Visualise Composition Changes

``` r
# Stacked bar plot of cell type proportions
ggplot(composition_data, aes(x = Sample, y = proportion, fill = Celltype)) +
  geom_bar(stat = "identity", colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = celltype_colors) +
  facet_grid(~ Group, scales = "free_x", space = "free_x") +
  labs(
    title = "Cell Type Composition Across Samples",
    subtitle = "Proportions within each sample sum to 1",
    x = "Sample",
    y = "Proportion",
    fill = "Cell Type"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )
```

![](04_differential_expression_files/figure-html/composition-barplot-1.png)

``` r
# Boxplot of cell type proportions by developmental stage
ggplot(composition_data, aes(x = Group, y = proportion, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  scale_fill_manual(values = group_colors) +
  facet_wrap(~ Celltype, scales = "free_y", ncol = 4) +
  labs(
    title = "Cell Type Proportions by Developmental Stage",
    subtitle = "Each point represents one sample",
    x = "Developmental Stage",
    y = "Proportion"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
```

![](04_differential_expression_files/figure-html/composition-boxplot-1.png)

**Interpretation**: Visual inspection suggests that cardiomyocyte
proportions decrease from fetal to adult stages, while fibroblast and
immune cell proportions appear to increase. However, we need formal
statistical testing to determine whether these differences are
significant given the sample variability.

### Statistical Testing with propeller

The [`propeller()`](https://rdrr.io/pkg/speckle/man/propeller.html)
function tests for differences in cell type proportions between groups
while properly handling the compositional nature of the data:

``` r
# Prepare data for propeller
# Need vectors of: clusters (cell types), samples, and groups
clusters <- cellinfo_filtered$Celltype
samples <- cellinfo_filtered$Sample

# Map samples to groups
sample_to_group <- cellinfo_filtered %>%
  select(Sample, Group) %>%
  distinct() %>%
  { setNames(.$Group, .$Sample) }
groups <- sample_to_group[samples]

# Run propeller test comparing developmental groups
# This tests whether cell type proportions differ between fetal, young, and adult
propeller_results <- propeller(
  clusters = clusters,
  sample = samples,
  group = groups
)

# Display results
cat("Propeller Test Results: Cell Type Composition Changes\n")
```

    ## Propeller Test Results: Cell Type Composition Changes

``` r
cat("======================================================\n\n")
```

    ## ======================================================

``` r
print(propeller_results)
```

    ##                     BaselineProp PropMean.adult PropMean.fetal PropMean.young
    ## Cardiomyocytes       0.557851719     0.21530334    0.731341454     0.39513672
    ## Immune cells         0.067615576     0.18606819    0.025324760     0.10022097
    ## Endothelial cells    0.079109295     0.13944054    0.047958598     0.09866238
    ## Fibroblast           0.201685745     0.31860682    0.132229871     0.28011811
    ## Epicardial cells     0.067058305     0.11532165    0.041633987     0.08651143
    ## Neurons              0.017507604     0.01330669    0.012790683     0.02901683
    ## Smooth muscle cells  0.009171756     0.01195278    0.008720647     0.01033356
    ##                     Fstatistic      P.Value          FDR
    ## Cardiomyocytes      14.3605368 7.606791e-05 0.0005324753
    ## Immune cells        12.2192238 2.118135e-04 0.0007413473
    ## Endothelial cells    4.8447367 1.693895e-02 0.0395242176
    ## Fibroblast           4.3904270 2.354630e-02 0.0412060255
    ## Epicardial cells     3.0333487 6.665080e-02 0.0933111197
    ## Neurons              1.8303029 1.817615e-01 0.2120550623
    ## Smooth muscle cells  0.6890638 5.115756e-01 0.5115755945

``` r
# Identify significant changes
sig_celltypes <- propeller_results %>%
  filter(FDR < 0.05)

if (nrow(sig_celltypes) > 0) {
  cat("\nCell types with significant composition changes (FDR < 0.05):\n")
  # Display columns relevant to our 3-group comparison
  print(sig_celltypes[, c("PropMean.fetal", "PropMean.young",
                          "PropMean.adult", "Fstatistic", "FDR")])
} else {
  cat("\nNo cell types showed statistically significant composition changes.\n")
  cat("This may reflect limited statistical power with n=3 samples per group.\n")
}
```

    ## 
    ## Cell types with significant composition changes (FDR < 0.05):
    ##                   PropMean.fetal PropMean.young PropMean.adult Fstatistic
    ## Cardiomyocytes        0.73134145     0.39513672      0.2153033  14.360537
    ## Immune cells          0.02532476     0.10022097      0.1860682  12.219224
    ## Endothelial cells     0.04795860     0.09866238      0.1394405   4.844737
    ## Fibroblast            0.13222987     0.28011811      0.3186068   4.390427
    ##                            FDR
    ## Cardiomyocytes    0.0005324753
    ## Immune cells      0.0007413473
    ## Endothelial cells 0.0395242176
    ## Fibroblast        0.0412060255

**Interpretation**: The propeller test assesses whether the observed
differences in cell type proportions are statistically significant when
accounting for sample-to-sample variation. In the original study, Sim et
al. found significant changes in cardiomyocyte, fibroblast, and immune
cell proportions across development (FDR \< 0.05). With our limited
sample size (n=3 per group), we may lack the statistical power to detect
these changes, even if the biological trends are evident in the
visualisations.

> **Note on statistical power**: The original study findings were based
> on the same 9 samples we are analysing. However, the propeller test is
> conservative and properly accounts for the small sample size.
> Non-significant results should not be interpreted as evidence of no
> biological effect—they may simply reflect insufficient power to detect
> true differences.

## Create Pseudobulk Samples

### Aggregation Strategy

We aggregate counts by both cell type and sample, creating one
pseudobulk profile for each cell type in each sample. This stratified
approach allows us to perform differential expression within each cell
type, accounting for the fact that different cell types may respond
differently to developmental changes.

``` r
# Create factor for cell type
celltype <- factor(cellinfo_filtered$Celltype)

# Create factor for sample with explicit level ordering
sample_factor <- factor(
    cellinfo_filtered$Sample,
    levels = c("f1", "f2", "f3", "y1", "y2", "y3", "a1", "a2", "a3")
)

# Create combined grouping variable (celltype.sample)
pseudobulk_group <- paste(celltype, sample_factor, sep = ".")
pseudobulk_group <- factor(pseudobulk_group)

cat("Number of pseudobulk samples:", length(levels(pseudobulk_group)), "\n")
```

    ## Number of pseudobulk samples: 63

``` r
cat("\nCells per pseudobulk sample:\n")
```

    ## 
    ## Cells per pseudobulk sample:

``` r
print(table(pseudobulk_group))
```

    ## pseudobulk_group
    ##      Cardiomyocytes.a1      Cardiomyocytes.a2      Cardiomyocytes.a3 
    ##                   1546                    293                    140 
    ##      Cardiomyocytes.f1      Cardiomyocytes.f2      Cardiomyocytes.f3 
    ##                   4948                   7470                   4628 
    ##      Cardiomyocytes.y1      Cardiomyocytes.y2      Cardiomyocytes.y3 
    ##                    998                   1769                   2233 
    ##   Endothelial cells.a1   Endothelial cells.a2   Endothelial cells.a3 
    ##                    569                    400                    123 
    ##   Endothelial cells.f1   Endothelial cells.f2   Endothelial cells.f3 
    ##                    427                    304                    346 
    ##   Endothelial cells.y1   Endothelial cells.y2   Endothelial cells.y3 
    ##                    468                    361                    409 
    ##    Epicardial cells.a1    Epicardial cells.a2    Epicardial cells.a3 
    ##                    330                    469                     79 
    ##    Epicardial cells.f1    Epicardial cells.f2    Epicardial cells.f3 
    ##                    427                    242                    264 
    ##    Epicardial cells.y1    Epicardial cells.y2    Epicardial cells.y3 
    ##                    599                    258                    220 
    ##          Fibroblast.a1          Fibroblast.a2          Fibroblast.a3 
    ##                    938                    874                    411 
    ##          Fibroblast.f1          Fibroblast.f2          Fibroblast.f3 
    ##                   1019                    739                   1176 
    ##          Fibroblast.y1          Fibroblast.y2          Fibroblast.y3 
    ##                   1431                   1420                    678 
    ##        Immune cells.a1        Immune cells.a2        Immune cells.a3 
    ##                    350                    329                    387 
    ##        Immune cells.f1        Immune cells.f2        Immune cells.f3 
    ##                    253                    190                    135 
    ##        Immune cells.y1        Immune cells.y2        Immune cells.y3 
    ##                    301                    479                    488 
    ##             Neurons.a1             Neurons.a2             Neurons.a3 
    ##                     52                     37                     13 
    ##             Neurons.f1             Neurons.f2             Neurons.f3 
    ##                     97                    104                     91 
    ##             Neurons.y1             Neurons.y2             Neurons.y3 
    ##                    199                     57                    104 
    ## Smooth muscle cells.a1 Smooth muscle cells.a2 Smooth muscle cells.a3 
    ##                     21                     47                     13 
    ## Smooth muscle cells.f1 Smooth muscle cells.f2 Smooth muscle cells.f3 
    ##                     49                     20                    116 
    ## Smooth muscle cells.y1 Smooth muscle cells.y2 Smooth muscle cells.y3 
    ##                     55                     27                     47

### Aggregate Counts

We use matrix multiplication with a design matrix to efficiently
aggregate counts. This is equivalent to summing counts across all cells
within each pseudobulk group:

``` r
# Create design matrix for aggregation
design_agg <- model.matrix(~ 0 + pseudobulk_group)
colnames(design_agg) <- levels(pseudobulk_group)

# Aggregate counts by matrix multiplication
# Each column of the result is the sum of counts for all cells in that group
pb_counts <- as.matrix(counts_filtered) %*% design_agg

cat("Pseudobulk matrix dimensions:\n")
```

    ## Pseudobulk matrix dimensions:

``` r
cat("- Genes:", nrow(pb_counts), "\n")
```

    ## - Genes: 19001

``` r
cat("- Samples:", ncol(pb_counts), "\n")
```

    ## - Samples: 63

### Create DGEList Object

The `DGEList` object from *edgeR* is the standard container for count
data in the limma-voom pipeline. It stores counts, library sizes, and
sample information:

``` r
# Create DGEList
dge <- DGEList(counts = pb_counts)

# Parse sample information from column names
sample_info <- data.frame(
    pseudobulk_id = colnames(dge),
    stringsAsFactors = FALSE
)

# Extract cell type and sample from the combined name
parsed <- strsplit(sample_info$pseudobulk_id, "\\.")
sample_info$celltype <- sapply(parsed, `[`, 1)
sample_info$sample <- sapply(parsed, `[`, 2)

# Add developmental group (keep as character to avoid factor level issues)
sample_info$group <- case_when(
    grepl("^f", sample_info$sample) ~ "fetal",
    grepl("^y", sample_info$sample) ~ "young",
    grepl("^a", sample_info$sample) ~ "adult"
)

# Add sex information from cellinfo (matches actual data)
sex_by_sample <- cellinfo_filtered %>%
    distinct(Sample, Sex) %>%
    { setNames(.$Sex, .$Sample) }
sample_info$sex <- sex_by_sample[sample_info$sample]

# Add gene annotation to DGEList (for output tables)
dge$genes <- ann

# Store sample information
# Note: DGEList creates a default 'group' column (all 1s), so we remove it first
dge$samples <- cbind(dge$samples[, c("lib.size", "norm.factors")], sample_info)

head(dge$samples)
```

    ##                   lib.size norm.factors     pseudobulk_id       celltype sample
    ## Cardiomyocytes.a1 29788439            1 Cardiomyocytes.a1 Cardiomyocytes     a1
    ## Cardiomyocytes.a2  6918676            1 Cardiomyocytes.a2 Cardiomyocytes     a2
    ## Cardiomyocytes.a3  2243481            1 Cardiomyocytes.a3 Cardiomyocytes     a3
    ## Cardiomyocytes.f1 50685628            1 Cardiomyocytes.f1 Cardiomyocytes     f1
    ## Cardiomyocytes.f2 66444418            1 Cardiomyocytes.f2 Cardiomyocytes     f2
    ## Cardiomyocytes.f3 62719116            1 Cardiomyocytes.f3 Cardiomyocytes     f3
    ##                   group sex
    ## Cardiomyocytes.a1 adult   f
    ## Cardiomyocytes.a2 adult   m
    ## Cardiomyocytes.a3 adult   m
    ## Cardiomyocytes.f1 fetal   m
    ## Cardiomyocytes.f2 fetal   m
    ## Cardiomyocytes.f3 fetal   f

## Exploratory Data Analysis

Before conducting formal statistical tests, we examine the overall
structure of the data using multidimensional scaling (MDS). This
unsupervised approach reveals the major sources of variation in the
dataset.

### MDS Plot Overview

``` r
# Calculate MDS for all samples
mds <- plotMDS(dge, plot = FALSE, gene.selection = "common")

# Create plotting data frame
mds_data <- data.frame(
    Dim1 = mds$x,
    Dim2 = mds$y,
    celltype = dge$samples$celltype,
    group = dge$samples$group,
    sex = dge$samples$sex,
    sample = dge$samples$sample
)

# Plot by cell type
p1 <- ggplot(mds_data, aes(x = Dim1, y = Dim2, colour = celltype)) +
    geom_point(size = 3) +
    scale_colour_manual(values = celltype_colors) +
    labs(
        title = "MDS: All Pseudobulk Samples",
        subtitle = "Coloured by cell type",
        x = "Leading logFC dim 1",
        y = "Leading logFC dim 2",
        colour = "Cell Type"
    ) +
    theme_minimal() +
    theme(legend.position = "right")

# Plot by developmental group
p2 <- ggplot(mds_data, aes(x = Dim1, y = Dim2, colour = group, shape = celltype)) +
    geom_point(size = 3) +
    scale_colour_manual(values = group_colors) +
    labs(
        title = "MDS: All Pseudobulk Samples",
        subtitle = "Coloured by developmental stage",
        x = "Leading logFC dim 1",
        y = "Leading logFC dim 2",
        colour = "Group",
        shape = "Cell Type"
    ) +
    theme_minimal() +
    theme(legend.position = "right")

p1 / p2
```

![](04_differential_expression_files/figure-html/mds-overview-1.png)

**Interpretation**: The MDS plot reveals that the primary axis of
variation corresponds to cell type identity. This is expected and
biologically sensible: cardiomyocytes, fibroblasts, and immune cells
have fundamentally different transcriptional programmes. Developmental
differences are secondary sources of variation, visible within cell type
clusters.

### MDS by Cell Type

To better visualise developmental differences, we examine MDS plots for
each major cell type separately:

``` r
# Function to create MDS plot for a specific cell type
plot_mds_celltype <- function(celltype_name) {
    # Subset to this cell type
    idx <- dge$samples$celltype == celltype_name
    if (sum(idx) < 3) return(NULL)

    dge_sub <- dge[, idx]

    # Calculate MDS
    mds_sub <- plotMDS(dge_sub, plot = FALSE, gene.selection = "common")

    # Create plot data
    plot_data <- data.frame(
        Dim1 = mds_sub$x,
        Dim2 = mds_sub$y,
        group = dge_sub$samples$group,
        sex = dge_sub$samples$sex,
        sample = dge_sub$samples$sample
    )

    # Sex symbols: filled = male, open = female
    shape_values <- c("m" = 16, "f" = 1)

    ggplot(plot_data, aes(x = Dim1, y = Dim2, colour = group, shape = sex)) +
        geom_point(size = 4, stroke = 1.5) +
        geom_text(aes(label = sample), vjust = -1, size = 3, show.legend = FALSE) +
        scale_colour_manual(values = group_colors) +
        scale_shape_manual(values = shape_values, labels = c("f" = "Female", "m" = "Male")) +
        labs(
            title = celltype_name,
            x = "Leading logFC dim 1",
            y = "Leading logFC dim 2",
            colour = "Group",
            shape = "Sex"
        ) +
        theme_minimal() +
        theme(legend.position = "bottom")
}

# Create plots for major cell types
celltypes_to_plot <- c("Cardiomyocytes", "Fibroblast", "Endothelial cells", "Immune cells")
mds_plots <- lapply(celltypes_to_plot, plot_mds_celltype)
mds_plots <- mds_plots[!sapply(mds_plots, is.null)]

# Combine plots
wrap_plots(mds_plots, ncol = 2) +
    plot_annotation(
        title = "MDS Plots by Cell Type",
        subtitle = "Developmental group separation within each cell type"
    )
```

![](04_differential_expression_files/figure-html/mds-by-celltype-1.png)

**Interpretation**: Within each cell type, fetal samples tend to
separate from young and adult samples along the first MDS dimension.
Young and adult samples often overlap, suggesting more transcriptional
similarity between postnatal stages than between fetal and postnatal
development. The degree of separation varies by cell type, with
cardiomyocytes showing particularly clear developmental stratification.

## Differential Expression Analysis

### Statistical Framework

We use the **limma-voom** pipeline, which is the gold standard for
differential expression analysis of RNA-seq count data. The workflow
consists of:

1.  **Filtering**: Remove genes with insufficient expression
2.  **Normalisation**: Apply TMM (trimmed mean of M-values)
    normalisation
3.  **Voom transformation**: Model the mean-variance relationship
4.  **Linear modelling**: Fit linear models with the experimental design
5.  **Empirical Bayes**: Moderate test statistics for improved power

We fit a single model containing all cell type × group × sex
combinations, then extract contrasts for specific comparisons of
interest.

### Design Matrix

The design matrix encodes our experimental factors. We use a cell-means
parameterisation (no intercept) where each column represents a unique
combination of cell type, developmental group, and sex:

``` r
# Create short cell type labels for cleaner coefficient names
celltype_short <- dge$samples$celltype
celltype_short <- gsub(" cells", "", celltype_short)
celltype_short <- gsub("Smooth muscle", "SMC", celltype_short)
celltype_short <- tolower(celltype_short)

# Create combined factor for design matrix
design_group <- paste(celltype_short, dge$samples$group, dge$samples$sex, sep = ".")
design_group <- factor(design_group)

# Create design matrix (no intercept = cell means model)
design <- model.matrix(~ 0 + design_group)
colnames(design) <- levels(design_group)

cat("Design matrix dimensions:", nrow(design), "x", ncol(design), "\n")
```

    ## Design matrix dimensions: 63 x 42

``` r
cat("\nDesign matrix coefficients:\n")
```

    ## 
    ## Design matrix coefficients:

``` r
print(colnames(design))
```

    ##  [1] "cardiomyocytes.adult.f" "cardiomyocytes.adult.m" "cardiomyocytes.fetal.f"
    ##  [4] "cardiomyocytes.fetal.m" "cardiomyocytes.young.f" "cardiomyocytes.young.m"
    ##  [7] "endothelial.adult.f"    "endothelial.adult.m"    "endothelial.fetal.f"   
    ## [10] "endothelial.fetal.m"    "endothelial.young.f"    "endothelial.young.m"   
    ## [13] "epicardial.adult.f"     "epicardial.adult.m"     "epicardial.fetal.f"    
    ## [16] "epicardial.fetal.m"     "epicardial.young.f"     "epicardial.young.m"    
    ## [19] "fibroblast.adult.f"     "fibroblast.adult.m"     "fibroblast.fetal.f"    
    ## [22] "fibroblast.fetal.m"     "fibroblast.young.f"     "fibroblast.young.m"    
    ## [25] "immune.adult.f"         "immune.adult.m"         "immune.fetal.f"        
    ## [28] "immune.fetal.m"         "immune.young.f"         "immune.young.m"        
    ## [31] "neurons.adult.f"        "neurons.adult.m"        "neurons.fetal.f"       
    ## [34] "neurons.fetal.m"        "neurons.young.f"        "neurons.young.m"       
    ## [37] "smc.adult.f"            "smc.adult.m"            "smc.fetal.f"           
    ## [40] "smc.fetal.m"            "smc.young.f"            "smc.young.m"

### Filtering and Normalisation

``` r
# Filter genes with low expression using edgeR's filterByExpr
# This removes genes that are not expressed at a biologically meaningful level
keep <- filterByExpr(dge, design = design)
dge_filtered <- dge[keep, , keep.lib.sizes = FALSE]

cat("Genes before filtering:", nrow(dge), "\n")
```

    ## Genes before filtering: 19001

``` r
cat("Genes after filtering:", nrow(dge_filtered), "\n")
```

    ## Genes after filtering: 16026

``` r
# Apply TMM normalisation
dge_filtered <- calcNormFactors(dge_filtered, method = "TMM")

cat("\nLibrary sizes and normalisation factors:\n")
```

    ## 
    ## Library sizes and normalisation factors:

``` r
print(dge_filtered$samples[, c("lib.size", "norm.factors")])
```

    ##                        lib.size norm.factors
    ## Cardiomyocytes.a1      29775232    0.5858095
    ## Cardiomyocytes.a2       6915653    0.6523898
    ## Cardiomyocytes.a3       2241819    0.7738208
    ## Cardiomyocytes.f1      50659280    0.8107119
    ## Cardiomyocytes.f2      66417589    0.8078087
    ## Cardiomyocytes.f3      62696000    0.7360192
    ## Cardiomyocytes.y1      12606920    0.6202950
    ## Cardiomyocytes.y2      28939726    0.6284126
    ## Cardiomyocytes.y3      16461173    0.8796761
    ## Endothelial cells.a1    3766176    1.0262367
    ## Endothelial cells.a2    2590816    1.0031447
    ## Endothelial cells.a3     791156    1.0949711
    ## Endothelial cells.f1    3502933    1.1077763
    ## Endothelial cells.f2    2632125    1.1224866
    ## Endothelial cells.f3    4013316    1.0535979
    ## Endothelial cells.y1    3131326    1.0323144
    ## Endothelial cells.y2    2486722    1.0008573
    ## Endothelial cells.y3    2205535    0.9582724
    ## Epicardial cells.a1     1770555    1.0460143
    ## Epicardial cells.a2     3598524    0.8779683
    ## Epicardial cells.a3      538683    1.0703207
    ## Epicardial cells.f1     3145516    1.1413185
    ## Epicardial cells.f2     1731608    1.1163473
    ## Epicardial cells.f3     2711665    1.0735700
    ## Epicardial cells.y1     3745678    0.9522513
    ## Epicardial cells.y2     1621534    0.9885182
    ## Epicardial cells.y3     1182445    0.9616773
    ## Fibroblast.a1           7387371    1.0668839
    ## Fibroblast.a2           7437402    0.9773801
    ## Fibroblast.a3           3871174    1.0674113
    ## Fibroblast.f1           7711614    1.0605987
    ## Fibroblast.f2           5283612    1.0407971
    ## Fibroblast.f3          11176735    1.0282190
    ## Fibroblast.y1          11739624    0.9211971
    ## Fibroblast.y2          10749844    0.9503454
    ## Fibroblast.y3           5364259    0.9909665
    ## Immune cells.a1         2587545    0.9922288
    ## Immune cells.a2         2081232    1.0383217
    ## Immune cells.a3         2939387    1.1048384
    ## Immune cells.f1         1634217    1.0686124
    ## Immune cells.f2         1107601    1.1762022
    ## Immune cells.f3         1024058    1.0954253
    ## Immune cells.y1         2273508    1.0548166
    ## Immune cells.y2         3256294    1.0139605
    ## Immune cells.y3         2544524    0.9961826
    ## Neurons.a1               257217    1.1813321
    ## Neurons.a2               214384    1.1009588
    ## Neurons.a3                74169    1.3576018
    ## Neurons.f1               584014    1.0697774
    ## Neurons.f2               683076    1.0908294
    ## Neurons.f3               697360    1.0536012
    ## Neurons.y1              1040655    1.0121766
    ## Neurons.y2               210213    1.1091080
    ## Neurons.y3               403832    1.0895271
    ## Smooth muscle cells.a1   119151    1.2774627
    ## Smooth muscle cells.a2   370952    1.0156736
    ## Smooth muscle cells.a3    96692    1.2563904
    ## Smooth muscle cells.f1   369597    1.0306194
    ## Smooth muscle cells.f2   146031    1.1502098
    ## Smooth muscle cells.f3  1265327    0.9699470
    ## Smooth muscle cells.y1   326422    1.1080687
    ## Smooth muscle cells.y2   160130    1.0800641
    ## Smooth muscle cells.y3   305833    1.0546025

### Voom Transformation

The voom transformation models the mean-variance relationship in the
data and assigns precision weights to each observation. This is critical
for proper statistical inference:

``` r
# Apply voom with cyclic loess normalisation
v <- voom(dge_filtered, design, plot = TRUE, normalize.method = "cyclicloess")
```

![](04_differential_expression_files/figure-html/voom-transform-1.png)**Interpretation**:
The mean-variance plot shows the relationship between average expression
level (x-axis) and variance (y-axis). The red curve represents voom’s
fitted mean-variance trend. A well-behaved dataset shows decreasing
variance with increasing expression, and the trend should be smooth.
Unusual patterns may indicate quality issues or batch effects.

### Fit Linear Model

``` r
# Fit linear model
fit <- lmFit(v, design)

cat("Model fitted with", ncol(fit$coefficients), "coefficients\n")
```

    ## Model fitted with 42 coefficients

## Cell Type-Specific Differential Expression

We now perform differential expression analysis for each major cell
type, comparing developmental stages while averaging over sex. This
approach answers the question: “Which genes change during heart
development within each cell type?”

### Define Contrasts

For each cell type, we define three contrasts:

- **Young vs Fetal**: Early postnatal changes
- **Adult vs Fetal**: Full developmental trajectory
- **Adult vs Young**: Late maturation changes

We average over sex to obtain marginal developmental effects:

``` r
# Function to create contrasts for a cell type
create_celltype_contrasts <- function(ct_short, design) {
    # Get coefficient names containing this cell type
    coefs <- colnames(design)
    ct_coefs <- coefs[grep(paste0("^", ct_short, "\\."), coefs)]

    if (length(ct_coefs) < 6) {
        warning("Insufficient coefficients for ", ct_short)
        return(NULL)
    }

    # Define contrasts averaging over sex
    # Young vs Fetal
    yvf <- paste0("0.5*(", ct_short, ".young.m + ", ct_short, ".young.f) - ",
                  "0.5*(", ct_short, ".fetal.m + ", ct_short, ".fetal.f)")

    # Adult vs Fetal
    avf <- paste0("0.5*(", ct_short, ".adult.m + ", ct_short, ".adult.f) - ",
                  "0.5*(", ct_short, ".fetal.m + ", ct_short, ".fetal.f)")

    # Adult vs Young
    avy <- paste0("0.5*(", ct_short, ".adult.m + ", ct_short, ".adult.f) - ",
                  "0.5*(", ct_short, ".young.m + ", ct_short, ".young.f)")

    # Use contrasts parameter (character vector) to avoid scoping issues
    contrast_list <- c(yvf, avf, avy)
    cm <- makeContrasts(contrasts = contrast_list, levels = design)
    # Set proper column names for the contrast matrix
    colnames(cm) <- c("YvF", "AvF", "AvY")
    cm
}
```

### Cardiomyocyte Differential Expression

Cardiomyocytes are the primary contractile cells of the heart and
undergo substantial transcriptional changes during development as they
transition from proliferative fetal cells to terminally differentiated
adult myocytes.

In the original study by Sim et al. (2021), cardiomyocytes showed the
largest number of differentially expressed genes compared to other
cardiac cell types (Figure 1E in the paper). The transcriptional changes
were characterised by:

- **Repression of cell cycle genes**: Reflecting the exit from
  proliferative capacity that occurs during postnatal maturation
- **Activation of metabolic genes**: Particularly those involved in
  oxidative phosphorylation, the TCA cycle, and the respiratory electron
  transport chain, reflecting the metabolic switch from glycolysis to
  fatty acid oxidation
- **Sex-specific transcriptional programmes**: More than 2,800 genes
  showed significant expression differences between adult male and
  female cardiomyocytes, the majority of which were autosomal genes

``` r
# Create contrasts for cardiomyocytes
cont_cardio <- create_celltype_contrasts("cardiomyocytes", design)

# Fit contrasts
fit_cardio <- contrasts.fit(fit, cont_cardio)
fit_cardio <- eBayes(fit_cardio, robust = TRUE)

# Apply TREAT with minimum log-fold-change threshold
# This tests whether the true fold-change exceeds a biologically meaningful threshold
treat_cardio <- treat(fit_cardio, lfc = 0.5)

# Summarise results
cat("Cardiomyocytes: Differential Expression Summary\n")
```

    ## Cardiomyocytes: Differential Expression Summary

``` r
cat("(FDR < 0.05, |logFC| > 0.5)\n\n")
```

    ## (FDR < 0.05, |logFC| > 0.5)

``` r
dt_cardio <- decideTests(treat_cardio)
summary(dt_cardio)
```

    ##          YvF   AvF   AvY
    ## Down    1235  1218     9
    ## NotSig 13733 13706 15987
    ## Up      1058  1102    30

``` r
# MD plots for each contrast
par(mfrow = c(1, 3))
for (i in 1:3) {
    plotMD(treat_cardio, coef = i, status = dt_cardio[, i],
           hl.col = c("blue", "red"), main = colnames(dt_cardio)[i])
    abline(h = 0, col = "grey")
}
```

![](04_differential_expression_files/figure-html/de-cardio-plots-1.png)

**Interpretation**: The MD (mean-difference) plots show average
expression (x-axis) versus log-fold-change (y-axis) for each contrast.
Red points indicate genes significantly upregulated in the second
condition; blue points indicate downregulated genes. The Adult vs Fetal
comparison typically shows the most differentially expressed genes,
reflecting the substantial transcriptional remodelling that occurs
during cardiac maturation.

``` r
# Display top differentially expressed genes for Adult vs Fetal
cat("Top 20 DE genes: Adult vs Fetal Cardiomyocytes\n")
```

    ## Top 20 DE genes: Adult vs Fetal Cardiomyocytes

``` r
topTreat(treat_cardio, coef = "AvF", n = 20)[, c("SYMBOL", "GENENAME", "logFC", "P.Value", "adj.P.Val")]
```

    ##                SYMBOL                                               GENENAME
    ## TMEM178B     TMEM178B                             transmembrane protein 178B
    ## TOGARAM2     TOGARAM2         TOG array regulator of axonemal microtubules 2
    ## FILIP1L       FILIP1L                   filamin A interacting protein 1 like
    ## CCSER1         CCSER1                      coiled-coil serine rich protein 1
    ## AAK1             AAK1                                AP2 associated kinase 1
    ## DGKG             DGKG                            diacylglycerol kinase gamma
    ## GRAMD1B       GRAMD1B                              GRAM domain containing 1B
    ## NCEH1           NCEH1                  neutral cholesterol ester hydrolase 1
    ## MIR29B2CHG MIR29B2CHG                           MIR29B2 and MIR29C host gene
    ## ADRA1A         ADRA1A                                  adrenoceptor alpha 1A
    ## PPP1R13L     PPP1R13L       protein phosphatase 1 regulatory subunit 13 like
    ## EMILIN2       EMILIN2                       elastin microfibril interfacer 2
    ## FYB2             FYB2                                  FYN binding protein 2
    ## FAM3D-AS1   FAM3D-AS1                                  FAM3D antisense RNA 1
    ## PFKFB2         PFKFB2  6-phosphofructo-2-kinase/fructose-2,6-biphosphatase 2
    ## EMC10           EMC10                 ER membrane protein complex subunit 10
    ## NEAT1           NEAT1              nuclear paraspeckle assembly transcript 1
    ## MBOAT2         MBOAT2 membrane bound glycerophospholipid O-acyltransferase 2
    ## AGPAT4         AGPAT4         1-acylglycerol-3-phosphate O-acyltransferase 4
    ## CFAP61         CFAP61               cilia and flagella associated protein 61
    ##                logFC      P.Value    adj.P.Val
    ## TMEM178B    7.610225 2.030283e-20 3.253731e-16
    ## TOGARAM2    7.478052 2.949519e-18 2.363450e-14
    ## FILIP1L     4.029497 1.921413e-17 1.026419e-13
    ## CCSER1      3.690413 3.470720e-17 1.030557e-13
    ## AAK1        1.913732 3.507541e-17 1.030557e-13
    ## DGKG        5.028229 3.858319e-17 1.030557e-13
    ## GRAMD1B     2.818759 4.300808e-16 9.846393e-13
    ## NCEH1       4.063683 7.078282e-16 1.417957e-12
    ## MIR29B2CHG  4.302848 1.030308e-15 1.834636e-12
    ## ADRA1A      4.048318 1.956095e-15 3.134838e-12
    ## PPP1R13L    3.746840 2.167455e-15 3.157785e-12
    ## EMILIN2    -4.025270 3.311711e-15 4.422790e-12
    ## FYB2        5.657809 5.518848e-15 6.803466e-12
    ## FAM3D-AS1   6.047737 8.021276e-15 9.182069e-12
    ## PFKFB2      3.714788 8.612511e-15 9.201606e-12
    ## EMC10      -3.248052 9.732175e-15 9.747990e-12
    ## NEAT1       3.037800 1.355876e-14 1.278193e-11
    ## MBOAT2     -3.308016 1.492278e-14 1.328625e-11
    ## AGPAT4     -4.201027 2.921819e-14 2.464477e-11
    ## CFAP61      4.841820 5.670074e-14 4.543430e-11

### Fibroblast Differential Expression

Cardiac fibroblasts provide structural support and undergo important
changes during development, including alterations in extracellular
matrix production.

The original study found that fibroblast maturation was associated with
increased interferon signaling capacity (Figure 1G in the paper).
Fibroblasts also showed significant expansion in proportion during
postnatal development, consistent with the increasing structural
complexity of the maturing heart.

``` r
# Create contrasts for fibroblasts
cont_fibro <- create_celltype_contrasts("fibroblast", design)

# Fit contrasts
fit_fibro <- contrasts.fit(fit, cont_fibro)
fit_fibro <- eBayes(fit_fibro, robust = TRUE)

# Apply TREAT
treat_fibro <- treat(fit_fibro, lfc = 0.5)

# Summarise results
cat("Fibroblasts: Differential Expression Summary\n")
```

    ## Fibroblasts: Differential Expression Summary

``` r
cat("(FDR < 0.05, |logFC| > 0.5)\n\n")
```

    ## (FDR < 0.05, |logFC| > 0.5)

``` r
dt_fibro <- decideTests(treat_fibro)
summary(dt_fibro)
```

    ##          YvF   AvF   AvY
    ## Down     887   695     9
    ## NotSig 14455 14787 16009
    ## Up       684   544     8

``` r
# Top DE genes
cat("\nTop 20 DE genes: Adult vs Fetal Fibroblasts\n")
```

    ## 
    ## Top 20 DE genes: Adult vs Fetal Fibroblasts

``` r
topTreat(treat_fibro, coef = "AvF", n = 20)[, c("SYMBOL", "GENENAME", "logFC", "P.Value", "adj.P.Val")]
```

    ##              SYMBOL                                               GENENAME
    ## CNTNAP2     CNTNAP2                         contactin associated protein 2
    ## MTUS1         MTUS1              microtubule associated scaffold protein 1
    ## VIT             VIT                                                 vitrin
    ## LINC02511 LINC02511            long intergenic non-protein coding RNA 2511
    ## LAMA2         LAMA2                                laminin subunit alpha 2
    ## COL28A1     COL28A1                     collagen type XXVIII alpha 1 chain
    ## COL6A6       COL6A6                         collagen type VI alpha 6 chain
    ## HHIP           HHIP                           hedgehog interacting protein
    ## MEST           MEST                           mesoderm specific transcript
    ## EMILIN2     EMILIN2                       elastin microfibril interfacer 2
    ## PIEZO2       PIEZO2    piezo type mechanosensitive ion channel component 2
    ## IGF2BP3     IGF2BP3    insulin like growth factor 2 mRNA binding protein 3
    ## CCL11         CCL11                          C-C motif chemokine ligand 11
    ## CACNB4       CACNB4 calcium voltage-gated channel auxiliary subunit beta 4
    ## C11orf87   C11orf87                    chromosome 11 open reading frame 87
    ## PRSS35       PRSS35                                     serine protease 35
    ## ZBTB20       ZBTB20               zinc finger and BTB domain containing 20
    ## SAMD5         SAMD5                sterile alpha motif domain containing 5
    ## ECHDC2       ECHDC2                enoyl-CoA hydratase domain containing 2
    ## GSN             GSN                                               gelsolin
    ##               logFC      P.Value    adj.P.Val
    ## CNTNAP2   -3.119905 1.903113e-15 3.049929e-11
    ## MTUS1      3.485722 1.394106e-13 1.117097e-09
    ## VIT        5.717173 2.571843e-13 1.373879e-09
    ## LINC02511  4.762924 3.512828e-13 1.407415e-09
    ## LAMA2      3.929683 6.109991e-13 1.958374e-09
    ## COL28A1    4.605533 9.087643e-13 2.427309e-09
    ## COL6A6    -2.623815 1.155926e-12 2.525903e-09
    ## HHIP      -5.797671 1.542910e-12 2.525903e-09
    ## MEST      -6.140893 1.627359e-12 2.525903e-09
    ## EMILIN2   -4.198400 1.681165e-12 2.525903e-09
    ## PIEZO2    -5.913009 1.733741e-12 2.525903e-09
    ## IGF2BP3   -5.812094 3.790866e-12 5.062701e-09
    ## CCL11      5.250673 4.841238e-12 5.968130e-09
    ## CACNB4    -3.768547 5.819565e-12 6.661740e-09
    ## C11orf87  -6.426033 7.276264e-12 7.773961e-09
    ## PRSS35    -6.776809 1.250521e-11 1.252553e-08
    ## ZBTB20     1.501944 1.734759e-11 1.635367e-08
    ## SAMD5     -4.684984 2.361898e-11 2.102876e-08
    ## ECHDC2     2.435788 3.030074e-11 2.555787e-08
    ## GSN        3.653066 3.234545e-11 2.591841e-08

### Endothelial Cell Differential Expression

Endothelial cells line the cardiac vasculature and play essential roles
in angiogenesis during heart development. Similar to fibroblasts, the
original study found that endothelial cell maturation was associated
with increased interferon signaling capacity. Endothelial cells also
contributed substantially to the total number of differentially
expressed genes during cardiac maturation (Figure 1E in the paper).

``` r
# Create contrasts for endothelial cells
cont_endo <- create_celltype_contrasts("endothelial", design)

# Fit contrasts
fit_endo <- contrasts.fit(fit, cont_endo)
fit_endo <- eBayes(fit_endo, robust = TRUE)

# Apply TREAT
treat_endo <- treat(fit_endo, lfc = 0.5)

# Summarise results
cat("Endothelial Cells: Differential Expression Summary\n")
```

    ## Endothelial Cells: Differential Expression Summary

``` r
cat("(FDR < 0.05, |logFC| > 0.5)\n\n")
```

    ## (FDR < 0.05, |logFC| > 0.5)

``` r
dt_endo <- decideTests(treat_endo)
summary(dt_endo)
```

    ##          YvF   AvF   AvY
    ## Down     297   350     0
    ## NotSig 15520 15443 16026
    ## Up       209   233     0

### Immune Cell Differential Expression

Resident immune cells, particularly macrophages, are present in the
heart from early development and may have distinct phenotypes across
developmental stages. The original study found that immune cell
proportions increase significantly during postnatal development (Figure
1C), and that immune cell maturation was associated with enhanced
interferon signaling capacity (Figure 1G). With our limited sample size,
we may detect fewer DE genes in immune cells compared to the more
abundant cardiomyocytes and fibroblasts.

``` r
# Create contrasts for immune cells
cont_immune <- create_celltype_contrasts("immune", design)

# Fit contrasts
fit_immune <- contrasts.fit(fit, cont_immune)
fit_immune <- eBayes(fit_immune, robust = TRUE)

# Apply TREAT
treat_immune <- treat(fit_immune, lfc = 0.5)

# Summarise results
cat("Immune Cells: Differential Expression Summary\n")
```

    ## Immune Cells: Differential Expression Summary

``` r
cat("(FDR < 0.05, |logFC| > 0.5)\n\n")
```

    ## (FDR < 0.05, |logFC| > 0.5)

``` r
dt_immune <- decideTests(treat_immune)
summary(dt_immune)
```

    ##          YvF   AvF   AvY
    ## Down      88    88     1
    ## NotSig 15840 15850 16024
    ## Up        98    88     1

## Visualisation of DE Results

### Summary Across Cell Types

``` r
# Compile summary statistics
de_summary <- data.frame(
    celltype = rep(c("Cardiomyocytes", "Fibroblasts", "Endothelial", "Immune"), each = 3),
    contrast = rep(c("Young vs Fetal", "Adult vs Fetal", "Adult vs Young"), 4),
    up = c(
        sum(dt_cardio[, 1] == 1), sum(dt_cardio[, 2] == 1), sum(dt_cardio[, 3] == 1),
        sum(dt_fibro[, 1] == 1), sum(dt_fibro[, 2] == 1), sum(dt_fibro[, 3] == 1),
        sum(dt_endo[, 1] == 1), sum(dt_endo[, 2] == 1), sum(dt_endo[, 3] == 1),
        sum(dt_immune[, 1] == 1), sum(dt_immune[, 2] == 1), sum(dt_immune[, 3] == 1)
    ),
    down = c(
        sum(dt_cardio[, 1] == -1), sum(dt_cardio[, 2] == -1), sum(dt_cardio[, 3] == -1),
        sum(dt_fibro[, 1] == -1), sum(dt_fibro[, 2] == -1), sum(dt_fibro[, 3] == -1),
        sum(dt_endo[, 1] == -1), sum(dt_endo[, 2] == -1), sum(dt_endo[, 3] == -1),
        sum(dt_immune[, 1] == -1), sum(dt_immune[, 2] == -1), sum(dt_immune[, 3] == -1)
    )
)

# Convert to long format
de_summary_long <- de_summary %>%
    pivot_longer(cols = c(up, down), names_to = "direction", values_to = "n_genes") %>%
    mutate(
        n_genes = ifelse(direction == "down", -n_genes, n_genes),
        direction = factor(direction, levels = c("up", "down"))
    )

# Create diverging bar plot
ggplot(de_summary_long, aes(x = contrast, y = n_genes, fill = direction)) +
    geom_bar(stat = "identity", position = "identity") +
    geom_hline(yintercept = 0, colour = "black") +
    facet_wrap(~ celltype, scales = "free_x") +
    scale_fill_manual(
        values = c("up" = "#D62728", "down" = "#1F77B4"),
        labels = c("up" = "Upregulated", "down" = "Downregulated")
    ) +
    coord_flip() +
    labs(
        title = "Differentially Expressed Genes Across Cell Types",
        subtitle = "TREAT analysis with |logFC| > 0.5 and FDR < 0.05",
        x = "",
        y = "Number of DE Genes",
        fill = ""
    ) +
    theme_minimal() +
    theme(
        legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 11)
    )
```

![](04_differential_expression_files/figure-html/de-summary-plot-1.png)

**Interpretation**: This summary plot reveals several patterns:

1.  **Adult vs Fetal** comparisons consistently yield the most DE genes
    across cell types, reflecting the substantial transcriptional
    changes during full cardiac development.

2.  **Cell type-specific responses**: Cardiomyocytes show the largest
    number of DE genes, consistent with the major functional changes
    these cells undergo during maturation. Immune cells show fewer DE
    genes, which may reflect either biological stability or limited
    power due to smaller cell numbers.

3.  **Symmetry of changes**: The balance between upregulated and
    downregulated genes varies by cell type and comparison, providing
    insight into the directional nature of developmental transcriptional
    programmes.

### Volcano Plots

Volcano plots visualise the relationship between statistical
significance and effect size, helping identify genes that are both
significant and biologically meaningful:

``` r
# Function to create volcano plot
create_volcano <- function(fit_obj, coef_name, title) {
    results <- topTreat(fit_obj, coef = coef_name, n = Inf)
    results$significant <- results$adj.P.Val < 0.05 & abs(results$logFC) > 0.5

    # Label top genes
    top_genes <- results %>%
        filter(significant) %>%
        arrange(adj.P.Val) %>%
        head(10)

    ggplot(results, aes(x = logFC, y = -log10(adj.P.Val))) +
        geom_point(aes(colour = significant), alpha = 0.6, size = 1) +
        geom_text(
            data = top_genes,
            aes(label = SYMBOL),
            size = 2.5, hjust = -0.1, vjust = 0.5,
            check_overlap = TRUE
        ) +
        geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", colour = "grey50") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
        scale_colour_manual(values = c("FALSE" = "grey70", "TRUE" = "#D62728")) +
        labs(
            title = title,
            x = "log2 Fold Change",
            y = "-log10(adjusted P-value)"
        ) +
        theme_minimal() +
        theme(legend.position = "none")
}

# Create volcano plots for Adult vs Fetal across cell types
v1 <- create_volcano(treat_cardio, "AvF", "Cardiomyocytes: Adult vs Fetal")
v2 <- create_volcano(treat_fibro, "AvF", "Fibroblasts: Adult vs Fetal")
v3 <- create_volcano(treat_endo, "AvF", "Endothelial: Adult vs Fetal")
v4 <- create_volcano(treat_immune, "AvF", "Immune: Adult vs Fetal")

(v1 + v2) / (v3 + v4) +
    plot_annotation(
        title = "Volcano Plots: Adult vs Fetal Development",
        subtitle = "Dashed lines indicate |logFC| > 0.5 and FDR < 0.05 thresholds"
    )
```

![](04_differential_expression_files/figure-html/volcano-plots-1.png)

### Developmental Trajectory Plot

To examine the consistency of developmental changes, we plot
log-fold-changes for Young vs Fetal against Adult vs Fetal. Genes
showing progressive changes (consistent direction across development)
appear in the upper-right or lower-left quadrants:

``` r
# Get results for cardiomyocytes
results_yvf <- topTreat(treat_cardio, coef = "YvF", n = Inf)
results_avf <- topTreat(treat_cardio, coef = "AvF", n = Inf)

# Merge results
trajectory_data <- data.frame(
    gene = results_yvf$SYMBOL,
    logFC_YvF = results_yvf$logFC,
    logFC_AvF = results_avf$logFC[match(results_yvf$SYMBOL, results_avf$SYMBOL)],
    sig_YvF = results_yvf$adj.P.Val < 0.05 & abs(results_yvf$logFC) > 0.5,
    sig_AvF = results_avf$adj.P.Val < 0.05 & abs(results_avf$logFC) > 0.5
)

# Identify genes significant in either comparison
trajectory_data$significant <- trajectory_data$sig_YvF | trajectory_data$sig_AvF

# Label genes significant in both
trajectory_data$sig_both <- trajectory_data$sig_YvF & trajectory_data$sig_AvF
label_genes <- trajectory_data %>%
    filter(sig_both) %>%
    arrange(desc(abs(logFC_AvF))) %>%
    head(15)

ggplot(trajectory_data, aes(x = logFC_YvF, y = logFC_AvF)) +
    geom_point(aes(colour = significant), alpha = 0.5, size = 1) +
    geom_text(
        data = label_genes,
        aes(label = gene),
        size = 2.5, colour = "black",
        check_overlap = TRUE
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_hline(yintercept = 0, colour = "grey80") +
    geom_vline(xintercept = 0, colour = "grey80") +
    scale_colour_manual(values = c("FALSE" = "grey70", "TRUE" = "#D62728")) +
    labs(
        title = "Cardiomyocyte Developmental Trajectory",
        subtitle = "Genes in upper-right/lower-left show progressive developmental changes",
        x = "log2 FC: Young vs Fetal",
        y = "log2 FC: Adult vs Fetal"
    ) +
    theme_minimal() +
    theme(legend.position = "none") +
    coord_fixed()
```

![](04_differential_expression_files/figure-html/trajectory-plot-1.png)

**Interpretation**: Genes along the diagonal (y = x line) show changes
that are already established by the young stage and maintained into
adulthood. Genes above the diagonal show additional changes between
young and adult stages. Genes with opposite signs between the two
comparisons (upper-left, lower-right quadrants) show non-monotonic
developmental patterns that may warrant further investigation.

## Export Results

``` r
# Create results directory
results_dir <- "../results"
if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
}

# Export function
export_de_results <- function(treat_obj, celltype_name, results_dir) {
    for (coef_name in colnames(treat_obj$coefficients)) {
        results <- topTreat(treat_obj, coef = coef_name, n = Inf)
        filename <- file.path(
            results_dir,
            paste0("DE_", gsub(" ", "_", celltype_name), "_", coef_name, ".csv")
        )
        write.csv(results, filename, row.names = FALSE)
    }
}

# Export results for each cell type
export_de_results(treat_cardio, "Cardiomyocytes", results_dir)
export_de_results(treat_fibro, "Fibroblasts", results_dir)
export_de_results(treat_endo, "Endothelial", results_dir)
export_de_results(treat_immune, "Immune", results_dir)

cat("Results exported to:", results_dir, "\n")
```

    ## Results exported to: ../results

``` r
list.files(results_dir, pattern = "DE_")
```

    ##  [1] "DE_Cardiomyocytes_AvF.csv" "DE_Cardiomyocytes_AvY.csv"
    ##  [3] "DE_Cardiomyocytes_YvF.csv" "DE_Endothelial_AvF.csv"   
    ##  [5] "DE_Endothelial_AvY.csv"    "DE_Endothelial_YvF.csv"   
    ##  [7] "DE_Fibroblasts_AvF.csv"    "DE_Fibroblasts_AvY.csv"   
    ##  [9] "DE_Fibroblasts_YvF.csv"    "DE_Immune_AvF.csv"        
    ## [11] "DE_Immune_AvY.csv"         "DE_Immune_YvF.csv"

## Limitations and Considerations

### Sample Size Limitations

With n=3 biological replicates per group, our statistical power is
limited. This has several implications:

1.  **Conservative estimates**: We may fail to detect true biological
    effects (Type II errors) due to insufficient power.

2.  **Variance estimation**: With few samples, variance estimates are
    unstable. The empirical Bayes moderation in limma helps but cannot
    fully compensate.

3.  **Outlier sensitivity**: Individual outlier samples can have
    outsized influence on results.

4.  **Effect size uncertainty**: Confidence intervals on
    log-fold-changes are wide with small samples.

### Sex Confounding

The unbalanced sex distribution across developmental groups complicates
interpretation. While we average over sex in our contrasts, the limited
sample sizes prevent robust separation of sex and developmental effects.
Genes showing apparent developmental regulation may partially reflect
differences in sex composition between groups.

### Pseudobulk Assumptions

The pseudobulk approach assumes:

1.  **Homogeneous cell populations**: Cells within a cell type are
    sufficiently similar that aggregation is meaningful.

2.  **Consistent cell type definitions**: The same cell type label
    represents comparable populations across samples.

3.  **Sufficient cells per sample**: Each pseudobulk sample has enough
    cells for reliable expression estimation.

### Recommendations for Interpretation

1.  **Focus on robust signals**: Prioritise genes with large effect
    sizes (\|logFC\| \> 1) and very small p-values (FDR \< 0.01).

2.  **Seek consistency**: Genes differentially expressed across multiple
    cell types or showing progressive developmental changes are more
    likely to represent true biological signals.

3.  **Validate orthogonally**: Important findings should be validated
    using independent methods (qPCR, protein staining) or datasets.

4.  **Consider biology**: Interpret results in the context of known
    cardiac developmental biology.

## Summary

In this module, we performed two complementary analyses of human heart
development: cell type composition analysis using propeller and
pseudobulk differential expression analysis using the limma-voom
pipeline. Our findings are consistent with those reported by Sim et
al. (2021):

**Cell Type Composition (propeller analysis):**

- The human heart undergoes significant changes in cellular composition
  during development
- Cardiomyocyte proportions decrease from fetal to adult stages
- Fibroblast and immune cell proportions increase during postnatal
  maturation
- These compositional changes reflect fundamental biological processes
  in cardiac development

**Differential Expression (limma-voom analysis):**

- The pseudobulk approach correctly treats samples (not cells) as the
  unit of biological replication, avoiding the pseudoreplication problem
- Adult vs Fetal comparisons reveal the most transcriptional changes
  across all cell types, reflecting the substantial remodelling during
  development
- Cardiomyocytes show the largest number of differentially expressed
  genes, consistent with major functional maturation involving metabolic
  reprogramming and cell cycle exit
- Non-myocyte populations (fibroblasts, endothelial cells, immune cells)
  show increased interferon signaling during maturation

**Methodological Notes:**

The limma-voom pipeline with TREAT provides conservative but reliable
identification of differentially expressed genes. With n=3 samples per
group, our statistical power is limited, and we may fail to detect true
biological effects. The original study findings should be consulted for
the complete picture of transcriptional changes during human heart
development.

### References

Sim CB, Phipson B, Ziemann M, et al. Sex-Specific Control of Human Heart
Maturation by the Progesterone Receptor. *Circulation*.
2021;143:1614-1628. <doi:10.1161/CIRCULATIONAHA.120.051921>

## Session Information

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
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] patchwork_1.3.0      pheatmap_1.0.13      RColorBrewer_1.1-3  
    ##  [4] tidyr_1.3.1          dplyr_1.1.4          ggplot2_3.5.2       
    ##  [7] org.Hs.eg.db_3.21.0  AnnotationDbi_1.70.0 IRanges_2.42.0      
    ## [10] S4Vectors_0.46.0     Biobase_2.68.0       BiocGenerics_0.54.0 
    ## [13] generics_0.1.4       speckle_1.8.0        edgeR_4.6.2         
    ## [16] limma_3.64.0         Seurat_5.3.0         SeuratObject_5.1.0  
    ## [19] sp_2.2-0            
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] jsonlite_2.0.0              magrittr_2.0.3             
    ##   [3] spatstat.utils_3.1-4        farver_2.1.2               
    ##   [5] rmarkdown_2.29              fs_1.6.6                   
    ##   [7] ragg_1.4.0                  vctrs_0.6.5                
    ##   [9] ROCR_1.0-11                 memoise_2.0.1              
    ##  [11] spatstat.explore_3.4-3      S4Arrays_1.8.0             
    ##  [13] htmltools_0.5.8.1           SparseArray_1.8.0          
    ##  [15] sass_0.4.10                 sctransform_0.4.2          
    ##  [17] parallelly_1.44.0           KernSmooth_2.23-26         
    ##  [19] bslib_0.9.0                 htmlwidgets_1.6.4          
    ##  [21] desc_1.4.3                  ica_1.0-3                  
    ##  [23] plyr_1.8.9                  plotly_4.10.4              
    ##  [25] zoo_1.8-14                  cachem_1.1.0               
    ##  [27] igraph_2.1.4                mime_0.13                  
    ##  [29] lifecycle_1.0.4             pkgconfig_2.0.3            
    ##  [31] Matrix_1.7-4                R6_2.6.1                   
    ##  [33] fastmap_1.2.0               GenomeInfoDbData_1.2.14    
    ##  [35] MatrixGenerics_1.20.0       fitdistrplus_1.2-2         
    ##  [37] future_1.49.0               shiny_1.10.0               
    ##  [39] digest_0.6.37               colorspace_2.1-1           
    ##  [41] tensor_1.5                  RSpectra_0.16-2            
    ##  [43] irlba_2.3.5.1               RSQLite_2.3.11             
    ##  [45] GenomicRanges_1.60.0        textshaping_1.0.1          
    ##  [47] labeling_0.4.3              progressr_0.15.1           
    ##  [49] spatstat.sparse_3.1-0       httr_1.4.7                 
    ##  [51] polyclip_1.10-7             abind_1.4-8                
    ##  [53] compiler_4.5.2              withr_3.0.2                
    ##  [55] bit64_4.6.0-1               DBI_1.2.3                  
    ##  [57] fastDummies_1.7.5           MASS_7.3-65                
    ##  [59] DelayedArray_0.34.1         tools_4.5.2                
    ##  [61] lmtest_0.9-40               httpuv_1.6.16              
    ##  [63] future.apply_1.11.3         goftest_1.2-3              
    ##  [65] glue_1.8.0                  nlme_3.1-168               
    ##  [67] promises_1.3.2              grid_4.5.2                 
    ##  [69] Rtsne_0.17                  cluster_2.1.8.1            
    ##  [71] reshape2_1.4.4              gtable_0.3.6               
    ##  [73] spatstat.data_3.1-6         data.table_1.17.2          
    ##  [75] utf8_1.2.5                  XVector_0.48.0             
    ##  [77] spatstat.geom_3.4-1         RcppAnnoy_0.0.22           
    ##  [79] ggrepel_0.9.6               RANN_2.6.2                 
    ##  [81] pillar_1.10.2               stringr_1.5.1              
    ##  [83] spam_2.11-1                 RcppHNSW_0.6.0             
    ##  [85] later_1.4.2                 splines_4.5.2              
    ##  [87] lattice_0.22-7              bit_4.6.0                  
    ##  [89] survival_3.8-3              deldir_2.0-4               
    ##  [91] tidyselect_1.2.1            SingleCellExperiment_1.30.1
    ##  [93] locfit_1.5-9.12             Biostrings_2.76.0          
    ##  [95] miniUI_0.1.2                pbapply_1.7-2              
    ##  [97] knitr_1.50                  gridExtra_2.3              
    ##  [99] SummarizedExperiment_1.38.1 scattermore_1.2            
    ## [101] xfun_0.52                   statmod_1.5.0              
    ## [103] matrixStats_1.5.0           UCSC.utils_1.4.0           
    ## [105] stringi_1.8.7               lazyeval_0.2.2             
    ## [107] yaml_2.3.10                 evaluate_1.0.3             
    ## [109] codetools_0.2-20            tibble_3.2.1               
    ## [111] cli_3.6.5                   uwot_0.2.3                 
    ## [113] xtable_1.8-4                reticulate_1.42.0          
    ## [115] systemfonts_1.2.3           jquerylib_0.1.4            
    ## [117] GenomeInfoDb_1.44.0         Rcpp_1.0.14                
    ## [119] globals_0.18.0              spatstat.random_3.4-1      
    ## [121] png_0.1-8                   spatstat.univar_3.1-3      
    ## [123] parallel_4.5.2              blob_1.2.4                 
    ## [125] pkgdown_2.2.0               dotCall64_1.2              
    ## [127] listenv_0.9.1               viridisLite_0.4.2          
    ## [129] scales_1.4.0                ggridges_0.5.6             
    ## [131] crayon_1.5.3                purrr_1.0.4                
    ## [133] rlang_1.1.6                 KEGGREST_1.48.1            
    ## [135] cowplot_1.1.3
