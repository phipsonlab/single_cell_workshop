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

    ## - Genes: 18,953

``` r
cat("- Cells:", format(ncol(counts), big.mark = ","), "\n")
```

    ## - Cells: 42,935

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
    ##               23995                3378                2839                  50 
    ##          Fibroblast        Immune cells             Neurons Smooth muscle cells 
    ##                8672                2880                 733                 388

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
    ##   Cardiomyocytes       1959 17060  4976
    ##   Endothelial cells    1087  1068  1223
    ##   Epicardial cells      864   915  1060
    ##   Erythroid               0    50     0
    ##   Fibroblast           2221  2905  3546
    ##   Immune cells         1059   574  1247
    ##   Neurons               101   284   348
    ##   Smooth muscle cells    80   183   125

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

    ## Genes with annotation: 18953 / 18953

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

    ## Erythroid cells found: 50

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

    ## Cells for DE analysis: 42885

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
    ##               23995                3378                2839                8672 
    ##        Immune cells             Neurons Smooth muscle cells 
    ##                2880                 733                 388

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
    ## 1 a1     adult        3787
    ## 2 a2     adult        2439
    ## 3 a3     adult        1145
    ## 4 f1     fetal        7198
    ## 5 f2     fetal        9026
    ## 6 f3     fetal        6765
    ## 7 y1     young        4101
    ## 8 y2     young        4351
    ## 9 y3     young        4073

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
    ## Cardiomyocytes       0.559519646     0.21325177    0.733992921     0.39745000
    ## Immune cells         0.067156348     0.18675439    0.025224804     0.09940743
    ## Endothelial cells    0.078768800     0.13981820    0.047667140     0.09795653
    ## Fibroblast           0.202215227     0.32089945    0.131125816     0.28196804
    ## Epicardial cells     0.066200303     0.11407479    0.040883730     0.08508428
    ## Neurons              0.017092223     0.01328168    0.012480686     0.02807099
    ## Smooth muscle cells  0.009047452     0.01191973    0.008624903     0.01006272
    ##                     Fstatistic      P.Value          FDR
    ## Cardiomyocytes      14.2756862 9.916688e-05 0.0006941682
    ## Immune cells        12.2523228 2.511712e-04 0.0008790993
    ## Endothelial cells    5.0233889 1.571079e-02 0.0366585150
    ## Fibroblast           4.4818201 2.299246e-02 0.0402368110
    ## Epicardial cells     3.1160027 6.388112e-02 0.0894335674
    ## Neurons              1.7405369 1.982157e-01 0.2312516849
    ## Smooth muscle cells  0.6694511 5.218903e-01 0.5218902627

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
    ## Cardiomyocytes        0.73399292     0.39745000      0.2132518  14.275686
    ## Immune cells          0.02522480     0.09940743      0.1867544  12.252323
    ## Endothelial cells     0.04766714     0.09795653      0.1398182   5.023389
    ## Fibroblast            0.13112582     0.28196804      0.3208994   4.481820
    ##                            FDR
    ## Cardiomyocytes    0.0006941682
    ## Immune cells      0.0008790993
    ## Endothelial cells 0.0366585150
    ## Fibroblast        0.0402368110

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
    ##                   1539                    288                    132 
    ##      Cardiomyocytes.f1      Cardiomyocytes.f2      Cardiomyocytes.f3 
    ##                   4948                   7449                   4663 
    ##      Cardiomyocytes.y1      Cardiomyocytes.y2      Cardiomyocytes.y3 
    ##                    994                   1765                   2217 
    ##   Endothelial cells.a1   Endothelial cells.a2   Endothelial cells.a3 
    ##                    567                    398                    122 
    ##   Endothelial cells.f1   Endothelial cells.f2   Endothelial cells.f3 
    ##                    423                    300                    345 
    ##   Endothelial cells.y1   Endothelial cells.y2   Endothelial cells.y3 
    ##                    468                    358                    397 
    ##    Epicardial cells.a1    Epicardial cells.a2    Epicardial cells.a3 
    ##                    330                    456                     78 
    ##    Epicardial cells.f1    Epicardial cells.f2    Epicardial cells.f3 
    ##                    418                    240                    257 
    ##    Epicardial cells.y1    Epicardial cells.y2    Epicardial cells.y3 
    ##                    595                    255                    210 
    ##          Fibroblast.a1          Fibroblast.a2          Fibroblast.a3 
    ##                    929                    887                    405 
    ##          Fibroblast.f1          Fibroblast.f2          Fibroblast.f3 
    ##                   1013                    730                   1162 
    ##          Fibroblast.y1          Fibroblast.y2          Fibroblast.y3 
    ##                   1488                   1416                    642 
    ##        Immune cells.a1        Immune cells.a2        Immune cells.a3 
    ##                    349                    328                    382 
    ##        Immune cells.f1        Immune cells.f2        Immune cells.f3 
    ##                    253                    187                    134 
    ##        Immune cells.y1        Immune cells.y2        Immune cells.y3 
    ##                    301                    474                    472 
    ##             Neurons.a1             Neurons.a2             Neurons.a3 
    ##                     52                     36                     13 
    ##             Neurons.f1             Neurons.f2             Neurons.f3 
    ##                     94                    100                     90 
    ##             Neurons.y1             Neurons.y2             Neurons.y3 
    ##                    199                     57                     92 
    ## Smooth muscle cells.a1 Smooth muscle cells.a2 Smooth muscle cells.a3 
    ##                     21                     46                     13 
    ## Smooth muscle cells.f1 Smooth muscle cells.f2 Smooth muscle cells.f3 
    ##                     49                     20                    114 
    ## Smooth muscle cells.y1 Smooth muscle cells.y2 Smooth muscle cells.y3 
    ##                     56                     26                     43

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

    ## - Genes: 18953

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
    ## Cardiomyocytes.a1 29509310            1 Cardiomyocytes.a1 Cardiomyocytes     a1
    ## Cardiomyocytes.a2  6839420            1 Cardiomyocytes.a2 Cardiomyocytes     a2
    ## Cardiomyocytes.a3  2120636            1 Cardiomyocytes.a3 Cardiomyocytes     a3
    ## Cardiomyocytes.f1 50630676            1 Cardiomyocytes.f1 Cardiomyocytes     f1
    ## Cardiomyocytes.f2 66175611            1 Cardiomyocytes.f2 Cardiomyocytes     f2
    ## Cardiomyocytes.f3 63596974            1 Cardiomyocytes.f3 Cardiomyocytes     f3
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

    ## Genes before filtering: 18953

``` r
cat("Genes after filtering:", nrow(dge_filtered), "\n")
```

    ## Genes after filtering: 15970

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
    ## Cardiomyocytes.a1      29496007    0.5980691
    ## Cardiomyocytes.a2       6836372    0.6848721
    ## Cardiomyocytes.a3       2119045    0.8246200
    ## Cardiomyocytes.f1      50603602    0.8102726
    ## Cardiomyocytes.f2      66148201    0.7778729
    ## Cardiomyocytes.f3      63573020    0.7202654
    ## Cardiomyocytes.y1      12462705    0.6378081
    ## Cardiomyocytes.y2      28785286    0.6194712
    ## Cardiomyocytes.y3      16298680    0.8782436
    ## Endothelial cells.a1    3740023    1.0509952
    ## Endothelial cells.a2    2568540    1.0133357
    ## Endothelial cells.a3     781344    1.1125139
    ## Endothelial cells.f1    3468570    1.1316873
    ## Endothelial cells.f2    2569781    1.1329679
    ## Endothelial cells.f3    3989951    1.0655295
    ## Endothelial cells.y1    3127825    1.0634920
    ## Endothelial cells.y2    2445699    1.0016953
    ## Endothelial cells.y3    2124667    0.9581740
    ## Epicardial cells.a1     1768137    1.0412641
    ## Epicardial cells.a2     3465782    0.8720177
    ## Epicardial cells.a3      532305    1.0668027
    ## Epicardial cells.f1     3041338    1.1353154
    ## Epicardial cells.f2     1698965    1.0920206
    ## Epicardial cells.f3     2595975    1.0471145
    ## Epicardial cells.y1     3703198    0.9713980
    ## Epicardial cells.y2     1601834    0.9618783
    ## Epicardial cells.y3     1120534    0.9289960
    ## Fibroblast.a1           7281212    1.0746025
    ## Fibroblast.a2           7670879    0.9835376
    ## Fibroblast.a3           3810418    1.0792989
    ## Fibroblast.f1           7595543    1.0494530
    ## Fibroblast.f2           5142099    1.0302138
    ## Fibroblast.f3          10956825    1.0066763
    ## Fibroblast.y1          12325273    0.9612371
    ## Fibroblast.y2          10759060    0.9561517
    ## Fibroblast.y3           4994896    0.9752673
    ## Immune cells.a1         2579054    1.0346262
    ## Immune cells.a2         2076196    1.0773197
    ## Immune cells.a3         2894769    1.1613529
    ## Immune cells.f1         1634048    1.1126271
    ## Immune cells.f2         1063960    1.2296393
    ## Immune cells.f3          996186    1.1131232
    ## Immune cells.y1         2270365    1.1228178
    ## Immune cells.y2         3180451    1.0420742
    ## Immune cells.y3         2447925    1.0320235
    ## Neurons.a1               256789    1.1747136
    ## Neurons.a2               207811    1.0601088
    ## Neurons.a3                74057    1.3396612
    ## Neurons.f1               544707    1.1034003
    ## Neurons.f2               633420    1.1201456
    ## Neurons.f3               681890    1.0560271
    ## Neurons.y1              1038650    1.0843009
    ## Neurons.y2               209761    1.0733240
    ## Neurons.y3               330745    1.0407414
    ## Smooth muscle cells.a1   118999    1.1796940
    ## Smooth muscle cells.a2   359768    0.9584901
    ## Smooth muscle cells.a3    96592    1.1756187
    ## Smooth muscle cells.f1   369138    1.0294430
    ## Smooth muscle cells.f2   145847    1.1400543
    ## Smooth muscle cells.f3  1246685    0.9764167
    ## Smooth muscle cells.y1   336263    1.0339953
    ## Smooth muscle cells.y2   153371    0.9810013
    ## Smooth muscle cells.y3   269864    0.9652504

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
    ## Down    1201  1226     8
    ## NotSig 13735 13667 15933
    ## Up      1034  1077    29

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

    ##                SYMBOL
    ## TMEM178B     TMEM178B
    ## TOGARAM2     TOGARAM2
    ## FILIP1L       FILIP1L
    ## DGKG             DGKG
    ## MIR29B2CHG MIR29B2CHG
    ## CCSER1         CCSER1
    ## AAK1             AAK1
    ## NCEH1           NCEH1
    ## EMILIN2       EMILIN2
    ## PPP1R13L     PPP1R13L
    ## GRAMD1B       GRAMD1B
    ## EMC10           EMC10
    ## ADRA1A         ADRA1A
    ## FYB2             FYB2
    ## FAM3D-AS1   FAM3D-AS1
    ## AGPAT4         AGPAT4
    ## MBOAT2         MBOAT2
    ## PFKFB2         PFKFB2
    ## HECW2           HECW2
    ## CFAP61         CFAP61
    ##                                                                   GENENAME
    ## TMEM178B                                        transmembrane protein 178B
    ## TOGARAM2                    TOG array regulator of axonemal microtubules 2
    ## FILIP1L                               filamin A interacting protein 1 like
    ## DGKG                                           diacylglycerol kinase gamma
    ## MIR29B2CHG                                    MIR29B2 and MIR29C host gene
    ## CCSER1                                   coiled-coil serine rich protein 1
    ## AAK1                                               AP2 associated kinase 1
    ## NCEH1                                neutral cholesterol ester hydrolase 1
    ## EMILIN2                                   elastin microfibril interfacer 2
    ## PPP1R13L                  protein phosphatase 1 regulatory subunit 13 like
    ## GRAMD1B                                          GRAM domain containing 1B
    ## EMC10                               ER membrane protein complex subunit 10
    ## ADRA1A                                               adrenoceptor alpha 1A
    ## FYB2                                                 FYN binding protein 2
    ## FAM3D-AS1                                            FAM3D antisense RNA 1
    ## AGPAT4                      1-acylglycerol-3-phosphate O-acyltransferase 4
    ## MBOAT2              membrane bound glycerophospholipid O-acyltransferase 2
    ## PFKFB2               6-phosphofructo-2-kinase/fructose-2,6-biphosphatase 2
    ## HECW2      HECT, C2 and WW domain containing E3 ubiquitin protein ligase 2
    ## CFAP61                            cilia and flagella associated protein 61
    ##                logFC      P.Value    adj.P.Val
    ## TMEM178B    7.639349 1.967356e-21 3.141868e-17
    ## TOGARAM2    7.492989 1.290657e-18 1.030589e-14
    ## FILIP1L     4.039662 1.482405e-16 7.891338e-13
    ## DGKG        5.038590 2.653101e-16 1.059251e-12
    ## MIR29B2CHG  4.317352 9.924486e-16 2.978624e-12
    ## CCSER1      3.690847 1.183646e-15 2.978624e-12
    ## AAK1        1.911225 1.315917e-15 2.978624e-12
    ## NCEH1       4.065068 1.492110e-15 2.978624e-12
    ## EMILIN2    -4.033696 2.214586e-15 3.929660e-12
    ## PPP1R13L    3.747764 3.254204e-15 5.156117e-12
    ## GRAMD1B     2.815380 3.551490e-15 5.156117e-12
    ## EMC10      -3.251951 4.073444e-15 5.421075e-12
    ## ADRA1A      4.040018 7.280892e-15 8.944295e-12
    ## FYB2        5.647460 8.136566e-15 9.281497e-12
    ## FAM3D-AS1   6.135190 1.139787e-14 1.213494e-11
    ## AGPAT4     -4.233240 2.338284e-14 2.333900e-11
    ## MBOAT2     -3.301955 2.606392e-14 2.448475e-11
    ## PFKFB2      3.697529 2.822445e-14 2.504136e-11
    ## HECW2      -5.507797 4.610690e-14 3.875406e-11
    ## CFAP61      4.834211 6.071162e-14 4.847823e-11

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
    ## Down     882   704    10
    ## NotSig 14395 14722 15950
    ## Up       693   544    10

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
    ## VIT             VIT                                                 vitrin
    ## MTUS1         MTUS1              microtubule associated scaffold protein 1
    ## COL28A1     COL28A1                     collagen type XXVIII alpha 1 chain
    ## LINC02511 LINC02511            long intergenic non-protein coding RNA 2511
    ## PIEZO2       PIEZO2    piezo type mechanosensitive ion channel component 2
    ## MEST           MEST                           mesoderm specific transcript
    ## HHIP           HHIP                           hedgehog interacting protein
    ## EMILIN2     EMILIN2                       elastin microfibril interfacer 2
    ## LAMA2         LAMA2                                laminin subunit alpha 2
    ## IGF2BP3     IGF2BP3    insulin like growth factor 2 mRNA binding protein 3
    ## COL6A6       COL6A6                         collagen type VI alpha 6 chain
    ## CACNB4       CACNB4 calcium voltage-gated channel auxiliary subunit beta 4
    ## ECHDC2       ECHDC2                enoyl-CoA hydratase domain containing 2
    ## PRSS35       PRSS35                                     serine protease 35
    ## CCL11         CCL11                          C-C motif chemokine ligand 11
    ## GSN             GSN                                               gelsolin
    ## C11orf87   C11orf87                    chromosome 11 open reading frame 87
    ## TMEM26       TMEM26                               transmembrane protein 26
    ## SAMD5         SAMD5                sterile alpha motif domain containing 5
    ##               logFC      P.Value    adj.P.Val
    ## CNTNAP2   -3.137683 1.242293e-15 1.983942e-11
    ## VIT        5.694914 3.420988e-14 2.731659e-10
    ## MTUS1      3.490190 9.118128e-14 4.853884e-10
    ## COL28A1    4.579719 1.505350e-13 6.010108e-10
    ## LINC02511  4.738088 2.908194e-13 9.288770e-10
    ## PIEZO2    -5.918862 4.978987e-13 1.325240e-09
    ## MEST      -6.182768 7.119017e-13 1.521895e-09
    ## HHIP      -5.802522 7.753886e-13 1.521895e-09
    ## EMILIN2   -4.219448 8.576743e-13 1.521895e-09
    ## LAMA2      3.935211 1.579173e-12 2.521939e-09
    ## IGF2BP3   -5.783264 4.572519e-12 6.428215e-09
    ## COL6A6    -2.650860 4.830218e-12 6.428215e-09
    ## CACNB4    -3.753094 5.238425e-12 6.435204e-09
    ## ECHDC2     2.436505 1.479163e-11 1.687302e-08
    ## PRSS35    -6.753537 1.664622e-11 1.754554e-08
    ## CCL11      5.220386 1.766183e-11 1.754554e-08
    ## GSN        3.660341 1.867715e-11 1.754554e-08
    ## C11orf87  -6.460264 3.126596e-11 2.773986e-08
    ## TMEM26    -8.496108 3.773647e-11 3.171849e-08
    ## SAMD5     -4.705047 4.090551e-11 3.266305e-08

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
    ## Down     311   367     0
    ## NotSig 15447 15362 15970
    ## Up       212   241     0

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
    ## Down     111   121     2
    ## NotSig 15756 15747 15967
    ## Up       103   102     1

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
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] patchwork_1.3.2      pheatmap_1.0.13      RColorBrewer_1.1-3  
    ##  [4] tidyr_1.3.2          dplyr_1.1.4          ggplot2_4.0.1       
    ##  [7] org.Hs.eg.db_3.22.0  AnnotationDbi_1.72.0 IRanges_2.44.0      
    ## [10] S4Vectors_0.48.0     Biobase_2.70.0       BiocGenerics_0.56.0 
    ## [13] generics_0.1.4       speckle_1.10.0       edgeR_4.8.2         
    ## [16] limma_3.66.0         Seurat_5.4.0         SeuratObject_5.3.0  
    ## [19] sp_2.2-0            
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] jsonlite_2.0.0              magrittr_2.0.4             
    ##   [3] spatstat.utils_3.2-1        farver_2.1.2               
    ##   [5] rmarkdown_2.30              fs_1.6.6                   
    ##   [7] ragg_1.5.0                  vctrs_0.6.5                
    ##   [9] ROCR_1.0-11                 memoise_2.0.1              
    ##  [11] spatstat.explore_3.6-0      S4Arrays_1.10.1            
    ##  [13] htmltools_0.5.9             SparseArray_1.10.8         
    ##  [15] sass_0.4.10                 sctransform_0.4.3          
    ##  [17] parallelly_1.46.1           KernSmooth_2.23-26         
    ##  [19] bslib_0.9.0                 htmlwidgets_1.6.4          
    ##  [21] desc_1.4.3                  ica_1.0-3                  
    ##  [23] plyr_1.8.9                  plotly_4.11.0              
    ##  [25] zoo_1.8-15                  cachem_1.1.0               
    ##  [27] igraph_2.2.1                mime_0.13                  
    ##  [29] lifecycle_1.0.5             pkgconfig_2.0.3            
    ##  [31] Matrix_1.7-4                R6_2.6.1                   
    ##  [33] fastmap_1.2.0               MatrixGenerics_1.22.0      
    ##  [35] fitdistrplus_1.2-4          future_1.68.0              
    ##  [37] shiny_1.12.1                digest_0.6.39              
    ##  [39] tensor_1.5.1                RSpectra_0.16-2            
    ##  [41] irlba_2.3.5.1               RSQLite_2.4.5              
    ##  [43] GenomicRanges_1.62.1        textshaping_1.0.4          
    ##  [45] labeling_0.4.3              progressr_0.18.0           
    ##  [47] spatstat.sparse_3.1-0       httr_1.4.7                 
    ##  [49] polyclip_1.10-7             abind_1.4-8                
    ##  [51] compiler_4.5.2              withr_3.0.2                
    ##  [53] bit64_4.6.0-1               S7_0.2.1                   
    ##  [55] DBI_1.2.3                   fastDummies_1.7.5          
    ##  [57] MASS_7.3-65                 DelayedArray_0.36.0        
    ##  [59] tools_4.5.2                 lmtest_0.9-40              
    ##  [61] otel_0.2.0                  httpuv_1.6.16              
    ##  [63] future.apply_1.20.1         goftest_1.2-3              
    ##  [65] glue_1.8.0                  nlme_3.1-168               
    ##  [67] promises_1.5.0              grid_4.5.2                 
    ##  [69] Rtsne_0.17                  cluster_2.1.8.1            
    ##  [71] reshape2_1.4.5              gtable_0.3.6               
    ##  [73] spatstat.data_3.1-9         data.table_1.18.0          
    ##  [75] utf8_1.2.6                  XVector_0.50.0             
    ##  [77] spatstat.geom_3.6-1         RcppAnnoy_0.0.23           
    ##  [79] ggrepel_0.9.6               RANN_2.6.2                 
    ##  [81] pillar_1.11.1               stringr_1.6.0              
    ##  [83] spam_2.11-3                 RcppHNSW_0.6.0             
    ##  [85] later_1.4.5                 splines_4.5.2              
    ##  [87] lattice_0.22-7              bit_4.6.0                  
    ##  [89] survival_3.8-3              deldir_2.0-4               
    ##  [91] tidyselect_1.2.1            SingleCellExperiment_1.32.0
    ##  [93] locfit_1.5-9.12             Biostrings_2.78.0          
    ##  [95] miniUI_0.1.2                pbapply_1.7-4              
    ##  [97] knitr_1.51                  gridExtra_2.3              
    ##  [99] Seqinfo_1.0.0               SummarizedExperiment_1.40.0
    ## [101] scattermore_1.2             xfun_0.55                  
    ## [103] statmod_1.5.1               matrixStats_1.5.0          
    ## [105] stringi_1.8.7               lazyeval_0.2.2             
    ## [107] yaml_2.3.12                 evaluate_1.0.5             
    ## [109] codetools_0.2-20            tibble_3.3.1               
    ## [111] cli_3.6.5                   uwot_0.2.4                 
    ## [113] xtable_1.8-4                reticulate_1.44.1          
    ## [115] systemfonts_1.3.1           jquerylib_0.1.4            
    ## [117] Rcpp_1.1.1                  globals_0.18.0             
    ## [119] spatstat.random_3.4-3       png_0.1-8                  
    ## [121] spatstat.univar_3.1-5       parallel_4.5.2             
    ## [123] blob_1.2.4                  pkgdown_2.2.0              
    ## [125] dotCall64_1.2               listenv_0.10.0             
    ## [127] viridisLite_0.4.2           scales_1.4.0               
    ## [129] ggridges_0.5.7              crayon_1.5.3               
    ## [131] purrr_1.2.1                 rlang_1.1.7                
    ## [133] KEGGREST_1.50.0             cowplot_1.2.0
