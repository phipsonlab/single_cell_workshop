# Module 3: Cell Type Annotation

## Introduction

In this module, we identify cell types by analysing marker genes that
distinguish each cluster. Cell type annotation is a critical step in
single cell analysis, as it transforms abstract cluster numbers into
biologically meaningful cell populations.

There are two main approaches to cell type annotation:

1.  **Marker-based annotation** (manual): Examine differentially
    expressed genes in each cluster and match them to known cell type
    markers from the literature. This approach requires domain knowledge
    but provides full control and interpretability.

2.  **Reference-based annotation** (automated): Use tools like Azimuth,
    SingleR, or scArches to transfer labels from a reference dataset.
    This is faster but may miss context-specific populations.

In this workshop, we use marker-based annotation to teach the
fundamental concepts. The human heart contains several major cell types:

- **Cardiomyocytes**: The contractile muscle cells, marked by TNNT2,
  TTN, MYH7
- **Fibroblasts**: Structural cells producing extracellular matrix (DCN,
  COL1A1)
- **Endothelial cells**: Line blood vessels (PECAM1, VWF, CDH5)
- **Pericytes/Smooth muscle**: Support vascular structures (ACTA2, RGS5)
- **Immune cells**: Macrophages and lymphocytes (PTPRC, CD68, CD163)
- **Epicardial cells**: Cover the heart surface (WT1, TBX18)
- **Neural cells**: Neurons and Schwann cells (NRXN1, PLP1)

## Load Libraries and Data

``` r
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(RColorBrewer)
library(pheatmap)
```

We load the integrated and clustered Seurat object from Module 2.

``` r
# Load the clustered data
seu <- readRDS("../data/processed/02_integrated_clustered.rds")

# Verify the data
cat("Loaded Seurat object:\n")
```

    ## Loaded Seurat object:

``` r
cat("- Cells:", ncol(seu), "\n")
```

    ## - Cells: 10000

``` r
cat("- Genes:", nrow(seu), "\n")
```

    ## - Genes: 17743

``` r
cat("- Clusters:", length(unique(seu$seurat_clusters)), "\n")
```

    ## - Clusters: 17

Let us first examine the cluster distribution to understand what we are
annotating.

``` r
# Cluster sizes
cluster_sizes <- table(seu$seurat_clusters)

cat("Cluster sizes:\n")
```

    ## Cluster sizes:

``` r
for (i in seq_along(cluster_sizes)) {
    cat(sprintf("  Cluster %s: %d cells\n", names(cluster_sizes)[i], cluster_sizes[i]))
}
```

    ##   Cluster 0: 2482 cells
    ##   Cluster 1: 1827 cells
    ##   Cluster 2: 1211 cells
    ##   Cluster 3: 661 cells
    ##   Cluster 4: 593 cells
    ##   Cluster 5: 559 cells
    ##   Cluster 6: 554 cells
    ##   Cluster 7: 363 cells
    ##   Cluster 8: 238 cells
    ##   Cluster 9: 224 cells
    ##   Cluster 10: 215 cells
    ##   Cluster 11: 213 cells
    ##   Cluster 12: 199 cells
    ##   Cluster 13: 190 cells
    ##   Cluster 14: 180 cells
    ##   Cluster 15: 177 cells
    ##   Cluster 16: 114 cells

## Define Colour Palette

We maintain consistent colours throughout our analysis.

``` r
# Group colours (developmental stage - lowercase to match data)
group_colors <- c(
    "fetal" = "#E64B35",
    "young" = "#4DBBD5",
    "adult" = "#3C5488"
)

# Generate a palette for clusters
n_clusters <- length(unique(seu$seurat_clusters))
cluster_colors <- colorRampPalette(brewer.pal(12, "Paired"))(n_clusters)
names(cluster_colors) <- levels(seu$seurat_clusters)
```

## Visualise Current Clusters

Before annotation, let us examine the clusters and their composition by
developmental stage.

``` r
# UMAP coloured by cluster and by group
p1 <- DimPlot(seu, group.by = "seurat_clusters", label = TRUE,
              label.size = 4, repel = TRUE, cols = cluster_colors) +
    ggtitle("Clusters") +
    theme(legend.position = "none")

p2 <- DimPlot(seu, group.by = "group", cols = group_colors) +
    ggtitle("Developmental Stage")

p1 + p2
```

![](03_cell_type_annotation_files/figure-html/umap-clusters-1.png)

Understanding the developmental composition of each cluster provides
important context for annotation.

``` r
# Calculate cluster composition by group
comp_data <- seu@meta.data %>%
    group_by(seurat_clusters, group) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(seurat_clusters) %>%
    mutate(pct = 100 * n / sum(n))

ggplot(comp_data, aes(x = seurat_clusters, y = pct, fill = group)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = group_colors) +
    labs(x = "Cluster", y = "Percentage", fill = "Stage",
         title = "Developmental Stage Composition per Cluster") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 10))
```

![](03_cell_type_annotation_files/figure-html/cluster-composition-1.png)

We observe that some clusters are highly enriched for specific
developmental stages. For example, clusters 0, 7, 8, 9, 12, 13, and 16
are nearly exclusively fetal, while cluster 2 contains predominantly
young and adult cells. This developmental segregation often reflects
biological differences in cell state, particularly among cardiomyocytes
where fetal cells have distinct transcriptional profiles from mature
cells.

## Find Cluster Markers

We use
[`FindAllMarkers()`](https://satijalab.org/seurat/reference/FindAllMarkers.html)
to identify genes that distinguish each cluster from all other cells.
This is the foundation of marker-based annotation.

``` r
# Set default assay to SCT
DefaultAssay(seu) <- "SCT"

# Prepare for marker finding with SCT assay
seu <- PrepSCTFindMarkers(seu)

# Find markers for all clusters
# We use only.pos = TRUE to find upregulated markers
# min.pct = 0.25 requires the gene to be expressed in at least 25% of cells
# logfc.threshold = 0.5 requires at least 0.5 log2 fold change
markers <- FindAllMarkers(
    seu,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.5,
    test.use = "wilcox"
)

cat("Found", nrow(markers), "marker genes across all clusters\n")
```

    ## Found 19280 marker genes across all clusters

Let us examine the top markers for each cluster.

``` r
# Get top 5 markers per cluster
top_markers <- markers %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 5) %>%
    select(cluster, gene, avg_log2FC, pct.1, pct.2, p_val_adj)

# Display top markers for each cluster
for (cl in levels(seu$seurat_clusters)) {
    cl_markers <- top_markers %>% filter(cluster == cl)
    cat("\n--- Cluster", cl, "---\n")
    for (i in seq_len(nrow(cl_markers))) {
        cat(sprintf("  %s: log2FC=%.2f, pct.1=%.0f%%, pct.2=%.0f%%\n",
                    cl_markers$gene[i],
                    cl_markers$avg_log2FC[i],
                    cl_markers$pct.1[i] * 100,
                    cl_markers$pct.2[i] * 100))
    }
}
```

    ## 
    ## --- Cluster 0 ---
    ##   LINC02008: log2FC=2.75, pct.1=39%, pct.2=8%
    ##   GPR39: log2FC=2.73, pct.1=48%, pct.2=8%
    ##   CUX2: log2FC=2.68, pct.1=57%, pct.2=11%
    ##   RNF175: log2FC=2.66, pct.1=44%, pct.2=8%
    ##   C1QTNF1-AS1: log2FC=2.51, pct.1=37%, pct.2=7%
    ## 
    ## --- Cluster 1 ---
    ##   VIT: log2FC=6.10, pct.1=41%, pct.2=2%
    ##   ADH1B: log2FC=5.88, pct.1=58%, pct.2=5%
    ##   SCARA5: log2FC=5.86, pct.1=45%, pct.2=2%
    ##   ACSM1: log2FC=5.64, pct.1=42%, pct.2=2%
    ##   MGST1: log2FC=5.19, pct.1=50%, pct.2=3%
    ## 
    ## --- Cluster 2 ---
    ##   LINC01880: log2FC=6.24, pct.1=33%, pct.2=1%
    ##   UGT2B4: log2FC=5.23, pct.1=43%, pct.2=1%
    ##   ADRB1: log2FC=5.12, pct.1=37%, pct.2=1%
    ##   TMEM178B: log2FC=5.12, pct.1=83%, pct.2=8%
    ##   MYOM3: log2FC=5.11, pct.1=88%, pct.2=5%
    ## 
    ## --- Cluster 3 ---
    ##   FHL5: log2FC=7.73, pct.1=30%, pct.2=0%
    ##   EGFLAM: log2FC=6.60, pct.1=79%, pct.2=5%
    ##   AGAP2: log2FC=6.35, pct.1=38%, pct.2=1%
    ##   LINC02237: log2FC=6.23, pct.1=40%, pct.2=2%
    ##   COX4I2: log2FC=6.15, pct.1=55%, pct.2=2%
    ## 
    ## --- Cluster 4 ---
    ##   BTNL9: log2FC=6.96, pct.1=71%, pct.2=2%
    ##   CA4: log2FC=6.95, pct.1=40%, pct.2=1%
    ##   NR5A2: log2FC=6.91, pct.1=47%, pct.2=1%
    ##   TPO: log2FC=6.82, pct.1=49%, pct.2=3%
    ##   NOTCH4: log2FC=6.66, pct.1=76%, pct.2=2%
    ## 
    ## --- Cluster 5 ---
    ##   S100A1: log2FC=3.93, pct.1=32%, pct.2=4%
    ##   ACTA1: log2FC=3.93, pct.1=79%, pct.2=16%
    ##   RASD1: log2FC=3.43, pct.1=27%, pct.2=4%
    ##   MYL2: log2FC=3.42, pct.1=100%, pct.2=61%
    ##   HLA-DPA1: log2FC=3.38, pct.1=42%, pct.2=5%
    ## 
    ## --- Cluster 6 ---
    ##   LILRB5: log2FC=7.82, pct.1=50%, pct.2=0%
    ##   MARCO: log2FC=7.66, pct.1=27%, pct.2=1%
    ##   F13A1: log2FC=7.65, pct.1=91%, pct.2=11%
    ##   MS4A14: log2FC=7.30, pct.1=39%, pct.2=0%
    ##   SIGLEC1: log2FC=7.30, pct.1=66%, pct.2=1%
    ## 
    ## --- Cluster 7 ---
    ##   CSMD1: log2FC=4.70, pct.1=65%, pct.2=18%
    ##   OPCML: log2FC=4.14, pct.1=65%, pct.2=12%
    ##   BRINP3: log2FC=3.31, pct.1=92%, pct.2=24%
    ##   ZMAT4: log2FC=2.76, pct.1=63%, pct.2=16%
    ##   BRINP2: log2FC=2.65, pct.1=39%, pct.2=8%
    ## 
    ## --- Cluster 8 ---
    ##   E2F1: log2FC=4.78, pct.1=58%, pct.2=3%
    ##   CDC45: log2FC=4.30, pct.1=38%, pct.2=2%
    ##   DTL: log2FC=4.27, pct.1=63%, pct.2=4%
    ##   GINS2: log2FC=4.03, pct.1=39%, pct.2=2%
    ##   MCM10: log2FC=3.93, pct.1=36%, pct.2=2%
    ## 
    ## --- Cluster 9 ---
    ##   PANCR: log2FC=7.49, pct.1=52%, pct.2=1%
    ##   VWDE: log2FC=7.03, pct.1=70%, pct.2=1%
    ##   KCNJ3: log2FC=6.85, pct.1=89%, pct.2=2%
    ##   ZNF385B: log2FC=6.78, pct.1=86%, pct.2=9%
    ##   KCNH7: log2FC=6.77, pct.1=92%, pct.2=9%
    ## 
    ## --- Cluster 10 ---
    ##   KIF20A: log2FC=7.71, pct.1=34%, pct.2=0%
    ##   KIF18B: log2FC=7.43, pct.1=90%, pct.2=1%
    ##   PBK: log2FC=7.31, pct.1=80%, pct.2=1%
    ##   AURKB: log2FC=7.24, pct.1=66%, pct.2=1%
    ##   IQGAP3: log2FC=7.13, pct.1=78%, pct.2=1%
    ## 
    ## --- Cluster 11 ---
    ##   PKHD1L1: log2FC=8.72, pct.1=88%, pct.2=2%
    ##   MMRN1: log2FC=7.60, pct.1=46%, pct.2=1%
    ##   SMOC1: log2FC=6.80, pct.1=71%, pct.2=2%
    ##   PCDH15: log2FC=6.40, pct.1=71%, pct.2=12%
    ##   INHBA: log2FC=6.23, pct.1=55%, pct.2=5%
    ## 
    ## --- Cluster 12 ---
    ##   NRG1: log2FC=5.84, pct.1=71%, pct.2=11%
    ##   DCHS2: log2FC=5.44, pct.1=53%, pct.2=2%
    ##   HPSE2: log2FC=5.17, pct.1=60%, pct.2=8%
    ##   NDST4: log2FC=5.16, pct.1=40%, pct.2=2%
    ##   KHDRBS2: log2FC=5.13, pct.1=57%, pct.2=6%
    ## 
    ## --- Cluster 13 ---
    ##   SRARP: log2FC=6.87, pct.1=35%, pct.2=0%
    ##   OTUD1: log2FC=5.27, pct.1=67%, pct.2=5%
    ##   ATF3: log2FC=4.60, pct.1=70%, pct.2=6%
    ##   XIRP1: log2FC=4.58, pct.1=80%, pct.2=12%
    ##   PNMT: log2FC=4.52, pct.1=40%, pct.2=3%
    ## 
    ## --- Cluster 14 ---
    ##   SKAP1: log2FC=7.87, pct.1=68%, pct.2=3%
    ##   THEMIS: log2FC=7.68, pct.1=35%, pct.2=2%
    ##   ITK: log2FC=7.31, pct.1=48%, pct.2=1%
    ##   SCML4: log2FC=7.12, pct.1=36%, pct.2=0%
    ##   CD247: log2FC=7.07, pct.1=48%, pct.2=3%
    ## 
    ## --- Cluster 15 ---
    ##   NRXN1: log2FC=8.39, pct.1=99%, pct.2=12%
    ##   INSC: log2FC=8.32, pct.1=70%, pct.2=0%
    ##   TFAP2A: log2FC=8.12, pct.1=40%, pct.2=0%
    ##   GRIK3: log2FC=8.08, pct.1=60%, pct.2=1%
    ##   XKR4: log2FC=8.05, pct.1=99%, pct.2=10%
    ## 
    ## --- Cluster 16 ---
    ##   DHFR: log2FC=2.97, pct.1=90%, pct.2=44%
    ##   MYH7: log2FC=2.94, pct.1=100%, pct.2=62%
    ##   TNNI1: log2FC=2.90, pct.1=99%, pct.2=34%
    ##   MTRNR2L12: log2FC=2.55, pct.1=91%, pct.2=35%
    ##   MTRNR2L8: log2FC=2.53, pct.1=68%, pct.2=18%

### Visualise Top Markers with DotPlot

The DotPlot below shows the top 5 differentially expressed genes for
each cluster. This provides an unbiased view of what distinguishes each
cluster before we examine the cell type markers from the original study.

``` r
# Get top 5 unique markers per cluster, ordered by cluster
top5_genes <- top_markers %>%
    arrange(cluster, desc(avg_log2FC)) %>%
    pull(gene) %>%
    unique()

DotPlot(seu, features = top5_genes, group.by = "seurat_clusters") +
    RotatedAxis() +
    scale_colour_gradient2(low = "blue", mid = "white", high = "red",
                           midpoint = 0) +
    labs(title = "Top 5 Marker Genes per Cluster (from FindAllMarkers)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
          axis.text.y = element_text(size = 10))
```

![](03_cell_type_annotation_files/figure-html/dotplot-top-markers-1.png)

### Marker Heatmap

A heatmap of these top markers provides a comprehensive view of the
expression patterns across clusters.

``` r
# Get top 3 markers per cluster for heatmap
top3_markers <- markers %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 3) %>%
    pull(gene) %>%
    unique()

# Downsample for visualisation (50 cells per cluster max)
set.seed(42)
cells_subset <- seu@meta.data %>%
    mutate(cell_id = rownames(seu@meta.data)) %>%
    group_by(seurat_clusters) %>%
    slice_sample(n = 50) %>%
    pull(cell_id)

# Create heatmap
DoHeatmap(subset(seu, cells = cells_subset),
          features = top3_markers,
          group.by = "seurat_clusters",
          size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0) +
  theme(axis.text.y = element_text(size = 6))
```

![](03_cell_type_annotation_files/figure-html/marker-heatmap-top-1.png)

## Cell Type Marker Analysis

While the differentially expressed genes provide a starting point, we
typically compare against known marker genes for the cell types expected
in our tissue. The human heart contains several well-characterised cell
populations.

### Cell Type Marker Genes

The marker genes below were identified by Sim et al. 2021 (Circulation)
from their analysis of the full human heart development dataset. These
markers distinguish the major cell populations across foetal, young, and
adult hearts.

#### Visualise with DotPlot

The DotPlot is a powerful visualisation that shows both the percentage
of cells expressing a gene (dot size) and the average expression level
(colour intensity).

``` r
# Define markers and cell types from Sim et al. 2021
marker_info <- data.frame(
    gene = c(
        # Cardiomyocytes
        "TNNT2", "ACTN2", "TNNI3", "RYR2", "MYH6", "TTN", "MYBPC3", "PLN", "MYL2",
        # Fibroblasts
        "COL1A1", "COL3A1", "DCN", "PDGFRA", "POSTN", "TCF21", "VIM",
        # Pericytes
        "KCNJ8", "VTN", "ABCC9", "HEYL",
        # Endothelial
        "PECAM1", "KDR", "CDH5", "EMCN", "FLT1",
        # Smooth Muscle
        "TAGLN", "MYH11", "MYLK", "LMOD1",
        # Macrophages
        "CD68", "CD163", "CD74", "F13A1", "LGALS3",
        # T Cells
        "CD3D", "CD3G", "CD3E", "IL7R",
        # Schwann/Neural
        "PLP1", "CNP", "S100B",
        # Epicardial
        "UPK3B", "WT1",
        # Proliferating
        "MKI67", "TOP2A"
    ),
    celltype = c(
        rep("Cardiomyocytes", 9),
        rep("Fibroblasts", 7),
        rep("Pericytes", 4),
        rep("Endothelial", 5),
        rep("Smooth Muscle", 4),
        rep("Macrophages", 5),
        rep("T Cells", 4),
        rep("Schwann", 3),
        rep("Epicardial", 2),
        rep("Proliferating", 2)
    )
)

# Define colors for each cell type
celltype_colors <- c(
    "Cardiomyocytes" = "#D62728",
    "Fibroblasts" = "#1F77B4",
    "Pericytes" = "#2CA02C",
    "Endothelial" = "#9467BD",
    "Smooth Muscle" = "#FF7F0E",
    "Macrophages" = "#8C564B",
    "T Cells" = "#E377C2",
    "Schwann" = "#7F7F7F",
    "Epicardial" = "#17BECF",
    "Proliferating" = "#333333"
)

# Filter to available genes
marker_info <- marker_info[marker_info$gene %in% rownames(seu), ]
marker_genes <- marker_info$gene
gene_colors <- celltype_colors[marker_info$celltype]
marker_info$celltype_f <- factor(marker_info$celltype, levels = names(celltype_colors))

# Order clusters by cell type identity for diagonal pattern
# CM (atrial then ventricular), Fibroblasts, Pericytes, Endothelial, Immune, Neural, Proliferating
cluster_order <- c("9", "2", "0", "7", "8", "13", "16", "5", "1", "12", "3", "4", "11", "6", "14", "15", "10")
seu$cluster_ordered <- factor(seu$seurat_clusters, levels = cluster_order)

# Calculate cell type label positions
celltype_positions <- marker_info %>%
    mutate(pos = row_number()) %>%
    group_by(celltype_f) %>%
    summarise(start = min(pos) - 0.5, end = max(pos) + 0.5,
              mid = mean(pos), .groups = "drop")

# Create DotPlot with cell type annotations
DotPlot(seu, features = marker_genes, group.by = "cluster_ordered") +
    RotatedAxis() +
    scale_colour_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                           midpoint = 0, name = "Average\nExpression") +
    scale_size_continuous(name = "Percent\nExpressed", range = c(0.5, 6)) +
    labs(x = NULL, y = "Cluster",
         title = "Cell Type Marker Expression (Sim et al. 2021)") +
    theme_bw() +
    theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10,
                                   color = gene_colors, face = "italic"),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "grey40", linewidth = 0.8),
        legend.position = "right",
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10),
        plot.margin = margin(10, 10, 45, 10)
    ) +
    # Add cell type color bars
    annotate("rect",
             xmin = celltype_positions$start, xmax = celltype_positions$end,
             ymin = -1.8, ymax = -1.2,
             fill = celltype_colors[as.character(celltype_positions$celltype_f)]) +
    # Add cell type labels inside bars
    annotate("text",
             x = celltype_positions$mid, y = -1.5,
             label = celltype_positions$celltype_f,
             size = 3.5, angle = 0, hjust = 0.5, fontface = "bold",
             color = "white") +
    coord_cartesian(clip = "off", ylim = c(0.5, 17.5))
```

![](03_cell_type_annotation_files/figure-html/dotplot-celltype-markers-1.png)

The DotPlot reveals a clear diagonal pattern where each cell type shows
specific marker expression. Gene names are coloured by their cell type,
matching the annotation bar at the bottom:

- **Cardiomyocytes** (red): TNNT2, ACTN2, TNNI3, RYR2, MYH6, TTN,
  MYBPC3, PLN, MYL2
- **Fibroblasts** (blue): COL1A1, COL3A1, DCN, PDGFRA, POSTN, TCF21, VIM
- **Pericytes** (green): KCNJ8, VTN, ABCC9, HEYL
- **Endothelial** (purple): PECAM1, KDR, CDH5, EMCN, FLT1
- **Smooth Muscle** (orange): TAGLN, MYH11, MYLK, LMOD1
- **Macrophages** (brown): CD68, CD163, CD74, F13A1, LGALS3
- **T Cells** (pink): CD3D, CD3G, CD3E, IL7R
- **Schwann** (grey): PLP1, CNP, S100B
- **Epicardial** (teal): UPK3B, WT1
- **Proliferating** (black): MKI67, TOP2A

### Feature Plots for Key Markers

FeaturePlots show the spatial distribution of marker expression on the
UMAP, helping us understand the relationship between clusters.

``` r
# Cardiomyocyte markers
FeaturePlot(seu,
            features = c("TNNT2", "TTN", "MYH7", "MYH6"),
            ncol = 2, order = TRUE) &
    scale_colour_gradient(low = "lightgrey", high = "darkred") &
    theme_minimal()
```

![](03_cell_type_annotation_files/figure-html/featureplot-cm-1.png)

The cardiomyocyte markers reveal an interesting pattern: MYH6 (atrial
isoform) is enriched in clusters 2 and 9, while MYH7 (ventricular
isoform) is more expressed in clusters 0 and 3. This suggests we can
distinguish atrial from ventricular cardiomyocytes.

``` r
# Stromal cell markers
FeaturePlot(seu,
            features = c("DCN", "PECAM1", "KCNJ8", "ACTA2"),
            ncol = 2, order = TRUE) &
    scale_colour_gradient(low = "lightgrey", high = "darkblue") &
    theme_minimal()
```

![](03_cell_type_annotation_files/figure-html/featureplot-stromal-1.png)

``` r
# Immune and other markers
FeaturePlot(seu,
            features = c("PTPRC", "F13A1", "NRXN1", "MKI67"),
            ncol = 2, order = TRUE) &
    scale_colour_gradient(low = "lightgrey", high = "darkgreen") &
    theme_minimal()
```

![](03_cell_type_annotation_files/figure-html/featureplot-immune-1.png)

## Assign Cell Type Annotations

Based on the marker analysis above, we can now assign cell type labels
to each cluster. We consider both the marker gene expression and the
developmental stage composition.

**Note:** These annotations are based on the 10,000 cell subset used for
this workshop. The full dataset (~54,000 cells) may reveal additional
rare populations or finer subclusters.

The table below summarises our annotations with the supporting marker
evidence:

| Cluster | Cell Type | Key Markers | Notes |
|----|----|----|----|
| 0 | Fetal CM | TNNT2+, MYH7+ | Fetal ventricular cardiomyocytes (100% fetal) |
| 1 | Fibroblasts | DCN+, COL1A1+ | Present across all developmental stages |
| 2 | Atrial CM | MYH6+ (atrial isoform) | Mature atrial CM (young/adult enriched) |
| 3 | Pericytes | KCNJ8+, RGS5+, ACTA2+ | Vascular support cells, some smooth muscle |
| 4 | Endothelial | PECAM1+, VWF+ | Main endothelial population |
| 5 | Ventricular CM | TNNT2+, MYH7+ | Mature ventricular CM (young enriched) |
| 6 | Macrophages | CD68+, F13A1+ | Tissue-resident macrophages (100% immune) |
| 7 | Fetal CM | TNNT2+, MYH7+ | Fetal ventricular cardiomyocytes |
| 8 | Fetal CM | TNNT2+, MYH7+ | Fetal ventricular cardiomyocytes |
| 9 | Fetal Atrial CM | MYH6+, TNNT2+ | Fetal atrial cardiomyocytes |
| 10 | Proliferating | MKI67+, TOP2A+, TNNT2+ | Proliferating cells (fetal enriched) |
| 11 | Endothelial | PECAM1+, EMCN+ | Endothelial cells (100% endothelial) |
| 12 | Fetal Fibroblasts | DCN+, COL1A1+ | Fetal-enriched fibroblasts (99% fetal) |
| 13 | Fetal CM | TNNT2+, MYH7+ | Fetal ventricular cardiomyocytes |
| 14 | Immune | PTPRC+ (CD45) | Mixed immune cells |
| 15 | Neurons | NRXN1+, CNP+ | Neural/Schwann cells (100% neurons) |
| 16 | Fetal CM | TNNT2+, MYH7+ | Fetal ventricular cardiomyocytes |

``` r
# Define cell type annotations for all 17 clusters
cluster_annotations <- c(
  "0" = "Fetal CM",
  "1" = "Fibroblasts",
  "2" = "Atrial CM",
  "3" = "Pericytes",
  "4" = "Endothelial",
  "5" = "Ventricular CM",
  "6" = "Macrophages",
  "7" = "Fetal CM",
  "8" = "Fetal CM",
  "9" = "Fetal Atrial CM",
  "10" = "Proliferating",
  "11" = "Endothelial",
  "12" = "Fetal Fibroblasts",
  "13" = "Fetal CM",
  "14" = "Immune",
  "15" = "Neurons",
  "16" = "Fetal CM"
)

cell_type_df <- data.frame(
  cell_type = cluster_annotations[as.character(seu$seurat_clusters)],
  row.names = colnames(seu)
)
seu <- AddMetaData(seu, metadata = cell_type_df)

# Broad categories for summary analyses
broad_annotations <- c(
  "0" = "Cardiomyocytes",
  "1" = "Fibroblasts",
  "2" = "Cardiomyocytes",
  "3" = "Pericytes",
  "4" = "Endothelial",
  "5" = "Cardiomyocytes",
  "6" = "Immune",
  "7" = "Cardiomyocytes",
  "8" = "Cardiomyocytes",
  "9" = "Cardiomyocytes",
  "10" = "Proliferating",
  "11" = "Endothelial",
  "12" = "Fibroblasts",
  "13" = "Cardiomyocytes",
  "14" = "Immune",
  "15" = "Neural",
  "16" = "Cardiomyocytes"
)

cell_type_broad_df <- data.frame(
  cell_type_broad = broad_annotations[as.character(seu$seurat_clusters)],
  row.names = colnames(seu)
)
seu <- AddMetaData(seu, metadata = cell_type_broad_df)
```

## Visualise Annotated Cell Types

With cell type annotations assigned, we now examine the results to
verify our annotations make biological sense and to understand the
cellular composition of the heart across developmental stages.

### Annotated UMAP

The UMAP visualisation with cell type labels provides the key output of
our annotation workflow. We show both detailed annotations
(distinguishing fetal from adult cardiomyocytes, atrial from
ventricular) and broad categories.

``` r
# Define colours for cell types (11 unique detailed types)
celltype_colors <- c(
    "Fetal CM" = "#E64B35",
    "Fibroblasts" = "#00A087",
    "Atrial CM" = "#3C5488",
    "Ventricular CM" = "#4DBBD5",
    "Pericytes" = "#F39B7F",
    "Endothelial" = "#8491B4",
    "Macrophages" = "#91D1C2",
    "Fetal Atrial CM" = "#B09C85",
    "Fetal Fibroblasts" = "#00CED1",
    "Neurons" = "#CD534C",
    "Immune" = "#868686",
    "Proliferating" = "#9370DB"
)

# Broad categories
broad_colors <- c(
    "Cardiomyocytes" = "#E64B35",
    "Fibroblasts" = "#00A087",
    "Endothelial" = "#8491B4",
    "Pericytes" = "#F39B7F",
    "Immune" = "#91D1C2",
    "Neural" = "#CD534C",
    "Proliferating" = "#9370DB"
)

p1 <- DimPlot(seu, group.by = "cell_type", label = TRUE,
              label.size = 3, repel = TRUE) +
    scale_color_manual(values = celltype_colors) +
    ggtitle("Detailed Cell Types") +
    theme(legend.position = "none")

p2 <- DimPlot(seu, group.by = "cell_type_broad", label = TRUE,
              label.size = 4, repel = TRUE) +
    scale_color_manual(values = broad_colors) +
    ggtitle("Broad Cell Types")

p1 + p2
```

![](03_cell_type_annotation_files/figure-html/umap-annotated-1.png)

The annotated UMAP reveals a well-organised cellular landscape:

- **Cardiomyocytes** form the dominant population, with fetal cells (red
  tones) clustering separately from mature atrial and ventricular
  populations
- **Stromal cells** (fibroblasts, pericytes, endothelial) occupy
  distinct regions of the UMAP space
- **Immune and neural cells** form smaller but clearly defined clusters

### Cluster Size Distribution

Before examining developmental patterns, we assess how many cells belong
to each cluster. This helps interpret whether observed patterns reflect
genuine biology or sampling effects.

``` r
# Cluster sizes with cell type labels
cluster_size_data <- seu@meta.data %>%
    group_by(seurat_clusters, cell_type) %>%
    summarise(n = n(), .groups = "drop")

ggplot(cluster_size_data, aes(x = seurat_clusters, y = n, fill = cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = celltype_colors) +
  labs(x = "Cluster", y = "Number of Cells", fill = "Cell Type",
       title = "Cell Count per Cluster") +
  theme_minimal() +
  theme(legend.position = "right",
        legend.text = element_text(size = 8))
```

![](03_cell_type_annotation_files/figure-html/cluster-sizes-1.png)

Cluster 0 (Fetal CM) is the largest population with ~2,500 cells,
consistent with fetal samples contributing the most cells. Clusters 16
and 10 are among the smallest, representing fetal cardiomyocytes and
proliferating cells respectively.

### Developmental Composition per Cluster

A critical question in developmental studies is: which cell populations
are stage-specific versus shared across development? This stacked bar
chart shows the fetal/young/adult breakdown within each cluster.

``` r
# Calculate developmental composition per cluster
cluster_dev_comp <- seu@meta.data %>%
    group_by(seurat_clusters, group) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(seurat_clusters) %>%
    mutate(pct = 100 * n / sum(n))

ggplot(cluster_dev_comp, aes(x = seurat_clusters, y = pct, fill = group)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = group_colors) +
  labs(x = "Cluster", y = "Percentage", fill = "Stage",
       title = "Developmental Stage Composition per Cluster") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 10))
```

![](03_cell_type_annotation_files/figure-html/dev-composition-per-cluster-1.png)

This figure reveals striking developmental patterns:

- **Fetal-specific clusters** (0, 7, 8, 9, 12, 13, 16): These
  cardiomyocyte and fibroblast populations are almost exclusively fetal,
  reflecting the distinct transcriptional state of immature heart cells
- **Mature CM clusters** (2, 5): Atrial and ventricular CM clusters
  contain primarily young and adult cells, representing mature
  cardiomyocyte subtypes
- **Shared populations** (1, 3, 4, 6, 11, 14, 15): Stromal, immune, and
  neural cells are present across all developmental stages, suggesting
  these supporting cell types maintain similar identities throughout
  heart development
- **Proliferating cells** (10): Nearly exclusively fetal, consistent
  with the high proliferative capacity of fetal cardiac cells

### Cell Type Composition by Developmental Stage

We can also view this relationship from the opposite perspective: what
is the cellular composition within each developmental stage?

``` r
# Calculate composition
stage_comp <- seu@meta.data %>%
    group_by(group, cell_type_broad) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(group) %>%
    mutate(pct = 100 * n / sum(n))

ggplot(stage_comp, aes(x = group, y = pct, fill = cell_type_broad)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = broad_colors) +
    labs(x = "Developmental Stage", y = "Percentage", fill = "Cell Type",
         title = "Cell Type Composition by Developmental Stage") +
    theme_minimal()
```

![](03_cell_type_annotation_files/figure-html/celltype-by-stage-1.png)

This complementary view shows:

- **Cardiomyocytes** dominate all stages but show distinct subtypes
  (fetal vs mature)
- **Proliferating cells** are predominantly fetal, consistent with the
  known proliferative capacity of fetal cardiomyocytes
- **Stromal populations** (fibroblasts, endothelial, pericytes) are
  relatively stable across development

### Split UMAP by Developmental Stage

Finally, we split the UMAP by developmental stage to directly visualise
which cell types are present at each time point. This view clearly shows
how the cellular landscape changes from fetal to adult heart.

``` r
DimPlot(seu, group.by = "cell_type", split.by = "group",
        label = FALSE, cols = celltype_colors) +
  ggtitle("Cell Types by Developmental Stage") +
  theme(legend.position = "right",
        legend.text = element_text(size = 8))
```

![](03_cell_type_annotation_files/figure-html/split-umap-1.png)

The split view highlights several key observations:

- **Fetal heart** is dominated by Fetal CM populations, with
  proliferating cells visible (cluster 10)
- **Young heart** shows a transition state with both mature
  atrial/ventricular CM and some residual fetal-like populations
- **Adult heart** shows predominantly mature cardiomyocyte subtypes
  (Atrial and Ventricular CM)
- **Stromal and immune populations** are consistent across all stages

## Save Annotated Object

``` r
# Save the annotated Seurat object
saveRDS(seu, "../data/processed/03_annotated.rds")

# Also save the markers table
write.csv(markers, "../results/cluster_markers.csv", row.names = FALSE)

message("Saved annotated object to: data/processed/03_annotated.rds")
message("Saved marker genes to: results/cluster_markers.csv")
```

## Summary

In this module, we performed cell type annotation on the clustered data.
We:

- Identified marker genes for each cluster using
  [`FindAllMarkers()`](https://satijalab.org/seurat/reference/FindAllMarkers.html)
- Compared cluster markers against known heart cell type markers from
  Sim et al. 2021
- Visualised marker expression using DotPlots and FeaturePlots
- Assigned cell type labels based on marker gene expression and
  developmental composition
- Created both detailed (12 types) and broad (7 types) cell type
  annotations

The annotated dataset is now ready for differential expression and
composition analysis in Module 4.

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
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] pheatmap_1.0.13    RColorBrewer_1.1-3 patchwork_1.3.0    tidyr_1.3.1       
    ## [5] dplyr_1.1.4        ggplot2_3.5.2      Seurat_5.3.0       SeuratObject_5.1.0
    ## [9] sp_2.2-0          
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] deldir_2.0-4           pbapply_1.7-2          gridExtra_2.3         
    ##   [4] rlang_1.1.6            magrittr_2.0.3         RcppAnnoy_0.0.22      
    ##   [7] spatstat.geom_3.4-1    matrixStats_1.5.0      ggridges_0.5.6        
    ##  [10] compiler_4.5.2         png_0.1-8              systemfonts_1.2.3     
    ##  [13] vctrs_0.6.5            reshape2_1.4.4         stringr_1.5.1         
    ##  [16] pkgconfig_2.0.3        fastmap_1.2.0          labeling_0.4.3        
    ##  [19] promises_1.3.2         rmarkdown_2.29         ragg_1.4.0            
    ##  [22] purrr_1.0.4            xfun_0.52              cachem_1.1.0          
    ##  [25] jsonlite_2.0.0         goftest_1.2-3          later_1.4.2           
    ##  [28] spatstat.utils_3.1-4   irlba_2.3.5.1          parallel_4.5.2        
    ##  [31] cluster_2.1.8.1        R6_2.6.1               ica_1.0-3             
    ##  [34] spatstat.data_3.1-6    bslib_0.9.0            stringi_1.8.7         
    ##  [37] limma_3.64.0           reticulate_1.42.0      spatstat.univar_3.1-3 
    ##  [40] parallelly_1.44.0      lmtest_0.9-40          jquerylib_0.1.4       
    ##  [43] scattermore_1.2        Rcpp_1.0.14            knitr_1.50            
    ##  [46] tensor_1.5             future.apply_1.11.3    zoo_1.8-14            
    ##  [49] sctransform_0.4.2      httpuv_1.6.16          Matrix_1.7-4          
    ##  [52] splines_4.5.2          igraph_2.1.4           tidyselect_1.2.1      
    ##  [55] abind_1.4-8            yaml_2.3.10            spatstat.random_3.4-1 
    ##  [58] spatstat.explore_3.4-3 codetools_0.2-20       miniUI_0.1.2          
    ##  [61] listenv_0.9.1          lattice_0.22-7         tibble_3.2.1          
    ##  [64] plyr_1.8.9             withr_3.0.2            shiny_1.10.0          
    ##  [67] ROCR_1.0-11            evaluate_1.0.3         Rtsne_0.17            
    ##  [70] future_1.49.0          fastDummies_1.7.5      desc_1.4.3            
    ##  [73] survival_3.8-3         polyclip_1.10-7        fitdistrplus_1.2-2    
    ##  [76] pillar_1.10.2          KernSmooth_2.23-26     plotly_4.10.4         
    ##  [79] generics_0.1.4         RcppHNSW_0.6.0         scales_1.4.0          
    ##  [82] globals_0.18.0         xtable_1.8-4           glue_1.8.0            
    ##  [85] lazyeval_0.2.2         tools_4.5.2            data.table_1.17.2     
    ##  [88] RSpectra_0.16-2        RANN_2.6.2             fs_1.6.6              
    ##  [91] dotCall64_1.2          cowplot_1.1.3          grid_4.5.2            
    ##  [94] colorspace_2.1-1       nlme_3.1-168           cli_3.6.5             
    ##  [97] spatstat.sparse_3.1-0  textshaping_1.0.1      spam_2.11-1           
    ## [100] viridisLite_0.4.2      uwot_0.2.3             gtable_0.3.6          
    ## [103] sass_0.4.10            digest_0.6.37          progressr_0.15.1      
    ## [106] ggrepel_0.9.6          htmlwidgets_1.6.4      farver_2.1.2          
    ## [109] htmltools_0.5.8.1      pkgdown_2.2.0          lifecycle_1.0.4       
    ## [112] httr_1.4.7             statmod_1.5.0          mime_0.13             
    ## [115] MASS_7.3-65
