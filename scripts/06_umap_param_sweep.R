# scripts/06_umap_param_sweep.R
#
# Sweep UMAP parameters (n_neighbors x min_dist) on the Module 6 CM subset,
# for both embeddings used by the pseudotime module:
#
#   - PCA (20 PCs on scran logcounts of HVGs)
#   - PhiSpace scores
#
# Purpose: pick a (n_neighbors, min_dist) pair per embedding that gives a
# well-connected maturation trajectory rather than fragmented islands.
# Not part of the rendered site.
#
# Run from the repo root:
#
#   Rscript scripts/06_umap_param_sweep.R
#
# Inputs:  results/06_cm_pca.rds    (produced by knitting Module 6 once)
# Outputs: results/umap_sweep/umap_{pca|phispace}_nn{nn}_md{md}.rds
#          results/umap_sweep/sweep_{pca|phispace}_by_{group|phicelltype}.png

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(scater)
    library(ggplot2)
    library(patchwork)
    library(RColorBrewer)
})

# ---- Paths & inputs --------------------------------------------------------

results_dir <- "results"
sweep_dir   <- file.path(results_dir, "umap_sweep")
cm_cache    <- file.path(results_dir, "06_cm_pca.rds")

if (!file.exists(cm_cache)) {
    stop("Missing ", cm_cache, ". Knit vignettes/06_pseudotime.Rmd once to ",
         "produce the CM object with PCA + PhiSpace reducedDims, then re-run ",
         "this script.")
}

dir.create(sweep_dir, showWarnings = FALSE, recursive = TRUE)

cm <- readRDS(cm_cache)
message("Loaded CM object: ", ncol(cm), " cells, reducedDims = ",
        paste(reducedDimNames(cm), collapse = ", "))

stopifnot("PCA"      %in% reducedDimNames(cm),
          "PhiSpace" %in% reducedDimNames(cm))

# ---- Sweep grid ------------------------------------------------------------

nn_vals <- c(10, 15, 30, 50, 100)
md_vals <- c(0.01, 0.1, 0.3, 0.5)
embeddings <- c("PCA", "PhiSpace")

# Short key for filenames / titles
embed_key <- function(dimred) tolower(dimred)

# ---- UMAP runner (cached per config) --------------------------------------

run_one_umap <- function(cm, dimred, nn, md) {
    fp <- file.path(
        sweep_dir,
        sprintf("umap_%s_nn%d_md%s.rds",
                embed_key(dimred), nn, sub("\\.", "p", format(md)))
    )
    if (file.exists(fp)) {
        return(readRDS(fp))
    }
    set.seed(42)
    cm2 <- scater::runUMAP(
        cm,
        dimred      = dimred,
        name        = "UMAP_tmp",
        n_neighbors = nn,
        min_dist    = md
    )
    mat <- reducedDim(cm2, "UMAP_tmp")
    colnames(mat) <- c("UMAP1", "UMAP2")
    saveRDS(mat, fp)
    mat
}

# ---- Plot helper -----------------------------------------------------------

plot_one <- function(mat, colour_vec, palette, title, legend = FALSE) {
    df <- as.data.frame(mat)
    df$colour <- colour_vec
    p <- ggplot(df, aes(UMAP1, UMAP2, colour = colour)) +
        geom_point(size = 0.25, alpha = 0.6) +
        theme_bw(base_size = 9) +
        theme(
            aspect.ratio       = 1,
            axis.text          = element_blank(),
            axis.ticks         = element_blank(),
            panel.grid.minor   = element_blank(),
            plot.title         = element_text(size = 9),
            legend.position    = if (legend) "right" else "none"
        ) +
        labs(x = NULL, y = NULL, title = title, colour = NULL) +
        palette
    p
}

# Build a 5 x 4 grid of panels for a given embedding and colouring.
build_grid <- function(dimred, colour_vec, palette, subtitle) {
    panels <- vector("list", length(nn_vals) * length(md_vals))
    i <- 1
    for (nn in nn_vals) {
        for (md in md_vals) {
            message("  UMAP: ", dimred, "  nn=", nn, "  min_dist=", md)
            mat <- run_one_umap(cm, dimred, nn, md)
            panels[[i]] <- plot_one(
                mat, colour_vec, palette,
                sprintf("nn=%d, md=%s", nn, format(md))
            )
            i <- i + 1
        }
    }
    # Collect a single legend from a throwaway panel
    legend_panel <- plot_one(
        run_one_umap(cm, dimred, nn_vals[1], md_vals[1]),
        colour_vec, palette, "", legend = TRUE
    )
    (wrap_plots(panels, ncol = length(md_vals)) |
         legend_panel) +
        plot_layout(widths = c(length(md_vals), 1)) +
        plot_annotation(
            title    = sprintf("UMAP sweep on %s embedding", dimred),
            subtitle = subtitle
        )
}

# ---- Palettes & colour vectors --------------------------------------------

group_vec <- factor(cm$group)
pal_group <- scale_colour_brewer(palette = "Set1", drop = FALSE)

phict_vec <- factor(cm$PhiCellType)
# PhiCellType has many levels; pool Set1 + Dark2 + Set2 for enough colours
pool <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Dark2"), brewer.pal(8, "Set2"))
pal_phict <- scale_colour_manual(
    values = pool[seq_along(levels(phict_vec))],
    drop   = FALSE
)

# ---- Build and save --------------------------------------------------------

out_files <- character()

for (dr in embeddings) {
    # by donor group
    g_grp <- build_grid(dr, group_vec, pal_group,
                        "coloured by donor stage (group)")
    fp <- file.path(sweep_dir,
                    sprintf("sweep_%s_by_group.png", embed_key(dr)))
    ggsave(fp, g_grp, width = 14, height = 16, dpi = 120)
    out_files <- c(out_files, fp)
    message("Wrote ", fp)

    # by PhiCellType
    g_ct <- build_grid(dr, phict_vec, pal_phict,
                       "coloured by PhiCellType (fine Gao subclusters)")
    fp <- file.path(sweep_dir,
                    sprintf("sweep_%s_by_phicelltype.png", embed_key(dr)))
    ggsave(fp, g_ct, width = 14, height = 16, dpi = 120)
    out_files <- c(out_files, fp)
    message("Wrote ", fp)
}

message("\nSweep complete. Files written:")
for (f in out_files) message("  ", f)
message("\nInspect the PNGs in ", sweep_dir, " and pick one (n_neighbors, ",
        "min_dist) per embedding.")
