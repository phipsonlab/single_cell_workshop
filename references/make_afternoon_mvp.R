library(Seurat)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)

data_dir   <- "~/Dropbox/scWorkshop/data"
output_dir <- "~/Dropbox/scWorkshop/processed"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Load ----
counts   <- readRDS(file.path(data_dir, "heart-counts.Rds"))
cellinfo <- readRDS(file.path(data_dir, "cellinfo_updated.Rds"))

# ---- 2. Seurat object + metadata ----
seu <- CreateSeuratObject(counts, project = "heart_workshop",
                          min.cells = 3, min.features = 200)
m <- cellinfo[match(colnames(seu), cellinfo$CellID), ]
seu$sample <- m$Sample; seu$group <- m$Group; seu$sex <- m$Sex
seu[["percent.mt"]]   <- PercentageFeatureSet(seu, pattern = "^MT-")
seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern = "^RP[SL]")

# ---- 3. Cell QC ----
seu <- subset(seu, subset =
                nFeature_RNA >= 500 &
                nCount_RNA   >= 2500 & nCount_RNA <= 40000 &
                percent.mt   <= 20)

# ---- 4. Gene filter (MT / ribo / no-Entrez / sex chr) ----
g   <- rownames(seu)
ann <- AnnotationDbi::select(org.Hs.eg.db, keys = g, keytype = "SYMBOL",
                             columns = c("ENTREZID","GENENAME","CHR"))
ann <- ann[!duplicated(ann$SYMBOL), ]
ann <- ann[match(g, ann$SYMBOL), ]

drop <- unique(c(
  grep("^MT-", g),
  grep("mitochondrial", ann$GENENAME, ignore.case = TRUE),
  grep("^RP[SL][0-9]", g),
  grep("ribosomal", ann$GENENAME, ignore.case = TRUE),
  which(is.na(ann$ENTREZID)),
  which(ann$CHR %in% c("X","Y"))
))
seu <- subset(seu, features = setdiff(g, g[drop]))

# ---- 5. Stratified downsample to 10k ----
set.seed(42)
n_target <- 10000
tbl <- seu@meta.data %>% count(sample, name = "n_cells") %>%
  mutate(n_sample = round(n_cells / sum(n_cells) * n_target))
tbl$n_sample[which.max(tbl$n_cells)] <-
  tbl$n_sample[which.max(tbl$n_cells)] + (n_target - sum(tbl$n_sample))

set.seed(42)
keep <- unlist(lapply(tbl$sample, function(s) {
  cells <- colnames(seu)[seu$sample == s]
  sample(cells, min(tbl$n_sample[tbl$sample == s], length(cells)))
}))
seu <- subset(seu, cells = keep)

# ---- 6. Save ----
saveRDS(seu, file.path(output_dir, "afternoon_mvp.rds"))
