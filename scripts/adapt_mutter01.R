#!/usr/bin/env Rscript
# Stage B adapter -- Mutter_01: raw counts parquet -> per-sample common-format Seurat.
# Sparse RNA assay on the common panel; shared metadata joined on cell_id; Yi's vendor
# columns sequestered to yi_reference.tsv; objects split by Condition.
# Args: <samplesheet.tsv> <common_genes.tsv> <DATADIR> <yi_reference_out.tsv>
suppressMessages({library(data.table); library(arrow); library(Seurat); library(Matrix); library(readr); library(dplyr)})
args <- commandArgs(trailingOnly = TRUE)
SS <- args[1]; GENES <- args[2]; DATADIR <- args[3]; YI_OUT <- args[4]

ss <- read_tsv(SS, show_col_types = FALSE) %>% filter(dataset == "Mutter_01")
common <- readLines(GENES)
counts_path <- unique(ss$counts_path); meta_path <- unique(ss$metadata_path)
stopifnot(length(counts_path) == 1, length(meta_path) == 1)

## ---- counts -> sparse genes x cells (common panel) ----
cdt <- as.data.table(read_parquet(counts_path))
stopifnot(!anyDuplicated(cdt$cell_id))
gene_cols <- intersect(setdiff(names(cdt), c("Slide", "fov", "cell_id")), common)
cat(sprintf("counts: %d cells x %d common genes\n", nrow(cdt), length(gene_cols)))
mat <- as.matrix(cdt[, ..gene_cols]); rownames(mat) <- cdt$cell_id   # cells x genes
cell_ids <- cdt$cell_id; rm(cdt); gc()
sp <- Matrix::t(Matrix::Matrix(mat, sparse = TRUE)); rm(mat); gc()   # genes x cells
stopifnot(identical(colnames(sp), cell_ids))

## ---- metadata: sequester Yi vendor cols; keep shared cols ----
mdt <- as.data.table(read_parquet(meta_path))
stopifnot(all(colnames(sp) %in% mdt$cell_id))                        # join completeness guard
setkey(mdt, cell_id)

yi_cols <- intersect(c("cell_id", "Condition", "Block", "qcFlagsCell",
  "insitutype_ImmuneAtlas_ImmGen.csv_clust", "ImmuneAtlas_ImmGen_Main_cell_Types",
  "ImmuneAtlas_ImmGen_cellFamily_Main_cell_Types", "RNA_snn_res.0.1", "umap_1", "umap_2",
  "TypeI_interferon", "TypeII_interferon", "STING", "DNA_Damage_Repair"), names(mdt))
dir.create(dirname(YI_OUT), recursive = TRUE, showWarnings = FALSE)
fwrite(mdt[, ..yi_cols], YI_OUT, sep = "\t")

shared <- intersect(c("nCount_RNA", "nFeature_RNA", "propNegative", "qcFlagsCell",
  "x_slide_mm", "y_slide_mm", "fov", "Area", "Mean.PanCK", "Mean.CD45",
  "Mean.CD298.B2M", "Mean.DAPI", "Condition"), names(mdt))
meta <- as.data.frame(mdt[, c("cell_id", ..shared)]); rownames(meta) <- meta$cell_id

## ---- split by Condition -> per-sample objects ----
stopifnot(all(ss$extract_key %in% meta$Condition))   # every sample's split key resolves
outdir <- file.path(DATADIR, "processing", "raw"); dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
for (i in seq_len(nrow(ss))) {
  s <- ss[i, ]
  cells <- meta$cell_id[meta$Condition == s$extract_key]
  stopifnot(length(cells) > 0)
  md <- meta[cells, , drop = FALSE]
  md$sample_id <- s$sample_id; md$dataset <- "Mutter_01"; md$treatment <- s$treatment
  md$condition <- s$condition; md$timepoint_h <- s$timepoint_h; md$model <- s$model
  obj <- CreateSeuratObject(counts = sp[, cells, drop = FALSE], meta.data = md, project = s$sample_id)
  saveRDS(obj, file.path(outdir, paste0(s$sample_id, ".raw.rds")))
  cat(sprintf("%s (%s): %d cells\n", s$sample_id, s$condition, ncol(obj)))
}
cat("adapt_mutter01 done\n")
