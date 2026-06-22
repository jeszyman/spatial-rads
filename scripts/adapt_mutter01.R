#!/usr/bin/env Rscript
# Stage B adapter -- Mutter_01: per-slide raw Seurat RDS -> per-sample common-format Seurat.
# Harmonized with adapt_mutter02 on the per-slide-RDS data structure (M01 re-delivered as raw RDS
# with negprobes+falsecode; Yi 2026-06-10: only the raw object exists). RNA counts now come from the
# 4 per-slide RDS (bit-identical to the legacy MC-SOLVE counts parquet); Condition + Yi vendor labels
# + shared morphology/QC metadata are joined from the metadata parquet by reconstructed global
# cell_id (the parquet is their only source). Yi vendor cols -> yi_reference.tsv; objects split by
# Condition. A setequal() guard asserts the RDS cohort reproduces the parquet cell set exactly, so
# per-sample output is numerically INVARIANT vs the prior parquet-counts adapter.
# Args: <samplesheet.tsv> <common_genes.tsv> <DATADIR> <yi_reference_out.tsv>
suppressMessages({library(data.table); library(arrow); library(Seurat); library(Matrix); library(readr); library(dplyr)})
args <- commandArgs(trailingOnly = TRUE)
SS <- args[1]; GENES <- args[2]; DATADIR <- args[3]; YI_OUT <- args[4]

ss <- read_tsv(SS, show_col_types = FALSE) %>% filter(dataset == "Mutter_01")
common    <- readLines(GENES)
meta_path <- unique(ss$metadata_path); stopifnot(length(meta_path) == 1)
rds_paths <- sort(unique(ss$raw_input_path))                         # 4 per-slide RDS (was the counts parquet)

## ---- metadata parquet: Condition + Yi vendor + shared cols live ONLY here, keyed by cell_id ----
mdt <- as.data.table(read_parquet(meta_path)); setkey(mdt, cell_id)

## ---- counts from the per-slide RDS -> sparse genes x cells on the common panel ----
## reconstruct the global parquet cell_id from each object-local barcode c_1_<fov>_<cell> (scheme
## proven in validate_m01_rds.R / control_sidecar.R); assert 100% reconstruction (no NA-fill).
mats <- vector("list", length(rds_paths))
for (i in seq_along(rds_paths)) {
  o   <- readRDS(rds_paths[i])
  rtn <- unique(as.character(o$Run_Tissue_name)); stopifnot(length(rtn) == 1)
  sx  <- paste0("S", as.integer(sub("_.*", "", rtn)))
  slide_full <- grep(paste0("_", sx, "$"), unique(mdt$Slide), value = TRUE); stopifnot(length(slide_full) == 1)
  p   <- tstrsplit(colnames(o), "_")
  gid <- paste0(slide_full, "_", p[[3]], "_", p[[4]])
  stopifnot(sum(gid %in% mdt$cell_id) == ncol(o))
  cm  <- LayerData(o, assay = "RNA", layer = "counts"); colnames(cm) <- gid
  mats[[i]] <- cm[intersect(common, rownames(cm)), , drop = FALSE]
  cat(sprintf("%-42s %s  %7d cells\n", basename(rds_paths[i]), slide_full, ncol(cm)))
  rm(o, cm); invisible(gc())
}
genes_common <- Reduce(intersect, lapply(mats, rownames))
sp <- do.call(cbind, lapply(mats, function(m) m[genes_common, , drop = FALSE])); rm(mats); invisible(gc())
stopifnot(setequal(colnames(sp), mdt$cell_id))                       # RDS cohort == parquet cohort -> invariance
cat(sprintf("counts (from RDS): %d cells x %d common genes\n", ncol(sp), nrow(sp)))

## ---- metadata: sequester Yi vendor cols; keep shared cols (from the parquet, as before) ----
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
