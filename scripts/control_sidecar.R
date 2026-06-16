#!/usr/bin/env Rscript
# Per-cell control sidecar (negprobe + falsecode) for BOTH datasets, from raw RDS.
# ADDITIVE: writes NEW files; never the aggregate's cell_neg.tsv. M01 cells keyed by the
# reconstructed parquet cell_id (joinable downstream); join must be 100% or stop().
# The per-cell table is ~4.2M rows -> written as PARQUET to the datadir (heavy intermediate,
# not git); only the small per-FOV QC summary is a committed TSV. Repo convention: per-cell
# tables are parquet (cf. full_labels.parquet), small summaries are TSV in results/.
# Plus a per-(dataset,slide,fov) falsecode QC flag (3-MAD outlier, report-only).
# Args: <m01_dir> <m01_meta.parquet> <m02_dir> <out_cell.parquet> <out_fov.tsv>
suppressMessages({library(Seurat); library(arrow); library(data.table); library(Matrix)})
a <- commandArgs(trailingOnly = TRUE)
M01_DIR <- a[1]; M01_META <- a[2]; M02_DIR <- a[3]; OUT_CELL <- a[4]; OUT_FOV <- a[5]
NNEG <- 10; NFALSE <- 184
stopifnot(basename(OUT_CELL) != "cell_neg.tsv")          # never collide with the aggregate input

cells_of <- function(o, ds, key) {
  nc_neg   <- if ("nCount_negprobes" %in% colnames(o@meta.data)) as.numeric(o$nCount_negprobes)
             else as.numeric(Matrix::colSums(LayerData(o, assay = "negprobes", layer = "counts")))
  nc_false <- if ("nCount_falsecode" %in% colnames(o@meta.data)) as.numeric(o$nCount_falsecode)
             else as.numeric(Matrix::colSums(LayerData(o, assay = "falsecode", layer = "counts")))
  nc_rna   <- as.numeric(o$nCount_RNA)
  data.table(cell_id = key, dataset = ds,
             slide = unique(as.character(o$Run_Tissue_name))[1], fov = as.integer(o$fov),
             neg = nc_neg / NNEG, falsecode = nc_false / NFALSE,
             false_frac = nc_false / (nc_rna + nc_neg + nc_false + 1e-9))
}

# ---- M01: reconstruct global cell_id, assert 100% join to parquet ----
meta <- as.data.table(read_parquet(M01_META, col_select = c("cell_id","Slide")))
m01_files <- sort(list.files(M01_DIR, pattern = "Mutter_01_CosMmR\\.RDS$", full.names = TRUE))
m01 <- vector("list", length(m01_files))
for (i in seq_along(m01_files)) {
  o <- readRDS(m01_files[i])
  sx <- paste0("S", as.integer(sub("_.*", "", unique(as.character(o$Run_Tissue_name)))))
  slide_full <- grep(paste0("_", sx, "$"), unique(meta$Slide), value = TRUE); stopifnot(length(slide_full)==1)
  p <- tstrsplit(colnames(o), "_")
  gid <- paste0(slide_full, "_", p[[3]], "_", p[[4]])
  n_in <- sum(gid %in% meta$cell_id)
  if (n_in != ncol(o)) stop(sprintf("%s: %d/%d cells fail parquet-id reconstruction (no NA-fill allowed)",
                                     basename(m01_files[i]), n_in, ncol(o)))
  m01[[i]] <- cells_of(o, "Mutter_01", gid); rm(o); invisible(gc())
}
# ---- M02: native object-local barcode (no parquet); slide-scoped unique key ----
m02_files <- sort(list.files(M02_DIR, pattern = "Mutter_02_CosMmR\\.RDS$", full.names = TRUE))
m02 <- vector("list", length(m02_files))
for (j in seq_along(m02_files)) {
  o <- UpdateSeuratObject(readRDS(m02_files[j]))
  key <- paste0(unique(as.character(o$Run_Tissue_name))[1], ":", colnames(o))   # slide-scoped unique
  m02[[j]] <- cells_of(o, "Mutter_02", key); rm(o); invisible(gc())
}
cell <- rbindlist(c(m01, m02))
stopifnot(!anyNA(cell$false_frac), !anyNA(cell$neg))
dir.create(dirname(OUT_CELL), recursive = TRUE, showWarnings = FALSE)
arrow::write_parquet(cell, OUT_CELL)                      # ~4.2M rows -> parquet, not TSV

# ---- per-FOV falsecode QC flag: 3-MAD outlier on mean false_frac, within dataset (report-only) ----
fov <- cell[, .(n_cells = .N, mean_false_frac = mean(false_frac)), by = .(dataset, slide, fov)]
fov[, `:=`(med = median(mean_false_frac), mad = mad(mean_false_frac)), by = dataset]
fov[, falsecode_flag := mean_false_frac > med + 3 * mad][, c("med","mad") := NULL]
dir.create(dirname(OUT_FOV), recursive = TRUE, showWarnings = FALSE)
fwrite(fov, OUT_FOV, sep = "\t")
cat(sprintf("wrote %d cells -> %s\n", nrow(cell), OUT_CELL))
cat(sprintf("per-FOV QC: %d FOVs, %d flagged (report-only) -> %s\n",
            nrow(fov), sum(fov$falsecode_flag), OUT_FOV))
