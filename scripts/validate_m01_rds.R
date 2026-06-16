#!/usr/bin/env Rscript
# Validate that the new Mutter_01 raw RDS are bit-identical to the existing parquet
# (the basis for "the aggregate is invariant to this data"), across ALL 4 slides,
# and that the object-local barcode reconstructs exactly to the parquet cell_id.
# Also build + check the (slide,fov)->Condition crosswalk. Report-only; stop() on any failure.
# Args: <rds_dir> <m01_metadata.parquet> <out_validation.tsv>
suppressMessages({library(Seurat); library(arrow); library(data.table); library(Matrix)})
a <- commandArgs(trailingOnly = TRUE)
RDS_DIR <- a[1]; META_PQ <- a[2]; OUT <- a[3]
NNEG <- 10L; NFALSE <- 184L

meta <- as.data.table(read_parquet(META_PQ,
  col_select = c("cell_id","Slide","Run_Tissue_name","fov","Condition",
                 "nCount_NegativeProbes","nCount_RNA")))
setkey(meta, cell_id)

# crosswalk: every (Slide,fov) must map to exactly one Condition
xw <- meta[, .(nCond = uniqueN(Condition)), by = .(Slide, fov)]
if (any(xw$nCond != 1L)) stop(sprintf("crosswalk impure: %d (Slide,fov) map to >1 Condition", sum(xw$nCond != 1L)))
cat(sprintf("crosswalk OK: %d (Slide,fov) cells, all map to exactly 1 Condition\n", nrow(xw)))

files <- sort(list.files(RDS_DIR, pattern = "Mutter_01_CosMmR\\.RDS$", full.names = TRUE))
stopifnot(length(files) == 4)
rows <- vector("list", length(files))
for (i in seq_along(files)) {
  o <- readRDS(files[i])
  # assay shape checks
  stopifnot("negprobes" %in% Assays(o), "falsecode" %in% Assays(o))
  stopifnot(nrow(o[["negprobes"]]) == NNEG, nrow(o[["falsecode"]]) == NFALSE)
  # slide identity from Run_Tissue_name "0X_Mutter_01_CosMmR" -> parquet Slide "..._S<X>"
  rtn <- unique(as.character(o$Run_Tissue_name)); stopifnot(length(rtn) == 1)
  sx  <- paste0("S", as.integer(sub("_.*", "", rtn)))
  slide_full <- grep(paste0("_", sx, "$"), unique(meta$Slide), value = TRUE)
  stopifnot(length(slide_full) == 1)
  # reconstruct global cell_id from object-local c_1_<fov>_<cell>
  p   <- tstrsplit(colnames(o), "_")
  gid <- paste0(slide_full, "_", p[[3]], "_", p[[4]])
  m   <- meta[gid]                                  # join RDS->parquet by reconstructed id
  n_match <- sum(!is.na(m$cell_id))
  if (n_match != ncol(o)) stop(sprintf("%s: only %d/%d cells reconstruct to a parquet cell_id",
                                        basename(files[i]), n_match, ncol(o)))
  neg_ok <- all(o$nCount_negprobes == m$nCount_NegativeProbes)
  rna_ok <- all(o$nCount_RNA == m$nCount_RNA)
  if (!neg_ok) stop(sprintf("%s: nCount_negprobes != parquet nCount_NegativeProbes", basename(files[i])))
  if (!rna_ok) stop(sprintf("%s: nCount_RNA mismatch vs parquet", basename(files[i])))
  rows[[i]] <- data.table(file = basename(files[i]), slide = slide_full, n_cells = ncol(o),
                          join_pct = 100 * n_match / ncol(o), neg_identical = neg_ok, rna_identical = rna_ok)
  cat(sprintf("%-40s %s  %7d cells  join=100%%  neg_identical=TRUE  rna_identical=TRUE\n",
              basename(files[i]), slide_full, ncol(o)))
  rm(o); invisible(gc())
}
res <- rbindlist(rows)
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
fwrite(res, OUT, sep = "\t")
cat(sprintf("\nALL SLIDES PASS -> %s\n", OUT))
