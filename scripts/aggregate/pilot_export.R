#!/usr/bin/env Rscript
# Integration bake-off pilot, step 1: subset the merged cohort to the two pilot
# slides (M01 sld0002 + M02 sld0005, ~27/71 dataset split mirroring the full
# cohort), attach the per-cell negprobe background, and emit both a Seurat object
# (Harmony arm) and an MTX + obs parquet handoff (scVI arm). Raw counts only --
# both arms do their own normalization so the only difference is the integrator.
# Args: <merged.rds> <cell_neg.tsv> <outdir>
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(arrow); library(data.table)
})

a        <- commandArgs(trailingOnly = TRUE)
merged   <- a[1]
neg_tsv  <- a[2]
outdir   <- a[3]
SLIDES   <- c("sld0002", "sld0005")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "mtx"), showWarnings = FALSE)

cat("reading merged object ...\n")
o <- readRDS(merged)
stopifnot(all(c("dataset", "slide_id", "cell_type", "sample_id") %in% colnames(o@meta.data)))

o <- subset(o, subset = slide_id %in% SLIDES)
o <- JoinLayers(o)                                   # consolidate split per-sample layers
cat(sprintf("subset to %s -> %d cells x %d genes\n",
            paste(SLIDES, collapse = "+"), ncol(o), nrow(o)))

# per-cell negprobe background (continuous covariate for scVI)
neg <- as.data.table(fread(neg_tsv))
nv  <- setNames(neg$neg, neg$cell)
o$neg <- unname(nv[colnames(o)])
miss <- sum(is.na(o$neg))
if (miss > 0) { med <- median(o$neg, na.rm = TRUE); o$neg[is.na(o$neg)] <- med
  cat(sprintf("filled %d missing neg with median %.3f\n", miss, med)) }

saveRDS(o, file.path(outdir, "pilot_full.rds"))

# scVI handoff: raw counts MTX + features/barcodes + obs parquet
cts <- LayerData(o, assay = "RNA", layer = "counts")
Matrix::writeMM(cts, file.path(outdir, "mtx", "counts.mtx"))
writeLines(rownames(cts), file.path(outdir, "mtx", "features.tsv"))
writeLines(colnames(cts), file.path(outdir, "mtx", "barcodes.tsv"))

obs <- data.table(
  cell        = colnames(o),
  sample_id   = as.character(o$sample_id),
  dataset     = as.character(o$dataset),
  slide_id    = as.character(o$slide_id),
  cell_type   = as.character(o$cell_type),
  condition   = as.character(o$condition),
  treatment   = as.character(o$treatment),
  timepoint_h = as.character(o$timepoint_h),
  neg         = as.numeric(o$neg))
arrow::write_parquet(obs, file.path(outdir, "obs.parquet"))

cat("\n--- pilot composition ---\n")
print(obs[, .N, by = .(dataset, slide_id)][order(dataset, slide_id)])
cat(sprintf("\nneg: mean %.3f median %.3f range [%.3f, %.3f]\n",
            mean(obs$neg), median(obs$neg), min(obs$neg), max(obs$neg)))
cat(sprintf("wrote: pilot_full.rds, mtx/counts.mtx (%d x %d), obs.parquet\n",
            nrow(cts), ncol(cts)))
