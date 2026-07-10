#!/usr/bin/env Rscript
# Full-cohort scVI handoff (pilot promoted to full launch, 2026-06-01). Same shape as
# pilot_export.R but on the whole flank cohort (3.27M cells, 7 slides; tongue already
# excluded from merged.rds). Attaches the per-cell negprobe background and emits raw
# counts MTX + features/barcodes + obs parquet for the scVI arm. No rds save -- the R
# annotation step re-reads merged.rds.
# Args: <merged.rds> <cell_neg.tsv> <outdir>
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(arrow); library(data.table)
})

a       <- commandArgs(trailingOnly = TRUE)
merged  <- a[1]; neg_tsv <- a[2]; outdir <- a[3]
dir.create(file.path(outdir, "mtx"), recursive = TRUE, showWarnings = FALSE)

cat("reading merged object ...\n")
o <- readRDS(merged)
stopifnot(all(c("dataset", "slide_id", "sample_id") %in% colnames(o@meta.data)))
o <- JoinLayers(o)
cat(sprintf("cohort: %d cells x %d genes\n", ncol(o), nrow(o)))

neg <- as.data.table(fread(neg_tsv))
nv  <- setNames(neg$neg, neg$cell)
o$neg <- unname(nv[colnames(o)])
miss <- sum(is.na(o$neg))
if (miss > 0) { med <- median(o$neg, na.rm = TRUE); o$neg[is.na(o$neg)] <- med
  cat(sprintf("filled %d missing neg with median %.3f\n", miss, med)) }

cts <- LayerData(o, assay = "RNA", layer = "counts")
Matrix::writeMM(cts, file.path(outdir, "mtx", "counts.mtx"))
writeLines(rownames(cts), file.path(outdir, "mtx", "features.tsv"))
writeLines(colnames(cts), file.path(outdir, "mtx", "barcodes.tsv"))

obs <- data.table(
  cell        = colnames(o),
  sample_id   = as.character(o$sample_id),
  dataset     = as.character(o$dataset),
  slide_id    = as.character(o$slide_id),
  cell_type   = NA_character_,  # per-sample typing retired; benchmark label re-based at merged scale
  condition   = as.character(o$condition),
  treatment   = as.character(o$treatment),
  timepoint_h = as.character(o$timepoint_h),
  neg         = as.numeric(o$neg))
arrow::write_parquet(obs, file.path(outdir, "obs.parquet"))

cat("\n--- cohort composition ---\n")
print(obs[, .N, by = .(dataset, slide_id)][order(dataset, slide_id)])
print(obs[, .(cells = .N, frac = round(.N / nrow(obs), 4)), by = dataset])
cat(sprintf("wrote: mtx/counts.mtx (%d x %d), obs.parquet, neg mean %.3f\n",
            nrow(cts), ncol(cts), mean(obs$neg)))
