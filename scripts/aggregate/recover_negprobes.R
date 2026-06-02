#!/usr/bin/env Rscript
# Rebuild the per-cell mean negative-probe count ("neg") that InSituType needs as
# its background model -- dropped from the 950-gene common panel during the
# adapter, so it must come from the truly-raw inputs:
#   M01: metadata parquet column nCount_NegativeProbes (keyed by cell_id, which is
#        identical to the M01 scored barcode).
#   M02: raw per-slide Seurat objects carry a `negprobes` assay (10 probes) and a
#        meta column nCount_negprobes (keyed by the scored barcode, but the c_1_*
#        barcode is object-local so we resolve each sample's slide by coverage).
# Output barcodes match merge.R's globally-unique key: {sample_id}_{raw_barcode},
# where sample_id is the scored-file stem (merge.R uses the same stem).
# neg = nCount_neg / NNEG (mean over the panel's negative probes); same panel for
# both datasets so one definition spans M01/M02.
# Args: <out_neg.tsv> <m01_metadata.parquet> <m02_rawdir> <scored.rds...>
suppressPackageStartupMessages({library(Seurat); library(Matrix); library(arrow); library(data.table)})
a          <- commandArgs(trailingOnly = TRUE)
out_tsv    <- a[1]
m01_pq     <- a[2]
m02_dir    <- a[3]
scored     <- a[-(1:3)]
NNEG       <- 10                      # NanoString mouse 1k panel: 10 negative probes

# ---- M01 neg lookup: cell_id -> neg ----
m01 <- as.data.table(arrow::read_parquet(m01_pq,
        col_select = c("cell_id", "nCount_NegativeProbes")))
m01[, neg := nCount_NegativeProbes / NNEG]
m01_neg <- setNames(m01$neg, m01$cell_id)
cat(sprintf("M01 parquet: %d cells, neg range [%.3f, %.3f]\n",
            length(m01_neg), min(m01_neg), max(m01_neg)))

# ---- M02 neg lookups: one named vector per raw slide object ----
m02_files <- sort(list.files(m02_dir, pattern = "Mutter_02_CosMmR\\.RDS$", full.names = TRUE))
slide_neg <- vector("list", length(m02_files))
for (i in seq_along(m02_files)) {
  o <- readRDS(m02_files[i])
  nc <- if ("nCount_negprobes" %in% colnames(o@meta.data)) o$nCount_negprobes
        else Matrix::colSums(LayerData(o, assay = "negprobes", layer = "counts"))
  slide_neg[[i]] <- setNames(as.numeric(nc) / NNEG, colnames(o))
  cat(sprintf("M02 %s: %d cells\n", basename(m02_files[i]), length(slide_neg[[i]])))
  rm(o); invisible(gc())
}

# ---- walk each scored sample, emit merged-barcode -> neg ----
out <- vector("list", length(scored))
for (k in seq_along(scored)) {
  o   <- readRDS(scored[k])
  sid <- sub("\\.scored\\.rds$", "", basename(scored[k]))     # merge.R's prefix
  ds  <- unique(as.character(o$dataset))[1]
  bc  <- colnames(o)
  if (grepl("01", ds)) {
    neg <- unname(m01_neg[bc])
    src <- "M01"
  } else {
    cov <- vapply(slide_neg, function(v) mean(bc %in% names(v)), numeric(1))
    j   <- which.max(cov)
    if (cov[j] < 0.99) stop(sprintf("%s: best M02 slide coverage only %.3f", sid, cov[j]))
    neg <- unname(slide_neg[[j]][bc])
    src <- sprintf("M02_slide%d", j)
  }
  miss <- sum(is.na(neg))
  out[[k]] <- data.table(cell = paste0(sid, "_", bc), neg = neg)
  cat(sprintf("%-10s %-12s %7d cells  neg NA=%d\n", sid, src, length(bc), miss))
}
res <- rbindlist(out)
# safe fallback for any unmatched cell: dataset-agnostic median (rare; reported above)
if (anyNA(res$neg)) { med <- median(res$neg, na.rm = TRUE); res[is.na(neg), neg := med]
  cat(sprintf("filled %d NA neg with global median %.3f\n", sum(is.na(res$neg)), med)) }
dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)
fwrite(res, out_tsv, sep = "\t")
cat(sprintf("\nwrote %d cells -> %s | neg mean %.3f median %.3f\n",
            nrow(res), out_tsv, mean(res$neg), median(res$neg)))
