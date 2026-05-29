#!/usr/bin/env Rscript
# Shared per-sample QC -- identical four-criteria gate on both datasets, all from data in hand:
#   counts (nCount_RNA), complexity (nFeature_RNA), negative-probe proportion (propNegative),
#   and cell-area outliers (robust MAD on log10 Area). This reproduces what NanoString's
#   composite qcFlagsCell bundles, transparently and identically across datasets.
# Args: <in.raw.rds> <out.qc.rds> <out.summary.tsv> <min_counts> <min_features> <max_propneg> [area_nmads=3]
suppressMessages({library(Seurat); library(readr)})
args <- commandArgs(trailingOnly = TRUE)
IN <- args[1]; OUT_RDS <- args[2]; OUT_SUM <- args[3]
min_counts <- as.numeric(args[4]); min_features <- as.numeric(args[5]); max_propneg <- as.numeric(args[6])
area_nmads <- if (length(args) >= 7) as.numeric(args[7]) else 3

obj <- readRDS(IN)
req <- c("nCount_RNA", "nFeature_RNA", "propNegative", "Area")
miss <- setdiff(req, colnames(obj@meta.data))
if (length(miss)) stop(sprintf("%s: missing required QC columns: %s", basename(IN), paste(miss, collapse = ", ")))

pre <- ncol(obj)
la  <- log10(obj$Area)                                   # robust area-outlier gate (identical method both datasets)
med <- median(la, na.rm = TRUE); md <- mad(la, na.rm = TRUE)
area_ok <- !is.na(la) & la > (med - area_nmads * md) & la < (med + area_nmads * md)
keep <- obj$nCount_RNA > min_counts & obj$nFeature_RNA > min_features &
        obj$propNegative < max_propneg & area_ok

obj$qcFlagsCell <- NULL                                  # drop vendor flag (M01-only) so qc objects share one schema
obj <- obj[, keep]
post <- ncol(obj)

sid <- if ("sample_id" %in% colnames(obj@meta.data)) as.character(obj$sample_id[1]) else sub("\\.raw\\.rds$", "", basename(IN))
summ <- data.frame(sample_id = sid, pre_filter = pre, post_filter = post,
                   pct_retained = round(100 * post / pre, 1), area_flagged = sum(!area_ok))
dir.create(dirname(OUT_RDS), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(OUT_SUM), recursive = TRUE, showWarnings = FALSE)
saveRDS(obj, OUT_RDS)
write_tsv(summ, OUT_SUM)
cat(sprintf("%s QC: %d -> %d (%.1f%%; %d area-flagged)\n", sid, pre, post, 100 * post / pre, sum(!area_ok)))
