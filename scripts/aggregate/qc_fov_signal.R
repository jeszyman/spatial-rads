#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table); library(arrow); library(Seurat)
})
a <- commandArgs(trailingOnly = TRUE)
rds_path <- a[1]; labels_path <- a[2]; obs_path <- a[3]; out_tsv <- a[4]
SIGNAL_LOSS_THRESHOLD <- 0.4

dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)

seu <- readRDS(rds_path)
meta <- as.data.table(seu@meta.data, keep.rownames = "cell")
lab <- as.data.table(read_parquet(labels_path))[, .(cell)]
ob  <- as.data.table(read_parquet(obs_path))[, .(cell, condition, slide_id, dataset)]
d <- meta[cell %in% lab$cell, .(cell, nCount_RNA, fov)]
d <- ob[d, on = "cell"]

per_fov <- d[, .(mean_count = mean(nCount_RNA),
                 median_count = median(nCount_RNA),
                 n_cells = .N),
             by = .(dataset, slide_id, fov, condition)]

per_fov[, slide_median_fov_mean := median(mean_count), by = .(dataset, slide_id)]
per_fov[, ratio_to_slide_median := mean_count / slide_median_fov_mean]
per_fov[, flagged := ratio_to_slide_median < SIGNAL_LOSS_THRESHOLD]

setorder(per_fov, dataset, slide_id, fov)
fwrite(per_fov, out_tsv, sep = "\t")

n_flagged <- sum(per_fov$flagged)
if (n_flagged > 0) {
  arm_dist <- per_fov[flagged == TRUE, .N, by = condition]
  cat(sprintf("fov_signal_qc: %d / %d FOVs flagged (>60%% signal loss) | arm distribution: %s\n",
              n_flagged, nrow(per_fov),
              paste(arm_dist$condition, arm_dist$N, sep = "=", collapse = ", ")))
} else {
  cat(sprintf("fov_signal_qc: 0 / %d FOVs flagged\n", nrow(per_fov)))
}
