#!/usr/bin/env Rscript
# Per-cell slide coordinates + a per-sample necrosis flag for the spatial tracks.
# Each sample_id is one contiguous tissue region (a physical slide carries several
# non-overlapping M02 regions), so the necrosis rule is applied within sample: a
# cell is necrosis_zone if its mean distance to its 20 nearest neighbours exceeds
# the 90th percentile of that distance within its sample (sparse neighbourhoods =
# necrotic/acellular). Re-keys the dev necrosis rule (dev/mbrt_vs_sbrt/
# 00_load_and_filter.R) to sample-level coordinate space. The global cell id
# matches merge.R: paste0(sample_id, "_", colnames(counts)).
# Args: <samples.tsv> <out_coords_necrosis.parquet> <scored.rds...>
suppressPackageStartupMessages({library(Seurat); library(data.table); library(arrow); library(RANN)})
a <- commandArgs(trailingOnly = TRUE)
out <- a[2]; rds <- a[-(1:2)]
res <- rbindlist(lapply(rds, function(p){
  s <- sub("\\.scored\\.rds$", "", basename(p)); o <- readRDS(p); m <- o@meta.data
  d <- data.table(cell = paste0(s, "_", rownames(m)), sample_id = s,
                  x_slide_mm = m$x_slide_mm, y_slide_mm = m$y_slide_mm)
  xy <- as.matrix(d[, .(x_slide_mm, y_slide_mm)])
  nn <- RANN::nn2(xy, k = min(21, nrow(xy)))
  d[, mean_knn_dist := rowMeans(nn$nn.dists[, -1, drop = FALSE])]
  d[, necrosis_zone := mean_knn_dist > quantile(mean_knn_dist, 0.90)]
  d
}))
write_parquet(res, out)
cat(sprintf("coords: %d cells, %.1f%% necrosis-flagged\n", nrow(res), 100*mean(res$necrosis_zone)))
