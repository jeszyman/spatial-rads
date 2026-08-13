#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table); library(arrow); library(Seurat); library(Matrix); library(InSituCor); library(RANN)
})
a <- commandArgs(trailingOnly = TRUE)
rds_path <- a[1]; labels_path <- a[2]; coords_path <- a[3]; obs_path <- a[4]
out_modules <- a[5]; out_summary <- a[6]
KNN_COR <- 100L; MAX_CELLS <- 100000L; SEED <- 42L

dir.create(dirname(out_modules), recursive = TRUE, showWarnings = FALSE)

seu <- readRDS(rds_path)
lab <- as.data.table(read_parquet(labels_path))[, .(cell, cell_subtype)]
co  <- as.data.table(read_parquet(coords_path))[, .(cell, sample_id, x_slide_mm, y_slide_mm)]
ob  <- as.data.table(read_parquet(obs_path))[, .(cell, timepoint_h, condition, dataset)]

meta <- lab[co, on = "cell"]
meta <- ob[meta, on = "cell"]
meta <- meta[!is.na(cell_subtype) & timepoint_h == 4]

cells_4h <- intersect(meta$cell, colnames(seu))
counts <- GetAssayData(seu, assay = "RNA", layer = "counts")[, cells_4h]

# normalize: linear scale (count / total, floor 20) per site recommendation
totalcounts <- Matrix::colSums(counts)
norm <- sweep(counts, 2, pmax(totalcounts, 20), "/")

meta <- meta[cell %in% cells_4h]
setkey(meta, cell)
meta <- meta[colnames(norm)]

# subsample if needed — report neighbor-radius change so spatial scale is documented
if (ncol(norm) > MAX_CELLS) {
  set.seed(SEED)
  keep <- sample.int(ncol(norm), MAX_CELLS)
  norm <- norm[, keep]
  meta <- meta[keep]
  cat(sprintf("insitucor: subsampled %d -> %d cells\n", ncol(norm) + length(setdiff(seq_len(ncol(norm)), keep)), MAX_CELLS))
}
# report the k=100 neighbor radius after any subsampling
coords_pre <- as.matrix(meta[, .(x_slide_mm, y_slide_mm)])
nn_check <- RANN::nn2(coords_pre, k = min(KNN_COR + 1L, nrow(coords_pre)))
k100_dists <- nn_check$nn.dists[, ncol(nn_check$nn.dists)]
cat(sprintf("insitucor: k=%d neighbor radius (mm): median=%.4f, 95th=%.4f, max=%.4f\n",
            KNN_COR, median(k100_dists), quantile(k100_dists, 0.95), max(k100_dists)))

# build kNN graph on spatial coordinates, split by sample
coords <- as.matrix(meta[, .(x_slide_mm, y_slide_mm)])
rownames(coords) <- meta$cell

# InSituCor: k=100, defaults from site
isc <- insitucor(
  counts    = norm,
  condvar   = meta$cell_subtype,
  neighbors = nearestNeighborGraph(coords, k = KNN_COR,
                                    splitby = meta$sample_id),
  k                   = KNN_COR,
  min_module_size     = 2,
  max_module_size     = 25,
  min_module_cor      = 0.1,
  gene_weighting_rule = "inverse_sqrt",
  roundcortozero      = 0.1
)

# extract modules
mods <- isc$modules
if (length(mods) == 0) {
  cat("insitucor: 0 modules found\n")
  fwrite(data.table(module = character(), gene = character(),
                    gene_weight = numeric(), module_size = integer()), out_modules, sep = "\t")
  fwrite(data.table(module = character(), n_genes = integer(),
                    top_genes = character()), out_summary, sep = "\t")
} else {
  mod_dt <- rbindlist(lapply(seq_along(mods), function(i) {
    data.table(module = paste0("M", i),
               gene = names(mods[[i]]),
               gene_weight = as.numeric(mods[[i]]),
               module_size = length(mods[[i]]))
  }))
  fwrite(mod_dt, out_modules, sep = "\t")

  summ <- mod_dt[, .(n_genes = .N,
                      top_genes = paste(head(gene[order(-abs(gene_weight))], 5), collapse = ", ")),
                 by = module]
  fwrite(summ, out_summary, sep = "\t")
  cat(sprintf("insitucor: %d modules, %d total genes | top modules: %s\n",
              nrow(summ), nrow(mod_dt),
              paste(head(summ$module), head(summ$n_genes), sep = "=", collapse = ", ")))
}
