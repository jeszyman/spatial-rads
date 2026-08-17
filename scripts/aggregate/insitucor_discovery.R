#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table); library(arrow); library(Seurat); library(Matrix); library(InSituCor); library(RANN)
})
a <- commandArgs(trailingOnly = TRUE)
rds_path <- a[1]; labels_path <- a[2]; coords_path <- a[3]; obs_path <- a[4]
out_modules <- a[5]; out_summary <- a[6]
out_ct_attr <- a[7]; out_gene_attr <- a[8]
out_scores_sc <- a[9]; out_scores_env <- a[10]; out_condcor <- a[11]
# MAX_CELLS bounds peak memory and is applied before the counts matrix is built or
# normalized (not after, as previously). At 100k cells, the ~950-gene sparse counts matrix
# (~6% detection), the k=100 sparse cell-cell neighbor graph, and InSituCor's own dense
# per-module score matrices are all well under 1 GB, comfortably inside the shared-machine
# ceiling.
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

# subsample before touching the counts matrix, so no full-cohort-scale intermediate is ever
# built (the 4h cohort spans multiple datasets and can exceed 1M cells)
if (length(cells_4h) > MAX_CELLS) {
  set.seed(SEED)
  n_total <- length(cells_4h)
  cells_4h <- cells_4h[sample.int(n_total, MAX_CELLS)]
  cat(sprintf("insitucor: subsampled %d -> %d cells\n", n_total, MAX_CELLS))
}

counts <- GetAssayData(seu, assay = "RNA", layer = "counts")[, cells_4h]

meta <- meta[cell %in% cells_4h]
setkey(meta, cell)
meta <- meta[cells_4h]

# normalize: linear scale (count / total, floor 20) per site recommendation. Sparse diagonal
# scaling (not sweep(), which builds a dense full-size STATS array and forces a sparse->dense
# coercion) keeps the matrix sparse, and transposes genes x cells to InSituCor's expected
# cells x genes orientation.
totalcounts <- Matrix::colSums(counts)
norm <- Matrix::Diagonal(x = 1 / pmax(totalcounts, 20)) %*% Matrix::t(counts)
rownames(norm) <- colnames(counts); colnames(norm) <- rownames(counts)

# report the k=100 neighbor radius (diagnostic only; InSituCor builds its own graph below)
coords_pre <- as.matrix(meta[, .(x_slide_mm, y_slide_mm)])
nn_check <- RANN::nn2(coords_pre, k = min(KNN_COR + 1L, nrow(coords_pre)))
k100_dists <- nn_check$nn.dists[, ncol(nn_check$nn.dists)]
cat(sprintf("insitucor: k=%d neighbor radius (mm): median=%.4f, 95th=%.4f, max=%.4f\n",
            KNN_COR, median(k100_dists), quantile(k100_dists, 0.95), max(k100_dists)))

coords <- as.matrix(meta[, .(x_slide_mm, y_slide_mm)])
rownames(coords) <- meta$cell
conditionon <- data.frame(cell_subtype = meta$cell_subtype, row.names = meta$cell)

# InSituCor: k=100, neighbor search split per sample via tissue=. neighbors is left NULL —
# insitucor() builds its own graph from xy/k/tissue internally; nearestNeighborGraph() is not
# exported by the package, so it cannot be called directly from this script.
isc <- insitucor(
  counts              = norm,
  conditionon         = conditionon,
  celltype            = meta$cell_subtype,
  xy                  = coords,
  k                   = KNN_COR,
  tissue              = meta$sample_id,
  min_module_size     = 2,
  max_module_size     = 25,
  min_module_cor      = 0.1,
  gene_weighting_rule = "inverse_sqrt",
  roundcortozero      = 0.1
)

# extract modules: isc$modules is a data.table with one row per gene (module, gene, weight)
mod_dt <- copy(isc$modules)
if (nrow(mod_dt) == 0) {
  cat("insitucor: 0 modules found\n")
  fwrite(data.table(module = character(), gene = character(),
                    gene_weight = numeric(), module_size = integer()), out_modules, sep = "\t")
  fwrite(data.table(module = character(), n_genes = integer(),
                    top_genes = character()), out_summary, sep = "\t")
  fwrite(data.table(module = character(), cell_type = character(),
                    attribution = numeric()), out_ct_attr, sep = "\t")
  fwrite(data.table(module = character(), gene = character(), cell_type = character(),
                    attribution = numeric()), out_gene_attr, sep = "\t")
  write_parquet(data.frame(cell = character()), out_scores_sc)
  write_parquet(data.frame(cell = character()), out_scores_env)
  fwrite(data.table(gene_a = character(), gene_b = character(),
                    cond_cor = numeric()), out_condcor, sep = "\t")
  quit(save = "no", status = 0)
}
setnames(mod_dt, "weight", "gene_weight")
mod_dt[, module_size := .N, by = module]
fwrite(mod_dt, out_modules, sep = "\t")

summ <- mod_dt[, .(n_genes = .N,
                    top_genes = paste(head(gene[order(-abs(gene_weight))], 5), collapse = ", ")),
               by = module]
fwrite(summ, out_summary, sep = "\t")
cat(sprintf("insitucor: %d modules, %d total genes | top modules: %s\n",
            nrow(summ), nrow(mod_dt),
            paste(head(summ$module), head(summ$n_genes), sep = "=", collapse = ", ")))

# module-level cell type attribution: celltypeinvolvement holds the max attribution score
# of each cell type over each module's genes. Orientation is asserted against the module
# names rather than assumed (cell type names and module names never collide).
cti <- as.matrix(isc$celltypeinvolvement)
if (all(rownames(cti) %in% mod_dt$module)) cti <- t(cti)
stopifnot(all(colnames(cti) %in% mod_dt$module))
ct_attr <- as.data.table(as.table(cti))
setnames(ct_attr, c("cell_type", "module", "attribution"))
setcolorder(ct_attr, c("module", "cell_type", "attribution"))
setorder(ct_attr, module, -attribution)
fwrite(ct_attr, out_ct_attr, sep = "\t")

# gene-level attribution: one cell_type x gene matrix per module
gene_attr <- rbindlist(lapply(names(isc$attributionmats), function(m) {
  am <- as.matrix(isc$attributionmats[[m]])
  if (all(rownames(am) %in% mod_dt$gene)) am <- t(am)
  dt <- as.data.table(as.table(am))
  setnames(dt, c("cell_type", "gene", "attribution"))
  dt[, module := m]
  dt
}))
setcolorder(gene_attr, c("module", "gene", "cell_type", "attribution"))
setorder(gene_attr, module, gene, -attribution)
fwrite(gene_attr, out_gene_attr, sep = "\t")

# per-cell and per-neighborhood module scores, keyed by cell id for downstream joins to
# labels and sample metadata. Large (n_cells x n_modules), so they go to the
# heavy-intermediates directory as parquet, not results/.
sc  <- as.matrix(isc$scores_sc)
env <- as.matrix(isc$scores_env)
if (is.null(rownames(sc)))  rownames(sc)  <- rownames(norm)
if (is.null(rownames(env))) rownames(env) <- rownames(norm)
stopifnot(nrow(sc) == length(cells_4h), nrow(env) == length(cells_4h))
write_parquet(as.data.table(sc,  keep.rownames = "cell"), out_scores_sc)
write_parquet(as.data.table(env, keep.rownames = "cell"), out_scores_env)

# conditional correlation network: sparse after roundcortozero, serialized as the nonzero
# upper-triangle edge list
cc <- as(isc$condcor, "TsparseMatrix")
ut <- cc@i < cc@j
edges <- data.table(gene_a = rownames(cc)[cc@i[ut] + 1L],
                    gene_b = colnames(cc)[cc@j[ut] + 1L],
                    cond_cor = cc@x[ut])
setorder(edges, -cond_cor)
fwrite(edges, out_condcor, sep = "\t")

cat(sprintf("insitucor: attribution %d module x cell_type rows, %d gene-level rows | scores %d cells x %d modules | condcor %d edges\n",
            nrow(ct_attr), nrow(gene_attr), nrow(sc), ncol(sc), nrow(edges)))
