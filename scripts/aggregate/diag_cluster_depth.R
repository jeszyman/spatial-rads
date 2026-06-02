#!/usr/bin/env Rscript
# Read-only diagnostic: per-cell sequencing depth (nCount = total transcripts,
# nFeature = genes detected) grouped by InSituType cluster, to decide whether the
# dominant low-signal cluster (i) is genuine quiet tumor or low-quality/junk cells.
# Uses merged_typed.rds metadata if nCount/nFeature already present; else derives
# from counts. No writes to pipeline outputs.
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(data.table) })
o <- readRDS("/mnt/data/projects/spatial-rads/aggregate/merged_typed.rds")
md <- as.data.table(o@meta.data, keep.rownames = "cell")
cat("meta columns:", paste(names(md), collapse = ", "), "\n\n")

if (!all(c("nCount_RNA", "nFeature_RNA") %in% names(md))) {
  cnt <- LayerData(o, assay = "RNA", layer = "counts")
  md[, nCount_RNA   := Matrix::colSums(cnt)]
  md[, nFeature_RNA := Matrix::colSums(cnt > 0)]
}
clcol <- if ("insitutype_clust" %in% names(md)) "insitutype_clust" else "cell_type_atlas"
md[, clust := md[[clcol]]]

summ <- md[, .(n = .N,
               med_nCount   = round(median(nCount_RNA)),
               med_nFeature = round(median(nFeature_RNA)),
               frac_lt50_genes  = round(mean(nFeature_RNA < 50), 3),
               frac_lt20_genes  = round(mean(nFeature_RNA < 20), 3)),
           by = clust][order(-n)]
print(summ, nrows = 100)

# whole-cohort reference percentiles
cat("\ncohort nFeature quantiles:\n")
print(round(quantile(md$nFeature_RNA, c(.01,.05,.25,.5,.75,.95))))
cat(sprintf("\ncluster i vs confident-tumor (a,l,h) vs immune (T,Macrophage) median genes:\n"))
for (g in c("i","a","l","h","T","Macrophage")) {
  v <- md[clust == g, nFeature_RNA]
  if (length(v)) cat(sprintf("  %-12s n=%8d  med_genes=%4d  med_counts=%5d\n",
                             g, length(v), round(median(v)), round(median(md[clust==g, nCount_RNA]))))
}
