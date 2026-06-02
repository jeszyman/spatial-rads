#!/usr/bin/env Rscript
# Crux question: is the de-anchored T=0% a fixable LABELING bug (a real T cluster got mislabeled) or a
# DATA CEILING (T cells smeared across clusters, never isolated)? Per cluster show: n, Cd3+%, corrected
# Yi-T%, raw-mean T UCell score + its RANK among lineages, and the assigned label. If some cluster is
# clearly T-dominated (high Cd3+/Yi-T, T score rank 1-2) but mislabeled -> fixable. If T% is flat ~15%
# everywhere with no cluster dominated -> ceiling. Args: <norm.rds> <markers.yaml> <yi.tsv> [res=2.0]
suppressMessages({library(Seurat); library(Matrix); library(UCell); library(yaml); library(readr)})
args <- commandArgs(trailingOnly = TRUE)
RES <- if (length(args) >= 4) as.numeric(args[4]) else 2.0
obj <- readRDS(args[1]); sets <- yaml::read_yaml(args[2])
sets <- lapply(sets, function(g) intersect(g, rownames(obj)))
sets <- sets[vapply(sets, length, integer(1)) >= 2]
if (length(VariableFeatures(obj)) < 50) obj <- FindVariableFeatures(obj, nfeatures = 500)
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 30, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
obj <- AddModuleScore_UCell(obj, features = sets, name = "_lin", maxRank = 100)
sc  <- as.matrix(obj[[paste0(names(sets), "_lin")]]); colnames(sc) <- names(sets)
obj <- FindClusters(obj, resolution = RES, verbose = FALSE)
cl  <- obj$seurat_clusters
cm  <- t(apply(sc, 2, function(x) tapply(x, cl, mean)))   # lineage x cluster
Mz  <- t(scale(t(cm)))
zlab <- apply(Mz, 2, function(v) names(v)[which.max(v)])   # current (z) label
rlab <- apply(cm, 2, function(v) names(v)[which.max(v)])   # raw-mean label
Trank <- apply(cm, 2, function(v) rank(-v)["T"])           # rank of T among lineages per cluster

cnt <- GetAssayData(obj, layer = "counts")
cd3 <- (cnt["Cd3e",]>0) | (cnt["Cd3d",]>0)
yi  <- read_tsv(args[3], show_col_types = FALSE)
ylab <- setNames(yi$ImmuneAtlas_ImmGen_Main_cell_Types, yi$cell_id)[colnames(obj)]
yiT <- grepl("gdT|NKT|CD8.T|Treg|Thymic|DN[0-9]|^ISP$|preT|Spleen.Naive.CD|Spleen.LN.Naive|CD4Act", ylab)

cat(sprintf("===== sam %s  res=%.1f  %d clusters =====\n", as.character(obj$sample_id[1]), RES, nlevels(cl)))
cat(sprintf("%-5s %6s %7s %7s %7s %6s  %-14s %-14s\n",
            "clust","n","Cd3+%","YiT%","Tscore","Trank","raw_label","z_label"))
ord <- order(tapply(yiT, cl, mean), decreasing = TRUE)
for (k in levels(cl)[ord]) {
  cat(sprintf("cl%-3s %6d %6.1f%% %6.1f%% %7.3f %6.0f  %-14s %-14s\n",
      k, sum(cl==k), 100*mean(cd3[cl==k]), 100*mean(yiT[cl==k], na.rm=TRUE),
      cm["T", k], Trank[k], rlab[k], zlab[k]))
}
cat(sprintf("\noverall: Cd3+ %.1f%%, Yi-T %.1f%%\n", 100*mean(cd3), 100*mean(yiT, na.rm=TRUE)))
