#!/usr/bin/env Rscript
# Root-cause check: do T cells form a separable cluster, or are they smeared by ambient signal?
# Cluster sam0007 (res=2.0), then per cluster show % detection of canonical single markers + Yi's own
# T-label fraction. If some cluster is clearly Cd3e-enriched -> labeling is the problem (fixable);
# if Cd3e is flat across clusters -> the data can't separate T (ambient dominance). Args: <norm.rds> <yi.tsv>
suppressMessages({library(Seurat); library(Matrix); library(readr)})
args <- commandArgs(trailingOnly = TRUE)
obj <- readRDS(args[1]); YI <- args[2]
if (length(VariableFeatures(obj)) < 50) obj <- FindVariableFeatures(obj, nfeatures = 500)
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 30, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
obj <- FindClusters(obj, resolution = 2.0, verbose = FALSE)
cl  <- obj$seurat_clusters
cnt <- GetAssayData(obj, layer = "counts")

mk <- c(Epcam="Epcam", Krt8="Krt8", Cd3e="Cd3e", Cd3d="Cd3d", Cd68="Cd68", C1qa="C1qa",
        Col1a1="Col1a1", Pecam1="Pecam1", Cdh5="Cdh5", Ms4a1="Ms4a1", Ncr1="Ncr1", Myh11="Myh11", Mzb1="Mzb1")
mk <- mk[mk %in% rownames(cnt)]
pct <- sapply(mk, function(g) tapply(cnt[g, ] > 0, cl, mean) * 100)   # clusters x markers, % detected

yi   <- read_tsv(YI, show_col_types = FALSE)
ylab <- setNames(yi$ImmuneAtlas_ImmGen_Main_cell_Types, yi$cell_id)[colnames(obj)]
yiT  <- grepl("^T\\.|^Tgd|NKT", ylab)
yiTpc <- tapply(yiT, cl, function(x) mean(x, na.rm = TRUE) * 100)

cat(sprintf("sam0007  %d clusters  | overall: Cd3e+ %.1f%%, Yi-T %.1f%%\n",
            nlevels(cl), 100*mean(cnt["Cd3e",]>0), 100*mean(yiT, na.rm=TRUE)))
cat(sprintf("%-5s %6s | %s | %6s\n", "clust", "n", paste(sprintf("%6s", names(mk)), collapse=" "), "Yi-T%"))
for (k in levels(cl)) {
  cat(sprintf("cl%-3s %6d | %s | %5.0f%%\n", k, sum(cl==k),
      paste(sprintf("%5.0f%%", pct[k, ]), collapse=" "), yiTpc[k]))
}
cat(sprintf("\nGold Cd3e+ cells: Yi-T concordance = %.1f%%\n",
            100*mean(yiT[cnt["Cd3e",]>0], na.rm=TRUE)))
