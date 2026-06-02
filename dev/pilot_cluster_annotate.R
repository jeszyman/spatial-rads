#!/usr/bin/env Rscript
# Tier-1.5 de-anchored typing: CLUSTER-then-annotate (good-practice workhorse; pools signal across
# cells to beat the 53-genes/cell sparsity that sank per-cell argmax). Per sample: ScaleData -> PCA ->
# graph clusters (Louvain), then label each CLUSTER by its mean per-cell UCell lineage score,
# Z-SCORED across clusters per lineage (so a T cluster is called T even if broadly-expressed
# macrophage genes give a higher raw mean), propagate to cells. De-anchored, validated by the SAME
# gold-marker recall (chance ~7%). Sweeps resolutions to find where immune lineages separate.
# Args: <norm.rds> <markers.yaml> [maxRank=100]
suppressMessages({library(Seurat); library(Matrix); library(UCell); library(yaml)})
args <- commandArgs(trailingOnly = TRUE)
MR  <- if (length(args) >= 3) as.numeric(args[3]) else 100
obj <- readRDS(args[1]); sets <- yaml::read_yaml(args[2])
sets <- lapply(sets, function(g) intersect(g, rownames(obj)))
sets <- sets[vapply(sets, length, integer(1)) >= 2]
sid <- as.character(obj$sample_id[1]); ds <- as.character(obj$dataset[1])

if (length(VariableFeatures(obj)) < 50) obj <- FindVariableFeatures(obj, nfeatures = 500)
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 30, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
obj <- AddModuleScore_UCell(obj, features = sets, name = "_lin", maxRank = MR)
sc  <- as.matrix(obj[[paste0(names(sets), "_lin")]]); colnames(sc) <- names(sets)

cnt <- GetAssayData(obj, layer = "counts"); det <- function(g){ g <- intersect(g, rownames(cnt)); Matrix::colSums(cnt[g,,drop=FALSE]>0) }
gold <- list(tumor_epithelial = det(c("Epcam","Krt8"))>=2, T = det("Cd3e")>=1 | det("Cd3d")>=1,
             B = det(c("Cd79a","Ms4a1"))>=2, Endothelial = det(c("Pecam1","Cdh5"))>=2,
             Fibroblast = det(c("Col1a1","Dcn"))>=2, Macrophage = det(c("Cd68","Csf1r"))>=1,
             NK = det("Ncr1")>=1, SmoothMuscle = det(c("Myh11","Acta2"))>=2)

label_clusters <- function(cl, zscore) {
  cm <- t(apply(sc, 2, function(x) tapply(x, cl, mean)))     # lineage x cluster
  M  <- if (zscore) t(scale(t(cm))) else cm                  # z across clusters per lineage
  lab <- apply(M, 2, function(v){ l <- names(v)[which.max(v)]; if (l=="Epithelial") "tumor_epithelial" else l })
  unname(lab[as.character(cl)])
}
report <- function(cluster_lab, tag) {
  fr <- sort(table(cluster_lab), decreasing = TRUE)
  cat(sprintf("   fractions: %s\n", paste(sprintf("%s %.0f%%", names(fr), 100*fr/length(cluster_lab)), collapse=", ")))
  cat(sprintf("   recall %s:\n", tag))
  for (k in names(gold)){ sub <- gold[[k]]; if(!sum(sub)) next
    m <- if(k=="tumor_epithelial") cluster_lab=="tumor_epithelial" else cluster_lab==k
    cat(sprintf("     %-16s n=%6d  recall=%5.1f%%\n", k, sum(sub), 100*mean(m[sub]))) }
}

for (RES in c(1.0, 2.0, 3.0)) {
  obj <- FindClusters(obj, resolution = RES, verbose = FALSE)
  cl  <- obj$seurat_clusters
  cat(sprintf("\n===== %s (%s) res=%.1f : %d clusters =====\n", sid, ds, RES, nlevels(cl)))
  report(label_clusters(cl, zscore = TRUE),  "[z-scored cluster argmax]")
}
