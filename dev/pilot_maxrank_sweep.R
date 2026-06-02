#!/usr/bin/env Rscript
# Pick UCell maxRank for de-anchored typing on the sparse CosMx panel (median ~53 genes/cell on 950).
# Default maxRank=1500 saturates -> near-random argmax (margin 0.029, held-out markers at background).
# Score sam0007 at several maxRanks; report median argmax margin + held-out-marker discrimination
# (in-assigned-lineage vs rest, using genes NOT in lineage_markers.yaml). Args: <norm.rds> <markers.yaml>
suppressMessages({library(Seurat); library(Matrix); library(UCell); library(yaml)})
args <- commandArgs(trailingOnly = TRUE)
obj  <- readRDS(args[1]); sets <- yaml::read_yaml(args[2])
sets <- lapply(sets, function(g) intersect(g, rownames(obj)))
sets <- sets[vapply(sets, length, integer(1)) >= 2]
cnt  <- GetAssayData(obj, layer = "counts")
held <- list(tumor_epithelial = c("Muc1","Sfn"), Macrophage = c("Tyrobp","Apoe"),
             Fibroblast = c("Col6a1","Thbs2"), Endothelial = c("Adgrl4","Ackr1"))
held <- lapply(held, function(g) intersect(g, rownames(cnt)))

for (mr in c(50, 100, 150, 300)) {
  o  <- AddModuleScore_UCell(obj, features = sets, name = "_lin", maxRank = mr)
  sc <- as.matrix(o[[paste0(names(sets), "_lin")]]); colnames(sc) <- names(sets)
  ti <- max.col(sc, ties.method = "first"); top1 <- sc[cbind(seq_len(nrow(sc)), ti)]
  sc2 <- sc; sc2[cbind(seq_len(nrow(sc)), ti)] <- -Inf
  top2 <- sc2[cbind(seq_len(nrow(sc)), max.col(sc2, ties.method = "first"))]
  lab <- colnames(sc)[ti]; lab[lab == "Epithelial"] <- "tumor_epithelial"
  cat(sprintf("\n=== maxRank=%d : median_score=%.3f median_margin=%.3f ===\n",
              mr, median(top1), median(top1 - top2)))
  for (lin in names(held)) {
    g <- held[[lin]]; if (!length(g)) next
    det <- Matrix::colSums(cnt[g, , drop = FALSE] > 0) > 0; inl <- lab == lin
    cat(sprintf("  %-16s [%s]  frac=%4.1f%%  held-out in %5.1f%% vs rest %5.1f%%  (ratio %.2f)\n",
        lin, paste(g, collapse="/"), 100*mean(inl),
        100*mean(det[inl]), 100*mean(det[!inl]), mean(det[inl])/max(mean(det[!inl]),1e-6)))
  }
}
