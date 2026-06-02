#!/usr/bin/env Rscript
# Decisive validation of de-anchored per-cell typing on sparse CosMx (median ~53 genes/cell).
# Score sam0007 with UCell (maxRank=100), then GOLD-MARKER RECALL: take cells unambiguously
# expressing a lineage's strong defining marker(s) and report what fraction the argmax typer assigns
# to the matching lineage -- under (a) force-argmax and (b) an unassigned threshold (top1 & margin).
# Chance is ~1/14=7%. If recall is high, tier-1 marker typing works; if ~chance, escalate to tier-2
# (cluster-level annotation). Args: <norm.rds> <markers.yaml>
suppressMessages({library(Seurat); library(Matrix); library(UCell); library(yaml)})
args <- commandArgs(trailingOnly = TRUE)
obj  <- readRDS(args[1]); sets <- yaml::read_yaml(args[2])
sets <- lapply(sets, function(g) intersect(g, rownames(obj)))
sets <- sets[vapply(sets, length, integer(1)) >= 2]
cnt  <- GetAssayData(obj, layer = "counts")
det  <- function(g) { g <- intersect(g, rownames(cnt)); Matrix::colSums(cnt[g,,drop=FALSE] > 0) }

o  <- AddModuleScore_UCell(obj, features = sets, name = "_lin", maxRank = 100)
sc <- as.matrix(o[[paste0(names(sets), "_lin")]]); colnames(sc) <- names(sets)
ti <- max.col(sc, ties.method = "first"); top1 <- sc[cbind(seq_len(nrow(sc)), ti)]
sc2 <- sc; sc2[cbind(seq_len(nrow(sc)), ti)] <- -Inf
top2 <- sc2[cbind(seq_len(nrow(sc)), max.col(sc2, ties.method = "first"))]
lab <- colnames(sc)[ti]; margin <- top1 - top2
lab_epi <- lab; lab_epi[lab_epi == "Epithelial"] <- "tumor_epithelial"

# gold subsets: cells clearly expressing a lineage's defining marker(s)
gold <- list(
  tumor_epithelial = det(c("Epcam","Krt8")) >= 2,   # both Epcam & Krt8
  T                = det("Cd3e") >= 1 | det("Cd3d") >= 1,
  B                = det(c("Cd79a","Ms4a1")) >= 2,
  Endothelial      = det(c("Pecam1","Cdh5")) >= 2,
  Fibroblast       = det(c("Col1a1","Dcn")) >= 2,
  Macrophage       = det(c("Cd68","Csf1r")) >= 1,
  NK               = det("Ncr1") >= 1,
  SmoothMuscle     = det(c("Myh11","Acta2")) >= 2)

report <- function(L, tag) {
  cat(sprintf("\n-- gold-marker recall (%s) : unassigned overall %.1f%% --\n",
              tag, 100*mean(L == "unassigned")))
  for (k in names(gold)) {
    sub <- gold[[k]]; if (sum(sub) == 0) next
    match <- if (k == "tumor_epithelial") L == "tumor_epithelial" else L == k
    cat(sprintf("  %-16s n=%6d  recall=%5.1f%%  (unassigned %4.1f%%)\n",
        k, sum(sub), 100*mean(match[sub]), 100*mean(L[sub] == "unassigned")))
  }
}
report(lab_epi, "force-argmax")
for (q in list(c(0.15,0.03), c(0.20,0.05))) {
  L <- lab_epi; L[top1 < q[1] | margin < q[2]] <- "unassigned"
  report(L, sprintf("min_score>=%.2f & margin>=%.2f", q[1], q[2]))
}
