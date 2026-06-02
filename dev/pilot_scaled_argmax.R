#!/usr/bin/env Rscript
# Does per-lineage normalization rescue immune recall? Naive raw-argmax fails on T (11.7%) / Mac (23.8%)
# because their defining markers are sparse vs broadly-expressed tumor/fibroblast/ambient genes. Test
# argmax on per-lineage-normalized UCell scores: (b) column z-score, (c) column percentile-rank.
# Report gold-marker recall (chance ~7%). Args: <norm.rds> <markers.yaml>
suppressMessages({library(Seurat); library(Matrix); library(UCell); library(yaml)})
args <- commandArgs(trailingOnly = TRUE)
obj  <- readRDS(args[1]); sets <- yaml::read_yaml(args[2])
sets <- lapply(sets, function(g) intersect(g, rownames(obj)))
sets <- sets[vapply(sets, length, integer(1)) >= 2]
cnt  <- GetAssayData(obj, layer = "counts")
det  <- function(g){ g <- intersect(g, rownames(cnt)); Matrix::colSums(cnt[g,,drop=FALSE] > 0) }

o  <- AddModuleScore_UCell(obj, features = sets, name = "_lin", maxRank = 100)
sc <- as.matrix(o[[paste0(names(sets), "_lin")]]); colnames(sc) <- names(sets)

gold <- list(tumor_epithelial = det(c("Epcam","Krt8")) >= 2, T = det("Cd3e")>=1 | det("Cd3d")>=1,
             B = det(c("Cd79a","Ms4a1"))>=2, Endothelial = det(c("Pecam1","Cdh5"))>=2,
             Fibroblast = det(c("Col1a1","Dcn"))>=2, Macrophage = det(c("Cd68","Csf1r"))>=1,
             NK = det("Ncr1")>=1, SmoothMuscle = det(c("Myh11","Acta2"))>=2)
assign_lab <- function(M){ ti <- max.col(M, ties.method="first"); l <- colnames(M)[ti]
                           l[l=="Epithelial"] <- "tumor_epithelial"; l }
report <- function(L, tag){ cat(sprintf("\n-- %s --  (marginal: %s)\n", tag,
    paste(sprintf("%s %.0f%%", names(sort(table(L),decreasing=TRUE))[1:4],
                  100*sort(prop.table(table(L)),decreasing=TRUE)[1:4]), collapse=", ")))
  for (k in names(gold)){ sub <- gold[[k]]; if(!sum(sub)) next
    m <- if(k=="tumor_epithelial") L=="tumor_epithelial" else L==k
    cat(sprintf("  %-16s n=%6d  recall=%5.1f%%\n", k, sum(sub), 100*mean(m[sub]))) } }

report(assign_lab(sc), "RAW argmax")
z <- scale(sc); report(assign_lab(z), "Z-SCORE per lineage, argmax")
r <- apply(sc, 2, function(x) rank(x)/length(x)); report(assign_lab(r), "PERCENTILE per lineage, argmax")
