#!/usr/bin/env Rscript
# Decisive, HONEST validation of de-anchored cluster-then-annotate. The gold-marker recall yardstick
# (Cd3+ -> should be T, etc.) is only meaningful against the bar the status-quo reference itself clears.
# So compute gold-marker recall for BOTH (a) the de-anchored cluster typer and (b) Yi's ImmGen reference
# labels mapped to the same taxonomy. If the de-anchored typer matches/beats Yi on each lineage, it is
# validated as "no worse than the M01-anchored status quo, and de-anchored." Args: <norm.rds> <markers.yaml> <yi.tsv> [res=2.0]
suppressMessages({library(Seurat); library(Matrix); library(UCell); library(yaml); library(readr)})
args <- commandArgs(trailingOnly = TRUE)
RES <- if (length(args) >= 4) as.numeric(args[4]) else 2.0
obj <- readRDS(args[1]); sets <- yaml::read_yaml(args[2])
sets <- lapply(sets, function(g) intersect(g, rownames(obj)))
sets <- sets[vapply(sets, length, integer(1)) >= 2]

# ---- de-anchored cluster-then-annotate (method A) ----
if (length(VariableFeatures(obj)) < 50) obj <- FindVariableFeatures(obj, nfeatures = 500)
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 30, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
obj <- AddModuleScore_UCell(obj, features = sets, name = "_lin", maxRank = 100)
sc  <- as.matrix(obj[[paste0(names(sets), "_lin")]]); colnames(sc) <- names(sets)
obj <- FindClusters(obj, resolution = RES, verbose = FALSE)
cl  <- obj$seurat_clusters
cm  <- t(apply(sc, 2, function(x) tapply(x, cl, mean)))      # lineage x cluster
M   <- t(scale(t(cm)))                                       # z across clusters per lineage
clab <- apply(M, 2, function(v) names(v)[which.max(v)])
method_lab <- clab[as.character(cl)]
method_lab[method_lab == "Epithelial"] <- "tumor_epithelial"

# ---- Yi reference labels mapped to the same taxonomy ----
yi   <- read_tsv(args[3], show_col_types = FALSE)
ylab <- setNames(yi$ImmuneAtlas_ImmGen_Main_cell_Types, yi$cell_id)[colnames(obj)]
map_yi <- function(x) {
  out <- rep(NA_character_, length(x))
  out[grepl("gdT|NKT|CD8.T|CD4|Treg|Thymic|DN[0-9]|^ISP$|preT|Spleen.Naive.CD|Spleen.LN.Naive|CD4Act", x)] <- "T"
  out[grepl("^B.cell|Memory.B|GC_centro|Spleen.CD19", x)] <- "B"
  out[grepl("Plasma", x)] <- "Plasma"
  out[grepl("^NK$|^ILC$", x)] <- "NK"
  out[grepl("Macrophage|macs|Microglia|monocyte", x)] <- "Macrophage"
  out[grepl("Neutrophil", x)] <- "Neutrophil"
  out[grepl("Dendritic", x)] <- "DC"
  out[grepl("endothelial", x)] <- "Endothelial"
  out[grepl("Fibroblastic", x)] <- "Fibroblast"
  out[grepl("^Pericyte$", x)] <- "Pericyte"
  out[grepl("^a$|Stem.Prog", x)] <- "tumor_epithelial"
  out
}
yi_lab <- map_yi(ylab)

# ---- gold subsets from strong, well-detected markers ----
cnt <- GetAssayData(obj, layer = "counts")
det <- function(g){ g <- intersect(g, rownames(cnt)); Matrix::colSums(cnt[g,,drop=FALSE]>0) }
gold <- list(tumor_epithelial = det(c("Epcam","Krt8"))>=2, T = det("Cd3e")>=1 | det("Cd3d")>=1,
             B = det(c("Cd79a","Ms4a1"))>=2, Endothelial = det(c("Pecam1","Cdh5"))>=2,
             Fibroblast = det(c("Col1a1","Dcn"))>=2, Macrophage = det(c("Cd68","Csf1r"))>=1,
             NK = det("Ncr1")>=1)

cat(sprintf("===== sam %s  res=%.1f  %d clusters =====\n",
            as.character(obj$sample_id[1]), RES, nlevels(cl)))
cat(sprintf("%-16s %8s | %18s | %18s\n", "lineage", "gold_n",
            "de-anchored recall", "Yi-ref recall"))
for (k in names(gold)) {
  sub <- gold[[k]]; if (!sum(sub)) next
  mm <- mean(method_lab[sub] == k) * 100
  yy <- mean(yi_lab[sub] == k, na.rm = TRUE) * 100
  cat(sprintf("%-16s %8d | %17.1f%% | %17.1f%%\n", k, sum(sub), mm, yy))
}
cat(sprintf("\nde-anchored fractions: %s\n",
    paste(sprintf("%s %.0f%%", names(sort(table(method_lab),decreasing=TRUE)),
                  100*sort(prop.table(table(method_lab)),decreasing=TRUE)), collapse=", ")))
