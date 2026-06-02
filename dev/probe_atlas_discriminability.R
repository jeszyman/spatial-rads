#!/usr/bin/env Rscript
# Pre-validate the aggregate annotation METHOD before the full run finishes.
# annotate.R labels clusters by Spearman-correlating per-cluster pseudobulk to the 93
# external atlas profiles, aggregating to a per-lineage mean, argmax = label. The open
# risk: on a 950-gene panel the atlas profiles may be too inter-correlated to separate
# lineages, so every group correlates ~equally with everything -> garbage labels.
# Cheap test: group ONE real sample by its known cell_type (proxy for clusters) and run
# the identical correlation rule. If the reliable structural types (tumor 'a', endothelial,
# fibroblast) land on the right lineage with a clear margin, the method discriminates.
# Args: <scored.rds> <ref_profiles.rds>
suppressPackageStartupMessages({library(Seurat); library(Matrix)})
args <- commandArgs(trailingOnly = TRUE)
obj <- readRDS(args[1]); ref <- readRDS(args[2])
DefaultAssay(obj) <- "RNA"; obj <- tryCatch(JoinLayers(obj), error = function(e) obj)

ct  <- as.character(obj$cell_type)
keep <- names(which(table(ct) >= 50))            # only types with enough cells to pseudobulk
ct[!ct %in% keep] <- NA
grp <- factor(ct[!is.na(ct)]); cells <- which(!is.na(ct))
dat <- LayerData(obj, layer = "data")[, cells, drop = FALSE]

# pseudobulk (mean log-norm) per cell_type, identical matmul to annotate.R
gi  <- as.integer(grp)
Ind <- sparseMatrix(i = seq_along(gi), j = gi, x = 1, dims = c(length(gi), nlevels(grp)))
pb  <- as.matrix(dat %*% Ind); pb <- sweep(pb, 2, Matrix::colSums(Ind), "/")
colnames(pb) <- levels(grp)

R   <- ref$mat[rownames(pb), , drop = FALSE]; lin <- ref$lineage
cor_cp <- matrix(NA_real_, ncol(pb), ncol(R), dimnames = list(colnames(pb), colnames(R)))
for (j in seq_len(ncol(R))) { ok <- !is.na(R[, j]); cor_cp[, j] <- cor(pb[ok, , drop = FALSE], R[ok, j], method = "spearman") }
lin_levels <- sort(unique(lin))
lin_score  <- sapply(lin_levels, function(L) rowMeans(cor_cp[, lin == L, drop = FALSE]))
rownames(lin_score) <- colnames(pb)

# expected lineage for the RELIABLE structural Yi types (memory: 'a'/Stem.Prog = tumor bucket;
# endothelial/fibroblast labels trustworthy; immune Yi labels unreliable -> reported, not scored)
expect <- c(a = "Epithelial", Stem.Prog = "Epithelial",
            Blood.endothelial = "Endothelial", Lymphatic.endothelial = "Endothelial",
            Fibroblastic.reticular = "Fibroblast")

cat(sprintf("sample %s: %d cells, %d cell_types pseudobulked\n",
            as.character(obj$sample_id[1]), length(cells), nlevels(grp)))
cat(sprintf("\n%-26s %6s | %-12s %6s  %-12s %6s  %6s  %-12s %s\n",
            "cell_type", "n", "top1_lin", "score", "top2_lin", "score", "margin", "expected", "HIT?"))
ord <- order(lin_score[, "Epithelial"], decreasing = TRUE)   # tumor-ish first, just for reading
for (k in rownames(lin_score)[ord]) {
  v <- lin_score[k, ]; o2 <- order(v, decreasing = TRUE)
  exp_l <- expect[k]; hit <- if (is.na(exp_l)) "(unscored)" else if (lin_levels[o2[1]] == exp_l) "HIT" else "MISS"
  cat(sprintf("%-26s %6d | %-12s %6.3f  %-12s %6.3f  %6.3f  %-12s %s\n",
      k, sum(grp == k), lin_levels[o2[1]], v[o2[1]], lin_levels[o2[2]], v[o2[2]],
      v[o2[1]] - v[o2[2]], ifelse(is.na(exp_l), "-", exp_l), hit))
}
# discriminability summary: how diagonal-dominant is the cell_type x lineage score matrix?
cat("\nfull cell_type x lineage Spearman (rows=Yi type, cols=atlas lineage):\n")
print(round(lin_score, 3))
