#!/usr/bin/env Rscript
# aggregate.smk Stage 1b -- de-anchored cluster annotation against an EXTERNAL atlas.
# Reads the integrated/clustered object (embed_celltype.R) and the external reference
# profile matrix (prepare_reference.R), then labels each joint Louvain cluster by
# correlation to the atlas -- NO M01 -> M02 transfer, so the M01 and M02 labels are
# defined by one identical rule and are structurally comparable (the whole point;
# see plan-processing-pipeline.md Tier-2 + memory project_yi_labels_unreliable_immune).
#
# Method (manual SingleR; celldex/SingleR absent from the env):
#   1. per-cluster pseudobulk = mean log-normalized `data` over the cluster's cells
#      (sparse data %*% cluster-indicator, then /n_k).
#   2. Spearman-correlate each cluster pseudobulk vs each of the 93 reference profiles
#      over that profile's non-NA genes. Spearman (rank) because the reference profiles
#      are linear-scale atlas means while pseudobulk is log-normalized -- only rank is
#      comparable.
#   3. per-lineage score = MEAN of its profiles' correlations (NOT max: max would favor
#      lineages with many profiles -- T has 22, Pericyte 1 -- by giving them more draws
#      at a spurious high; the mean is count-unbiased). argmax lineage = cluster label.
#   4. 4T1 overlay: the tumor is mammary carcinoma, so an atlas "Epithelial" cluster IS
#      tumor -- relabel Epithelial -> tumor_epithelial when the cluster actually expresses
#      Epcam/Krt8 (guards against a non-tumor epithelial-correlating cluster).
#   5. propagate the cluster label to every cell as `cell_type_atlas` (old per-sample
#      cell_type kept for cross-reference only).
#
# Validation: gold-marker recall per lineage (fraction of strong-marker+ cells that got
# the matching label), split M01/M02. This is the honest yardstick from the pilot work.
#
# Args: <merged_celltype.rds> <ref_profiles.rds> <out_typed.rds> <out_clusters.tsv>
#       <out_summary.tsv> <out_validation.tsv> <out_cell_labels.tsv>
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
})
set.seed(1)

args        <- commandArgs(trailingOnly = TRUE)
in_rds      <- args[1]
ref_rds     <- args[2]
out_typed   <- args[3]
out_clusters<- args[4]
out_summary <- args[5]
out_valid   <- args[6]
out_labels  <- args[7]

EPI_GATE <- 0.25   # min cluster fraction Epcam+/Krt8+ to relabel Epithelial -> tumor

for (d in unique(c(dirname(out_typed), dirname(out_clusters), dirname(out_labels))))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- load object + reference -------------------------------------------------
o   <- readRDS(in_rds)
DefaultAssay(o) <- "RNA"
o   <- tryCatch(JoinLayers(o), error = function(e) o)
ref <- readRDS(ref_rds)
stopifnot(!is.null(o$seurat_clusters))
cl   <- factor(o$seurat_clusters)
Klev <- levels(cl)
cat(sprintf("loaded %d cells, %d clusters; reference %d profiles x %d genes, %d lineages\n",
            ncol(o), nlevels(cl), ncol(ref$mat), nrow(ref$mat), length(unique(ref$lineage))))

# ---- per-cluster pseudobulk (mean of log-normalized data) --------------------
dat <- LayerData(o, layer = "data")                 # 950 x N sparse
ci  <- as.integer(cl)
Ind <- sparseMatrix(i = seq_along(ci), j = ci, x = 1,
                    dims = c(length(ci), nlevels(cl)))
nk  <- Matrix::colSums(Ind)
pb  <- as.matrix(dat %*% Ind)                        # 950 x K (cluster sums)
pb  <- sweep(pb, 2, nk, "/")                         # -> means
colnames(pb) <- Klev

# ---- cluster x profile Spearman, then lineage means --------------------------
R   <- ref$mat[rownames(pb), , drop = FALSE]         # align genes (both = 950 panel)
lin <- ref$lineage
cor_cp <- matrix(NA_real_, ncol(pb), ncol(R),
                 dimnames = list(Klev, colnames(R)))
for (j in seq_len(ncol(R))) {
  ok <- !is.na(R[, j])
  cor_cp[, j] <- cor(pb[ok, , drop = FALSE], R[ok, j], method = "spearman")
}
lin_levels <- sort(unique(lin))
lin_score  <- sapply(lin_levels, function(L)
  rowMeans(cor_cp[, lin == L, drop = FALSE]))        # K x nLineage
rownames(lin_score) <- Klev

ord_mat    <- t(apply(lin_score, 1, order, decreasing = TRUE))   # K x nLineage indices
rows       <- seq_len(nrow(lin_score))
top1_idx   <- ord_mat[, 1]; top2_idx <- ord_mat[, 2]
top1_lin   <- lin_levels[top1_idx]
top2_lin   <- lin_levels[top2_idx]
top1_score <- lin_score[cbind(rows, top1_idx)]
top2_score <- lin_score[cbind(rows, top2_idx)]
margin     <- top1_score - top2_score
lab        <- top1_lin
best_prof  <- colnames(cor_cp)[apply(cor_cp, 1, which.max)]
best_cor   <- apply(cor_cp, 1, max)

# ---- 4T1 overlay: atlas Epithelial that expresses Epcam/Krt8 -> tumor ---------
cnt <- LayerData(o, layer = "counts")
epi_genes <- intersect(c("Epcam", "Krt8"), rownames(cnt))
epi_det   <- if (length(epi_genes)) {
  Matrix::colSums(cnt[epi_genes, , drop = FALSE] > 0) > 0
} else rep(FALSE, ncol(o))
epi_frac  <- tapply(epi_det, cl, mean)[Klev]
clab      <- lab
is_epi    <- clab == "Epithelial" & !is.na(epi_frac) & epi_frac >= EPI_GATE
clab[is_epi] <- "tumor_epithelial"
names(clab) <- Klev

# ---- propagate to cells ------------------------------------------------------
o$cell_type_atlas <- unname(clab[as.character(cl)])
atl <- o$cell_type_atlas
ds  <- o$dataset

# ---- per-cluster diagnostics -------------------------------------------------
old_modal <- function(k) {
  t <- sort(table(o$cell_type[cl == k]), decreasing = TRUE)
  c(name = names(t)[1], frac = round(as.numeric(t[1]) / sum(t), 3)) }
om <- t(sapply(Klev, old_modal))
clusters <- data.table(
  cluster        = Klev,
  n_cells        = as.integer(nk),
  n_M01          = as.integer(tapply(ds == "Mutter_01", cl, sum)[Klev]),
  n_M02          = as.integer(tapply(ds == "Mutter_02", cl, sum)[Klev]),
  atlas_label    = clab,
  top1_lineage   = top1_lin,
  top1_score     = round(top1_score, 3),
  top2_lineage   = top2_lin,
  top2_score     = round(top2_score, 3),
  margin         = round(margin, 3),
  best_profile   = best_prof,
  best_profile_cor = round(best_cor, 3),
  epcam_krt8_frac  = round(as.numeric(epi_frac), 3),
  modal_old_celltype = om[, "name"],
  modal_old_frac     = as.numeric(om[, "frac"]))
setorder(clusters, -n_cells)
fwrite(clusters, out_clusters, sep = "\t")

# ---- composition summary: cell_type_atlas x dataset (comparability check) -----
tab <- table(atl, ds)
m01 <- if ("Mutter_01" %in% colnames(tab)) tab[, "Mutter_01"] else 0
m02 <- if ("Mutter_02" %in% colnames(tab)) tab[, "Mutter_02"] else 0
summary <- data.table(
  cell_type   = rownames(tab),
  n_M01       = as.integer(m01),
  n_M02       = as.integer(m02),
  frac_M01    = round(as.numeric(m01) / sum(m01), 4),
  frac_M02    = round(as.numeric(m02) / sum(m02), 4),
  frac_overall= round(rowSums(tab) / sum(tab), 4))
setorder(summary, -frac_overall)
fwrite(summary, out_summary, sep = "\t")

# ---- gold-marker recall validation -------------------------------------------
det <- function(g) {
  g <- intersect(g, rownames(cnt))
  if (!length(g)) return(rep(0L, ncol(cnt)))
  Matrix::colSums(cnt[g, , drop = FALSE] > 0) }
gold <- list(
  tumor_epithelial = det(c("Epcam", "Krt8"))     >= 2,
  T                = det(c("Cd3e", "Cd3d"))       >= 1,
  B                = det(c("Cd79a", "Ms4a1"))     >= 2,
  Plasma           = det(c("Jchain", "Sdc1", "Xbp1")) >= 2,
  NK               = det("Ncr1")                  >= 1,
  Macrophage       = det(c("Cd68", "Csf1r"))      >= 1,
  DC               = det(c("Itgax", "Flt3", "Cd209a")) >= 2,
  Neutrophil       = det(c("S100a8", "S100a9"))   >= 2,
  Endothelial      = det(c("Pecam1", "Cdh5"))     >= 2,
  Fibroblast       = det(c("Col1a1", "Dcn"))      >= 2,
  SmoothMuscle     = det(c("Myh11", "Acta2"))     >= 2,
  Pericyte         = det(c("Rgs5", "Pdgfrb", "Kcnj8")) >= 2)
m01c <- ds == "Mutter_01"
valid <- rbindlist(lapply(names(gold), function(k) {
  sub <- gold[[k]]
  if (!sum(sub, na.rm = TRUE)) return(NULL)
  data.table(
    lineage        = k,
    gold_n         = sum(sub),
    gold_n_M01     = sum(sub & m01c),
    gold_n_M02     = sum(sub & !m01c),
    recall_pct     = round(100 * mean(atl[sub] == k), 1),
    recall_M01_pct = round(100 * mean(atl[sub & m01c] == k), 1),
    recall_M02_pct = round(100 * mean(atl[sub & !m01c] == k), 1)) }))
fwrite(valid, out_valid, sep = "\t")

# ---- slim per-cell labels (joinable downstream without loading 3 GB object) --
fwrite(data.table(cell_id = colnames(o), dataset = ds,
                  cluster = as.character(cl), cell_type_atlas = atl),
       out_labels, sep = "\t")

saveRDS(o, out_typed)

cat("=== atlas cluster labels ===\n"); print(clusters[, .(cluster, n_cells, atlas_label, top1_score, margin, modal_old_celltype)])
cat("\n=== composition (M01 vs M02 fractions) ===\n"); print(summary)
cat("\n=== gold-marker recall ===\n"); print(valid)
cat(sprintf("\nannotate: %d cells labeled across %d atlas types; typed.rds written.\n",
            ncol(o), uniqueN(atl)))
