#!/usr/bin/env Rscript
# Read-only per-stage typing scorecard (no pipeline outputs written). Answers:
#  (1) coarse InSituType gold-marker recall BEFORE vs AFTER the Epcam/Krt8 tumor overlay;
#  (2) precision of the cells DIRECTLY anchored to a lineage (the supervised hits) vs the
#      de-novo cluster -> nearest-lineage fallback (the bulk);
#  (3) composition pre- vs post-overlay (where the mass went);
#  (4) per-cluster sequencing depth (is the dominant low-signal cluster quiet tumor or junk).
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(data.table) })

AGG <- "/mnt/data/projects/spatial-rads/aggregate"
o   <- readRDS(file.path(AGG, "merged_typed.rds"))
res <- readRDS(file.path(AGG, "merged_typed.insitutype_res.rds"))
ref <- readRDS(file.path(AGG, "ref_profiles.rds"))
cnt <- LayerData(o, assay = "RNA", layer = "counts")            # genes x cells
md  <- as.data.table(o@meta.data, keep.rownames = "cell")

# --- reconstruct de-novo -> nearest-lineage exactly as typing_insitutype.R ---
lin <- ref$lineage; lin_levels <- sort(unique(lin))
prof <- sapply(lin_levels, function(L) rowMeans(ref$mat[, lin == L, drop = FALSE], na.rm = TRUE))
prof <- prof[rownames(cnt), , drop = FALSE]
gmean <- rowMeans(prof, na.rm = TRUE); gmean[is.nan(gmean)] <- 0
na_ix <- which(is.na(prof), arr.ind = TRUE); if (nrow(na_ix)) prof[na_ix] <- gmean[na_ix[, 1]]
prof[!is.finite(prof)] <- 0
ref_types <- colnames(prof)

clust  <- res$clust; names(clust) <- colnames(o)
denovo <- setdiff(unique(clust), ref_types)
updated <- res$profiles
gg     <- intersect(rownames(updated), rownames(prof))
profA  <- prof[gg, , drop = FALSE]
cos1 <- function(u, V) { u <- u / sqrt(sum(u^2)); as.numeric(crossprod(sweep(V, 2, sqrt(colSums(V^2)), "/"), u)) }
denovo_lin <- setNames(rep(NA_character_, length(denovo)), denovo)
for (d in denovo) if (d %in% colnames(updated)) denovo_lin[d] <- colnames(profA)[which.max(cos1(updated[gg, d], profA))]

md[, clust       := clust]
md[, pre_overlay := ifelse(clust %in% ref_types, clust, denovo_lin[clust])]   # general atlas typing, NO tumor overlay
md[, post        := cell_type_atlas]                                          # final (with overlay)
md[, anchored    := clust %in% ref_types]

# --- gold markers (same panel as the validation step) ---
det <- function(g) { g <- intersect(g, rownames(cnt)); if (!length(g)) return(rep(0L, ncol(cnt))); Matrix::colSums(cnt[g, , drop = FALSE] > 0) }
gold <- list(
  T = det(c("Cd3e","Cd3d")) >= 1, B = det(c("Cd79a","Ms4a1")) >= 2, NK = det("Ncr1") >= 1,
  Macrophage = det(c("Cd68","Csf1r")) >= 1, DC = det(c("Itgax","Flt3","Cd209a")) >= 2,
  Neutrophil = det(c("S100a8","S100a9")) >= 2, Endothelial = det(c("Pecam1","Cdh5")) >= 2,
  Fibroblast = det(c("Col1a1","Dcn")) >= 2, Pericyte = det(c("Rgs5","Pdgfrb","Kcnj8")) >= 2,
  tumor_epithelial = det(c("Epcam","Krt8")) >= 2)

cat("=== (1) STAGE 1 coarse typing: gold-marker recall, BEFORE vs AFTER tumor overlay ===\n")
sb <- rbindlist(lapply(names(gold), function(k) { m <- gold[[k]]; if (!sum(m)) return(NULL)
  data.table(lineage = k, gold_n = sum(m),
    recall_pre_overlay  = round(100 * mean(md$pre_overlay[m] == k), 1),
    recall_post_overlay = round(100 * mean(md$post[m] == k), 1)) }))
print(sb)

cat("\n=== (2) STAGE 1 precision of DIRECT anchor hits (of cells anchored to X, %% carrying X gold marker) ===\n")
prec <- rbindlist(lapply(names(gold), function(k) { idx <- md$anchored & md$post == k; if (!sum(idx)) return(NULL)
  data.table(lineage = k, n_anchored = sum(idx), precision_pct = round(100 * mean(gold[[k]][idx]), 1)) }))
print(prec)
cat(sprintf("\nanchored cells: %d (%.1f%%) | de-novo cells: %d (%.1f%%)\n",
            sum(md$anchored), 100*mean(md$anchored), sum(!md$anchored), 100*mean(!md$anchored)))

cat("\n=== (3) composition: pre-overlay (general atlas) vs post-overlay (with tumor overlay) ===\n")
cmp <- merge(md[, .(n_pre = .N), by = .(label = pre_overlay)],
             md[, .(n_post = .N), by = .(label = post)], by = "label", all = TRUE)
for (j in c("n_pre","n_post")) set(cmp, which(is.na(cmp[[j]])), j, 0L)
cmp[, `:=`(pct_pre = round(100*n_pre/nrow(md),1), pct_post = round(100*n_post/nrow(md),1))]
print(cmp[order(-n_post)])

cat("\n=== (4) sequencing depth by InSituType cluster (junk vs real) ===\n")
if (!all(c("nCount_RNA","nFeature_RNA") %in% names(md))) {
  md[, nCount_RNA := Matrix::colSums(cnt)]; md[, nFeature_RNA := Matrix::colSums(cnt > 0)] }
depth <- md[, .(n = .N, med_genes = round(median(nFeature_RNA)), med_counts = round(median(nCount_RNA)),
                frac_lt20_genes = round(mean(nFeature_RNA < 20), 3)), by = clust][order(-n)]
print(depth, nrows = 100)
cat("\ncohort median genes/cell:", round(median(md$nFeature_RNA)),
    "| median counts/cell:", round(median(md$nCount_RNA)), "\n")
