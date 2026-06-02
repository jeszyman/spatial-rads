#!/usr/bin/env Rscript
# Pilot diagnostic for de-anchored typing (plan-processing-pipeline.md v2.0). Reads a typed.rds from
# the new celltype.R and reports: lineage fractions, unassigned %, score/margin distributions,
# HELD-OUT marker concordance (markers NOT in lineage_markers.yaml -> independent check), and (M01
# only) a cross-tab of new labels vs Yi's ImmGen labels (does de-anchoring resolve the 48.5% `a` blob?).
# Args: <typed.rds> <lineage_markers.yaml> [yi_reference.tsv]
suppressMessages({library(Seurat); library(Matrix); library(yaml); library(readr)})
args <- commandArgs(trailingOnly = TRUE)
TYPED <- args[1]; MARKERS <- args[2]; YI <- if (length(args) >= 3) args[3] else NA

obj <- readRDS(TYPED)
ct  <- obj$cell_type; n <- length(ct)
sid <- as.character(obj$sample_id[1]); ds <- as.character(obj$dataset[1])
cat(sprintf("\n===== %s (%s) : %d cells =====\n", sid, ds, n))

## (1) lineage fractions
tab <- sort(table(ct), decreasing = TRUE)
cat("\n-- lineage fractions --\n")
for (k in names(tab)) cat(sprintf("  %-16s %7d  %5.1f%%\n", k, tab[k], 100 * tab[k] / n))

## (2) score / margin distribution
cat(sprintf("\n-- scores --  median_score=%.3f  median_margin=%.3f  frac_score<0.10=%.1f%%\n",
            median(obj$celltype_score), median(obj$celltype_margin),
            100 * mean(obj$celltype_score < 0.10)))

## (3) HELD-OUT marker concordance (genes deliberately absent from lineage_markers.yaml)
cnt <- GetAssayData(obj, layer = "counts")
held <- list(tumor_epithelial = intersect(c("Muc1", "Sfn"), rownames(cnt)),
             Endothelial      = intersect(c("Adgrl4", "Ackr1"), rownames(cnt)),
             Fibroblast       = intersect(c("Col6a1", "Thbs2"), rownames(cnt)),
             T                = intersect(c("Ms4a4b"), rownames(cnt)),
             Macrophage       = intersect(c("Tyrobp", "Apoe"), rownames(cnt)))
cat("\n-- held-out marker detection (% cells with >0 count) : assigned-lineage vs rest --\n")
for (lin in names(held)) {
  g <- held[[lin]]; if (!length(g)) next
  det <- Matrix::colSums(cnt[g, , drop = FALSE] > 0) > 0
  inlin <- ct == lin
  cat(sprintf("  %-16s [%s]  in-lineage %5.1f%%   rest %5.1f%%\n",
              lin, paste(g, collapse = "/"), 100 * mean(det[inlin]), 100 * mean(det[!inlin])))
}

## (4) M01 only: new label vs Yi ImmGen label (does `a` blob resolve?)
if (!is.na(YI) && ds == "Mutter_01" && file.exists(YI)) {
  yi  <- read_tsv(YI, show_col_types = FALSE)
  ylab <- setNames(yi$ImmuneAtlas_ImmGen_Main_cell_Types, yi$cell_id)[colnames(obj)]
  cat("\n-- new label composition WITHIN Yi's `a` bucket (the 48.5% blob) --\n")
  a_new <- sort(table(ct[which(ylab == "a")]), decreasing = TRUE)
  for (k in names(a_new)) cat(sprintf("  a -> %-16s %5.1f%%\n", k, 100 * a_new[k] / sum(a_new)))
}
cat("\n")
