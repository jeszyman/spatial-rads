#!/usr/bin/env Rscript
# Build the four remaining lm_engine smoke inputs from the current producer outputs, so the
# engine can be verified to reproduce each producer's inline limma test before that test is
# removed. matrix path = long (sample_id, feature_id, value, condition, slide_id);
# proportion path = per-cell (cell, sample_id, label, condition, slide_id). M02 day-2 only.
# Args: <obs.parquet> <full_labels.parquet> <mixing_per_sample.tsv> <myeloid_scores.tsv>
#       <niche_per_cell.parquet> <fibroblast_substate.parquet> <out_dir>
suppressPackageStartupMessages({library(data.table); library(arrow)})
a <- commandArgs(trailingOnly = TRUE)
ob  <- as.data.table(read_parquet(a[1]))[, .(cell, sample_id, dataset, condition, slide_id, timepoint_h)]
lab <- as.data.table(read_parquet(a[2]))[, .(cell, cell_subtype)]
outd <- a[7]; dir.create(outd, recursive = TRUE, showWarnings = FALSE)

melt_metrics <- function(dt, metrics) {
  d <- dt[dataset == "Mutter_02" & timepoint_h == 48L]
  melt(d[, c("sample_id", "condition", "slide_id", metrics), with = FALSE],
       id.vars = c("sample_id", "condition", "slide_id"),
       variable.name = "feature_id", value.name = "value")
}
# mixing: 4 global metrics
mix <- fread(a[3])
fwrite(melt_metrics(mix, c("mixing_score","mean_immune_frac","immune_per_tumor","enrichment_over_random")),
       file.path(outd, "mixing_input.tsv"), sep = "\t")
# myeloid: 3 macrophage metrics
mye <- fread(a[4])
fwrite(melt_metrics(mye, c("M1_mean","M2_mean","M2_M1_ratio")),
       file.path(outd, "myeloid_input.tsv"), sep = "\t")
# niche: per-cell niche label, M02
ni <- as.data.table(read_parquet(a[5]))[, .(cell, sample_id, niche)]
ni <- merge(ni, ob[dataset == "Mutter_02", .(cell, condition, slide_id)], by = "cell")
fwrite(ni[, .(cell, sample_id, label = niche, condition, slide_id)],
       file.path(outd, "niche_input.tsv"), sep = "\t")
# substate: per-cell cell_subtype with Fibroblast split, M02
sub <- as.data.table(read_parquet(a[6]))                       # cell, substate
m <- merge(ob[dataset == "Mutter_02"], lab, by = "cell")[!is.na(cell_subtype)]
m <- merge(m, sub, by = "cell", all.x = TRUE)
m[cell_subtype == "Fibroblast", cell_subtype := paste0("Fibroblast_", substate)]
fwrite(m[, .(cell, sample_id, label = cell_subtype, condition, slide_id)],
       file.path(outd, "substate_input.tsv"), sep = "\t")
cat("lm inputs built: mixing, myeloid, niche, substate\n")
