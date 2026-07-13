#!/usr/bin/env Rscript
# Sub-state composition adjudication (plan-differential-robustness Fix 3b). Propeller on
# the M02 day-2 labels with Fibroblast split into resting/activated. If the activated-
# fibroblast proportion shifts with SBRT and the collagen/Acta2 fibroblast DE tracks it,
# the fibrosis program is (partly) composition -- more activated cells -- not within-state
# regulation. Mirrors composition.R's propeller block on the relabelled clusters.
# Args: <obs.parquet> <full_labels.parquet> <fibroblast_substate.parquet> <out_test.tsv>
suppressPackageStartupMessages({ library(data.table); library(arrow) })
a   <- commandArgs(trailingOnly = TRUE)
m   <- as.data.table(read_parquet(a[1]))
lab <- as.data.table(read_parquet(a[2]))[, .(cell, cell_type = cell_subtype)]
sub <- as.data.table(read_parquet(a[3]))                 # cell, substate
m[, cell_type := NULL]
m <- merge(m, lab, by = "cell", all.x = TRUE)[!is.na(cell_type)]
m <- merge(m, sub, by = "cell", all.x = TRUE)
m[cell_type == "Fibroblast", cell_type := paste0("Fibroblast_", substate)]

# engine input: per-cell labels (cell, sample_id, label, condition, slide_id), M02 day2, with
# Fibroblast split into resting/activated. The propeller arm test on these sub-state proportions
# now runs in the shared lm_engine (proportion/logit path, robust eBayes).
m2 <- m[dataset == "Mutter_02"]
lm_input <- m2[, .(cell, sample_id, label = cell_type, condition, slide_id)]
fwrite(lm_input, a[4], sep = "\t")
cat(sprintf("substate composition input: %d M02 cells, %d sub-state labels, %d samples\n",
            nrow(lm_input), uniqueN(lm_input$label), uniqueN(lm_input$sample_id)))
