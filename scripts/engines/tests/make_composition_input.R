#!/usr/bin/env Rscript
# Build the per-cell label input the lm_engine proportion path consumes, for the M02
# composition smoke test. Mirrors composition.R: join obs + unified labels, keep M02,
# label = cell_subtype (unassigned retained, as in the primary propeller test).
# Args: <obs.parquet> <full_labels.parquet> <out.tsv>
suppressPackageStartupMessages({library(data.table); library(arrow)})
a <- commandArgs(trailingOnly = TRUE)
ob  <- as.data.table(read_parquet(a[1]))[, .(cell, sample_id, dataset, condition, slide_id)]
lab <- as.data.table(read_parquet(a[2]))[, .(cell, label = cell_subtype)]
m <- merge(ob[dataset == "Mutter_02"], lab, by = "cell")[!is.na(label)]
fwrite(m[, .(cell, sample_id, label, condition, slide_id)], a[3], sep = "\t")
cat(sprintf("composition input: %d M02 cells, %d subtypes, %d samples\n",
            nrow(m), uniqueN(m$label), uniqueN(m$sample_id)))
