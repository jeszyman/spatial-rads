#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# verify_arm_tables.R
# Step-0 gate: assert the arm-level backing tables reflect the locked 15-subtype
# roster before any figure renders on top of them. Checks that the cell-type
# columns of composition_by_sample / composition_test / results_master are a
# subset of the full_labels roster, that no retired label names leak through,
# and that composition cell counts reconcile to the locked total. Prints PASS or
# a FAIL with the offending values; exits nonzero on failure.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({library(tidyverse); library(arrow)})

roster  <- read_parquet("results/aggregate/full_labels.parquet") %>% pull(cell_subtype) %>% unique()
retired <- c("tumor_epithelial", "Pericyte", "a", "b")
fail <- character(0)

chk_subset <- function(vals, label) {
  extra <- setdiff(unique(vals), roster)
  if (length(extra)) fail <<- c(fail, sprintf("%s has non-roster labels: %s",
                                              label, paste(extra, collapse = ", ")))
}
chk_retired <- function(vals, label) {
  bad <- intersect(unique(vals), retired)
  if (length(bad)) fail <<- c(fail, sprintf("%s carries retired labels: %s",
                                            label, paste(bad, collapse = ", ")))
}

cbs <- read_tsv("results/aggregate/composition_by_sample.tsv", show_col_types = FALSE)
cts <- read_tsv("results/aggregate/composition_test_m02day2.tsv", show_col_types = FALSE)
mas <- read_tsv("results/aggregate/results_master.tsv", show_col_types = FALSE)
de_units <- mas %>% filter(readout_class == "DE") %>% pull(unit)

chk_subset(cbs$cell_type, "composition_by_sample");  chk_retired(cbs$cell_type, "composition_by_sample")
chk_subset(cts$cell_type, "composition_test");       chk_retired(cts$cell_type, "composition_test")
chk_subset(de_units, "results_master DE");           chk_retired(de_units, "results_master DE")

# Composition cells reconcile to the locked total.
locked_n <- nrow(read_parquet("results/aggregate/full_labels.parquet"))
comp_n   <- sum(cbs$n_cells)
if (comp_n != locked_n)
  fail <- c(fail, sprintf("composition cell total %d != locked labels %d", comp_n, locked_n))

if (length(fail)) {
  cat("FAIL:\n"); cat(paste0("  - ", fail, collapse = "\n"), "\n"); quit(status = 1)
}
cat(sprintf("PASS: arm tables on locked roster (%d subtypes, %d cells reconciled)\n",
            length(roster), comp_n))
