#!/usr/bin/env Rscript
# Task 2 TDD: the lm_engine proportion path must reproduce the CURRENT composition.R
# propeller estimate/p per (cell_type x contrast). Reference = composition_test_m02day2.tsv
# on disk (present pre-refactor). Live source-compare, not a persisted golden baseline.
suppressPackageStartupMessages(library(data.table))
ref <- fread("results/aggregate/composition_test_m02day2.tsv")   # current script output
eng <- fread("results/aggregate/_lm_smoke/composition_lm.tsv")   # produced by the smoke run
m <- merge(ref[, .(feature_id = cell_type, contrast, ref_est = log2FC_logit, ref_p = pvalue)],
           eng[,  .(feature_id, contrast, estimate, p)],
           by = c("feature_id", "contrast"))
stopifnot(nrow(m) > 0)
bad_est <- m[abs(ref_est - estimate) / pmax(abs(ref_est), 1e-9) >= 1e-6]
bad_p   <- m[abs(ref_p - p) >= 1e-9]
if (nrow(bad_est) || nrow(bad_p)) {
  cat(sprintf("FAIL: %d estimate mismatches, %d p mismatches\n", nrow(bad_est), nrow(bad_p)))
  print(head(bad_est)); print(head(bad_p)); quit(status = 1L)
}
cat(sprintf("PASS test_lm_engine: %d (cell_type x contrast) rows match current composition.R\n", nrow(m)))
