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
retired <- c("tumor_epithelial", "Pericyte", "a", "b", "Epithelial cells")
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

# Inference rows must carry the per-arm replicate count their comparison actually has,
# not total-samples. Each cohort has its own arm size (day-2 n=2/arm, combined_4h n=2,
# combined_4h_treated n=3, pooled variants n=2/n=3), so the expected value is read per
# (comparison, contrast) from the registry rather than asserted as one study-wide number.
reg <- read_tsv("results/data_model/comparisons.tsv", show_col_types = FALSE) %>%
  filter(kind == "sample", !is.na(contrast_num_level)) %>%
  transmute(comparison = cohort, contrast = name,
            n_expected = pmin(n_group1, n_group2)) %>%
  distinct()
if (any(duplicated(reg[c("comparison","contrast")])))
  fail <- c(fail, "registry has conflicting n per (comparison, contrast)")

missing <- setdiff(c("comparison", "contrast", "n_per_arm", "readout_class"), names(mas))
if (length(missing)) {
  fail <- c(fail, sprintf("results_master lacks %s -- regenerate it with assemble_results.R",
                          paste(missing, collapse = ", ")))
} else {
  bad_n <- mas %>%
    filter(readout_class %in% c("DE","composition","pathway")) %>%
    left_join(reg, by = c("comparison", "contrast")) %>%
    filter(is.na(n_expected) | is.na(n_per_arm) | n_per_arm != n_expected)
  if (nrow(bad_n)) {
    detail <- bad_n %>% count(comparison, contrast, n_per_arm, n_expected) %>%
      mutate(txt = sprintf("%s/%s: master n_per_arm=%s registry=%s (%d rows)",
                           comparison, contrast, n_per_arm, n_expected, n)) %>% pull(txt)
    fail <- c(fail, sprintf("results_master n_per_arm disagrees with the registry on %d inference rows: %s",
                            nrow(bad_n), paste(detail, collapse = "; ")))
  }
}

if (length(fail)) {
  cat("FAIL:\n"); cat(paste0("  - ", fail, collapse = "\n"), "\n"); quit(status = 1)
}
cat(sprintf("PASS: arm tables on locked roster (%d subtypes, %d cells reconciled)\n",
            length(roster), comp_n))
