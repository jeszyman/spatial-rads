# Validates the Task 6 schema fixes on the re-emitted master (run after Task 8 emit).
library(testthat); library(data.table)
m <- fread("results/aggregate/results_master.tsv")

test_that("n_per_arm is the true per-arm count, n_samples_used is separate", {
  expect_true("n_samples_used" %in% names(m))
  expect_true(all(m[readout_class %in% c("DE","composition","pathway"), n_per_arm] == 4L))
})

test_that("padj_exploratory dropped, gate + IHW columns present", {
  expect_false("padj_exploratory" %in% names(m))
  expect_true(all(c("gate_program","gate_padj","padj_gated","padj_ihw") %in% names(m)))
})

test_that("tumor compartment uses one roster name (no Epithelial cells alias)", {
  expect_false("Epithelial cells" %in% m$unit)
})
