# Confirms the detection/level decomposition is fully retired.
library(testthat); library(data.table)

test_that("decomposition columns are gone from the master", {
  m <- fread("results/aggregate/results_master.tsv")
  for (c in c("detection_padj","level_padj","call_class","mean_among_expr_max"))
    expect_false(c %in% names(m))
})

test_that("detection_test.R is deleted", {
  expect_false(file.exists("scripts/aggregate/detection_test.R"))
})
