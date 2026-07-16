# Validates muscat differential-detection output and its ingestion into the master.
library(testthat); library(data.table)

test_that("DD table is detection-only and well-formed", {
  d <- fread("results/aggregate/differential_detection.tsv")
  expect_true(all(c("cell_type","contrast","gene","dd_log2fc","dd_padj") %in% names(d)))
  expect_true(all(d$dd_padj >= 0 & d$dd_padj <= 1, na.rm = TRUE))
  expect_gt(nrow(d), 0)
})

test_that("detection rows reach the master, always exploratory", {
  m <- fread("results/aggregate/results_master.tsv")
  det <- m[readout_class == "detection"]
  expect_gt(nrow(det), 0)
  # detection is never confirmatory: no hypothesis tag (NA/"" after TSV round-trip)
  expect_true(all(is.na(det$hypothesis) | det$hypothesis == ""))
  expect_true(all(det$tier == "exploratory"))
})
