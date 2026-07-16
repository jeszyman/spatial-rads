# Validates the CPU-rebuilt merged cohort object matches the locked labels exactly.
# Run on beast (the object is 1.2GB, spatial-rads env only).
library(testthat)
suppressPackageStartupMessages(library(Seurat))
m <- readRDS("/mnt/data/projects/spatial-rads/aggregate/full/merged.rds")

test_that("merged reconciles to the locked label cohort", {
  expect_equal(ncol(m), 3277090L)                     # locked full_labels rowcount
  expect_false(any(is.na(m$cell_subtype)))            # every cell labeled (subset to labeled)
  expect_equal(length(unique(m$cell_subtype)), 15L)   # locked roster
  expect_equal(nrow(m), 950L)                         # common panel
})

test_that("dropped-before-typing samples are excluded", {
  # sam0009/0010/0011 were removed before typing; not in the locked cohort.
  expect_equal(length(unique(m$sample_id)), 20L)
  expect_false(any(c("sam0009","sam0010","sam0011") %in% m$sample_id))
})
