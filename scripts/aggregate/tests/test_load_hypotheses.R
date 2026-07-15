# scripts/aggregate/tests/test_load_hypotheses.R
library(testthat); library(data.table)
source("scripts/aggregate/load_hypotheses.R")

test_that("registry loads to long tuple form", {
  dt <- load_hypotheses("config/hypotheses.yaml")
  expect_true(all(c("hypothesis","cell_type","contrast","program","gene") %in% names(dt)))
  expect_equal(sort(unique(dt$hypothesis)),
    sort(c("immune_activation","vascular_hypoxia","stromal_fibrosis",
           "myofibroblast_expansion","normal_vessel_sparing","myeloid_M2")))
  # vascular_hypoxia spans two cell types
  expect_setequal(unique(dt[hypothesis=="vascular_hypoxia", cell_type]), c("Endothelial","Tumor"))
  # normal_vessel_sparing is single-contrast
  expect_equal(unique(dt[hypothesis=="normal_vessel_sparing", contrast]), "MBRT_vs_SBRT")
})

test_that("programs are distinct-listable", {
  hp <- hypotheses_programs("config/hypotheses.yaml")
  expect_true(all(c("hypothesis","cell_type","contrast","program") %in% names(hp)))
  expect_true(nrow(hp) > 0)
})
