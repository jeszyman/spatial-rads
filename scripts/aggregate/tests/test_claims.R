# scripts/aggregate/tests/test_claims.R  (run from repo root, not test_file — see Global Constraints)
library(testthat); library(data.table)
source("scripts/aggregate/load_hypotheses.R"); source("scripts/aggregate/load_claims.R")

claims <- load_claims("config/confirmatory_claims.yaml", "config/hypotheses.yaml")

test_that("claims expand to one row per readout x unit x contrast x feature", {
  expect_true(all(c("hypothesis","readout_class","unit","contrast","feature") %in% names(claims)))
  # a stromal DE gene claim exists
  expect_true(claims[hypothesis=="stromal_fibrosis" & readout_class=="DE" &
                     unit=="Fibroblast" & feature=="Col1a1", .N] > 0)
  # myeloid is pathway-only (no DE claim)
  expect_equal(claims[hypothesis=="myeloid_M2" & readout_class=="DE", .N], 0)
  expect_true(claims[hypothesis=="myeloid_M2" & readout_class=="pathway", .N] > 0)
})

test_that("no result is claimed by more than one hypothesis", {
  expect_silent(assert_no_multiclaim(claims))
})

test_that("apply_claims joins hypothesis onto matching results only", {
  m <- data.table(readout_class="DE", unit="Fibroblast", feature="Col1a1",
                  contrast="SBRT_vs_Ctrl", hypothesis=NA_character_)
  out <- apply_claims(m, claims)
  expect_equal(out$hypothesis, "stromal_fibrosis")
  # a non-claimed row stays NA (exploratory)
  m2 <- data.table(readout_class="DE", unit="Tumor", feature="Foxp3",
                   contrast="SBRT_vs_Ctrl", hypothesis=NA_character_)
  expect_true(is.na(apply_claims(m2, claims)$hypothesis))
})
