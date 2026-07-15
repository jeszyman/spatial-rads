# scripts/aggregate/tests/test_ihw.R
library(testthat); library(data.table)
source("scripts/aggregate/fdr_helpers.R")   # exposes add_ihw
test_that("IHW column present and valid for DE rows", {
  set.seed(1)
  m <- data.table(readout_class="DE", unit="Fibroblast",
                  feature=paste0("g",1:500), contrast="SBRT_vs_Ctrl",
                  pvalue=runif(500), baseMean=10^runif(500,0,3))
  out <- add_ihw(m)
  expect_true("padj_ihw" %in% names(out))
  expect_true(all(out$padj_ihw >= 0 & out$padj_ihw <= 1, na.rm=TRUE))
})
