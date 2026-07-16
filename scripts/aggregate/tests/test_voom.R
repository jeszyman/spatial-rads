# Validates the limma-voom / DESeq2 concordance cross-check output.
library(testthat); library(data.table)
test_that("concordance table is well-formed and methods agree", {
  c <- fread("results/aggregate/deseq2_voom_concordance.tsv")
  expect_true(all(c("cell_type","contrast","spearman_rho","n_genes") %in% names(c)))
  expect_gt(nrow(c), 0)
  expect_gt(median(c$spearman_rho, na.rm = TRUE), 0.8)   # DESeq2 and voom broadly agree
})
