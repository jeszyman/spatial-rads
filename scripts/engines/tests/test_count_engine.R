#!/usr/bin/env Rscript
# Task 3 TDD: count_engine must reproduce the CURRENT deg_pseudobulk.R output exactly --
# apeglm log2FC + lfcSE and the unshrunk Wald stat/p, per (cell_type x contrast x gene).
# Reference = degs_pseudobulk_m02day2.tsv on disk (present pre-refactor). Live source-compare.
suppressPackageStartupMessages(library(data.table))
ref <- fread("results/aggregate/degs_pseudobulk_m02day2.tsv")   # current script output
eng <- fread("results/aggregate/_count_smoke/degs.tsv")
eng[, cell_type := unit]; eng[, gene := feature_id]
m <- merge(ref[, .(cell_type, contrast, gene, r_lfc = log2FC, r_se = lfcSE, r_stat = stat, r_p = pvalue)],
           eng[,  .(cell_type, contrast, gene, estimate, se, stat, p)],
           by = c("cell_type", "contrast", "gene"))
stopifnot(nrow(m) > 100)
chk <- function(x, y, tol, nm) {
  bad <- sum(abs(x - y) / pmax(abs(x), 1e-9) >= tol & !(is.na(x) & is.na(y)), na.rm = TRUE) +
         sum(is.na(x) != is.na(y))
  if (bad > 0) cat(sprintf("FAIL %s: %d mismatches\n", nm, bad)); bad
}
fail <- chk(m$r_lfc, m$estimate, 1e-6, "log2FC") +
        chk(m$r_se,  m$se,       1e-6, "lfcSE") +
        chk(m$r_stat, m$stat,    1e-6, "stat") +
        chk(m$r_p,    m$p,       1e-6, "pvalue")
if (fail > 0) quit(status = 1L)
cat(sprintf("PASS test_count_engine: %d (cell_type x contrast x gene) rows match current deg_pseudobulk.R\n", nrow(m)))
