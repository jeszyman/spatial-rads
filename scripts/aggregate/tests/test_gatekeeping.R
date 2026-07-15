# scripts/aggregate/tests/test_gatekeeping.R
library(testthat); library(data.table)
source("scripts/aggregate/load_hypotheses.R"); source("scripts/aggregate/load_claims.R")
source("scripts/aggregate/fdr_helpers.R")   # exposes add_gatekeeping

test_that("gate opens for a claimed program and BHs its claimed DE children", {
  hyp <- load_hypotheses("config/hypotheses.yaml")
  claims <- load_claims("config/confirmatory_claims.yaml", "config/hypotheses.yaml")
  # gsea parent (Fibrosis_remodeling in Fibroblast, SBRT_vs_Ctrl, enriched) + a claimed DE child gene
  child_gene <- claims[hypothesis=="stromal_fibrosis" & readout_class=="DE" &
                       unit=="Fibroblast", feature][1]
  m <- rbindlist(list(
    data.table(readout_class="gsea", unit="Fibroblast", feature="Fibrosis_remodeling",
               contrast="SBRT_vs_Ctrl", pvalue=1e-4),
    data.table(readout_class="DE", unit="Fibroblast", feature=child_gene,
               contrast="SBRT_vs_Ctrl", pvalue=1e-3)
  ), fill=TRUE)
  out <- add_gatekeeping(m, claims, hyp)
  expect_true(all(c("gate_program","gate_padj","padj_gated") %in% names(out)))
  expect_false(is.na(out[readout_class=="DE", padj_gated][1]))
  expect_equal(out[readout_class=="DE", gate_program][1], "Fibrosis_remodeling")
})

test_that("an unclaimed program does NOT create a gate", {
  hyp <- load_hypotheses("config/hypotheses.yaml")
  claims <- load_claims("config/confirmatory_claims.yaml", "config/hypotheses.yaml")
  # Fibrosis_remodeling in Fibroblast is claimed by stromal_fibrosis DE, so use a genuinely
  # unclaimed case: myofibroblast_expansion claims ONLY substate — its program has no DE claim
  # under a cell_type without a stromal_fibrosis DE claim. Tumor has no DE program claim.
  m <- rbindlist(list(
    data.table(readout_class="gsea", unit="Tumor", feature="Fibrosis_remodeling",
               contrast="SBRT_vs_Ctrl", pvalue=1e-6),
    data.table(readout_class="DE", unit="Tumor", feature="Col1a1",
               contrast="SBRT_vs_Ctrl", pvalue=1e-6)
  ), fill=TRUE)
  out <- add_gatekeeping(m, claims, hyp)
  expect_true(is.na(out[readout_class=="DE", padj_gated][1]))   # Tumor DE not a claimed child
})
