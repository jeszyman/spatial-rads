#!/usr/bin/env Rscript
# Which pseudobulk engine's test statistic is the REPORTED inference for a cohort.
#
# count_engine.R (DESeq2 Wald) and voom_engine.R (limma-voom moderated t) fit the same
# registry design on the same pseudobulk SE with the same filters, so the choice between
# them is a per-cohort property of the realized residual df, not of the readout. The
# workflow declares that choice by running voom_engine for the low-df cohorts only, and
# every consumer reads the declaration off one place: whether the engine directory holds a
# voom_engine_<cohort>.tsv. assemble_results.R reports those cohorts' moderated t and
# carries the Wald p as a secondary column; gsea.R ranks their genes on the same moderated
# t. A sourceable module so neither script carries its own cohort list.
MODERATED_T_LABEL <- "limma_voom_moderated_t"
WALD_LABEL        <- "DESeq2_Wald"

# `_skipped.tsv` siblings are dropped: they name cell types that produced no fit, not a
# cohort whose inference was swapped.
moderated_t_cohorts <- function(engine_dir) {
  f <- list.files(engine_dir, pattern = "^voom_engine_.*\\.tsv$")
  f <- f[!grepl("_skipped\\.tsv$", f)]
  sub("\\.tsv$", "", sub("^voom_engine_", "", f))
}

voom_engine_file <- function(engine_dir, cohort) {
  f <- file.path(engine_dir, sprintf("voom_engine_%s.tsv", cohort))
  if (file.exists(f)) f else NA_character_
}

reported_inference <- function(engine_dir, cohort) {
  if (cohort %in% moderated_t_cohorts(engine_dir)) MODERATED_T_LABEL else WALD_LABEL
}
