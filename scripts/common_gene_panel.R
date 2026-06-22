#!/usr/bin/env Rscript
# Common gene panel = Mutter_01 panel INT Mutter_02 panel, both read from the RNA-assay
# rownames of one per-slide raw RDS (both datasets are format=rds; the counts parquet is
# retired). Invariance guard: the freshly-built list must setequal the currently-committed
# common_genes.tsv before overwrite -- this panel is frozen because the locked aggregate
# was typed on it, so identity (not just cardinality) must hold.
# Args: <samplesheet.tsv> <out.tsv>
suppressMessages({library(readr); library(dplyr); library(Seurat)})
args <- commandArgs(trailingOnly = TRUE); SS <- args[1]; OUT <- args[2]
ss <- read_tsv(SS, show_col_types = FALSE)

# One per-slide RDS per dataset suffices for a panel (rownames) lookup; the RNA assay is
# the gene panel (negprobes/falsecode are separate assays in the same object).
rds_genes <- function(ds) {
  p <- ss %>% filter(name == ds) %>% pull(raw_input_path) %>% unique()
  stopifnot(length(p) >= 1, file.exists(p[1]))
  o <- readRDS(p[1]); g <- rownames(o[["RNA"]]); rm(o); gc(); g
}
m01 <- rds_genes("Mutter_01"); m02 <- rds_genes("Mutter_02")
common <- sort(intersect(m01, m02))
cat(sprintf("Mutter_01: %d | Mutter_02: %d | common: %d\n", length(m01), length(m02), length(common)))
stopifnot(length(common) > 500)  # sanity: panels should overlap heavily

## ---- invariance guard vs committed snapshot (identity, not count) ----
if (file.exists(OUT)) {
  old <- readLines(OUT)
  if (!setequal(common, old))
    stop(sprintf("panel invariance FAILED: %d added, %d removed vs committed %s",
                 length(setdiff(common, old)), length(setdiff(old, common)), OUT))
}
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
writeLines(common, OUT)
