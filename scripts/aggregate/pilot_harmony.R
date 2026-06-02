#!/usr/bin/env Rscript
# Integration bake-off pilot, Harmony arm. From pilot_full.rds: normalize -> HVG ->
# scale -> PCA(30) -> Harmony(slide_id, theta=2) -- slide_id nests dataset so a single
# covariate corrects both axes (the live r-harmony 2.0.2 LAPACK bug blocks >=2
# covariates; see project_harmony_lapack). Exports the 30-dim Harmony latent AND the
# pre-integration PCA latent (unintegrated baseline for the compare step). 30 dims to
# match scVI n_latent.
# Args: <pilot_full.rds> <outdir>
suppressPackageStartupMessages({
  library(Seurat); library(harmony); library(arrow); library(data.table); library(future)
})
set.seed(1)
plan("sequential"); options(future.globals.maxSize = 64 * 1024^3)

a      <- commandArgs(trailingOnly = TRUE)
infile <- a[1]
outdir <- a[2]
NPCS   <- 30

o <- readRDS(infile)
o <- NormalizeData(o, verbose = FALSE)
o <- FindVariableFeatures(o, nfeatures = 2000, verbose = FALSE)
o <- ScaleData(o, verbose = FALSE)
o <- RunPCA(o, npcs = NPCS, verbose = FALSE)
o <- RunHarmony(o, group.by.vars = "slide_id", theta = 2,
                reduction.use = "pca", dims.use = 1:NPCS,
                reduction.save = "harmony", verbose = FALSE)

dump_emb <- function(emb, path) {
  dt <- as.data.table(emb[, 1:NPCS], keep.rownames = "cell")
  arrow::write_parquet(dt, path)
}
dump_emb(Embeddings(o, "harmony"), file.path(outdir, "harmony_latent.parquet"))
dump_emb(Embeddings(o, "pca"),     file.path(outdir, "pca_latent.parquet"))

cat(sprintf("harmony arm done: %d cells x %d dims -> harmony_latent.parquet, pca_latent.parquet\n",
            ncol(o), NPCS))
