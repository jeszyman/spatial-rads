#!/usr/bin/env Rscript
# Build the cell-typing reference for Seurat anchor label transfer, on the common gene panel:
#   pooled Mutter_01 cells labeled with cell_type = Yi ImmGen label, with confident-epithelial
#   a/b/NA cells relabeled tumor_epithelial (identical rule to celltype.R's Mutter_01 branch, so
#   the reference labels match how M01 cells are actually typed). Per-label subsampled for balance
#   and anchor-finding speed; LogNormalized. -> reference Seurat object (reference.rds).
# Args: <yi_reference.tsv> <common_genes.tsv> <out_reference.rds> <m01_norm_rds...>
suppressMessages({library(Seurat); library(Matrix); library(readr)})
args <- commandArgs(trailingOnly = TRUE)
YI <- args[1]; GENES <- args[2]; OUT <- args[3]; M01_RDS <- args[-(1:3)]
MIN_CELLS <- 100      # drop Yi label classes rarer than this in the pool
N_TOTAL   <- 20000    # uniform random subsample of the pool (preserves natural abundance; bounds anchor cost / per-job memory: 40k OOM-killed sam0023 transfer at --cores 1)
EPI <- c("Epcam", "Krt8", "Krt18", "Krt19")
set.seed(1)

common <- readLines(GENES)
yi <- read_tsv(YI, show_col_types = FALSE)
yi_lab <- setNames(yi$ImmuneAtlas_ImmGen_Main_cell_Types, yi$cell_id)

cat(sprintf("Pooling %d Mutter_01 objects...\n", length(M01_RDS)))
mats   <- lapply(M01_RDS, function(f) GetAssayData(readRDS(f), layer = "counts"))
genes  <- Reduce(intersect, c(list(common), lapply(mats, rownames)))
counts <- do.call(cbind, lapply(mats, function(m) m[genes, , drop = FALSE]))
rm(mats); gc()

## cell_type = Yi label, with confident-epithelial a/b/NA -> tumor_epithelial (== celltype.R M01 logic)
ct  <- unname(yi_lab[colnames(counts)])
epi <- intersect(EPI, genes); stopifnot(length(epi) >= 2)
epi_detect <- Matrix::colSums(counts[epi, , drop = FALSE] > 0)
ct[(ct %in% c("a", "b") | is.na(ct)) & epi_detect >= 2] <- "tumor_epithelial"
cat(sprintf("Pooled: %d genes x %d cells; %d candidate types\n",
            nrow(counts), ncol(counts), length(unique(na.omit(ct)))))

## keep labels with >= MIN_CELLS; uniform random subsample PRESERVING natural abundance.
## Per-label balancing over-represents rare ImmGen subsets, flattening class priors so transfer
## scatters epithelial cells into spurious thymocyte/spleen types; natural abundance keeps the
## tumor/immune composition that anchors the transfer (smoke test sam0018: balancing put only
## 8% of epithelial-marker+ cells in tumor_epithelial; natural abundance restores it).
tc   <- table(ct)
keep <- names(tc[tc >= MIN_CELLS])
pool <- which(ct %in% keep)
if (length(pool) > N_TOTAL) pool <- sample(pool, N_TOTAL)
counts <- counts[, pool]; ct <- ct[pool]
cat(sprintf("Reference: %d cells across %d types (natural abundance; >= %d cells/type in pool; dropped %d rare/NA types)\n",
            length(pool), length(unique(ct)), MIN_CELLS, length(tc) - length(keep)))

ref <- CreateSeuratObject(counts = counts)
ref$cell_type <- ct
ref <- NormalizeData(ref, verbose = FALSE)
ref <- FindVariableFeatures(ref, nfeatures = length(genes), verbose = FALSE)
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
saveRDS(ref, OUT)
cat("wrote", OUT, "-", ncol(ref), "cells,", length(keep), "types\n")
