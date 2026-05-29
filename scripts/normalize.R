#!/usr/bin/env Rscript
# Shared per-sample normalization: LogNormalize + variable features (full common panel
# is effectively all-variable on a ~1000-gene targeted panel). Args: <in.qc.rds> <out.norm.rds> <scale_factor> <nfeatures>
suppressMessages(library(Seurat))
args <- commandArgs(trailingOnly = TRUE)
IN <- args[1]; OUT <- args[2]
sf <- as.numeric(args[3]); nfeat <- as.numeric(args[4])

obj <- readRDS(IN)
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = sf, verbose = FALSE)
nfeat <- min(nfeat, nrow(obj))
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = nfeat, verbose = FALSE)
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
saveRDS(obj, OUT)
cat(sprintf("normalized %s: %d cells, %d/%d variable features\n",
            sub("\\.qc\\.rds$", "", basename(IN)), ncol(obj), length(VariableFeatures(obj)), nrow(obj)))
