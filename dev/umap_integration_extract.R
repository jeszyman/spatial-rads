#!/usr/bin/env Rscript
# Extract a plot-ready UMAP from the aggregate Harmony embedding checkpoint.
# The pipeline's own UMAP (embed_celltype.R) runs only AFTER the LISI gate, which
# currently fails -- so embed_checkpoint.rds carries the harmony reduction but no
# umap. Here we compute UMAP on a 150k-cell subsample of the harmony coords (uwot,
# Seurat RunUMAP defaults) and write a slim TSV the plot script reads. Heavy step
# (loads the 11 GB checkpoint) is done once; figure iteration reads only the TSV.
suppressPackageStartupMessages({library(Seurat); library(uwot); library(data.table)})
set.seed(1)

CKPT <- "/mnt/data/projects/spatial-rads/aggregate/embed_checkpoint.rds"
OUT  <- "results/aggregate/umap_integration.tsv"
NPCS <- 30
NSUB <- 150000
stopifnot(file.exists(CKPT))

o    <- readRDS(CKPT)
harm <- Embeddings(o, "harmony")[, 1:NPCS]
md   <- o@meta.data
stopifnot(all(c("dataset", "slide_id", "cell_type") %in% colnames(md)))
cat(sprintf("checkpoint: %d cells | datasets %s\n",
            nrow(harm), paste(sprintf("%s=%d", names(table(md$dataset)), table(md$dataset)), collapse = " ")))

pidx <- sample(nrow(harm), min(NSUB, nrow(harm)))   # random subsample (honest proportions)
um <- uwot::umap(harm[pidx, , drop = FALSE], n_neighbors = 30, min_dist = 0.3,
                 metric = "cosine", n_threads = 8, verbose = TRUE)

d <- data.table(UMAP1 = um[, 1], UMAP2 = um[, 2],
                dataset = md$dataset[pidx], slide_id = md$slide_id[pidx],
                cell_type = md$cell_type[pidx])
d <- d[sample(.N)]   # shuffle draw order so neither dataset systematically overplots
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
fwrite(d, OUT, sep = "\t")
cat(sprintf("wrote %s: %d cells, %d cell types\n", OUT, nrow(d), uniqueN(d$cell_type)))
