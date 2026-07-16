#!/usr/bin/env Rscript
# CPU rebuild of the merged cohort object from the 23 per-sample norm.rds plus the
# locked cell labels. NOT a re-typing run: full_labels.parquet is read as-is. The
# norm.rds colnames are prefixed with the sample id (add.cell.ids) to align to the
# parquet `cell` key (sam000N_<orig>), which is otherwise a zero-overlap join.
suppressPackageStartupMessages({library(Seurat); library(arrow); library(data.table)})

NORM <- "/mnt/data/projects/spatial-rads/processing/norm"
LAB  <- "results/aggregate/full_labels.parquet"
OUT  <- "/mnt/data/projects/spatial-rads/aggregate/full/merged.rds"

labs  <- as.data.table(read_parquet(LAB))
setkey(labs, cell)
files <- list.files(NORM, pattern="\\.norm\\.rds$", full.names = TRUE)
sids  <- sub("\\.norm\\.rds$", "", basename(files))

# Subset each per-sample object to its labeled cells BEFORE merging. Subsetting the
# 20 labeled samples is far cheaper than one subset() of the 3.45M-cell merged object
# (Seurat v5 subset rebuilds layers). Samples absent from the locked label set
# (sam0009/0010/0011 were dropped before typing) are skipped entirely; the locked
# cohort is the analysis unit.
labeled_sids <- unique(sub("_.*$", "", labs$cell))
keep_i <- sids %in% labeled_sids
files <- files[keep_i]; sids <- sids[keep_i]
cat(sprintf("samples: %d labeled of %d on disk (skipped: %s)\n",
            length(sids), sum(keep_i) + sum(!keep_i),
            paste(setdiff(sub("\\.norm\\.rds$","",basename(list.files(NORM,pattern="norm.rds"))), sids), collapse=",")))
n_pre <- 0L
objs <- Map(function(f, sid) {
  o <- readRDS(f)
  o <- RenameCells(o, add.cell.id = sid)          # sam000N_<orig>, matches parquet key
  n_pre <<- n_pre + ncol(o)
  keep <- colnames(o)[colnames(o) %in% labs$cell]
  subset(o, cells = keep)
}, files, sids)

merged <- merge(objs[[1]], objs[-1])
merged <- JoinLayers(merged)

md <- labs[data.table(cell = colnames(merged))]   # keyed join
merged$compartment  <- md$compartment
merged$cell_subtype <- md$cell_subtype
merged$sample_id    <- sub("_.*$", "", colnames(merged))

dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
saveRDS(merged, OUT)
cat(sprintf("merged: %d cells kept (%d dropped unlabeled), %d subtypes\n",
            ncol(merged), n_pre - ncol(merged),
            length(unique(na.omit(merged$cell_subtype)))))
