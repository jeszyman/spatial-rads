#!/usr/bin/env Rscript
# Tier-2 immune subtyping, step 2/2. Cluster-level SingleR against the celldex ImmGen
# reference (label.main: B/T/NK/Macrophage/DC/Monocyte/Granulocyte/...). Annotates each
# immune subcluster from step 1 (one label per subcluster, not per cell) and maps the
# label back to every immune cell. Cluster-level keeps it conservative at panel sparsity
# and is the SingleR-recommended mode for noisy single-cell data.
# Args: <immune_mtx_dir> <immune_subclusters.parquet> <outdir>
suppressPackageStartupMessages({
  library(Matrix); library(arrow); library(data.table)
  library(SingleR); library(celldex)
})

a <- commandArgs(trailingOnly = TRUE)
mtxdir <- a[1]; sub_pq <- a[2]; outdir <- a[3]

cts <- as(Matrix::readMM(file.path(mtxdir, "counts.mtx")), "CsparseMatrix")
feats <- readLines(file.path(mtxdir, "features.tsv"))
bcs   <- readLines(file.path(mtxdir, "barcodes.tsv"))
rownames(cts) <- feats; colnames(cts) <- bcs
cat(sprintf("immune counts: %d genes x %d cells\n", nrow(cts), ncol(cts)))

sub <- as.data.table(read_parquet(sub_pq))
stopifnot(all(colnames(cts) == sub$cell))   # MTX cols and parquet rows share order

# library-size log2 normalization == scuttle::logNormCounts; sparse-preserving
# (transform applies to nonzeros only; structural zeros map to log2(1)=0).
sf <- Matrix::colSums(cts); sf <- sf / mean(sf)
logc <- cts %*% Matrix::Diagonal(x = 1 / sf)
logc@x <- log2(logc@x + 1)

ref <- celldex::ImmGenData()
cat(sprintf("ImmGen ref: %d genes, %d label.main classes\n",
            nrow(ref), length(unique(ref$label.main))))

pred <- SingleR(test = logc, ref = ref, labels = ref$label.main,
                clusters = sub$immune_subcluster)
map <- data.table(immune_subcluster = rownames(pred),
                  immune_subtype = pred$labels,
                  delta = pred$delta.next)
setorder(map, immune_subcluster)
fwrite(map, file.path(outdir, "immune_subtype_singler.tsv"), sep = "\t")
print(map)

sub[map, immune_subtype := i.immune_subtype, on = "immune_subcluster"]
write_parquet(sub[, .(cell, immune_subcluster, immune_subtype, dataset, slide_id)],
              file.path(outdir, "immune_subtypes.parquet"))

cat("\n--- immune subtype composition ---\n")
print(sub[, .N, by = immune_subtype][order(-N)])
cat("\n--- subtype x dataset ---\n")
print(dcast(sub[, .N, by = .(immune_subtype, dataset)],
            immune_subtype ~ dataset, value.var = "N", fill = 0))
cat(sprintf("\nwrote immune_subtypes.parquet (%d cells) + immune_subtype_singler.tsv\n",
            nrow(sub)))
