#!/usr/bin/env Rscript
# Per-cohort expression cache for the per-cell differential steps.
#
# The merged object carries all 20 samples, so a step that analyses a single cohort
# still pays the whole-cohort read before it subsets: peak memory is tens of GB per
# process and only two or three such processes fit on the box at once, which idles
# most of the cores. This step pays that read once per cohort and writes the cohort's
# expression matrices on their own, so the consumers' peak is their own cohort.
#
# Both RNA layers are kept: smide_de.R fits on `counts`, cellchat_communication.R
# scores on the log-normalised `data` layer. Cells are subset, genes never are --
# cellchat_communication.R measures ligand-receptor coverage against the whole panel.
# Nothing else from the Seurat object is carried; the consumers take cell metadata
# from full_labels.parquet / coords_necrosis.parquet / obs.parquet, not from the object.
#
# Cohort membership is the cohort_samples.tsv whitelist, the same cross-dataset-leakage
# guard the engines use. It is deliberately the looser of the two filters the consumers
# apply: they additionally restrict to the registry's contrast condition levels, which
# is a subset of the whitelist, so the cache is a superset of what each consumer keeps.
#
# Args: <merged.rds> <full_labels.parquet> <cohort_samples.tsv> <cohort> <out.rds>
suppressPackageStartupMessages({
  library(data.table); library(arrow); library(Matrix); library(SeuratObject)
})

# Consumers detect a cache by this class rather than by path, so either the merged
# object or a cache can be passed to them as the same argument.
CACHE_CLASS <- "cohort_counts_cache"
GB <- 1024^3

a <- commandArgs(trailingOnly = TRUE)
if (length(a) != 5) {
  stop("usage: build_cohort_counts.R <merged.rds> <full_labels.parquet> ",
       "<cohort_samples.tsv> <cohort> <out.rds>")
}
rds_path <- a[1]; labels_path <- a[2]; cohort_tsv <- a[3]
cohort_name <- a[4]; out_rds <- a[5]

cs <- fread(cohort_tsv)
stopifnot(all(c("cohort", "sample_id") %in% names(cs)))
sids <- unique(cs[cohort == cohort_name, sample_id])
if (!length(sids))
  stop(sprintf("cohort '%s' has no samples in %s; known cohorts: %s", cohort_name,
               basename(cohort_tsv), paste(sort(unique(cs$cohort)), collapse = ", ")),
       call. = FALSE)

lab_cells <- as.data.table(read_parquet(labels_path, col_select = "cell"))$cell

seu <- readRDS(rds_path)
if (!inherits(seu, "Seurat"))
  stop(sprintf("%s is not a Seurat object (class: %s)", rds_path,
               paste(class(seu), collapse = "/")), call. = FALSE)
cnt <- SeuratObject::LayerData(seu, assay = "RNA", layer = "counts")
dat <- SeuratObject::LayerData(seu, assay = "RNA", layer = "data")
if (is.null(cnt) || !nrow(cnt)) stop("merged object has no RNA 'counts' layer", call. = FALSE)
if (is.null(dat) || !nrow(dat)) stop("merged object has no RNA 'data' layer", call. = FALSE)
stopifnot(identical(dim(cnt), dim(dat)), identical(colnames(cnt), colnames(dat)))

n_in <- ncol(cnt)
full_gb <- as.numeric(object.size(seu)) / GB

# build_merged.R writes sample_id as object metadata and encodes it in the cell key
# (sam000N_<orig>); prefer the column and fall back to the key so the cache can also
# be built from an object that predates the column.
cell_sample <- if ("sample_id" %in% names(seu[[]])) as.character(seu$sample_id) else
  sub("_.*$", "", colnames(cnt))
keep <- which(cell_sample %in% sids & colnames(cnt) %in% lab_cells)
if (!length(keep))
  stop(sprintf("cohort '%s' matched 0 of the %d cells in %s", cohort_name, n_in,
               basename(rds_path)), call. = FALSE)

cnt <- cnt[, keep, drop = FALSE]
dat <- dat[, keep, drop = FALSE]
kept_samples <- sort(unique(cell_sample[keep]))
# Dropped here rather than at the end of the script: holding the merged object and the
# subset at the same time is the peak of this job, and everything after this point is
# bounded by the subset alone.
rm(seu); invisible(gc())

cache <- list(
  cohort  = cohort_name,
  samples = kept_samples,
  counts  = cnt,
  data    = dat,
  source  = normalizePath(rds_path),
  built   = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"))
class(cache) <- CACHE_CLASS

dir.create(dirname(out_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(cache, out_rds)

cache_gb <- as.numeric(object.size(cache)) / GB
disk_gb  <- file.size(out_rds) / GB
cat(sprintf("build_cohort_counts: cohort=%s | %d cells in, %d retained, %d genes, %d samples (%s)\n",
            cohort_name, n_in, ncol(cnt), nrow(cnt), length(kept_samples),
            paste(kept_samples, collapse = ", ")))
cat(sprintf("build_cohort_counts: merged object %.2f GB in memory -> cache %.2f GB in memory, %.2f GB on disk (%s)\n",
            full_gb, cache_gb, disk_gb, out_rds))
cat(sprintf("build_cohort_counts: R max-used heap %.2f GB\n", sum(gc()[, 6]) / 1024))
