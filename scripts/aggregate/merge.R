#!/usr/bin/env Rscript
# Stage 0 single-pass merge for aggregate.smk. All 20 flank scored.rds share the
# identical 950-gene common panel, so we bypass Seurat v5 merge() (pathologically
# slow renaming 3.3M barcodes across split layers) and instead column-bind the sparse
# counts directly, then build ONE clean Seurat v5 object with joined layers and
# re-normalize. slide_id (absent from per-cell metadata but needed by
# Harmony/propeller/DESeq2) is joined from the master samplesheet. Metadata schemas
# differ slightly between M01/M02 -> unioned via rbindlist(fill=TRUE).
# Args: <samplesheet.tsv> <scale_factor> <out_merged.rds> <out_summary.tsv> <scored.rds...>
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
})

args         <- commandArgs(trailingOnly = TRUE)
ss_path      <- args[1]
scale_factor <- as.numeric(args[2])
out_rds      <- args[3]
out_tsv      <- args[4]
out_meta     <- args[5]                          # slim per-cell metadata cache
rds_paths    <- args[-(1:5)]

peak_rss_gb <- function() {
  hwm <- grep("^VmHWM:", readLines("/proc/self/status"), value = TRUE)
  as.numeric(gsub("[^0-9]", "", hwm)) / 1024 / 1024
}

ss        <- fread(ss_path)
sids      <- sub("\\.norm\\.rds$", "", basename(rds_paths))
slide_map <- setNames(ss$slide_id, ss$sample_id)
stopifnot(all(sids %in% ss$sample_id))

counts_list <- vector("list", length(sids))
md_list     <- vector("list", length(sids))
summ        <- vector("list", length(sids))
ref_genes   <- NULL

for (i in seq_along(sids)) {
  s <- sids[i]
  o <- readRDS(rds_paths[i])
  cnt <- LayerData(o, assay = "RNA", layer = "counts")          # 950 x n, sparse

  if (is.null(ref_genes)) {
    ref_genes <- rownames(cnt)
  } else {
    stopifnot(setequal(rownames(cnt), ref_genes))               # same panel
    if (!identical(rownames(cnt), ref_genes)) cnt <- cnt[ref_genes, , drop = FALSE]
  }
  colnames(cnt) <- paste0(s, "_", colnames(cnt))                # globally-unique barcodes

  md <- o@meta.data
  md$nCount_RNA   <- NULL                                       # recomputed by CreateSeuratObject
  md$nFeature_RNA <- NULL
  md$orig.ident   <- NULL
  md$slide_id     <- slide_map[[s]]                             # join sample-level field
  rownames(md)    <- colnames(cnt)

  summ[[i]] <- data.frame(
    sample_id            = s,
    n_cells              = ncol(cnt),
    n_genes_detected     = sum(Matrix::rowSums(cnt) > 0),
    mean_counts          = round(mean(Matrix::colSums(cnt)), 1),
    max_counts           = max(Matrix::colSums(cnt)),
    frac_negative_probes = round(mean(o$propNegative), 4),
    ram_peak_gb          = NA_real_,
    stringsAsFactors     = FALSE
  )
  counts_list[[i]] <- cnt
  md_list[[i]]     <- md
  rm(o, cnt, md)
}

cell_names <- unlist(lapply(counts_list, colnames), use.names = FALSE)
big_counts <- do.call(cbind, counts_list); rm(counts_list)
big_md     <- as.data.frame(rbindlist(md_list, fill = TRUE)); rm(md_list)
rownames(big_md) <- cell_names
invisible(gc())

merged <- CreateSeuratObject(counts = big_counts, meta.data = big_md,
                             project = "spatial_rads_flank")
rm(big_counts, big_md); invisible(gc())
merged <- NormalizeData(merged, normalization.method = "LogNormalize",
                        scale.factor = scale_factor, verbose = FALSE)

# --- FK / integrity validation ---
md <- merged@meta.data
stopifnot(
  all(md$sample_id %in% ss$sample_id),    # every cell maps to a known sample
  !anyNA(md$slide_id),                    # slide_id populated for all cells
  identical(rownames(merged), ref_genes), # gene panel preserved + ordered
  ncol(merged) == sum(vapply(summ, function(x) x$n_cells, numeric(1)))
)

saveRDS(merged, out_rds)

# --- slim per-cell metadata cache for downstream metadata-only tracks ---
# (composition / niche_composition / pseudobulk colData read this 30 MB table
# instead of re-loading the ~17 GB merged object).
meta_cols <- c("sample_id", "condition", "treatment",
               "timepoint_h", "dataset", "slide_id", "x_slide_mm", "y_slide_mm")
meta_out  <- data.table(cell = rownames(md), md[, meta_cols])
fwrite(meta_out, out_meta, sep = "\t")

# --- per-sample summary + aggregate total row ---
tab   <- rbindlist(summ)
total <- data.table(
  sample_id            = "TOTAL",
  n_cells              = sum(tab$n_cells),
  n_genes_detected     = nrow(merged),
  mean_counts          = round(mean(merged$nCount_RNA), 1),
  max_counts           = max(merged$nCount_RNA),
  frac_negative_probes = round(mean(merged$propNegative), 4),
  ram_peak_gb          = round(peak_rss_gb(), 2)
)
out <- rbind(tab, total)
dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)
fwrite(out, out_tsv, sep = "\t")

cat(sprintf("merged %d samples -> %d cells x %d genes | peak RSS %.2f GB\n",
            length(sids), ncol(merged), nrow(merged), total$ram_peak_gb))
