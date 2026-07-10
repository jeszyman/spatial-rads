#!/usr/bin/env Rscript
# Per-sample technical-QC metrics (report-only): published SpatialQM/Spatial Touchstone metrics
# (Plummer/Segato-Dezem et al., Nat Biotechnol 2025; vendored in scripts/aggregate/spatialqm_metrics.R)
# computed once per sample. One row per sample; the cross-arm driver (qc_arm_balance.R) later compares
# these across Control/MBRT/SBRT -- the confound check on the day-2 composition (fraction-shift) result.
#   tpc            transcripts per cell (mean), on the analyzed common panel  -- SENSITIVITY
#   med_nFeature   genes detected per cell (median), on the analyzed panel    -- SENSITIVITY
#   sparsity       fraction of zero entries in the gene x cell matrix         -- SENSITIVITY (twin)
#   snr            log10 signal-to-noise vs the negprobe background           -- SpatialQM SigNoiseRatio
#   specificity_fdr global false-discovery rate from the negprobe background  -- SpatialQM specificityFDR
# Signal quantities come from the in-memory common-panel counts (self-consistent, = what we analyze);
# the negprobe background is recovered EXACTLY from the per-cell control fraction:
# propNegative = nCount_negprobes/(nCount_RNA+nCount_negprobes) => nCount_RNA*p/(1-p) = vendor negprobe
# total per cell (no raw-RDS re-read). NEG_N (number of negprobes) from the sample sheet. Report-only.
# Args: <qc_dir> <spatialqm_metrics.R> <samples.tsv> <out.tsv>
suppressMessages({library(Seurat); library(data.table); library(Matrix)})
a <- commandArgs(trailingOnly = TRUE)
QC_DIR <- a[1]; SQM <- a[2]; SAMPLES <- a[3]; OUT <- a[4]
source(SQM)                                            # sqm_sparsity, sqm_snr, sqm_specificity_fdr

ss     <- fread(SAMPLES)
neg_n  <- setNames(as.integer(ss$n_negprobe), as.character(ss$sample_id))   # sample_id -> NEG_N

files <- sort(list.files(QC_DIR, pattern = "\\.qc\\.rds$", full.names = TRUE)); stopifnot(length(files) > 0)
rows <- list()
for (f in files) {
  o  <- readRDS(f); md <- o@meta.data
  cm <- LayerData(o, assay = "RNA", layer = "counts")
  sid <- as.character(md$sample_id[1]); ds <- as.character(md$dataset[1])
  NEG_N <- neg_n[[sid]]; stopifnot(is.finite(NEG_N), NEG_N > 0)

  ## negprobe totals recovered exactly from the control fraction (guard p->1, absent post-QC)
  p <- pmin(md$propNegative, 0.999)
  neg_cell   <- md$nCount_RNA * p / (1 - p)             # = vendor nCount_negprobes per cell
  gene_means <- Matrix::rowMeans(cm)
  gene_total <- sum(cm)

  rows[[f]] <- data.table(
    sample_id       = sid,
    dataset         = ds,
    n_cells         = ncol(cm),
    tpc             = mean(Matrix::colSums(cm)),
    med_nFeature    = as.numeric(median(Matrix::colSums(cm > 0))),
    sparsity        = sqm_sparsity(cm),
    snr             = sqm_snr(gene_means, neg_bg = mean(neg_cell) / NEG_N),
    specificity_fdr = sqm_specificity_fdr(gene_total, exp_neg = sum(neg_cell),
                                          n_genes = nrow(cm), n_neg = NEG_N))
  cat(sprintf("%-28s %-10s n=%7d  tpc=%6.1f  sparsity=%.3f  snr=%.3f  fdr=%.4f\n",
              sid, ds, ncol(cm), rows[[f]]$tpc, rows[[f]]$sparsity, rows[[f]]$snr, rows[[f]]$specificity_fdr))
  rm(o, cm); invisible(gc())
}
tab <- rbindlist(rows)
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
fwrite(tab, OUT, sep = "\t")
cat(sprintf("sample technical metrics: %d samples -> %s\n", nrow(tab), OUT))
