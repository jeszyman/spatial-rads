#!/usr/bin/env Rscript
# Published imaging-spatial QC metrics from SpatialQM / Spatial Touchstone, vendored
# (algorithm-unchanged): sample-level segmentation contamination (MECR), immune-annotation
# cross-check (pseudobulk correlation), and the per-sample technical-QC metrics (sparsity,
# signal-to-noise, global specificity FDR) that feed the cross-arm balance check.
#
# Source (pinned for divergence detection):
#   SpatialQM -- Plummer, Segato Dezem et al., "Standardized metrics for assessment and
#     reproducibility of imaging-based spatial transcriptomics datasets," Nature Biotechnology
#     2025, doi:10.1038/s41587-025-02811-9 (peer-reviewed, VOR 2025-12-03).
#   Repo github.com/Center-for-Spatial-OMICs/SpatialQM, R/utils_update_final.R
#   pinned at commit 36e1a59d4ca3 (2026-05-06). getMECR ~L1701, getCorrelation ~L2331,
#   getSparsity ~L2037, getMeanSignalRatio ~L1334, getGlobalFDR ~L743.
#   getMECR's algorithm is Hartman & Satija, bioRxiv 2024 (per the SpatialQM source).
#
# Why vendor the algorithms (not install the package):
#   - SpatialQM is GitHub-only (no CRAN/conda) so installing breaks the conda pin, and its
#     Imports pull 6 heavy Bioconductor deps (Voyager, SpatialFeatureExperiment, scater,
#     spacexr, bluster, BioQC) used only by unrelated functions; the metrics below need none.
#   - The source functions wrap platform-specific file loaders (Xenium/CosMx/Merscope raw
#     output dirs). Our objects are already adapted + in memory, so only the core algorithm is
#     vendored; the per-sample / reference plumbing lives in the calling driver.
#   - getMeanSignalRatio is vendored here (sqm_snr) for the per-SAMPLE technical-QC role. It is
#     the same formula our tier-1/2 detectability gate applies at CLUSTER grain (`signal_to_bg`
#     in full_cluster.py / tier2_detectability.py) -- a different application, not a duplicate.
#   - getEntropy/getComplexity are deliberately NOT vendored: getEntropy needs BioQC (uninstalled)
#     and getComplexity's source cumsum is gene-order-dependent (correct only on a pre-sorted
#     matrix); both are low-decision-weight descriptors here, so they are omitted rather than
#     shipped subtly wrong.

suppressPackageStartupMessages({ library(Matrix) })

# --- SpatialQM getMECR (Hartman & Satija 2024) -- SAMPLE-LEVEL segmentation-contamination QC ----
# Mutually Exclusive Co-expression Rate: over all marker pairs from DIFFERENT cell types, the
# Jaccard co-expression sum(g1>0 & g2>0)/sum(g1>0 | g2>0), averaged. LOW = clean, separable
# populations (markers of different types rarely co-occur in one cell); HIGH = markers bleed
# across cells (segmentation contamination). In the SpatialQM source getMECR returns ONE value
# per SAMPLE -- it is a whole-sample QC, not a per-lineage gate -- so the driver calls this once
# per slide/sample and compares MECR comparatively across slides. Verbatim logic incl. the
# source's 25-marker cap. marker_df (gene, cell_type) is the published "custom marker table"
# input -- here our compartment lineage markers.
sqm_mecr <- function(counts, marker_df, max_markers = 25, seed = 0) {
  rownames(marker_df) <- marker_df$gene
  genes <- intersect(rownames(counts), marker_df$gene)
  if (length(genes) > max_markers) { set.seed(seed); genes <- sample(genes, max_markers) }
  mtx <- as.matrix(counts[genes, , drop = FALSE])
  rates <- c()
  for (g1 in genes) for (g2 in genes) {
    if (g1 != g2 && g1 > g2 && marker_df[g1, "cell_type"] != marker_df[g2, "cell_type"]) {
      c1 <- mtx[g1, ]; c2 <- mtx[g2, ]
      rates <- c(rates, sum(c1 > 0 & c2 > 0) / sum(c1 > 0 | c2 > 0))
    }
  }
  round(mean(rates), 3)
}

# --- SpatialQM getCorrelation -- IMMUNE annotation cross-check vs a reference -----------------
# Pseudobulk Spearman correlation over genes common to query and reference. This is the
# annotation-accuracy metric in the BMC Xenium benchmark (Cheng et al. 2025; >0.7 = high).
# Used ONLY where a non-circular reference exists: the immune SingleR/ImmGen labels cross-checked
# against the SAME ImmGen reference. NOT used for stroma (no non-circular bulk reference -- the
# reason MouseRNAseqData was dropped) and NOT part of the retain/drop detectability gate.
sqm_pseudobulk_cor <- function(counts_query, counts_ref) {
  g <- intersect(rownames(counts_query), rownames(counts_ref))
  cor(Matrix::rowMeans(counts_query[g, , drop = FALSE]),
      Matrix::rowMeans(counts_ref[g, , drop = FALSE]), method = "spearman")
}

# ============================ PER-SAMPLE TECHNICAL-QC METRICS =================================
# One scalar per sample; the calling driver (scripts/sample_metrics.R) loops the per-sample QC
# objects and the cross-arm driver compares these across Control/MBRT/SBRT. Report-only.

# --- SpatialQM getSparsity (coop::sparsity, source ~L2037) -- fraction of zero entries in the -----
# gene x cell counts matrix. HIGH = sparser (fewer genes detected per cell); the structural twin of
# per-cell sensitivity. coop::sparsity(m) = 1 - nnz/length, reimplemented directly on the sparse
# matrix (coop is a SpatialQM dep, uninstalled). Verbatim definition, exact.
sqm_sparsity <- function(counts) {
  1 - Matrix::nnzero(counts) / prod(dim(counts))
}

# --- SpatialQM SigNoiseRatio (getMeanSignalRatio, source ~L1394) -- per-SAMPLE signal/background ---
# mean over genes of log10(per-gene mean + 0.1) minus log10(negprobe background + 0.1): a log10
# signal-to-noise. gene_means = per-gene mean counts on the analyzed panel; neg_bg = mean per-
# negprobe-per-cell background. The source takes per-negprobe row means then means log10; post-adapter
# the negprobe assay is collapsed to a per-cell control fraction, so negprobes are pooled to one
# background (they are exchangeable by design) -- the one intentional deviation from the source.
sqm_snr <- function(gene_means, neg_bg) {
  mean(log10(gene_means + 0.1)) - log10(neg_bg + 0.1)
}

# --- SpatialQM specificityFDR (getGlobalFDR, source ~L816) -- per-SAMPLE global false-discovery -----
# rate. Verbatim formula: (exp_neg/(gene_total+exp_neg)) * (n_genes/n_neg) / 100, where exp_neg /
# gene_total are total negprobe / gene transcript counts and n_genes / n_neg the feature counts. The
# source derives exp_neg/gene_total from a transcript table via `.N by target`, which equals per-
# feature count sums, so this reproduces EXACTLY from the count matrix (no transcript table needed).
# ~0.05 => ~5% of per-gene transcripts expected to be background false positives.
sqm_specificity_fdr <- function(gene_total, exp_neg, n_genes, n_neg) {
  (exp_neg / (gene_total + exp_neg)) * (n_genes / n_neg) / 100
}
