#!/usr/bin/env Rscript
# Two published imaging-spatial QC metrics from SpatialQM / Spatial Touchstone, vendored
# (algorithm-unchanged) for the two roles where they are NOT redundant with machinery we
# already run.
#
# Source (pinned for divergence detection):
#   SpatialQM -- Plummer, Segato Dezem et al., "Standardized metrics for assessment and
#     reproducibility of imaging-based spatial transcriptomics datasets," Nature Biotechnology
#     2025, doi:10.1038/s41587-025-02811-9 (peer-reviewed, VOR 2025-12-03).
#   Repo github.com/Center-for-Spatial-OMICs/SpatialQM, R/utils_update_final.R
#   pinned at commit 36e1a59d4ca3 (2026-05-06). getMECR ~L1701, getCorrelation ~L2331.
#   getMECR's algorithm is Hartman & Satija, bioRxiv 2024 (per the SpatialQM source).
#
# Why vendor only these two (not the package, not getMeanSignalRatio):
#   - getMeanSignalRatio (marker-set signal over the negprobe background) is NOT vendored here:
#     it is the metric our tier-1/tier-2 detectability gate already computes as `signal_to_bg`
#     (full_cluster.py / tier2_detectability.py). SpatialQM is the field-standard CITATION for
#     that gate, not a second implementation of it.
#   - The SpatialQM package Imports pull 6 heavy Bioconductor deps (Voyager,
#     SpatialFeatureExperiment, scater, spacexr, bluster, BioQC) used only by other functions;
#     the two below need none, and SpatialQM is GitHub-only (no CRAN/conda), so vendor + pin.
#   - The source getMECR/getCorrelation wrap platform-specific file loaders (Xenium/CosMx/
#     Merscope). Our merged object is already in memory, so only the core algorithm is vendored;
#     the per-sample / reference plumbing lives in the calling driver.

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
