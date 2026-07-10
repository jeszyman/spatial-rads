# Plan — published imaging-ST QC metrics (SpatialQM) as pipeline steps

**Status: BUILT 2026-07-10.** Four QC signals from Plummer et al. (SpatialQM / Spatial Touchstone,
*Nat Biotechnol* 2025, doi:10.1038/s41587-025-02811-9) wired as report-only steps in the existing
workflows — a confound check sitting next to the day-2 composition claim it guards.

## Why

The day-2 headline is a **cell-state composition shift** (more cells expressing collagen/Acta2;
per-expresser level n.s.). That detection-fraction claim is most directly threatened by three
technical variables the pipeline did not gate on: per-cell **sensitivity** (the axis that already
broke per-cell typing), the **negprobe/false-positive floor** the detection-vs-level test measures
fractions against, and **segmentation contamination** (MECR — computed, but never compared across
arms). All CosMx here, so the paper's dominant variance driver (platform) is constant; the residual
to police is per-sample/slide, and the decision-relevant output is **balance across
Control/MBRT/SBRT**.

## Metrics built (the endorsed set)

| metric | what | role |
|---|---|---|
| TPC, median features/cell | transcripts & genes per cell | sensitivity |
| sparsity | fraction of zero entries | sensitivity (twin) |
| SNR (SigNoiseRatio) | log10 signal vs negprobe background | negprobe floor |
| specificity FDR (getGlobalFDR) | global false-discovery rate | negprobe floor |
| MECR | cross-lineage marker bleed | contamination (already run) |
| replicate reproducibility | per-arm pseudobulk concordance + metric-PCA | outlier-slide check |

**Dropped** (Tier-C, low decision weight, and reproduction hazards): entropy (needs uninstalled
BioQC), complexity (source `cumsum` is gene-order-dependent — correct only on a pre-sorted matrix),
dynamic range. **Not computable** (flag, don't attempt): TPN and FTC need off-cell / nucleus-level
transcripts we don't have (no transcript coordinates delivered — same reason FastReseg was
infeasible). TPC therefore carries the sensitivity load, with the caveat that it is
segmentation-coupled (unlike the segmentation-free TPN).

## Decisions

**1. Vendor the algorithms, do NOT install SpatialQM.** The package is GitHub-only (breaks the conda
pin), Imports 6 heavy Bioconductor deps (Voyager, SpatialFeatureExperiment, scater, spacexr,
bluster, BioQC) used only by unrelated functions, and its metric entry points are platform file
loaders (raw Xenium/CosMx output dirs) that don't fit our adapted in-memory Seurat objects. The
metric formulas were transcribed verbatim from the pinned source (commit `36e1a59d4ca3`,
`R/utils_update_final.R`: getSparsity ~L2037, getMeanSignalRatio ~L1394, getGlobalFDR ~L816) into
`scripts/aggregate/spatialqm_metrics.R` (`sqm_sparsity`/`sqm_snr`/`sqm_specificity_fdr`), joining the
already-vendored `sqm_mecr`/`sqm_pseudobulk_cor`. This is the same decision made for MECR — we use
SpatialQM's cited *definitions*, not its dependency tree. The FDR reproduces exactly: the source's
transcript-table `.N by target` equals per-feature count sums, all derivable from the count matrix +
control fraction. The one intentional deviation: SNR's negprobe background is pooled across the
10 negprobes (exchangeable by design), since the adapter collapses the negprobe assay to a per-cell
control fraction.

**2. Compute per-sample → preprocessing; compare across arms → aggregate.** Per-sample metric compute
is a deterministic function of one QC object (report-only, alongside probe/control/contamination
QC). The arm-balanced comparison and replicate reproducibility need the arm design + all samples
together — QC *of the differential result* — so they live in the aggregate workflow next to the
differential layer.

## Wiring

- **`preproc_sample_metrics`** (org `preprocessing.smk`, tangles `scripts/sample_metrics.R`) → one
  row per sample in `results/processing/sample_tech_metrics.tsv`. Negprobe background recovered
  exactly from `propNegative` + `nCount_RNA` (= vendor negprobe total per cell); NEG_N from the
  sample sheet. No raw-RDS re-read.
- **`agg_qc_arm_balance`** (`scripts/aggregate/qc_arm_balance.R`) → `qc_arm_balance.tsv` +
  `plots/qc_arm_balance.png`. Joins the per-sample metrics + MECR to the arm design, M02 day-2
  (n=4/arm); effect-size-forward (per-arm mean, between-arm gap in within-arm sd units), no p-values.
- **`agg_qc_reproducibility`** (`scripts/aggregate/qc_reproducibility.R`) → `qc_reproducibility.tsv` +
  `plots/qc_reproducibility.png`. Per-arm pseudobulk Spearman (cell types collapsed from
  `pseudobulk_se.rds`) + PCA of the scaled technical metrics.

## Result (M02 day-2)

- **5/6 metrics balanced across arms.** SBRT has neither higher sensitivity (TPC/SNR) nor higher
  contamination (MECR) than Control — the fraction-shift signal is **not** a QC artifact. The one
  flagged metric (specificity FDR, standardized gap 1.26) is immaterial: all arms ~7e-4 and SBRT is
  the *cleanest*.
- **Replicate concordance high** (per-arm mean Spearman: Control 0.94 / MBRT 0.94 / SBRT 0.97 — SBRT
  tightest, matching its strong consistent fibrosis signal), **0 outlier slides**; PCA shows no arm
  separates on a batch axis (PC1 = sensitivity gradient, both arms span it).

Net: the day-2 composition result survives the QC confound check. Extends the robustness pass in
`plan-mbrt-vs-sbrt.md`.
