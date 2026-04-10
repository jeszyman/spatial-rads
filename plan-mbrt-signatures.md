# Autonomous Run: MBRT Peak/Valley Signature Persistence and Replication

## Scientific Intention

> Do transcriptomic signatures derived from spatially validated MBRT peak and valley zones at 4h remain detectable — as scoreable gene programs — at later timepoints and in independent biological replicates (Mutter_02), and are these spatial-fractionation signatures absent in uniform-dose SBRT?

### Interpretive boundary (MANDATORY framing for all results)

- **At 4h MBRT**: peak/valley labels are validated against gamma-H2AX IHC ground truth. Signature genes are causally linked to spatial dose zones.
- **At all other timepoints**: we SCORE cells for the 4h-derived program but CANNOT confirm those cells were in peaks or valleys. High signature scores may reflect: (a) original peak/valley cells retaining the program, (b) infiltrating cells acquiring it secondarily, (c) microenvironment shaped by original damage. Frame as "signature persistence" NOT "peak/valley identity."
- **Spatial maps of scores**: if striped patterns of signature scores emerge at later timepoints, report this as converging evidence but not proof of spatial memory. Absence of striping does NOT mean absence of the program (cell migration, tissue remodeling).

## Data Map

### Mutter_01 (processed, Seurat object)

- **Clustered object**: `$SCRATCH/seurat_clustered.rds` (1.7 GB, staged from GCS)
- **Peak/valley classification**: `$REPO/dev/peak_valley_analysis/data/mbrt4h_peak_valley.tsv`
- **Peak/valley DEGs (signature source)**: `$REPO/dev/peak_valley_analysis/data/pv_degs_bulk.tsv`
- **Per-cell-type DEGs**: `$REPO/dev/peak_valley_analysis/data/pv_degs_by_celltype.tsv`
- **Stripe model**: `$REPO/dev/peak_valley_analysis/data/stripe_model.rds`
- **Conditions**: Control, MBRT_1h, MBRT_4h, SBRT_4h, MBRT_day2, SBRT_day2, MBRT_day6, SBRT_day6, Tongue_MBRT_day8, Tongue_MBRT_day10, Tongue_SBRT_day10
- **Pathway scores**: TypeI_interferon, TypeII_interferon, STING, DNA_Damage_Repair (pre-computed by Yi, in metadata)

### Mutter_02 (unprocessed, raw Seurat objects)

Staged on VM local disk at `$SCRATCH/mutter02/`:
- `seuratObject_01_Mutter_02_CosMmR.RDS` (393 MB, ~999K cells, 462 FOVs)
- `seuratObject_02_Mutter_02_CosMmR.RDS` (210 MB, ~649K cells, 366 FOVs)
- `seuratObject_03_Mutter_02_CosMmR.RDS` (143 MB, ~352K cells, 281 FOVs)
- `seuratObject_04_Mutter_02_CosMmR.RDS` (315 MB, ~793K cells, 477 FOVs)

**Gene panel**: 950 main panel (Mouse UCC) + 21 custom add-on genes. Same panel as Mutter_01.

**Expected conditions (per Jenn Fazzari Mar 2 2026 email):**
- Blocks 1-3: 0 Gy controls
- Blocks 16, 18: 20 Gy 4h (= CRT/SBRT uniform dose)
- Blocks 19, 20: MBRT 4h
- Blocks 28, 30: 20 Gy 2d (= CRT/SBRT uniform dose)
- Blocks 34, 35: MBRT 2d

**Slides 1-2**: 4T1 model (same as Mutter_01) — primary replication target.
**Slides 3-4**: New tumor model pilot — analyze separately, do NOT pool with 4T1.

**FOV-to-condition mapping**: Provided by Rob/Jenn (file: `$SCRATCH/mutter02/fov_condition_map.tsv`). Yi confirmed she does not have this mapping. The executor MUST use this file to label cells. If the file is missing, STOP.

**Important**: These are RAW Seurat objects — no QC, no cell type annotation, no pathway scores. The executor must compute everything from scratch. Pathway scores should be computed using NanoString standard gene modules (see Pathway Gene Lists below).

### NanoString Pathway Gene Lists (RESOLVED)

Source: `$REPO/config/CosMx-Mouse-Universal-Cell-Characterization-Gene-List-(1).XLSX`, Annotations sheet. These are NanoString's standard panel modules — confirmed as the source of Yi's pre-computed Mutter_01 pathway scores.

**Type I Interferon Signaling (19 genes):**
Bst2, Ifi27, Ifit1, Ifit3/b, Ifitm1, Ifitm3, Ifna, Ifnar1, Irf3, Irf4, Isg15, Jak1, Mx1, Oas1a/g, Oas2, Oas3, Oasl1, Rigi, Tyk2

**Type II Interferon Signaling (17 genes):**
B2m, Cd44, Ciita, Icam1, Ifng, Ifngr1, Ifngr2, Irf3, Irf4, Jak1, Jak2, Ncam1, Oas1a/g, Oas2, Oas3, Oasl1, Vcam1

**DNA Damage Repair (18 genes):**
Abl1, Atm, Atr, Brca1, Cdkn1a, Chek1, Chek2, Crip1, H2az1, Hmgb2, Isg15, Mlh1, Msh2, Parp1, Pcna, Sod1, Trp53, Uba52

**STING**: Not a named NanoString module. Yi may have defined this custom or used a subset of Pattern Recognition Receptors (30 genes). For Mutter_02, compute STING score using cGAS-STING pathway genes available in the panel: Irf3, Ifna, Isg15, Cxcl10, plus any Tmem173/Sting1 if present. Document the gene list used.

**Additional modules relevant to radiation biology:**
- p53 Signaling (28 genes)
- Apoptosis (52 genes)
- Cellular Stress (59 genes)
- Pattern Recognition Receptors (30 genes)

Use `UCell::ScoreSignatures()` (rank-based, no background gene sampling) as the primary scoring method for all signatures and pathway scores. UCell is preferred over AddModuleScore for targeted panels (~1000 genes) because AddModuleScore's random background pool is thin and may include biologically meaningful genes. Also compute AddModuleScore (seed=42) in parallel for comparison, but UCell scores are primary for all downstream analyses.

### Circularity rules

The stripe model (script 05) used **Cdkn1a (p21) expression** as the optimization target to find beam positions. The DDR pathway score was also visualized during stripe fitting.

**Exclude from peak/valley signatures**: All 18 DNA Damage Repair genes listed above. These genes were directly or indirectly involved in defining zone boundaries. Including them would be circular.

**Safe to include**: All other genes in the panel, including Type I/II IFN, STING, apoptosis, etc. These were not used for stripe detection.

### Stripe model geometry

The Mutter_01 stripe model parameters:
- **Tilt**: 15 degrees — this is the tissue MOUNTING ANGLE on the slide, not beam tilt. The beam itself is straight.
- **Spacing**: 1.02mm center-to-center (determined by collimator, should be the same across runs)
- **Peaks**: 4 beam centers at corrected-y 1.14, 2.16, 3.18, 4.20 mm

For Mutter_02 stripe validation: spacing should be conserved (~1.02mm) but tilt and offset MUST be re-fit per slide (different tissue mounting angle). Re-fit using the same p21 optimization approach as script 05. Do NOT transfer Mutter_01 tilt/offset directly.

## Disk Budget

| Item | Size | Notes |
|------|------|-------|
| Mutter_02 RDS files | ~1.1 GB | cp from Box to VM local disk, NOT symlink |
| Mutter_01 Seurat object | 1.7 GB | cp from GCS FUSE to VM local disk |
| Mutter_01 analysis data | ~50 MB | TSVs + RDS from dev/peak_valley_analysis/ |
| Working Seurat objects | ~4 GB | Mutter_02 processed + merged objects |
| Output figures + tables | ~200 MB | |
| **Total needed** | **~7 GB** | |

Confirm >=15 GB free on VM scratch disk before proceeding.

## Step 2 Setup Prerequisites (interactive, BEFORE autonomous run)

These must be completed interactively on the VM before launching `run-mbrt-signatures.sh`:

- [ ] Clone repo on VM, checkout feature branch
- [ ] Conda env `spatial-rads` created from `config/spatial-rads-conda-env.yaml`
- [ ] Install additional packages: `r-harmony`, `bioconductor-deseq2` (or `r-limma`), `r-ucell`
- [ ] Mutter_02 RDS files copied from Box to `$SCRATCH/mutter02/` (cp, NOT symlink)
- [ ] Mutter_01 `seurat_clustered.rds` copied from GCS FUSE to `$SCRATCH/`
- [ ] FOV-to-condition mapping from Rob/Jenn saved as `$SCRATCH/mutter02/fov_condition_map.tsv`
- [ ] Gene list XLSX copied to VM: `$REPO/config/CosMx-Mouse-Universal-Cell-Characterization-Gene-List-(1).XLSX`
- [ ] Load each Mutter_02 RDS, verify: cell counts match QC email (~999K, 649K, 352K, 793K), spatial coordinates present, gene panel size ~971 genes
- [ ] Verify gene panel overlap with Mutter_01 (expect >95%)
- [ ] Verify FOV-condition mapping covers all FOVs in data
- [ ] Disk budget confirmed (>=15 GB free)
- [ ] `~/.claude/settings.local.json` with `{"hooks": {}}` (VM has no Emacs)
- [ ] H2AX IHC status confirmed (from Rob) — document whether available for Mutter_02

## Executor Tasks

### Task 1: Process Mutter_02 (QC, normalize, annotate)

Process each Mutter_02 slide:

1. Load raw Seurat objects from `$SCRATCH/mutter02/`.
2. Apply FOV-to-condition mapping from `fov_condition_map.tsv`. Parse treatment (MBRT/SBRT/NT) and timepoint_h. If mapping is ambiguous, document assumptions.
3. **QC filter** (same criteria as Mutter_01): nCount_RNA > 20, nFeature_RNA > 10. If `qcFlagsCell` column exists, require "Pass". If `propNegative` exists, require < 0.5. Adapt if columns differ and document.
4. **Normalize**: LogNormalize, scale factor 1e4.
5. **Variable features**: VST, 2000 requested.
6. **PCA**: 20 dims.
7. **Clustering**: Louvain, resolution 0.4.
8. **UMAP**: 30 neighbors, 0.3 min_dist.
9. **Cell type annotation**: Use Seurat `TransferData` from Mutter_01 reference to predict labels. Validate with canonical markers (Epcam, Cd3e, Cd4, Cd8a, Cd19, Ptprc/CD45, etc.).
10. **Compute pathway scores**: Use `AddModuleScore()` with the NanoString gene lists (Type I IFN, Type II IFN, DDR, STING, p53, Apoptosis, Cellular Stress). Use seed=42.
11. **Separate 4T1 from new tumor model**: Slides 1-2 = 4T1 (primary). Slides 3-4 = new model (exploratory only, Task 7).

**Outputs:**
- `results/signatures/data/mutter02_qc_summary.tsv`
- `results/signatures/plots/mutter02_umap_landscape.png`
- `results/signatures/plots/mutter02_spatial_celltype.png`
- Processed Mutter_02 Seurat object saved to `$SCRATCH/`

### Task 2: Batch assessment

Assess batch effects before integrating datasets:

1. Merge Mutter_01 and Mutter_02 4T1 objects (no integration, just concatenate).
2. Run PCA + UMAP on merged object.
3. Color by: (a) dataset (Mutter_01 vs Mutter_02), (b) condition, (c) cell type.
4. Quantify mixing: compute silhouette scores by dataset within each cell type. If dataset dominates over biology, flag for integration.
5. If batch effect is strong: apply Harmony integration on the merged object, re-cluster, re-UMAP. Document whether integration changed the landscape.
6. If batch effect is minimal: proceed without integration (preferred — simpler, fewer assumptions).

**Outputs:**
- `results/signatures/plots/batch_umap_merged.png` (2x3 panel)
- `results/signatures/data/batch_assessment.tsv`
- Decision documented in SUMMARY.md

### Task 3: Define peak and valley gene signatures

Extract compact, validated signatures from Mutter_01:

1. Load `pv_degs_bulk.tsv`.
2. **Exclude all 18 DDR genes** (listed in Circularity Rules above) from both signatures.
3. **Per-cell-type signatures (PRIMARY)**: from `pv_degs_by_celltype.tsv`, define separate peak and valley signatures for the top 2-3 most abundant cell types (at least 10 non-DDR genes passing avg_log2FC > 0.25 / < -0.25 and pct threshold > 0.1 per type). Cell-type-stratified signatures avoid the composition confound: bulk peak/valley DEGs may partly reflect cell type proportion differences near beam centers rather than within-cell-type transcriptional changes. These are the primary signatures for all replication claims.
4. **Bulk signatures (SECONDARY/EXPLORATORY)**: all remaining genes with avg_log2FC > 0.25 AND pct.1 > 0.1 (peak) or avg_log2FC < -0.25 AND pct.2 > 0.1 (valley). If this yields <10 or >50 genes, report and use top 20 by effect size as fallback. Bulk signatures are reported for completeness but should not be the basis for replication conclusions.
5. Verify all signature genes exist in Mutter_02 gene panel. Note any missing.

**Outputs:**
- `results/signatures/data/peak_signature_genes.tsv` (gene, avg_log2FC, pct_peak, pct_valley)
- `results/signatures/data/valley_signature_genes.tsv`
- `results/signatures/data/celltype_signatures.tsv` (if applicable)
- `results/signatures/data/excluded_ddr_genes.tsv` (the 18 excluded genes with rationale)

### Task 4: Score signatures across all Mutter_01 timepoints

Signature persistence analysis within Mutter_01:

1. Score every Mutter_01 cell for peak and valley signatures (both cell-type-stratified and bulk) using **UCell::ScoreSignatures()** as primary. Also compute `AddModuleScore()` (seed=42) in parallel for comparison. Report correlation between UCell and AddModuleScore scores.
2. Compute mean signature scores per condition x cell type. Report as table. Use UCell scores for all downstream analyses.
3. **Kinetic plots**: signature score (y-axis) vs timepoint (x-axis), one line per treatment (MBRT vs SBRT vs Control), faceted by cell type. Error bars = cell-level SD (NOT SEM — n=1 per condition, no biological replication within Mutter_01). Plot cell-type-stratified signatures for the corresponding cell type; plot bulk signatures separately.
4. **Spatial signature maps**: for each flank MBRT condition (MBRT_1h, MBRT_4h, MBRT_day2, MBRT_day6), plot cells colored by peak_signature score on spatial coordinates. Same color scale across all timepoints. Do the same for valley_signature.
5. **Striping analysis**: for MBRT conditions only, compute rolling mean of signature score along the corrected-Y axis (same axis used for stripe detection at 4h). If periodic patterns emerge, report the spatial frequency and compare to the known 1.02mm beam spacing. Additionally, compute **Moran's I** spatial autocorrelation of signature scores using a distance weight matrix at the beam spacing scale (~1.02mm) to formally test whether scores are spatially structured. Frame as "consistent with" or "no evidence of" persistent spatial patterning — NEVER claim spatial identity.
6. **SBRT comparison (descriptive, Mutter_01 only)**: score SBRT_4h cells with both signatures. Report as distribution comparison (violin/box plots) against MBRT_4h peak and valley cells. Does SBRT resemble peaks, valleys, or neither?

**Outputs:**
- `results/signatures/data/signature_scores_by_condition.tsv`
- `results/signatures/data/signature_kinetics.tsv`
- `results/signatures/plots/signature_kinetics_peak.png`
- `results/signatures/plots/signature_kinetics_valley.png`
- `results/signatures/plots/spatial_peak_sig_all_timepoints.png` (multi-panel)
- `results/signatures/plots/spatial_valley_sig_all_timepoints.png`
- `results/signatures/plots/striping_analysis.png`
- `results/signatures/plots/sbrt_vs_peak_valley_signatures.png`

### Task 5: Score signatures in Mutter_02

Apply the Mutter_01-derived signatures to the independent dataset:

1. Score all Mutter_02 4T1 cells with **UCell::ScoreSignatures()** (primary) and `AddModuleScore()` (seed=42, secondary) for peak and valley signatures (both cell-type-stratified and bulk).
2. **Spatial maps**: plot signature scores on spatial coordinates for each Mutter_02 condition.
3. **Comparison**: violin/box plots of signature scores in Mutter_02 MBRT vs Control vs CRT/SBRT (if present), faceted by cell type.
4. **Pseudobulk statistical testing**: Aggregate counts by FOV (technical replicates within a slide — acknowledge this limitation). Use DESeq2 or limma-voom for MBRT vs Control. If CRT/SBRT replicates exist, also test MBRT vs CRT. Report adjusted p-values alongside Mutter_01 effect sizes. **Caution**: FOVs from the same slide are NOT independent biological replicates. Report this caveat prominently.
5. **Concordance**: for genes tested in both Mutter_01 (effect sizes) and Mutter_02 (pseudobulk), report direction concordance and effect size correlation. Scatter plot of Mutter_01 log2FC vs Mutter_02 log2FC.

**Outputs:**
- `results/signatures/data/mutter02_signature_scores.tsv`
- `results/signatures/data/mutter02_pseudobulk_degs.tsv`
- `results/signatures/data/concordance_m01_m02.tsv`
- `results/signatures/plots/mutter02_spatial_signatures.png`
- `results/signatures/plots/mutter02_signature_violins.png`
- `results/signatures/plots/concordance_scatter.png`

### Task 6: Stripe model validation on Mutter_02 4h MBRT (if available)

If Mutter_02 contains an MBRT 4h condition with spatial coordinates:

1. Re-fit the stripe model for this slide: conserve ~1.02mm spacing but optimize tilt angle and offset fresh (different tissue mounting angle).
2. Use p21/Cdkn1a as the optimization signal (same as Mutter_01 script 05).
3. If a clear periodic pattern emerges: classify cells as peak/valley using the new fit. Check whether peak/valley DEGs match Mutter_01 direction.
4. If H2AX IHC is available for this slide (check `$SCRATCH/mutter02/` for PPTX or image files): validate against external ground truth.
5. If no clear periodic pattern or if p21 contrast is weak: report as negative finding. Do NOT force-fit.

**Outputs (conditional):**
- `results/signatures/plots/mutter02_stripe_validation.png`
- `results/signatures/data/mutter02_stripe_model.rds`
- `results/signatures/data/mutter02_4h_peak_valley.tsv`

### Task 7: New tumor model exploratory analysis

For Mutter_02 slides 3-4 (new tumor model, separate from 4T1):

1. Basic QC, normalize, cluster, annotate (same pipeline as Task 1).
2. Score with peak/valley signatures.
3. Compare signature scores: new model MBRT vs Control.
4. Report as exploratory — do NOT pool with 4T1 results.

**Outputs:**
- `results/signatures/data/new_model_signature_scores.tsv`
- `results/signatures/plots/new_model_signatures.png`

### Task 8: Write SUMMARY.md

Write `results/signatures/SUMMARY.md` addressing:

1. **Data inventory**: what Mutter_02 actually contained after loading (conditions, cell counts, gene overlap).
2. **Batch assessment**: was integration needed? What was done?
3. **Signature definition**: which genes, how many, rationale for inclusion/exclusion, DDR genes excluded.
4. **Signature persistence (Mutter_01)**: do peak/valley signatures decay, persist, or evolve over time? Spatial striping patterns? Frame correctly per interpretive boundary.
5. **Replication (Mutter_02)**: do signatures score differently in MBRT vs Control? Concordance with Mutter_01?
6. **SBRT distinction**: is SBRT distinguishable from MBRT peak/valley biology? (Descriptive, Mutter_01 only.)
7. **Stripe model validation** (if applicable).
8. **New tumor model** (exploratory).
9. **Answer to scientific intention**: one clear paragraph.
10. **Limitations**: n=1 within Mutter_01, targeted panel blind spots, interpretive boundary for >4h results, FOV pseudobulk is not true biological replication, STING gene list uncertainty.
11. **Recommended next steps**.

Reference all figures and tables by path. Commit everything to git.

## Reviewer Checklist

### Statistical issues
- Pseudobulk: are the independent units truly independent? FOVs from the same slide are NOT independent biological replicates.
- Multiple testing: FDR correction applied for all gene-level tests?
- UCell used as primary scoring method? AddModuleScore comparison reported?
- AddModuleScore: was seed set? Is the background gene set appropriate?
- Concordance analysis: is the correlation inflated by shared noise?
- n=1 per condition within Mutter_01: were any formal p-values reported where only effect sizes are justified?

### Biological issues
- Circularity: were ALL 18 DDR genes excluded from signatures? Check the excluded list against the actual signature genes.
- Cell type composition confound: are cell-type-stratified signatures used as primary (not bulk)? Could apparent signature persistence be driven by changing cell type proportions rather than within-cell-type transcriptional changes?
- STING gene list: was it documented? Is it defensible?
- New tumor model: were results appropriately separated from 4T1?

### Methodological issues
- Batch correction: if Harmony was applied, did it remove biological signal along with batch?
- Stripe model re-fit: for Mutter_02 4h, was tilt re-fit (not transferred from Mutter_01)?
- Signature size: was sensitivity to signature size assessed?
- Spatial analysis: were coordinate systems consistent? Was Moran's I computed for spatial autocorrelation?

### Interpretive issues
- Was the interpretive boundary (scoring vs classification at >4h) maintained throughout?
- Were SBRT comparisons framed as descriptive (Mutter_01 only, n=1)?
- Were spatial striping results at later timepoints appropriately hedged?
- Did SUMMARY.md directly answer the scientific intention?

### Reproducibility
- Random seeds set for all stochastic steps (seed=42)?
- All parameters documented?
- Could results be driven by 1-2 outlier cells or FOVs?

## Reviser Protocol

For each numbered critique in REVIEW.md:
- **Valid and fixable**: Fix it. Re-run analysis, update figures and SUMMARY.md.
- **Valid but unfixable**: Acknowledge in SUMMARY.md under "Limitations" with specific impact assessment.
- **Reviewer is wrong**: Explain why under "Reviewer Response" with evidence.

Add a "Revision Notes" section to SUMMARY.md documenting all changes. Run validation gates after changes. Final commit.

## Storage Rules — MANDATORY

- `/mnt/gcs/` is READ-ONLY via FUSE mount. **NEVER** write to any `/mnt/gcs/` path. **NEVER** use `gsutil` to write.
- `/mnt/rclone/box/` is READ-ONLY. **NEVER** write to Box.
- Copy data FROM mounts to VM local disk. Use `cp`, NOT symlink.
- Write pipeline intermediates to `$SCRATCH` (VM local disk).
- Write all deliverables to `~/repos/spatial-rads/results/signatures/` (git-tracked).
- Do not modify any existing files in `dev/peak_valley_analysis/`.

## Environment

- Conda env: `spatial-rads` (from `config/spatial-rads-conda-env.yaml`)
- Additional packages needed: `r-harmony` (batch integration), `bioconductor-deseq2` (pseudobulk testing), `r-ucell` (rank-based signature scoring for targeted panels)
- All R runs via: `conda run -n spatial-rads Rscript <script>`
