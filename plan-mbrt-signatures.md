# Autonomous Run: MBRT Peak/Valley Signature Persistence and Replication

## v2.0 — use external-atlas typing, not M01 TransferData (2026-05-31)

The "Cell type annotation" step below currently specifies TransferData from a
Mutter_01 reference. That is the v1 scheme and is being **replaced** — see
`plan-processing-pipeline.md` Atlas stem. When this autonomous run executes, type
cells with the **v2.0 external-atlas profile-native** approach (NanoString
CellProfileLibrary mouse mammary atlas + ImmGen + marker malignant, scored
directly), not by transferring M01 labels. This matters because the
**per-cell-type signatures are the PRIMARY readout**: if the cell-type vocabulary
changes, the stratified signature definitions and all per-type claims change too.

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

**Confirmed sample layout (Jenn Fazzari 2026-04-09):** Each of the 4 slides contains 3 samples: 0 Gy control, MBRT 2 days, SBRT 20 Gy 2 days. **All Mutter_02 samples are at 2-day timepoint — no 4h data.**

| Slide | Sample ID | Treatment         |
|-------|-----------|-------------------|
| 1     | 1         | 0 Gy              |
| 1     | 34        | MBRT 2 days       |
| 1     | 28        | SBRT 20 Gy 2 days |
| 2     | 2         | 0 Gy              |
| 2     | 35        | MBRT 2 days       |
| 2     | 30        | SBRT 20 Gy 2 days |
| 3     | 3         | 0 Gy              |
| 3     | 19        | MBRT 2 days       |
| 3     | 16        | SBRT 20 Gy 2 days |
| 4     | 1         | 0 Gy              |
| 4     | 20        | MBRT 2 days       |
| 4     | 18        | SBRT 20 Gy 2 days |

**No FOV-level map**: Jenn confirmed no FOV maps exist for this run. The executor MUST infer sample-to-FOV assignment from the raw Seurat object metadata (Sample ID field, tissue region annotations, or spatial clustering of FOVs on each slide). Save the inferred mapping as `$RESULTS/data/mutter02_fov_sample_map.tsv` and document the inference method in SUMMARY.md.

**H2AX IHC for Mutter_02**: Not available in this run but adjacent tissue sections were cut during original prep. Jenn received them back from Florida 2026-04-09 and is arranging staining at Mayo. Not blocking this run — stripe validation on Mutter_02 must use p21/Cdkn1a as the optimization signal like Mutter_01, but H2AX validation can be added in a follow-up run once staining is complete.

**Important**: These are RAW Seurat objects — no QC, no cell type annotation, no pathway scores. The executor must compute everything from scratch. Pathway scores should be computed using NanoString standard gene modules (see Pathway Gene Lists below).

### NanoString Pathway Gene Lists (RESOLVED)

Source: `$REPO/data/sources/2026-06-22-bruker-mouse-ucc-gene-list.xlsx`, Annotations sheet. These are NanoString's standard panel modules — confirmed as the source of Yi's pre-computed Mutter_01 pathway scores.

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
- [ ] Sample layout TSV (from Jenn Fazzari 2026-04-09 email) saved as `$SCRATCH/mutter02/sample_layout.tsv` (3 columns: slide, sample_id, treatment). FOV-to-sample inference is done by the executor from Seurat metadata — no FOV map exists.
- [ ] Gene list XLSX copied to VM: `$REPO/data/sources/2026-06-22-bruker-mouse-ucc-gene-list.xlsx`
- [ ] Load each Mutter_02 RDS, verify: cell counts match QC email (~999K, 649K, 352K, 793K), spatial coordinates present, gene panel size ~971 genes
- [ ] Verify gene panel overlap with Mutter_01 (expect >95%)
- [ ] Verify raw Seurat object metadata contains Sample ID or equivalent field that can be mapped to the sample layout
- [ ] Disk budget confirmed (>=15 GB free)
- [ ] `~/.claude/settings.local.json` with `{"hooks": {}}` (VM has no Emacs)
- [x] H2AX IHC status confirmed: NOT available for Mutter_02 in this run (adjacent sections being stained at Mayo, Jenn 2026-04-09). Stripe validation on Mutter_02 uses p21 only.

## Executor Tasks

### Task 1: Process Mutter_02 (QC, normalize, annotate)

Process each Mutter_02 slide:

1. Load raw Seurat objects from `$SCRATCH/mutter02/`.
2. **Infer FOV-to-sample assignment from metadata.** No FOV map exists. Inspect each raw Seurat object for a Sample ID field (check fields like `Sample`, `sample_id`, `tissue`, `slide_sample`, or integer codes). Cross-reference with `sample_layout.tsv` (Jenn 2026-04-09) which gives Slide × Sample ID × Treatment for all 12 samples. If no Sample ID field exists, fall back to spatial clustering of FOVs into 3 groups per slide and assign sample IDs by tissue position. Save the inferred map as `$RESULTS/data/mutter02_fov_sample_map.tsv`. Document inference method. All Mutter_02 samples are at the 2-day timepoint.
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

### Task 6: Exploratory stripe detection on Mutter_02 2d MBRT (p21-only, no H2AX)

**Important constraints:** (a) Mutter_02 has no 4h timepoint — all MBRT samples are at 2 days. (b) No H2AX IHC is available for Mutter_02 in this run (adjacent sections being stained at Mayo, deferred to a follow-up). This task is therefore exploratory and interpretively weaker than the Mutter_01 4h stripe model.

For each Mutter_02 MBRT 2d sample (one per slide, 4 total):

1. Extract the cells belonging to that sample using the inferred FOV-to-sample map.
2. Re-fit the stripe model: conserve ~1.02mm spacing (same collimator) but optimize tilt angle and offset fresh per sample (different tissue mounting angle).
3. Use p21/Cdkn1a as the optimization signal (same as Mutter_01 script 05). **p21 is induced early and decays**, so at 2d the contrast will likely be weaker than at 4h. If p21 contrast is not significantly above a random-stripe null, report as "no detectable stripe pattern at 2d" — this is an expected and interpretable negative finding, NOT a failure.
4. If a stripe pattern IS detected at 2d: classify cells as peak/valley using the new fit. Check whether the peak/valley signatures (from Task 3) are enriched in the corresponding classes. Frame as "consistent with persistent spatial memory of original beam geometry" — do NOT claim this validates the 4h classification, because there is no independent ground truth (no H2AX) at 2d.
5. Do NOT force-fit. If 2/4 samples show stripes and 2/4 do not, report honestly.
6. Random-null baseline: permute y-coordinates within each sample and re-run the stripe fit 100 times. Report the observed p21 peak/valley contrast against this null distribution.

**Outputs:**
- `results/signatures/plots/mutter02_stripe_exploration.png` (4-panel, one per MBRT 2d sample)
- `results/signatures/data/mutter02_stripe_fits.tsv` (tilt, offset, p21 contrast, null percentile for each sample)
- `results/signatures/data/mutter02_2d_peak_valley.tsv` (only populated for samples with detected stripes)

### Task 7: Write SUMMARY.md

Write `results/signatures/SUMMARY.md` addressing:

1. **Data inventory**: what Mutter_02 actually contained after loading (samples, cell counts, gene overlap, FOV-to-sample inference method).
2. **Batch assessment**: was integration needed? What was done?
3. **Signature definition**: which genes, how many, rationale for inclusion/exclusion, DDR genes excluded. Cell-type-stratified signatures primary.
4. **Signature persistence (Mutter_01)**: do peak/valley signatures decay, persist, or evolve over time? Spatial striping patterns? Frame correctly per interpretive boundary.
5. **Replication (Mutter_02 at 2d only)**: do signatures score differently in MBRT vs Control vs SBRT at 2 days? Concordance with Mutter_01 2d effect sizes.
6. **SBRT distinction**: is SBRT distinguishable from MBRT peak/valley biology at 2d, both within Mutter_01 (descriptive, n=1) and Mutter_02 (pseudobulk with caveats)?
7. **Mutter_02 stripe exploration (2d, p21-only)**: did p21-based stripe detection find periodic patterns at 2d? Frame honestly — absence of stripes is an informative negative finding consistent with p21 decay between 4h and 2d.
8. **Answer to scientific intention**: one clear paragraph.
9. **Limitations**:
   - n=1 within Mutter_01
   - Targeted ~1000-gene panel blind spots
   - Interpretive boundary for >4h results (scoring, not classification)
   - FOV pseudobulk is not true biological replication
   - STING gene list uncertainty
   - **No H2AX validation for Mutter_02** — stripe fits at 2d rely solely on p21, which has decayed from its 4h peak; negative stripe detection at 2d cannot distinguish "no spatial memory" from "p21 too low for detection". Adjacent-section H2AX IHC is in process at Mayo (Jenn Fazzari, 2026-04-09) and will enable follow-up validation.
   - **No 4h Mutter_02 data** — cannot directly replicate the 4h peak/valley classification; replication is at 2d only
10. **Recommended next steps**: (a) re-run stripe validation on Mutter_02 once H2AX staining is complete, (b) consider non-DDR contrast signals for 2d stripe detection if p21 is too weak, (c) if 2d replication is strong, design a Mutter_03 experiment with 4h + 2d timepoints AND H2AX on Cos slides.

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
- Mutter_02 H2AX absence: was the lack of ground-truth peak/valley validation clearly stated? Were Mutter_02 stripe results explicitly framed as exploratory / pattern-detection only, never as validated peak-vs-valley biology?
- p21 decay at 2d: was the expectation of weaker stripe signal acknowledged? Was a null-distribution baseline (random angle/spacing permutations) used to judge whether any detected stripes exceed chance?

### Methodological issues
- Batch correction: if Harmony was applied, did it remove biological signal along with batch?
- Stripe model re-fit: for Mutter_02 2d, was tilt re-fit per-slide (not transferred from Mutter_01)? Was spacing held fixed at 1.02mm (collimator-conserved)?
- FOV-to-sample inference: was the method for recovering 3-samples-per-slide structure from metadata documented and sanity-checked (expected ~equal FOV counts per sample)?
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
