# M01 peak/valley signature: define, rebuild on current labels, characterize (H2AX-anchored)

> **Spec brainstormed 2026-07-13, adjudicated by a 3-lens independent panel**
> (spatial-omics methods, radiotherapy biology, biostatistics; concordant, not split).
> Splits the retired `plan-mbrt-signatures.md` at the cohort seam: **this plan is M01-only**
> (the H2AX-anchored 4h signature and its within-M01 characterization). All M02 replication and
> any future stripe re-detection live in a **separate companion plan**
> (`plan-m02-signature-replication.md`, to be written), because M02 is the replication/inference
> cohort and folds into the differential layer. Spec location follows project convention
> (`plan-*.md` at root).

## What this plan is

Define the MBRT peak-versus-valley transcriptional signature at the one timepoint with ground truth
(M01 4h flank block, Block_21, gamma-H2AX-marked beam tracks), **rebuilt on the current QC'd object
and merged-scale atlas labels**, and characterize it within M01. The peak/valley contrast is the
only internally-controlled comparison in the project: within-slide, within-animal, so it does not
inherit the cross-animal confounding that limits everything else. That is exactly why it is worth
doing rigorously and worth freezing as a reusable asset.

This supersedes the retired signature-projection plan, whose core defects the panel identified: it
was built on now-retired vendor ImmGen labels, treated an n=1 cross-animal series as a kinetic
trajectory, and proposed de-novo stripe re-detection that aliases with the CosMx FOV acquisition
grid. This plan keeps only what survives that critique.

## The question (correctly bounded)

At 4h, does the beam-peak transcriptional state differ from the valley state, and what programs
define that difference? This is answerable with ground truth **only here**. Persistence over time,
spatial memory, and any peak/valley identity off this block are **not** in scope: they are gated on
staining that M01's later blocks and M02 do not have, and belong to the companion plan as a gated
future test.

## Panel triage carried into this plan

- **4h peak/valley signature: DERIVE fresh, tumor-restricted, on honest replication.** There is no
  trustworthy prior signature to reproduce: the old per-cell Wilcoxon lists
  (`peak_signature_bulk.tsv` etc.) were pseudoreplicated (149k cells treated as independent
  replicates), are DELETED, and are not a concordance target. Derive the signature de novo on the
  current QC'd 4h cells with H2AX zone labels, **within the tumor compartment** (peak/valley is a
  within-tumor-tissue contrast; comparing across a whole slide confounds zone with cell-type mix),
  using FOV-level pseudobulk so the replicate unit is the FOV, not the cell. Whatever clears (or does
  not clear) at FOV-level FDR is the signature; report it honestly with no external yardstick.
- **Cell-type-stratified signature: RETIRE the artifact, REBUILD at coarse resolution.** The old
  per-ImmGen-subtype lists (`celltype_signatures.tsv`) are dead: retired labels, known-bad for
  immune in flank, and fine-subtype-by-zone splitting fragments one n=1 block into uninterpretable
  bins. The composition-versus-intrinsic question it was meant to answer is now *more* central (the
  robustness pass showed day-2 biology is composition shift), so rebuild it **only** at the current
  coarse tumor/stroma/immune compartments, never fine subtypes, never the retired field.
- **Cross-timepoint module scoring on M01: keep the scoreboard, RETIRE the kinetics framing.**
  Recompute the module scores across the M01 conditions as a **descriptive per-condition
  scoreboard**. Every M01 timepoint is a different mouse and tumor block (n=1), so time is aliased
  with animal: the scores may be shown across timepoints as **cross-animal snapshots**, but never
  connected as a decay/persistence trajectory or given a temporal-causal reading.
- **De-novo stripe re-detection: RETIRE** (out of scope here; gated to the companion plan when
  later-block H2AX exists). ~1 mm beam periodicity aliases with the FOV acquisition pitch; the
  y-shuffle null cannot exclude it; p21 is circular (defines and scores the stripe) and decayed.
- **The "persistence" claim: RETIRE as a result** (positive-only, sunk by null-ambiguity); preserve
  as the gated question for a future 4h-plus-2d-plus-H2AX experiment.

## Analysis design: nested grains on FOV pseudobulk

The signature lives at the **gene** grain; the pathway tables are built on top of it, not parallel
to it. All grains run on the same within-block peak-versus-valley contrast. A fourth,
panel-independent readout (morphology) corroborates the transcriptional result through segmentation
geometry.

### Step 0: zone-differential QC attrition audit (the rebuild-and-verify crux)

The panel's central worry is that the current QC preferentially drops peak cells, and every
count- or quality-based cutoff is mechanistically anti-correlated with radiation damage: heavily
irradiated 4h peak cells in mitotic catastrophe / early apoptosis show transcriptional shutdown, RNA
degradation, and abnormal nuclei, so a counts floor, a complexity floor, an area filter, or a
negprobe-proportion cap each fail those cells preferentially. Damage is what defines a peak, so QC
and signal are coupled, not independent. Before rebuilding the signature, audit this directly on the
4h MBRT block (`sam0003`, Block_21):

- Join the raw per-cell QC metrics (`nCount_RNA`, `nFeature_RNA`, `propNegative`, `Area`) to the H2AX
  zone labels via the `<slide>_<fov>_<local>` cell-id crosswalk.
- Report per-zone retention at each of the four current criteria (counts > 20, features > 10,
  propNegative < 0.5, area 3-MAD) and jointly; expose whether peak attrition exceeds valley.
- Report the peak/valley effect-size for the signature genes **before versus after** QC, so
  attenuation is measured, not assumed.
- Deciding tension: loosening QC to retain damaged peak cells also readmits genuine technical junk,
  and the two look alike (both low-count, high-negprobe). Do not silently loosen; if attrition is
  strongly zone-biased, report a sensitivity analysis at a relaxed floor alongside the default, and
  flag which signature genes depend on the retained low-count cells.

This step gates interpretation of the whole signature: if peak attrition is modest and the
effect-size is stable across QC, the rebuild is trustworthy; if not, the signature is partly a
survival-biased artifact and must be reported as such.

**Result (run 2026-07-13 on beast, sam0003 raw object, 149,111 zoned cells).** The feared
zone-selective attrition is empirically absent. Applying the current four-criteria gate
(counts > 20, features > 10, propNegative < 0.5, area 3-MAD) to the zoned cells drops only 0.06% of
peak and 0.05% of valley cells (a 0.01-point gap); zero cells fail the counts, features, or negprobe
floors, and only the area-MAD filter removes anything. The mechanistic premise also fails: peak and
valley cells have matched library depth (median 84 counts / 70 features versus 84 / 69), so 4h peak
cells are not depth-depressed. Scope: zone labels exist only for cells that passed the original QC,
so this measures whether the current gate perturbs the exact population the signature was built on
(it does not, ~99.94% invariant, not zone-biased); it does not speak to cells dropped before zoning,
which were never in the signature. Conclusion: the derivation runs on essentially the original cell
set, so the fresh tumor-restricted FOV-level DE is not survival-biased; it stands on its own FDR, not
on agreement with any prior list.

### Prerequisite gate: FOV footprint versus beam spacing

Before choosing the design, measure the CosMx FOV footprint (from `x_slide_mm`/`y_slide_mm` extents
per `fov`) against the 1.02 mm collimator spacing. This decides the unit:

**Result (run 2026-07-13):** CosMx FOVs are 0.506 mm square (median long-axis 0.507, tight IQR),
half the 1.02 mm beam period; the unaligned FOV grid means 74/80 FOVs (92%) contain both peak and
valley cells, only 6 single-zone. Zones near-balanced (73,558 peak / 75,553 valley). **Paired
`~ FOV + zone` design adopted** (well-powered regime).

- **FOVs smaller than the spacing and straddling beam boundaries (expected, ~0.5-0.9 mm CosMx):**
  use the **paired FOV-by-zone design** below. The beam stripes cut across FOVs at the mounting
  tilt, so a typical FOV contains both peak and valley cells, giving within-FOV zone contrast to
  block on. The design's validity does **not** assume FOV/beam alignment; it works *because* the two
  grids are unaligned. Its power scales with how many FOVs straddle a boundary.
- **FOVs larger than the spacing or sitting entirely within one stripe:** those FOVs hold one zone
  and drop out of the paired contrast; fall back to the unpaired FOV-mean design (each FOV assigned
  to its majority zone). Report how many FOVs this discards.

### 1. Per-gene DE (grain = gene): the signature

FOV-pseudobulk counts, cells zoned per-cell by the H2AX stripe geometry. Preferred model
`~ FOV + zone` (paired: pseudobulk within FOV separately for peak-cells and valley-cells, block on
FOV); this absorbs the spatial-region confound (peak versus valley FOVs sitting in different tumor
territory) that an unpaired FOV-mean cannot. limma-voom or DESeq2 on the ~14 peak / ~13 valley
pseudobulk units (unpaired) or the paired FOV-by-zone units.

Output columns: `gene, mean_peak, mean_valley, log2FC, se, stat, df, p, padj, pct_expr_peak,
pct_expr_valley, direction`. Thresholding yields the signature gene list; the ranked statistic feeds
camera/GSEA.

**Guard, DDR-gene exclusion (circularity).** The stripe geometry was fit by maximizing p21/DDR
contrast, so DDR genes are guaranteed to separate peak/valley by construction. Exclude the DDR gene
set (the existing 18-gene exclusion list) from the reported signature; keep it documented as the
fitting basis, not a finding.

### 2. Gene-set enrichment (grain = gene set): the module readout

Two complementary tests on the same contrast, kept together because each catches what the other
misses:
- **camera** (preferred over a FOV-mean t-test): controls for inter-gene correlation; row grain one
  per gene set. Output: `pathway, n_genes, direction, camera_stat, p, padj`.
- **FOV-mean UCell** module scores: a set can move via a coherent sub-threshold shift the per-gene
  FDR misses, and vice-versa. Output: `pathway, n_peak_fov, n_valley_fov, mean_peak, mean_valley,
  delta, ci_low, ci_high, stat, df, p, padj, direction`.

The FOV-level score matrix (`fov, zone, n_cells, compartment, pathway, fov_mean_score`) is the tidy
intermediate the UCell test consumes; carry `n_cells` to weight or filter thin FOVs.

### 3. Optional coarse-compartment stratification (grain = gene-by-compartment)

Same contrast, re-run within each current coarse compartment (tumor / stroma / immune from the
atlas), to separate a tumor-intrinsic peak/valley shift from an immune-driven one. Leading
`compartment` key on the grain-2 output. **Coarse only:** not fine subtypes, not the retired field.
At n=1, expect only tumor / fibroblast / macrophage-abundant compartments to be interpretable;
report cell counts per compartment-by-zone and drop thin bins.

### 4. Morphology (grain = cell): an orthogonal, panel-independent damage phenotype

Segmentation-derived shape metrics (cell area, nuclear area, circularity, solidity, aspect ratio;
present in the CosMx metadata parquet) contrast peak versus valley at 4h as a physical readout that
does not touch panel sparsity, corroborating or contradicting the transcriptional signature through
an independent channel. The completed dev script `17_morphology_analysis.R` already computes these,
but must be re-fit to this plan's discipline: within-block peak/valley on the current zone labels,
coarse atlas compartments (not the retired ImmGen `a`-to-tumor relabel), and per-condition snapshots
(not the kinetic line plots it currently draws).

**Guard, segmentation-artifact check.** Morphology comes from segmentation, and segmentation quality
can itself vary by zone (dense or damaged peak tissue segments differently than valley). Before any
peak/valley shape difference is called biology, check it against a segmentation-quality confound
(per-zone QC-flag rate, transcript-density or count covariate); report a shape shift only if it
survives that adjustment.

## Preserve-set (freeze as reusable assets)

Regardless of outcome, freeze with provenance so the gated persistence question reopens cleanly if
later-block or M02 H2AX ever arrives:
1. **The 4h H2AX per-cell peak/valley labels** (`dev/peak_valley_analysis/data/mbrt4h_peak_valley.tsv`):
   the project's only ground-truth spatial dose-zone partition, irreplaceable.
2. **The re-derived, concordance-checked, DDR-excluded signature gene list**, registered as a
   first-class named module alongside `config/pathway_gene_lists.yaml`, with provenance (QC object
   version, cell set, thresholds).
3. **The DDR-exclusion list and the fitted stripe geometry** (1.02 mm spacing, mounting tilt/offset).
4. **The gated question statement plus the M03 design requirement** (4h plus 2d plus H2AX on the
   same cohort) that a true persistence/spatial-memory test demands.

## Interpretive boundary (non-negotiable)

- Within-block peak/valley at 4h is causal and internally controlled; **off this block there is no
  peak/valley geometry**, so score, never classify.
- The M01 cross-condition scoreboard is cross-animal snapshots, **not** kinetics; no decay/persistence
  trajectory.
- No dose or "sparing" language anywhere (MBRT mean dose unrecorded).
- Results are n=1 descriptive for the effect itself; the inferential test of whether the program
  replicates is the companion M02 plan's job, and even there is positive-only (null-ambiguity: a
  diluted whole-compartment null cannot distinguish loss from a still-active peak-restricted
  program).

## Not in this plan (explicitly deferred)

- **M02 day-2 replication:** sample-level (not FOV) positive-only UCell readout folded into
  `aggregate_differential.smk`; companion plan.
- **Stripe re-detection at later timepoints / M02:** gated on adjacent-section H2AX; companion plan.
- **VM/autonomous execution:** the retired plan's failed autonomous run is a compute post-mortem
  (`retro-mbrt-signatures.md`); this analysis runs locally against the current object, not a merged
  monolith on a VM.
