# MBRT vs SBRT (vs Control) downstream differential analysis — design & plan

The cross-dataset differential analysis that consumes the unified per-cell
label table (`results/aggregate/full_labels.parquet`, 3,277,090 cells, M01+M02
jointly typed). Goal: characterize how microbeam radiation therapy (**MBRT** —
dose delivered as an array of narrow ~1 mm-spaced parallel beams, producing
alternating high-dose "peak" and low-dose "valley" tissue) differs from
stereotactic body radiation therapy (**SBRT** — here, a uniform broad-beam
comparator delivering the same prescription everywhere) and from unirradiated
Control, in the 4T1 mouse mammary flank tumor and its surrounding normal
tissue.

Brainstormed and full-panel-reviewed 2026-06-02. Companion to
`plan-aggregate.md` (the upstream cell-typing this consumes) and
`plan-processing-pipeline.md` (the per-sample pipeline beneath that).

> **FRAMING SUPERSEDED 2026-06-22 by `plan-mbrt-vs-sbrt-reanalysis.md`, then its PREMISE RETRACTED
> 2026-07-13.** The reanalysis reframed the "MBRT essentially null" reading below as spatial
> *dilution* of a real M01 4h peak signature (37-gene MHC-I/IFN/damage from the rotation null). That
> peak signature was subsequently RESOLVED to an artifact (`plan-rotation-null-reconcile.md`, Steps
> 0-2): it does not survive current atlas tumor labels + QC or FOV-level replication, and a phase
> null is geometrically futile here. So there is no established peak signal to have been diluted. The
> still-valid part: whole-compartment averaging is the wrong scale to detect any peak-restricted
> effect, so the arm-level MBRT contrasts **bound, they do not resolve**. Read the arm-level MBRT
> result as **underpowered / unresolved at whole-compartment scale**, and a dose-resolved 4h MBRT
> effect as an **open, currently-unsupported** question. This doc remains the run-record of the
> executed layer.
>
> **STATUS — EXECUTED (2026-06-03).** The full layer ran and is committed on
> `main` (the 15-task build; commits `ce660ca`..`b8be903`. The spent
> implementation checklist `plan-mbrt-vs-sbrt-impl.md` was removed once executed —
> recoverable from git history). Single canonical output `results/aggregate/results_master.tsv`
> (32,961 tier-tagged rows; 621 confirmatory / 32,340 exploratory), produced by
> `workflows/aggregate_differential.smk`. **Realized result:** the dominant day-2 signal was
> SBRT-driven stromal fibrosis (collagens + Acta2; H3 = 20 of 27 SBRT_vs_Ctrl
> confirmatory hits, padj_confirmatory to ~2e-11); MBRT was essentially null at
> whole-compartment scale (MBRT_vs_Ctrl 2/207 confirmatory hits, all 207 below
> the n=4 MDE); MBRT_vs_SBRT confirmed SBRT > MBRT for the same stromal genes
> (19/20 hits H3). Critically, 578/621 (93%) confirmatory rows fell below their
> n=4 MDE, so most MBRT non-results are underpowered, not established nulls.
> See the Limitations section below for realized magnitudes and the
> dose-mismatch + spatial-dilution caveats.
>
> **QC CONFOUND CHECK (2026-07-10).** The day-2 `fraction_shift`
> reading (more cells expressing collagen/Acta2) was checked against SpatialQM per-sample
> technical metrics compared across arms: 5/6 balanced (SBRT neither more sensitive nor more
> contaminated than Control; the FDR flag immaterial), replicate concordance 0.94–0.97 with 0
> outlier slides. The composition shift is **not** a sensitivity/contamination/outlier artifact.

## Governing principle: look first

This is an initial, agnostic spatial-omics survey, not a confirmatory test of a
single pre-chosen mechanism. We test literature-grounded a priori hypotheses
**and** scan the data for new structure. The manuscript leads with a biological
idea only if the data surface a compelling one — we do not assume in advance
that a given mechanism (or any mechanism) is present or absent. The known prior
from this project's own earlier per-sample work is weak-to-null (whole-sample
pseudobulk MBRT-vs-Control collapsed toward zero, immune-activation markers
silent — `project_mbrt_mechanism_status.md`), but that prior was generated
before any per-cell-type or modality-direct (MBRT-vs-SBRT) contrast had been
run, so it constrains expectations, not the decision to look.

## Scope

**In scope.** Whole-tumor and whole-normal-tissue differences between the three
arms (MBRT, SBRT, Control), resolved by cell type, at the merged cross-dataset
scale. Composition shifts, cell-type-resolved differential expression, pathway/
program activity, and tissue spatial organization (cellular niches, tumor-immune
mixing).

**Out of scope (deferred to a separate reanalysis).** The within-tumor
peak-versus-valley contrast — comparing the high-dose beam tracks against the
low-dose gaps inside a single MBRT sample. That analysis needs the FOV-
extrapolation peak/valley assignment (`feedback_pv_method.md`) and the H2AX
ground truth, and is the natural follow-up once the arm-level picture is in
hand. All `dev/peak_valley_analysis/` and `dev/03*_classify.R` work belongs to
that deferred track.

## Experimental design and the two cohorts

Two CosMx datasets, jointly typed, analyzed with cohort-appropriate statistics:

**Mutter_02 (M02) — the inference cohort.** 4 slides (sld0005–0008), each
carrying one Control + one MBRT + one SBRT region, all harvested at 48 h (day
2). This is a **randomized complete block design**: slide is the block, arm is
the treatment, n = 4 per arm. Because each slide contains all three arms, the
slide (batch) effect is removed in the contrast. This cohort carries the formal
statistics.

**Mutter_01 (M01) — the descriptive time-course cohort.** n = 1 per
arm×timepoint: Control at t = 0; MBRT at 1 h / 4 h / 48 h / 144 h; SBRT at 4 h /
48 h / 144 h (flank model; the tongue samples sam0009–0011 are a separate model,
excluded here). With n = 1 there is no within-arm variance, so M01 yields
**effect sizes and kinetic trajectories only — no p-values** (per the project
stats convention and `feedback_pv_method.md`). M01's role is to show whether an
M02 day-2 finding sits on a coherent time course and to provide a qualitative
cross-cohort concordance check.

**The dose-mismatch confound (open input, documented, not a blocker).** In the
metadata, M02 SBRT = 20 Gy uniform; M02 Control = 0 Gy; M02 MBRT dose is blank
because a microbeam array has no single dose number (peak and valley doses
differ by an order of magnitude). The mean/integral dose of the MBRT arm is not
recorded in the data model and must come from the beam physics (Fazzari/Mutter).
Until then, an MBRT-vs-SBRT difference cannot be cleanly attributed to dose
*pattern* versus dose *magnitude*. This is recorded as an open lab input and a
stated limitation, not a reason to defer the analysis.

## A priori hypotheses (literature-grounded)

Three hypotheses with established support in the spatially-fractionated /
microbeam radiotherapy literature. Each names the readout and the expected
direction. Citations were verified against Scopus (peer-reviewed, by
citation count) on 2026-06-02.

### H1 — Immune activation, greater under MBRT than SBRT

Spatially fractionated and microbeam radiotherapy have been reported to drive
stronger immune engagement (interferon signaling, cytokine induction,
T-cell infiltration, abscopal responses) than uniform delivery at matched dose.

- **Readout:** Type I interferon program (Ifit/Isg/Mx/Oas family), Type II
  interferon program (Gbp/Irf1/Cxcl9-10-11/MHC-II), cGAS-STING; immune-cell
  compartment fraction and infiltration.
- **Expected direction:** MBRT > SBRT > Control.
- **Support:** Kanagavelu 2014, *Radiation Research* (doi:10.1667/RR3819.1);
  Billena & Khan 2019, *IJROBP* (doi:10.1016/j.ijrobp.2019.01.073); Yan 2020,
  *Clin Transl Radiat Oncol* (doi:10.1016/j.ctro.2019.10.004).

### H2 — Vascular/endothelial remodeling and tumor oxygenation

Microbeam radiotherapy preferentially damages immature tumor vasculature while
sparing mature normal vessels; downstream this shifts tumor oxygenation
(hypoxia/reoxygenation is the readout one layer below the vascular change, not a
separate mechanism).

- **Readout:** endothelial-cell fraction and state; angiogenesis/vascular-
  maturation program; hypoxia program (Hif targets, Vegfa, Car9, Slc2a1) as the
  oxygenation proxy.
- **Expected direction:** tumor-vessel disruption + hypoxia signature greater
  under MBRT; normal-vessel sparing (less endothelial perturbation in
  normal-tissue regions) greater under MBRT than SBRT.
- **Support:** Bouchet 2010, *IJROBP* (doi:10.1016/j.ijrobp.2010.06.021);
  Bouchet 2012, *Radiation Research* (doi:10.1667/RR2784.1, vascular
  architecture + oxygenation); Bouchet 2013, *Radiother Oncol*
  (doi:10.1016/j.radonc.2013.05.013, MRT induces hypoxia in gliosarcoma but not
  normal brain); Sabatasso 2015, *Physica Medica*
  (doi:10.1016/j.ejmp.2015.04.014); Serduc 2006, *IJROBP*
  (doi:10.1016/j.ijrobp.2005.11.047); Slatkin 1992, *Med Phys*
  (doi:10.1118/1.596771).

### H3 — Normal-tissue (stroma) sparing

The defining preclinical claim for MBRT is that normal tissue tolerates it
better than uniform radiation at matched dose: less fibrosis, less stromal
stress/damage response in non-tumor regions.

- **Readout:** stromal compartment differential expression in normal-tissue
  regions; fibrosis/remodeling program (collagens, Acta2, Tgfb targets);
  stress/senescence and damage-response programs in fibroblasts.
- **Expected direction:** less stromal damage/fibrotic response under MBRT than
  SBRT in normal tissue.
- **Support:** Dilmanian 2006, *PNAS* (doi:10.1073/pnas.0603567103); Crosbie
  2010, *IJROBP* (doi:10.1016/j.ijrobp.2010.01.035, tumor vs normal differ);
  Bazyar 2018, *Sci Rep* (doi:10.1038/s41598-018-30543-1, MRT vs conventional
  toxicity).

### Moved to exploratory: DNA-damage / cell-cycle arrest / senescence

A DNA-damage-and-arrest hypothesis was considered as a fourth a priori line but
its spatial-fractionation literature is dominated by within-sample peak/valley
and bystander effects (the deferred track), not arm-level MBRT-vs-SBRT
differences. Scopus did not surface a compelling arm-level a priori basis, so
per the 2026-06-02 decision it is analyzed in the **exploratory** phase (the
DDR program is already scored), not pre-registered as confirmatory. No fourth
standalone a priori hypothesis is added (hypoxia folds into H2).

## Statistical design and principles

**Confirmatory family (M02, day 2, n = 4/arm).** Pre-registered: the H1/H2/H3
readouts above, tested with the three arm contrasts. This family carries its own
multiple-testing correction, kept separate from the exploratory grid.

**Three contrasts, every confirmatory readout:**
1. MBRT vs Control — does microbeam change the tissue at all?
2. SBRT vs Control — does uniform radiation change the tissue at all?
3. **MBRT vs SBRT — the modality-specific test.** This direct contrast is the
   analysis that isolates dose *pattern*. The "shared-vs-modality-specific"
   scatter of the two vs-Control contrasts is illustration only: both axes share
   the Control arm, which manufactures apparent agreement, so it is never the
   inferential test.

**Methods (all field-standard, established tools):**
- **Composition:** propeller (speckle package) on cell-type proportions,
  block on slide. Reports proportion shifts per cell type per contrast.
- **Cell-type-resolved differential expression:** pseudobulk per
  sample×cell-type, DESeq2 with design `~ slide_id + condition`, apeglm
  log-fold-change shrinkage, the three contrasts above. Pseudobulk (not
  per-cell) is the established way to get calibrated p-values with n = 4
  biological replicates.
- **Program / pathway activity:** UCell (rank-based, primary) + AddModuleScore
  (seed 42, secondary) per cell, summarized to per-sample×cell-type means, then
  limma on those means with the slide block. Inference is on sample-level means,
  not cells, so n = 4 is the real n.

**Statistical principles (from the review-panel consensus):**
- **Tiered FDR.** The confirmatory H1/H2/H3 family is corrected independently of
  the exploratory genome-/program-wide scans. Findings are labeled by tier.
- **Power / minimum detectable effect, precomputed.** At n = 4/arm in a blocked
  design there are 6 residual degrees of freedom. Before interpreting any null,
  compute and report the minimum detectable effect size (composition shift,
  log2FC, program-score delta) at 80% power. A null is only meaningful against a
  stated detectable effect.
- **Effect sizes + 95% CIs always**, alongside (not instead of) p-values, for
  every contrast — and effect sizes only for the n = 1 M01 cohort.
- **Confounds modeled or stated.** Tumor size / total cell count differs across
  arms and inflates apparent composition differences; account for it
  (proportions, and cell-count as a covariate where appropriate). Necrotic-zone
  cells are excluded before aggregation, following the convention already used in
  `dev/mbrt_vs_sbrt/01/02/05` (low-density necrosis gated out at day 2/6).
- **Panel-coverage verification.** Every gene in every readout/program is
  checked for presence and adequate detection in the common 950-gene panel
  before a result rests on it. Genes absent or near-zero are dropped and noted;
  a "null" on a gene the panel cannot see is not a biological null. Missing
  programs (angiogenesis/vascular-maturation, hypoxia, fibrosis/remodeling,
  stromal stress/senescence) are **built gene-by-gene against the panel** —
  `config/pathway_gene_lists.yaml` currently holds only TypeI_interferon,
  TypeII_interferon, DNA_Damage_Repair, STING.
- **Broaden thin readouts.** The earlier stromal readout leaned on a single
  marker (p21/Cdkn1a); broaden to multi-gene stromal programs. Pool endothelial
  cells across samples to lift sparse immature-vessel markers above the
  detection floor before testing H2's vascular readout.

## Confirmatory analyses (the pre-registered plan)

For each of H1, H2, H3, run all three contrasts on the inference cohort (M02,
day 2):
1. **Composition** — propeller, per cell type, slide-blocked.
2. **Cell-type-resolved DE** — pseudobulk DESeq2 within the compartment the
   hypothesis names (immune for H1; endothelial/stroma for H2; stroma in
   normal-tissue regions for H3).
3. **Program scoring** — UCell + AddModuleScore on the hypothesis's gene
   set(s), limma on per-sample×cell-type means.

Report each as: effect size + 95% CI + tiered-FDR p-value, with the precomputed
minimum detectable effect for that readout, and the panel-coverage status of its
genes.

## Exploratory analyses (hypothesis generation — the point of an omics survey)

- **Unbiased cell-type-resolved DE** — no gene list; full pseudobulk DESeq2 per
  cell type, all three contrasts, exploratory-tier FDR.
- **Broad program scan** — the ~50 MSigDB mouse Hallmark sets (extending the
  6-set UCell scan already in `dev/mbrt_vs_sbrt/02_pathway_kinetics.R`), scored
  per cell type, screened for arm differences.
- **DNA-damage / arrest / senescence** — the DDR program (already scored) plus
  senescence, examined here rather than as confirmatory.
- **Spatial organization — promoted to a first-class, early track.** Cellular
  niches via per-slide RANN k-nearest-neighbor cell-type composition → k-means
  niches (the proven lighter alternative to BANKSY already built in
  `dev/mbrt_vs_sbrt/07_niche_clustering.R`; BANKSY remains the heavier option if
  niches prove central). Plus the Keren 2018 tumor-immune mixing score
  (`dev/.../12_mixing_score.R`) and the k = 20 immune-neighbor fraction
  (`dev/.../05_spatial_nn.R`). Spatial structure is where a microbeam dose
  pattern is most likely to leave a signature that whole-sample pseudobulk
  averages away, so it runs early, not last.
- **Myeloid M1/M2 polarization** — the prior-work hypothesis that MBRT day-2
  myeloid cells skew M2 (Cd163/Mrc1 up), from
  `dev/mbrt_vs_sbrt/11_m1_m2_polarization.R`.

## Cross-cutting analyses

- **MBRT-vs-SBRT direct contrast** — run for every readout above; it is the
  modality test, surfaced explicitly in every results table.
- **M01 descriptive time course** — for any M02 day-2 finding, plot the M01
  1 h→144 h trajectory (effect sizes only) to show whether it sits on a coherent
  kinetic curve. Reuses `dev/mbrt_vs_sbrt/01/02/03/05`.
- **Cross-cohort concordance gate (qualitative, not validation).** M01-vs-M02
  day-2 effect-size scatter (`dev/mbrt_vs_sbrt/08_set2_validation.R`,
  Spearman). Agreement raises confidence; it is a qualitative gate, not formal
  external validation (different n, different cohort).

## Reuse from `dev/` (do not rebuild)

Already built, lift into the pipeline:
- Composition trajectories — `dev/mbrt_vs_sbrt/03_composition_trajectories.R`
- Pathway kinetics (UCell IFN/STING/DDR/hypoxia/senescence; extend 6 → ~50
  Hallmark) — `dev/mbrt_vs_sbrt/02_pathway_kinetics.R`
- M01 time course — `dev/mbrt_vs_sbrt/01,02,03,05`
- M01↔M02 concordance scatter — `dev/mbrt_vs_sbrt/08_set2_validation.R`
- Per-slide kNN → k-means niches — `dev/mbrt_vs_sbrt/07_niche_clustering.R`
- Keren tumor-immune mixing — `dev/mbrt_vs_sbrt/12_mixing_score.R`
- k = 20 immune-neighbor fraction — `dev/mbrt_vs_sbrt/05_spatial_nn.R`
- M1/M2 polarization — `dev/mbrt_vs_sbrt/11_m1_m2_polarization.R`
- Pseudobulk DE — `scripts/aggregate/deg_pseudobulk.R` already exists (output
  `deg_summary_m02day2.tsv`)
- **Cell-type-label QC** — adapt `dev/.../03_cell_type_validation.R` marker
  DotPlot as a first step, to confirm the new unified labels behave before any
  contrast rests on them.

Out of scope (peak/valley deferred track): all `dev/peak_valley_analysis/05-19`,
all `dev/03*_classify.R`, all `dev/if_channel_check/`.

## Pipeline prerequisites (fix before running)

These are defects in the aggregate scripts/workflow that block correct
execution; fix them before the first confirmatory run:
- **Pseudobulk label re-point.** `scripts/aggregate/pseudobulk_build.R:23` reads
  `md$cell_type_atlas` — the rejected per-cell InSituType labels (84.6%
  mistyped tumor). Re-point to the unified `full_labels.parquet`.
  (`deg_pseudobulk.R` itself is fine; it only consumes the SummarizedExperiment.)
- **Retain spatial coordinates.** `scripts/aggregate/merge.R:100-101` drops
  `x_slide_mm`/`y_slide_mm`, so `merged.rds` has no coordinates and the spatial/
  niche track has nothing to run on. Add coords back (cheap merge re-run); they
  exist in the 20 per-sample scored.rds.
- **Snakemake arg desync.** `aggregate_differential.smk` (dated 2026-06-01, before typing was
  locked) passes wrong/off-by-one args to the rewritten scripts — the
  `pathway_summary` rule (smk:271-274) and `composition` rule (smk:192,204) pass
  a TSV where a parquet is expected. Reconcile rule args with the current script
  signatures before running via snakemake. (The currently-running pathway job was
  launched standalone with correct args, so it is unaffected; the desync only
  bites the snakemake path.) BLAS hygiene in the workflow is already correct
  (shell.prefix BLAS-pin + threads:1).

## Limitations (the spatial-dilution and power magnitudes written after looking)

- **Spatial dilution of the MBRT signal (realized magnitude).** Whole-compartment
  aggregation averages MBRT peak and valley together, and the run bears this out
  starkly. In the confirmatory family (207 rows per contrast), `MBRT_vs_Ctrl`
  produced 2 hits at `padj_confirmatory < 0.05` with a median |effect| of 0 and a
  max |effect| of 1.23, whereas `SBRT_vs_Ctrl` produced 27 hits (median |effect|
  0.057, max 1.63) and the modality-direct `MBRT_vs_SBRT` produced 20 hits, 19 of
  them H3 stromal genes all in the SBRT > MBRT direction. At whole-compartment
  scale MBRT therefore shows essentially no day-2 transcriptional response while
  uniform SBRT does. This reproduces the earlier whole-sample pseudobulk pattern
  (`project_mbrt_mechanism_status.md`). Whether MBRT's near-zero arm-level effect
  is genuine tissue sparing or a peak-restricted signal cancelled by averaging
  across the ~50% valley volume cannot be separated here; the per-cell-type and
  spatial-niche tracks did not rescue an MBRT effect, and the full mitigation
  remains the deferred peak/valley reanalysis (needs the FOV-condition map).
- **Dose mismatch (still open).** SBRT = 20 Gy uniform; the M02 MBRT mean/integral
  dose is unrecorded (metadata field blank — a microbeam array has no single dose
  number). The stronger SBRT response may therefore reflect higher mean dose, not
  uniform *pattern*, and the `MBRT_vs_SBRT` contrast cannot cleanly attribute its
  20 stromal hits to dose pattern until the beam physics arrive from Fazzari/Mutter.
  Asked of the lab; not returned as of 2026-06-02.
- **Power at n = 4 (realized).** 6 residual df. 578 / 621 (93%) of confirmatory
  rows fell below their minimum detectable effect; by contrast 207/207 (100%) for
  `MBRT_vs_Ctrl`, 185/207 (89%) for `SBRT_vs_Ctrl`, 186/207 (90%) for
  `MBRT_vs_SBRT`. The only confirmatory signal that clears its MDE is SBRT-driven
  H3 stromal activation (collagens + `Acta2`) plus a few endothelial/ISG rows.
  Every MBRT confirmatory non-result is underpowered (100% below MDE), so the
  absence of an MBRT immune-activation signal is NOT established as a true null —
  it is jointly underpowered and spatially diluted.
- **Panel blind spots.** ~950 genes; a negative on an off-panel or near-zero
  gene is not a biological negative. The four new stromal/vascular sets carry
  uneven panel coverage (`results/aggregate/gene_set_panel_coverage.tsv`):
  Angiogenesis 14/19, Fibrosis_remodeling 13/18, Stromal_stress_senescence
  8/12, but **Hypoxia only 5/13** (Hif1a, Vegfa, Slc2a1, Ldha, Ndrg1; Car9,
  Pgk1, Bnip3, Pdk1, Eno1, Aldoa, Hk2, Egln3 all off-panel). The Hypoxia score
  rests on five genes, so a null on it is panel-limited, not biological.
- **M01 is n = 1.** Descriptive only; its time course supports, never confirms.

## Status and next steps

- **Now running (standalone):** `scripts/aggregate/pathway_summary.R` over the
  merged object with the unified labels (8 UCell workers active as of
  2026-06-02 ~14:40). On completion: read `pathway_test_m02day2.tsv` and view
  the three PNGs (m02 heatmap, m01 timecourse, ucell/ams scatter).
- **Before the confirmatory run:** apply the three pipeline prerequisite fixes;
  build the missing gene sets (angiogenesis/vascular-maturation, hypoxia,
  fibrosis/remodeling, stromal stress/senescence) gene-by-gene against the
  panel; precompute the power / minimum-detectable-effect table.
- **Then:** run the confirmatory H1/H2/H3 plan (composition + cell-type DE +
  program scoring × three contrasts), the exploratory scans, and the spatial
  track; assemble the single narrative under "look first."
