# Comparison registry + marker provenance in the data-model stage — design

> **Brainstormed 2026-07-12.** Shaped by a structural survey of the 14 `plan-*.md`
> docs, a survey of how contrasts are currently defined in code, and an independent
> spatial-methods expert review of the comparison set and practices (verdict: comparisons
> and statistics are sound; three technical guards folded in as method notes, not new
> scope). Spec location per the superpowers convention; a `plan-comparison-registry.md`
> at project root follows from writing-plans.

## Job to be done

Make the study's **experimental comparison design** a first-class, machine-readable,
validated artifact of the data-model stage — the same treatment pathway gene sets already
receive — and give **cell-lineage markers** the provenance + panel-coverage audit that
pathways have but markers lack. Consolidate the plan-doc sprawl in the same pass.

Three asks, in priority order:
1. **Comparison design** → a curated registry the data-model workflow resolves against the
   samplesheet and emits as a validated table; downstream reads it instead of re-declaring
   contrasts in five scripts.
2. **Cell-marker sources/validation** → provenance TSV + panel-coverage artifact, mirroring
   the pathway pattern exactly.
3. **Pathway sources/validation** → already shipped; confirm as canonical, no rebuild.

## What already exists (do not rebuild)

The data-model workflow already emits, from a curated YAML source plus a provenance TSV, a
validated tiered gene-set table (`pathway_sets.tsv`) with a per-set panel-coverage table
(`gene_set_panel_coverage.tsv`): tier-1 curated (frozen, git-HEAD freeze check),
tier-2 pinned MSigDB Hallmark, per-set `n_total / n_panel / usable / thin`. Eight downstream
scripts read `pathway_sets.tsv` as the single source of truth. **This is the exact idiom the
two new artifacts imitate.** Pillar 3 is therefore a no-op: confirm `pathway_sets.tsv` stays
canonical and align the marker artifact's shape to it.

The samplesheet (`results/data_model/samples.tsv`) already carries every field a comparison
resolves against: `treatment`, `condition`, `timepoint_h`, `model`, `slide_id`, `mouse_id`,
`dataset_id`, `dose_gy`. Arm membership is already data-driven from it; only the *level order*
and the *contrast list* are hardcoded downstream.

## The gap being closed

- **Comparisons are nowhere machine-readable.** The three day-2 contrasts and the
  `~0+condition+slide_id` formula are copy-pasted across `deg_pseudobulk.R`,
  `pathway_scores.R`, `composition.R`, `niches.R`, `spatial_mixing.R`. The
  confirmatory/exploratory tiering and the H1/H2/H3 family live as literals in
  `assemble_results.R`. There is no single artifact that states what comparisons the sample
  set admits, at what replication, with what status.
- **Markers carry no provenance and no coverage audit.** `lineage_markers.yaml` and
  `substate_markers.yaml` are consumed by ~8 typing/QC scripts. Their sources and their
  known weak separations (Pericyte↔SmoothMuscle, DC↔Macrophage) live in prose comments; the
  "every marker verified on-panel" claim is a comment, not a check.

---

## Artifact A — comparison registry

### A.1 Source: `config/comparisons.yaml`

One entry per comparison the sample set admits. The entry **declares intent**; the builder
**computes realized facts** against the samplesheet. Group membership is expressed as
predicates on samplesheet columns, so arm assignment stays data-driven (no embedded sample
lists).

```yaml
# config/comparisons.yaml
# Curated registry of the experimental comparisons this sample set admits.
# Group predicates filter results/data_model/samples.tsv; the builder resolves realized n,
# model rank, residual df, and inference-capability from the resolved sample subset.
# Pre-registered: tier + hypotheses are frozen (git-HEAD freeze check, like pathway sets).

comparisons:
  # ---- sample-level, day-2 inference cohort (n=4, formal) ----
  - name: SBRT_vs_Ctrl
    kind: sample
    resolution: [whole, compartment]          # scales at which the contrast is run
    cohort: mutter02_day2
    group1: {treatment: SBRT, timepoint_h: 48, model: flank}
    group2: {treatment: NT,   timepoint_h: 48, model: flank}
    unit: mouse
    formula: "~0 + condition + slide_id"       # slide = randomized complete block
    tier: confirmatory
    hypotheses: [H1_immune, H2_vascular, H3_stromal]
    dose_confounded: false

  - name: MBRT_vs_Ctrl
    kind: sample
    resolution: [whole, compartment]
    cohort: mutter02_day2
    group1: {treatment: MBRT, timepoint_h: 48, model: flank}
    group2: {treatment: NT,   timepoint_h: 48, model: flank}
    unit: mouse
    formula: "~0 + condition + slide_id"
    tier: confirmatory
    hypotheses: [H1_immune, H2_vascular, H3_stromal]
    dose_confounded: false
    notes: >
      Whole-compartment averaging cancels a peak-restricted MBRT signal (dilution).
      Bounds, does not resolve, the MBRT question — see region-level MBRT_peak_vs_valley.

  - name: MBRT_vs_SBRT
    kind: sample
    resolution: [whole, compartment]
    cohort: mutter02_day2
    group1: {treatment: MBRT, timepoint_h: 48, model: flank}
    group2: {treatment: SBRT, timepoint_h: 48, model: flank}
    unit: mouse
    formula: "~0 + condition + slide_id"        # estimated as a single-model contrast
    tier: confirmatory
    hypotheses: [H1_immune, H2_vascular, H3_stromal]
    dose_confounded: true                       # SBRT 20 Gy; MBRT mean/integral dose unrecorded
    view:
      differential_response_scatter: true       # per-gene MBRT_vs_Ctrl effect vs SBRT_vs_Ctrl effect;
                                                 # algebraically = this contrast at the linear level.
                                                 # Visualization only; NOT counted as independent evidence.

  # ---- region-level, within-sample ----
  - name: MBRT_peak_vs_valley
    kind: region
    axis: dist_to_peak                           # continuous, modeled nonlinearly (spline)
    test: periodicity_rotation_null              # Lomb-Scargle/Fourier at ~1 mm beam freq + rotation null
    cohorts: [mutter01_mbrt_4h, mutter02_day2]
    unit: mouse                                  # per-mouse peak-vs-valley delta, then across mice
    pairing: within_mouse                        # dose-clean: peak vs valley same mouse/same delivery
    requires: h2ax_registration                  # M01 4h (1 block); M02 day2 slides 2/3/4 have OME-TIFF
    tier: confirmatory                           # M02 n>1; M01 4h is descriptive (n=1)
    method_notes: [fov_alias_guard, dose_dropout_check, boundary_buffer, nonlinear_distance]

  - name: MBRT_valley_vs_SBRT
    kind: region
    cohorts: [mutter02_day2]
    unit: mouse
    requires: h2ax_registration
    tier: exploratory
    dose_confounded: true

  - name: niche_DA_by_arm
    kind: region
    test: kmeans_knn_composition                 # K swept; robustness of arm effect to K reported
    cohort: mutter02_day2
    unit: mouse                                   # niche frequency collapsed to one value per mouse
    formula: "~0 + condition + slide_id"
    tier: exploratory

  - name: tumor_immune_mixing_by_arm
    kind: region
    test: keren_mixing_score
    cohort: mutter02_day2
    unit: mouse                                   # per-mouse mixing score
    formula: "~0 + condition + slide_id"
    tier: exploratory

  - name: substate_composition_by_arm
    kind: region
    test: substate_fraction                       # e.g. resting vs activated fibroblast
    cohort: mutter02_day2
    unit: mouse
    formula: "~0 + condition + slide_id"
    tier: exploratory

  # ---- descriptive, time-course cohort (n=1, effect sizes only) ----
  - name: treated_vs_baseline_timecourse
    kind: sample
    cohort: mutter01_flank
    group1: {treatment: [MBRT, SBRT], model: flank}
    group2: {treatment: NT, timepoint_h: 0, model: flank}   # single shared 0h baseline (NOT time-matched)
    unit: sample
    tier: descriptive
    notes: single 0h control conflates time-since-implant with treatment; descriptive only.

  - name: MBRT_vs_SBRT_matched_timepoint
    kind: sample
    cohort: mutter01_flank
    strata: [4, 48, 144]                          # hours; one contrast per matched timepoint
    unit: sample
    tier: descriptive

  - name: MBRT_vs_SBRT_tongue_day10
    kind: sample
    cohort: mutter01_tongue
    unit: sample
    tier: descriptive
    notes: different tumor + site; timepoint-confounded cross-model context only.

  # ---- context / replication ----
  - name: cohort_concordance
    kind: context
    test: effect_size_concordance                 # M01 day2 vs M02 day2 per-feature
    tier: exploratory
    notes: M01 n=1 anti-correlates with M02 on the SBRT signal; contextual, cannot corroborate.
```

### A.2 Builder: `scripts/build_comparisons.R`

Reads `comparisons.yaml` + the samplesheet; for every `sample`/`region` entry with group
predicates it resolves the matching samples and **computes** the facts that validate the
design:

- `n_group1`, `n_group2`, `n_total` — resolved sample counts.
- `n_mouse_group1/2` — distinct mice (the true replication unit).
- `model_rank`, `resid_df` — from `model.matrix(formula)` on the resolved subset
  (self-checks estimability: e.g. `~0+condition+slide_id` on the 12 day-2 samples → rank 6,
  resid df 6).
- `inference_capable` — `TRUE` iff `unit`-level replication ≥ 2 per group **and** tier ≠
  descriptive.
- `gate` — for `region` entries, records `requires`; if the required input path is reachable
  from the build host, stamps `available`, else `unchecked` (keeps Stage A metadata-pure — it
  never hard-depends on heavy data).

Emits `results/data_model/comparisons.tsv` (one row per resolved comparison × resolution).
Frozen like pathway sets: a git-HEAD freeze check on `tier`/`hypotheses` hard-errors on silent
change to the confirmatory family. Hard-errors if a group predicate resolves to zero samples
(a typo'd contrast).

### A.3 Downstream hookup (phased, `aggregate_differential.smk`-touching — sequenced, not required for A)

The five scripts that hardcode the contrast list, and `assemble_results.R` (confirmatory
family + tiered FDR), all live in `aggregate_differential.smk` (the differential half of the
2026-07-12 split of the monolithic `aggregate.smk` into `aggregate_typing.smk` →
`full_labels.parquet` and `aggregate_differential.smk` → `results_master.tsv`). They would read
`comparisons.tsv` instead of embedded literals. **This edit touches the differential
workflow's code**, so it rides the next `results_master.tsv` regeneration rather than triggering
a separate one. Artifact A builds and validates with zero dependence on that rerun; the hookup is
a follow-on phase.

---

## Artifact B — cell-marker provenance + coverage

Keep both marker YAMLs as the source. Add a provenance TSV and a builder that mirrors
`build_gene_sets.R`.

### B.1 Source: `config/marker_sets_provenance.tsv`

```
set                 tier      source                                        weak_separation
T                   immune    ImmGen canonical T-lineage                     
B                   immune    ImmGen canonical B-lineage                     
Macrophage          immune    ImmGen canonical myeloid                       
DC                  immune    ImmGen canonical DC                            shared_with_Macrophage (Itgax/Ciita)
Pericyte            stroma    canonical perivascular                         shared_with_SmoothMuscle (3 panel markers; Pdgfrb fibroblast-shared)
SmoothMuscle        stroma    canonical SMC                                  
Fibroblast          stroma    canonical fibroblast                          
...
fibroblast_resting  substate  pre-registered (Zhang 2019 double-dip guard)  
fibroblast_activated substate pre-registered (Zhang 2019 double-dip guard)  
```

### B.2 Builder: `scripts/build_marker_sets.R`

Reads `lineage_markers.yaml` + `substate_markers.yaml` + the provenance TSV + `common_genes.tsv`.
Emits `results/data_model/marker_panel_coverage.tsv`: `set, tier, source, n_total, n_panel,
frac_panel, thin, weak_separation`. Enforces the two integrity claims the YAMLs currently only
assert in comments:
- **On-panel completeness:** the coarse-lineage YAML claims every marker is on-panel — the
  builder *checks* it and errors (or warns with an explicit off-panel list) on violation,
  instead of trusting the comment.
- **Freeze check:** markers are pre-registered; a working-tree vs git-HEAD divergence in the
  marker YAMLs hard-errors, same as tier-1 pathway sets.

Missing provenance for any set is a hard error (mirrors `build_gene_sets.R`).

---

## Snakemake wiring (`data_model.smk`)

Two new rules, both `threads: 1`, R via the existing `RSCRIPT` `conda run` pattern; scripts
declared as `input:` deps so code edits trigger reruns (per the script-tracking convention).
Both are Stage-A metadata-only — no heavy per-cell objects.

```
rule build_comparisons:
    input:
        script = "scripts/build_comparisons.R",
        yaml   = "config/comparisons.yaml",
        ss     = config["samplesheet"],
    output: "results/data_model/comparisons.tsv"

rule build_marker_sets:
    input:
        script = "scripts/build_marker_sets.R",
        lineage  = config["...lineage_markers.yaml"],
        substate = "config/substate_markers.yaml",
        prov     = "config/marker_sets_provenance.tsv",
        panel    = "results/data_model/common_genes.tsv",
    output:
        coverage = "results/data_model/marker_panel_coverage.tsv"
```

Add both outputs to `rule all` (the manifest-of-all-deliverables convention). One figure,
`marker_panel_coverage.png`, mirrors the existing `geneset_coverage.png` rule. A rendered
comparison-status table (comparison × realized n × resid df × inference_capable × gate) is a
low-cost addition mirroring the design-grid figure; include it.

---

## Method notes folded from the independent review (documentation, not new scope)

Carried as `method_notes` on the relevant registry entries and as prose in the eventual plan;
these harden *existing* methods, they do not add comparisons:

1. **FOV-grid aliasing** (`fov_alias_guard`): CosMx FOV tiles (~0.5–1 mm) sit near the ~1 mm
   beam period. The periodicity test must confirm the detected period is not phase-locked to
   the FOV grid and covary FOV position / total counts.
2. **Dose-correlated QC dropout** (`dose_dropout_check`): low-RNA peak/necrotic cells are
   preferentially filtered, which can mask or manufacture a peak signature. Check cell
   retention vs distance-to-peak; apply the necrosis mask; covary total counts.
3. **Per-mouse collapse** for every region-level between-arm test (niche DA, mixing,
   substate): one summary statistic per mouse (n=4) before testing; extend the rotation null
   beyond peak/valley to these.
4. **Nonlinear distance-to-peak** (`nonlinear_distance`): model with a spline, since valley
   priming can be maximal at max distance (opposite sign to a monotonic dose term).
5. **Segmentation boundary buffer** (`boundary_buffer`): restrict peak/valley calls to cells
   interior to a stratum by the segmentation-uncertainty margin (transcript-coord
   repartitioning is infeasible — no coords).

**Map correction (not an opinion):** M02 gamma-H2AX OME-TIFF is present on disk for slides
2/3/4 (each covering all three arms); slide 1 is on the erroring GCS mount. So within-tumor
peak/valley on the replicate cohort is a **data-available** comparison pending scan-to-cell
registration, not gated on staining. The registry records this via `requires:
h2ax_registration` + the builder's `gate` stamp. Promoting/depriotitizing any comparison is a
scientific decision left to the user and is **out of scope** for this spec — the registry
captures the set verbatim and status-stamps it.

---

## Plan-doc consolidation

**Eliminate (shipped or superseded; git history is the archive):**
- `plan-pathway-genesets.md` + `plan-pathway-genesets-impl.md` — Phase 1 shipped, Phase 2
  absorbed by the preprocessing split.
- `plan-mbrt-vs-sbrt.md` — superseded run-record; verify its content is captured in
  `plan-mbrt-vs-sbrt-reanalysis.md` + `CLAUDE.md` before deleting.

**Fold:**
- `plan-qc-metrics.md` → into `plan-preprocessing.md` (both complete, QC metrics live inside
  preprocessing).

**Correct (stale "v2.0 external-atlas" banners describing an approach never taken):**
- `plan-aggregate.md`, `plan-mbrt-signatures.md`, `plan-outcomes-de.md` — the actual typing
  is scVI cluster-then-annotate, not an external atlas.

**Keep, merge later (blocked on the pending `results_master.tsv` rerun — cannot delete a
doc for unshipped work):**
- `plan-aggregate-refactor.md`, `plan-differential-robustness{,-impl}.md`,
  `plan-mbrt-vs-sbrt-reanalysis{,-impl}.md`.

**Keep (forward-looking, now more live):**
- `plan-svg-v2.md`, `plan-mbrt-signatures.md`, `plan-outcomes-de.md` — the peak/valley +
  outcomes future work; `plan-svg-v2.md` is now runnable on M02 given the confirmed H2AX.

**New:** `plan-comparison-registry.md` (from writing-plans) becomes the implementation plan and
the narrative home for the comparison-space description currently scattered across the
reanalysis / outcomes / svg docs. Update the `spatial-rads.org` `** Plans` index + status
ledger to match all of the above.

---

## Verification

1. `data_model.smk` dry-run clean; `build_comparisons` + `build_marker_sets` run to green.
2. `comparisons.tsv`: every entry resolves to ≥1 sample; the three day-2 contrasts show
   `n_mouse=4/arm`, `resid_df=6`, `inference_capable=TRUE`; descriptive entries
   `inference_capable=FALSE`; region entries carry the correct `gate`.
3. `marker_panel_coverage.tsv`: coarse-lineage on-panel completeness check passes (or the
   off-panel list is surfaced); Pericyte/DC weak-separation flags present; freeze check active.
4. `pathway_sets.tsv` unchanged (confirm canonical, no rebuild drift).
5. Both coverage figures render and are inspected (axis ranges match, nothing clipped).
6. Plan consolidation: deleted docs' content confirmed captured elsewhere; org Plans index +
   ledger updated; no dangling cross-references.

## Scope

- **In:** `config/comparisons.yaml` + builder + `comparisons.tsv`; marker provenance +
  coverage builder + `marker_panel_coverage.tsv`; two `data_model.smk` rules + figures;
  plan-doc consolidation; method notes as documentation.
- **Phased (sequenced after the pending rerun):** rewiring the five `aggregate_differential.smk`
  scripts + `assemble_results.R` to read `comparisons.tsv`.
- **Out:** any change to which comparisons are run or how they are prioritized; new hypothesis
  families; the BANKSY typing benchmark; H2AX scan-to-cell registration itself; the
  `results_master.tsv` rerun.
```