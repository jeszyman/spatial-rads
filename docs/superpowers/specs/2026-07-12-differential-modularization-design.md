# Modularizing the differential layer — design

> **Status: IMPLEMENTED + validated 2026-07-13.** Built per
> `docs/superpowers/plans/2026-07-12-differential-modularization.md`; the convergent run passed
> run-once scientific validation and the results master is byte-identical to pre-refactor. The
> architecture below is now the live code (engines in `scripts/engines/`, correction consolidated
> in `assemble_results.R`). Implementation notes and gotchas: memory `project_differential_modularization`.

## Goal

Split the bespoke downstream differential workflow (`workflows/aggregate_differential.smk`)
into (1) reusable, experiment-agnostic **test engines**, (2) **feature-producers** that emit
standard-shaped inputs to those engines, and (3) a study-specific **results layer** that frames
the scientific questions. The workflow currently hard-codes this one experiment: the literal
`m02day2` appears in ~30 output filenames, and the three treatment contrasts (`MBRT_vs_Ctrl`,
`SBRT_vs_Ctrl`, `MBRT_vs_SBRT`) plus the design formula `~0 + condition + slide_id` are
copy-pasted as literals across five scripts. The target makes the study an *analysis that calls
tools*, not a monolith.

Scope of reuse: **within spatial-rads only.** The engines are their own workflow(s) in this
repo, not a separate package.

## Architecture

Three tiers, decoupled at two well-defined table schemas.

### Tier 1 — Two test engines (modular, experiment-agnostic)

Kept as **two separate engines** (not one family-switched engine): their input data types and
model families differ enough to keep distinct.

- **Count engine** — negative-binomial pseudobulk differential expression (DESeq2/edgeR).
  Input: a `sample × cell_type × gene` count matrix + a design formula + a contrast (all from
  the comparison registry). Output: per-gene sufficient statistics. This is the field-standard
  method for pseudobulk-per-cell-type DE.
- **Linear-model engine** — `limma` on a `sample × feature` matrix, with a transform parameter:
  identity for continuous features, logit-via-`speckle::getTransformedProps` for proportions.
  Input: a `sample × feature` matrix + design + contrast. Output: per-feature sufficient
  statistics. Serves cell-type composition, spatial niches, substate composition (proportions),
  and pathway scores, tumor-immune mixing, myeloid M1/M2 (continuous).

Both engines read their contrast, formula, and unit from `results/data_model/comparisons.tsv`
(the comparison registry built 2026-07-12) — never from embedded literals.

**Engine output schema (both engines, identical):**
`comparison, contrast, unit, feature_type, feature_id, estimate, se, stat, p`.
Engines emit **sufficient statistics — estimate, SE, df — NOT adjusted verdicts.** The results
tier needs these to re-pool corrections across a family; an engine that emitted only
per-comparison `padj` would force the results tier to reach into engine internals. Hard
requirement.

Also modular (per-comparison or cohort-level checks with clean tabular output; NOT results-tier):
technical QC (arm balance, replicate reproducibility), power/MDE, panel coverage.

### Tier 2 — Feature-producers (live in the study layer; emit standard-shaped tables)

A producer computes the thing a test consumes and writes it in one of two standard schemas:

- **Per-cell label parquet:** `cell, sample_id, label` (+ coords if the producer needs them) —
  consumed by the linear-model engine's proportion path. Producers: cell-type (upstream), niche
  clustering, fibroblast substate.
- **Sample × feature matrix:** `sample_id, feature_id, value` (continuous) or
  `sample_id, cell_type, gene, count` (counts) — consumed by the respective engine.

**Producers live OUTSIDE the engines, in the study layer.** They are experiment-specific
recipes (kNN k-means spatial niches; Keren tumor-immune mixing score; myeloid M1/M2 UCell
scoring on macrophages; fibroblast resting/activated substate gate; UCell/AddModuleScore pathway
scoring; plus generic pseudobulk count aggregation). The engine never knows how a feature was
derived — only its standard-shaped output crosses the boundary. This keeps the reusable engines
free of study biology. A producer being "bespoke" does not make it non-modular at its *output*:
a k-means niche label is shaped exactly like a cell-type label, so the proportion test is blind
to its origin.

### Tier 3 — Results / questions layer (study-specific, stays bespoke)

Consumes engine outputs (sufficient statistics) and produces the study's answers:
- Pre-registered hypothesis assignment (H1 immune / H2 vascular / H3 stromal).
- The cross-family multiple-testing correction (currently BH-pooled → `padj_confirmatory`). This
  runs *after* all per-comparison statistics exist, so it is structurally a results-tier step.
  **The correction method itself is an open methods/best-practice question, out of scope here** —
  the architecture only fixes *where* it runs.
- Detectability table + MBRT spatial-dilution narrative; figures + framing; M01↔M02 concordance
  read; and selection of which comparisons are worth making.

## Decisions (settled)

| Decision | Choice |
|---|---|
| Reuse scope | This repo only (engines are workflow(s) in spatial-rads, not a package) |
| Count vs. linear-model | **Two separate engines** |
| Feature-producer placement | **Study layer**, outside the engines; only standard-shaped output crosses |
| Engine output | Sufficient statistics (estimate/SE/df), never adjusted p |
| Parameter source | `results/data_model/comparisons.tsv` (registry), never embedded literals |
| Multiple-testing method | Deferred to a separate methods question; architecture fixes only its tier |
| Build order | Implementer's discretion |

## Inputs that already exist

- **Comparison registry** (`config/comparisons.yaml` → `results/data_model/comparisons.tsv`,
  built by `scripts/build_comparisons.R`, wired in `data_model.smk`, 2026-07-12): one row per
  comparison with groups, replication `unit`, blocking, model `formula`, `tier`
  (confirmatory/exploratory/descriptive), computed `inference_capable` / realized n / resid_df /
  `gate`. This is the engines' parameter source; nothing consumes it yet — wiring it is job one.
  12 comparisons: M02 day-2 arm contrasts; M01 n=1 descriptive time-course; gated within-tumor
  peak/valley region contrasts (`requires: h2ax_registration`, unwired).
- **Provenance/coverage artifacts** (built 2026-07-12): `pathway_sets.tsv`,
  `marker_panel_coverage.tsv` — the pathway/marker producers read these as gene-set inputs.
- **Migration inventory:** the current `aggregate_differential.smk` — 14 rules, 43 outputs
  (`rule all`).

## Design risk to respect

The generalized multi-comparison engine must not silently assume the M02 day-2 shape
(group = arm, unit = mouse, blocking = slide). The within-tumor peak/valley comparison has a
different unit (regions within a tumor) and grouping (pseudobulk aggregated within region within
sample), so an engine that bakes in the arm/mouse shape will need `if comparison == …` branches
the moment peak/valley is wired — a leaked abstraction. The `unit` and `formula` fields in the
registry exist precisely so the engine reads structure rather than assuming it. (Peak/valley
itself is out of scope here — separate dev-stage track, H2AX-gated: `plan-svg-v2.md`,
`plan-mbrt-signatures.md`.)

## Out of scope

- The cell-typing workflow (`aggregate_typing.smk`) — upstream, separate, untouched.
- Peak/valley within-tumor spatial analysis — separate track; it is the second comparison shape
  the engine must not preclude, not part of this build.
- Re-running the typing rebuild or refreshing `results_master.tsv`.
- The choice of multiple-testing correction method.

## Starting pointers

- Refactor target: `workflows/aggregate_differential.smk` (14 rules, 43 outputs).
- Count engine, from: `scripts/aggregate/deg_pseudobulk.R` + `pseudobulk_build.R`.
- Linear-model engine, from: `composition.R`, `pathway_scores.R`, `niches.R`, `spatial_mixing.R`,
  `myeloid_polarization.R`, `substate_composition.R` (all call `lmFit` / `getTransformedProps`).
- Feature-producers (→ study layer): `niches.R` (kNN+kmeans), `spatial_mixing.R` (Keren),
  `myeloid_polarization.R` (UCell M1/M2), `substate_split.R` (UCell + specificity-anchor gate).
- Results tier, from: `assemble_results.R` (hypothesis tagging + tiered FDR → `results_master.tsv`).
- Parameter source to wire in: `results/data_model/comparisons.tsv`.
- Compute: jeff-beast, `spatial-rads` conda env. SSH from jeff-pad needs
  `source ~/miniconda3/etc/profile.d/conda.sh` first (memory `reference_pad_beast_compute`).
