# Planning brief — spatially-variable-gene analysis v2

> **Spatial axis of the deferred MBRT track** (see `plan-mbrt-vs-sbrt-reanalysis.md`). The continuous
> `dist_to_peak` / beam-frequency method needs peak/valley localization — available for M01 4h,
> **absent for the M02 inference cohort** (the staining gate). Runs on M01 now; M02 awaits staining.

**Status: plan-of-plan (handoff seed, NOT a committed plan).** Hand off to a
fresh session to brainstorm + devils-advocate before writing a real
`plan-*.md`. Captures the 2026-05-31 audit of the v1 one-shot so planning can
start cold.

**v2.0 note — low Atlas dependence (2026-05-31).** This plan's core method
(continuous `dist_to_peak`, beam-frequency periodicity) is **not** keyed on
cell-type labels, so the v2.0 Atlas re-typing (`plan-processing-pipeline.md`)
does not block it. Only contact point: the *tumor-compartment subset* used as
the worked example — re-pull it from re-typed labels when running.

This is the spatial axis of the peak/valley question; the condition-level /
whole-compartment arm axis is owned by `plan-mbrt-vs-sbrt-reanalysis.md`. Scope: M01 flank Block_21 MBRT
4h is the worked example (peak FOVs known from H2AX IHC); generalizing to
other samples is an open question below.

---

## Job to be done (de-anchored)

Not "run MERINGUE again better." The job is: **rigorously test whether MBRT
produces dose-driven, beam-periodic transcriptional zonation** — signal at
the known ~1 mm beam frequency that exceeds what tissue architecture alone
would produce — and make a null result interpretable (no transcriptional
periodicity at the beam frequency), not ambiguous.

Key reframe from v1: treat dose position as a **continuous spatial
covariate** (`dist_to_peak`, the phase within the ~1 mm beam period), not a
binary peak/valley label. The binary split mis-assigns boundary cells,
discards graded signal, and a group-mean DE on n=1 has no legitimate null.

---

## What v1 did (build on, audit-aware)

All in `dev/if_channel_check/` and `dev/peak_valley_analysis/`:

- **`07_peak_valley_classify.R`** — computes continuous per-cell
  `dist_to_peak` (`:21`) from a tilt-corrected stripe model, then **discards
  it** by thresholding to peak/valley at spacing/4 (`:22`). Geometric stripe
  model is DEPRECATED (memory `feedback_pv_method`: prefer FOV
  extrapolation).
- **`08_pv_degs.R`** — binary `FindMarkers(peak vs valley)` on n=1 MBRT 4h
  (`:16`); recovered little. The motivating failure.
- **`meringue_full.R` / `meringue_tumor.R`** — proper MERINGUE pipeline
  (`getSpatialNeighbors` → `getSpatialPatterns` → `filterSpatialPatterns`
  BH α=0.05 → `spatialCrossCorMatrix` → `groupSigSpatialPatterns`). Tool
  invoked correctly.
- **`meringue_moran.R`** — hand-rolled Moran's I on a fixed-k kNN graph
  (`:40-41`) but **no significance test** (`:62` just ranks I).
- **`pv_gam_spline.R`** — continuous GAM on position (the right modeling
  shape; review how it was parameterized).
- **`pv_rotation_null.R`** — architecture control: fits stripe geometry to 14
  manual peak FOVs (`:11-12,38-39`), extrapolates labels to all tumor cells
  (`:86-100`), rotates to 200 random angles ≥30° off real (`:105-116`), keeps
  a gene only if `rot_null_p<0.05 AND |delta|>2*null_sd` (`:164`). **Recorded
  result: 37/670 genes survive at 4h, null sd 0.027 vs SBRT 0.016** (org
  bullet; `/tmp` artifacts deleted, so re-run to re-derive).

---

## Audit findings — why v1 was a one-shot screen, not a test (2026-05-31)

1. **Heavy subsampling + fixed-radius graph = degraded neighbor density.**
   `meringue_tumor.R:9` subsamples ~100k tumor cells to `N_SUB=2500` (~40×),
   then builds neighbors at fixed `filterDist=0.05 mm` (`:49`). At 40× lower
   density a 50 µm radius catches <1 neighbor — guts Moran's I power and
   confounds it with density. The robust-graph version (`meringue_moran.R`,
   fixed-k) drops the significance test. So stats and robustness never
   co-occurred.
2. **Moran's I detects ANY spatial structure, not beam-periodic structure.**
   A high-I gene may track tumor-core/margin or a single patch, not the
   stripes. The **periodicity-specific test at the known ~1 mm frequency
   (Lomb-Scargle / Fourier on `dist_to_peak`) was never run** — this is the
   real gap.
3. **No architecture null in the MERINGUE path.** `pv_rotation_null.R` is the
   only dose-specificity control and it was **never crossed** with the
   MERINGUE/Moran SVG hits.

---

## Candidate v2 components (the menu to brainstorm)

- **Primary test = periodicity at the known beam frequency**: Lomb-Scargle or
  Fourier power of each gene's expression vs `dist_to_peak` (phase within the
  ~1 mm period). A within-sample **coordinate-permutation null** sidesteps
  the n=1 biological-replication problem and makes a null result
  interpretable.
- **Continuous-covariate regression (GAM on `dist_to_peak`) as the model**,
  replacing binary peak/valley (`pv_gam_spline.R` is the seed).
- **Consistent neighbor graph**: fixed-k kNN throughout, or drop/greatly
  reduce subsampling if a fixed-radius graph is kept — don't pair heavy
  subsampling with a density-sensitive radius.
- **Intersect spatial significance with dose specificity**: keep genes that
  are MERINGUE/periodicity-significant AND survive the rotation null.
- **Model-free cross-checks**: Moran's I / SpatialDE / SPARK as sanity, not
  primary.

---

## Open questions for the planning session

- **Geometry for samples without manual peak FOVs.** Block_21 has H2AX-picked
  peak FOVs; M02 (day2, n=4, no H2AX ground truth) does not. Can beam phase
  be recovered de novo (e.g., spectral peak in cell density or in a damage
  marker), or is the periodicity analysis limited to M01 4h?
- **Combining across samples**: M02 replicates could give a real biological
  null — but only if beam geometry is recoverable per sample.
- **Multiple testing** across ~670–950 genes; FDR on periodicity p-values.
- **Which marker anchors phase** when H2AX has faded (memory
  `project_alternative_staining`: TUNEL/Ki67/H&E on adjacent sections).

---

## Housekeeping

- **Revise the org SVG bullet.** The current "Key results" bullet in
  `spatial-rads.org` (after the dilution caveat) frames SVG as *proposed /
  untried*. It overstates: MERINGUE (`meringue_*.R`), Moran's I
  (`meringue_moran.R`), GAM (`pv_gam_spline.R`), and a rotation-null
  architecture control (`pv_rotation_null.R`) already ran. Rewrite to: v1
  done as exploratory screen → 3 audit weaknesses → v2 = periodicity test +
  rotation-null intersection. (org-edit + sci-write skills apply.)
- **Rotation-null is recorded-only.** `/tmp` inputs/outputs deleted; re-run
  `pv_rotation_null.R` to re-derive the 37-gene figure before citing it as
  fresh.

---

## Cross-references

- Continuous covariate: `dev/peak_valley_analysis/07_peak_valley_classify.R:21-22`
- Binary DE failure: `dev/peak_valley_analysis/08_pv_degs.R:16`
- MERINGUE: `dev/if_channel_check/meringue_full.R`, `meringue_tumor.R`, `meringue_moran.R`
- GAM: `dev/if_channel_check/pv_gam_spline.R`
- Architecture null: `dev/if_channel_check/pv_rotation_null.R`
- Stripe model conventions: 15° tilt, 1.02 mm spacing, 4 peaks (CLAUDE.md);
  FOV-extrapolation preferred over geometric model (memory `feedback_pv_method`).
