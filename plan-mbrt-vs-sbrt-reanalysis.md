# MBRT vs SBRT vs Control — arm-level reanalysis (consolidate · reframe · look honestly)

> **Spec — brainstormed 2026-06-22, revised after a devils-advocate panel (4/4
> solid-with-caveats).** Reframes and hardens the *executed* `plan-mbrt-vs-sbrt.md` after the
> reproducible typing rebuild (`plan-aggregate-refactor.md`, Steps 1–5, bit-identical) and a
> 2026-06-22 data audit. Companion to `plan-mbrt-vs-sbrt.md` (original design + run record).
> Spec location follows the project convention (`plan-*.md` at root).

## What this re-plan is

Look at all three arm contrasts — **MBRT vs Control, MBRT vs SBRT, SBRT vs Control** — at the
merged cross-dataset scale, on the now-reproducible labels, and report them **honestly by effect
size and trend, not by significance gating**. It (1) **consolidates** the executed analysis onto
the reproducible labels + the consolidated tiered gene-sets; (2) **reframes** the MBRT result
correctly (real but spatially diluted, not an established null); (3) **strengthens** the read with
effect-size + coverage transparency. The within-tumor peak/valley spatial analysis is **gated on
staining the M02 cohort does not yet have** (see Gating inputs) — so the n=4 inference cohort can
only be analyzed at arm-level; that is a constraint, not a choice.

## The correct framing of the MBRT result (this fixes the panel's circularity critique)

MBRT delivers dose as ~1 mm-spaced peaks/valleys. The dev `dev/peak_valley_analysis/` work on
**M01 4h (with H2AX ground truth)** already found a **real, peak-localized MBRT signature** — 37
genes surviving a within-block rotation null, dominated by **MHC-I (H2-K1/H2-D1/B2m), IFN-response
(Ifitm2), and damage/stress**. So MBRT *has* a response; it is **spatially restricted to the beam
peaks**. Whole-compartment averaging across the ~50% valley volume therefore **cancels** it. The
arm-level MBRT result is that cancellation — **not** a biological absence. We report the arm-level
MBRT contrasts as **bounded / diluted**, never as a "coherence-tested null," and we cite the M01 4h
spatial result as why.

## What the data show (2026-06-22 audit of `results_master.tsv`, effect-size view)

- **SBRT vs Control — the strong whole-tissue signal.** Stromal fibrosis (collagens + `Acta2`),
  ~27 confirmatory / 706 exploratory hits; the only contrast with a robust whole-compartment
  effect. **MBRT vs SBRT** recovers SBRT > MBRT for the same stromal genes.
- **MBRT vs Control — weak and diluted at arm-level.** 2 confirmatory hits (`Ifitm1` ↑ in ILC;
  `Ciita` ↓), and program-level interferon scores flat/slightly-down across cell types. Consistent
  with a peak-restricted signal averaged out (above), **not** evidence of no MBRT effect.
- **Dose caveat (carried in the headline, not a footnote).** SBRT = 20 Gy uniform; the M02 MBRT
  mean/integral dose is **unrecorded**. So SBRT > MBRT differences cannot be attributed to dose
  *pattern* vs *magnitude*. We label every MBRT-vs-SBRT difference **dose-confounded** and never
  use the word "sparing" until the beam physics arrive.

## Design

### 1. Three contrasts, every readout — effect-size-forward
For each readout below, report **effect size + 95% CI + n=4 MDE**, for all three contrasts
(MBRT_vs_Ctrl, MBRT_vs_SBRT, SBRT_vs_Ctrl). Significance is shown but **does not gate** what is
looked at — sub-threshold MBRT effects are displayed (the user's "look at the trend" ask). Trend
ranking and direction views are labeled **exploratory** (no confirmatory multiplicity claim on
them).

### 2. A priori program panel — symmetric coverage
The 8 curated primary programs (`results/data_model/pathway_sets.tsv`), UCell + AddModuleScore, per
cell type × 3 contrasts. **Every set carries its panel-coverage fraction explicitly** — including
**TypeI_interferon 7/19** (~37%, the same sparsity as the flagged Hypoxia 5/13 and STING 2/3); no
set is silently protected. A null on a thin set is read against a **higher effect-size bar** scaled
to set size (the gate stringency tracks coverage), so a 2-gene-set null never reads as confident.

### 3. Detectability table (replaces the over-claimed "decision-gate")
Per readout × contrast, a transparent table: **effect, MDE, panel coverage, and `effect ≥ MDE?`**
— plus a magnitude floor: a direction is only called a "trend" when |effect| clears a stated
fraction of the MDE (direction at sub-MDE magnitude is noise, not signal). No binary
"coherent-null vs underpowered-null" label — that was one-off judgment; we report the components
and let the reader see them. The honest summary per readout: *detectable & present / detectable &
absent / not detectable at this n.*

### 4. The SBRT-positive story (whole-tissue, dose-gated)
SBRT-driven stromal fibrosis as the robust arm-level finding: stromal pseudobulk DE +
Fibrosis/Stromal-stress programs, via MBRT_vs_SBRT and SBRT_vs_Ctrl — carried with the
dose-confound label throughout. No volcano plots (sparse panel).

### 5. M01 role — descriptive context only (demoted)
M01 (n=1) provides **descriptive trajectories and the 4h spatial reframing evidence**, never a
falsification/coherence arm (it is n=1 by the project convention, and `concordance_m01_m02.tsv`
shows M01 anti-correlates with M02 for the SBRT signal, so it cannot corroborate). The M01 4h
peak/valley result is **cited** (from dev) as the reframe; it is not re-run here.

### 6. Supporting (kept, not leading)
Composition (propeller, slide-blocked), the ~50-Hallmark exploratory scan, k-means kNN niches,
Keren tumor-immune mixing, myeloid M1/M2.

### 7. Engineering / consolidation
- Differential layer reads the reproducible `full_labels.parquet` (done) + the `data_model.smk`
  tiered gene-sets (`pathway_summary.R`/`gsea.R` source `pathway_sets.tsv`, no live `msigdbr`).
- **Step 6 of `plan-aggregate-refactor.md`** regenerates `results_master.tsv` reproducibly on the
  bit-identical labels — result-invariant, the substrate for this re-plan.
- **Verify the old prereq fixes are closed before Step 6:** pseudobulk label re-point (confirmed —
  reads `full_labels.parquet`); spatial coords retained (`coords_necrosis`); snakemake arg desync
  (clean dry-run of the differential rules against the current scripts).

## Gating inputs (the two things that would unlock the MBRT question)

1. **M02 peak/valley staining** (H2AX or peak-FOV localization). The M02 cohort has none, so the
   dose-confound-immune within-slide spatial contrast — the analysis that can actually resolve
   MBRT — **cannot run on the inference cohort yet.** This is why the n=4 cohort is arm-level only.
2. **MBRT mean/integral dose** (Fazzari/Mutter beam physics) — gates pattern-vs-magnitude
   interpretation of every MBRT-vs-SBRT difference.

Until these arrive, the arm-level read + the M01 4h spatial citation are what we have.

## Scope

- **In:** all three arm contrasts at merged scale; M02 day-2 inference (n=4, slide-blocked) + M01
  n=1 descriptive; consolidation + reframe + effect-size/coverage transparency.
- **Out / gated:** M02 within-tumor peak/valley spatial (blocked on staining); MBRT dose
  interpretation (blocked on physics). M01 4h spatial is cited, not re-run.

## Limitations (panel-informed)

- Whole-compartment averaging is the **wrong scale for MBRT** — it cannot see a peak-restricted
  signal, so the arm-level MBRT contrasts **bound, they do not resolve** the MBRT question.
- Dose mismatch unresolved → MBRT-vs-SBRT is dose-confounded inference, labeled as such.
- n=4 (6 residual df); panel blind spots flagged symmetrically (TypeI IFN 7/19, Hypoxia 5/13,
  STING 2/3).
- M01 n=1 — descriptive only.

## Verification

1. Differential rules dry-run clean against `full_labels.parquet` + `data_model` gene-sets.
2. `results_master.tsv` regenerated reproducibly; SBRT fibrosis + the diluted MBRT contrasts
   reproduce.
3. Detectability table produced for every readout × contrast (effect / MDE / coverage / clears-MDE),
   with the magnitude floor applied; coverage shown for all sets.
4. A priori program-panel figure renders and is inspected (symmetric coverage annotation, not
   clipped/overplotted).

## Provenance

Brainstormed 2026-06-22; reframed after a blinded devils-advocate panel (4/4 solid-with-caveats,
run `wf_8f641118`) — folded fixes: dropped the over-claimed "coherence-tested null," gated the dose
language, made panel-coverage flagging symmetric + coverage-scaled, demoted M01 from a falsification
arm, replaced the binary decision-gate with a transparent detectability table + magnitude floor,
and reframed the MBRT arm-level result as spatial dilution of the known M01 4h peak signature. The
panel's "elevate spatial to primary" recommendation is answered by the staining constraint (M02 has
none). Consumes the reproducible typing from `plan-aggregate-refactor.md`. Companion to the executed
`plan-mbrt-vs-sbrt.md`.
