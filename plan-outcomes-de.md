# Planning brief — aggregate DE → clinically relevant outcomes

> **Harmonization note (2026-06-22):** the DE findings recorded below are on PRE-atlas labels and
> pre-reframe. The current canonical MBRT/SBRT framing is `plan-mbrt-vs-sbrt-reanalysis.md` (MBRT
> arm-level result is spatial *dilution*, not a null; all 3 contrasts effect-size-forward; spatial
> gated on M02 staining). This brief is a still-open *outcomes/clinical-relevance* direction —
> distinct from the mechanism reanalysis; re-pull any numbers against the reproducible labels.

**Status: plan-of-plan (handoff seed, NOT a committed plan).** Hand off to a
fresh session to brainstorm + devils-advocate before writing a real
`plan-*.md`. Captures context from the 2026-05-31 discussion so planning can
start cold without re-deriving.

Companion to `plan-aggregate.md` (the workflow that produces the DE tables
this plan consumes). Where `plan-aggregate.md` is about *building the
machine*, this brief is about *pointing it at outcome-relevant biology*.

## v2.0 — findings below are on PRE-atlas labels (2026-05-31)

The recorded findings here (`deg_summary_m02day2.tsv` counts, the tumor `a` /
`tumor_epithelial` / Pericyte / Stem.Prog cuts, propeller composition) were
computed on v1 ImmGen + de-novo-`a` labels. Cell typing has since been **replaced**
by merged-scale re-typing: scVI integration then Leiden cluster-then-annotate, plus
tier-2 SingleR/UCell and marker rescue (canonical in `CLAUDE.md`; the external-atlas
approach named in the 2026-05-31 draft never shipped). The `a` bucket and
`tumor_epithelial` conventions changed. **Re-verify every
cell-type-keyed number against re-typed outputs when promoting this brief to a
committed plan.** The dilution and MBRT-vs-SBRT/p21 findings (compartment-level,
not label-fragile) should hold; the precise per-cell-type counts will not.

---

## Job to be done (de-anchored)

Not "add more DE stages." The job is: **translate the aggregate cohort's
expression differences into readouts a radiation oncologist would recognize
as outcome-relevant** — tumor control vs immune activation vs normal-tissue
sparing — and state, honestly, what the current data can and cannot support.

The originating scientific question (CLAUDE.md): do MBRT *peaks* show direct
damage signatures while *valleys* show immune priming / bystander effects,
and how does the whole-tumor response compare to uniform SBRT? Plan 2
(`plan-svg-v2.md`) owns the peak/valley spatial axis; THIS plan owns the
**cell-type-resolved, condition-level, outcome-mapped** axis.

---

## What already exists (build on, don't rebuild)

`aggregate_differential.smk` DE machinery (per `plan-aggregate.md`), stages built +
verified this session: merge, composition (propeller), pseudobulk_build,
deg_pseudobulk, gsea, pathway_summary. Remaining aggregate stages
(deg_percell, concordance_m01_m02, niche tracks) are still pending and may
feed this.

Outputs already on disk in `results/aggregate/`:
- `deg_summary_m02day2.tsv` — per-cell-type DE counts + top genes, three
  contrasts (MBRT_vs_Ctrl, SBRT_vs_Ctrl, MBRT_vs_SBRT), n=4 replicated.
- `composition_test_m02day2.tsv` — propeller cell-type proportion shifts.
- `pathway_scores_summary.tsv` (+ test, concordance) — UCell/AddModuleScore
  Hallmark + curated (IFN-I/II, DDR, STING) per cell type, M02 limma test.
- `dev/if_channel_check/cdkn1a_kinetics_summary.tsv` — M01 p21 timecourse
  by cell type (exploratory, not in aggregate yet).

---

## Recorded findings to build on (re-verify against cited files when planning)

1. **Whole-sample dilution is the central obstacle.** M02 day2 pseudobulk:
   MBRT-vs-Control ≈ 5 DE genes in tumor "a" vs 127 for SBRT-vs-Control
   (similar ~20–30× gap in tumor_epithelial 0/99, Pericyte 2/130,
   Stem.Prog 3/97). Averaging directly-irradiated peak cells with spared
   valley cells collapses the MBRT signal. (`deg_summary_m02day2.tsv`)
2. **MBRT and SBRT are NOT identical where they differ.** MBRT-vs-SBRT
   direct contrast: **Cdkn1a/p21 up in MBRT** across ~19 cell types;
   proliferation (Top2a/Mki67/Hmgb2) and collagens (Col1a1/Col3a1) up in
   SBRT. (`deg_summary_m02day2.tsv`)
3. **Dynamics differ in shape, not just magnitude.** p21 timecourse:
   **MBRT = sharp spike that resolves; SBRT = lower plateau that persists.**
   Myeloid MBRT p21 exceeds SBRT (macrophage 2d 1.31 vs 0.83); by day6 MBRT
   has largely resolved while SBRT stays elevated in tumor/T/stroma.
   (`cdkn1a_kinetics_summary.tsv`, verified 2026-05-31)
4. **Composition shifts are nominal-only.** propeller n=4: SBRT-vs-Ctrl shows
   Plasma/GC-centroblast/B/CD8-T depletion and FRC/DC enrichment at nominal
   p, but **nothing survives global BH** (min padj ≈ 0.11). MBRT far fewer
   nominal hits. (`composition_test_m02day2.tsv`)

---

## Candidate directions (the menu to brainstorm)

- **Map DE to outcome-relevant programs**, not gene counts: proliferation
  arrest / senescence (p21, Cdkn1a), immunogenic cell death, antigen
  presentation (MHC-I, B2m, H2-K1/D1), IFN-I/II response, DDR resolution,
  fibrosis/normal-tissue toxicity (collagens — note these rose in SBRT).
- **Tumor control vs immune priming vs normal-tissue sparing as three named
  axes**, scored per condition/timepoint, so the MBRT "sparing" hypothesis
  becomes testable rather than rhetorical.
- **Cross to published RT-response / outcome signatures** (radiosensitivity
  index, IFN/T-cell-inflamed, senescence-associated secretory phenotype) —
  caution: 1000-gene CosMx panel has blind spots; check coverage first.
- **Per-cell-type DE (deg_percell)** as the un-diluted complement to
  pseudobulk: the dilution that kills whole-sample MBRT may be partly
  recoverable by conditioning on cell type (myeloid p21 already shows this).
- **Tumor vs stroma vs immune compartment decomposition** before contrasts.
- **Decide the inference boundary explicitly:** M02 day2 (n=4) is the only
  formal-inference stratum; M01 timecourse (n=1, Control only at t=0) is
  descriptive and CANNOT test "does MBRT change mirror Control over time."

---

## Deferred optimization (carry into planning)

`pathway_summary.R` spends a long single-threaded `AddModuleScore` pass to
reproduce UCell, which it parallelizes — smoke showed **UCell vs AMS
concordance r ≈ 0.97–0.98** per pathway. Once verified end-to-end, drop AMS
(or run it on a subsample purely as a concordance check) to cut this rule's
wall time. Low-risk; do when folding the stage into `spatial-rads.org`.

---

## Open questions for the planning session

- What counts as "clinically relevant" for a preclinical 4T1 flank model —
  surrogate outcome signatures only, or is there a path to actual tumor
  control / survival data from the Mutter lab?
- Pseudobulk vs per-cell-type DE as the primary inference unit?
- How hard to lean on M01 timecourse given n=1 and no untreated time series?
- Multiple-testing strategy across cell types × contrasts × pathways.

---

## Inputs / data available

- Aggregate DE tables in `results/aggregate/` (above).
- 23 per-sample `scored.rds` at `/mnt/data/projects/spatial-rads/processing/scored/`.
- Merged object `/mnt/data/projects/spatial-rads/aggregate/merged.rds`.
- Curated pathway sets `config/pathway_gene_lists.yaml`.
- Cell-type conventions: tumor compartment = `a` + `tumor_epithelial`
  (CLAUDE.md).
