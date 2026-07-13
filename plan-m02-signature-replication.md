# M02 peak/valley signature replication (skeleton)

> **Skeleton spec, 2026-07-13.** Companion to `plan-m01-peak-valley-signature.md`
> (which defines and freezes the 4h H2AX-anchored signature). This plan covers the M02 day-2
> replication cohort: the sample-level module test plus, now that H2AX exists for M02, the
> within-slide peak/valley re-analysis that was previously gated. Skeleton only: structure and
> decisions fixed, step detail to be filled when the FOV overlay is done. Spec location follows
> project convention (`plan-*.md` at root).

## Status of the gating input

**H2AX figures now exist for M02** (adjacent-section staining, previously deferred at Mayo). The
remaining human step is the **FOV overlay**: manually register the H2AX beam tracks onto each M02
MBRT slide's FOV grid to assign per-cell peak/valley zones. This uses manual FOV extrapolation (the
established method), not the geometric stripe model. **This overlay is a prerequisite the analyst
performs; every downstream step is gated on it.**

## What this plan is

Two readouts on the M02 day-2 cohort (n=4/arm, Control/MBRT/SBRT), consuming the frozen M01
signature from the companion plan:

1. **Whole-compartment module replication (no overlay needed).** Score every M02 cell for the frozen
   DDR-excluded M01 peak signature (UCell), aggregate **to sample** (not FOV: FOVs within a slide
   are pseudoreplication), and test MBRT vs Control and MBRT vs SBRT through the existing
   `aggregate_differential.smk` limma engine with the confirmatory-family FDR. **Positive-only
   interpretation** (null-ambiguity: a diluted whole-compartment null cannot separate loss from a
   still-active peak-restricted program). Lead with M01-4h-vs-M02-2d effect-size concordance, not a
   lone arm score. This readout can run now and does not depend on the overlay.

2. **Within-slide peak/valley re-analysis (gated on the FOV overlay).** Once zones are assigned, run
   the M01 plan's FOV-pseudobulk design on the M02 MBRT slides: paired `~ FOV + zone` peak-vs-valley
   DE, camera/UCell enrichment, coarse-compartment stratification, orthogonal morphology. This is
   the genuine spatial test the M01 signature was built to enable, now with independent-cohort H2AX
   ground truth at a later timepoint (2d, not 4h).

## Design carried from the M01 plan

- **FOV footprint vs 1.02 mm spacing gate** decides paired FOV-by-zone vs unpaired FOV-mean (re-check
  per M02 slide; CosMx FOV size may differ from M01).
- **DDR-gene exclusion** (circularity guard) in any newly derived M02 peak/valley gene list.
- **Coarse tumor/stroma/immune compartments only** for stratification (current atlas labels).
- **Morphology segmentation-artifact guard** (per-zone QC-flag/density adjustment before calling a
  shape shift biology).
- **No dose/sparing language** (M02 MBRT dose unrecorded).

## What replication at 2d can and cannot claim

- **Can:** whether the 4h-defined program is detectable in an independent cohort at a later
  timepoint (positive module elevation), and whether a within-slide peak/valley contrast re-emerges
  at 2d with its own H2AX ground truth.
- **Cannot:** persistence as a within-cell trajectory (no 4h anchor in M02; 4h and 2d are different
  cohorts and animals), or anything from a flat/null whole-compartment score.
- **Note the p21 caveat:** at 2d p21/DDR has decayed from its 4h peak, so a 2d stripe is expected to
  be weaker than 4h; a weak or absent within-slide contrast at 2d is an interpretable negative, not
  a failure, and must not be read as program loss.

## Retro evaluation and cleanup (close the loop)

The retired autonomous run left a post-mortem (`retro-mbrt-signatures.md`) whose portable lessons
were never harvested. As the final step of this plan, **after the new M02 run completes**:

1. **Evaluate** the retro's action items against how the new run actually went: which portable
   lessons (checkpoint discipline, parallelization preflight, api_error retry-burst detection in the
   nanny filter, elapsed-vs-timeout dashboard, sample-not-FOV aggregation) proved out or are now
   moot because this run is local rather than a merged-monolith VM job.
2. **Harvest** anything still portable into its correct upstream home (virtual-scientist skill, the
   global jsonl filter, a `feedback_vm_runs` memory note), not back into a project markdown.
3. **Delete** `retro-mbrt-signatures.md` once its content is either harvested upstream or made moot
   by the new run. The retro exists to be spent, not kept.

## To fill in when the overlay is done

- Per-slide FOV overlay method + zone-assignment QC (peak/valley cell counts per MBRT slide).
- Exact rule tables/scripts for the two readouts (reuse M01 plan's producers).
- Whether the whole-compartment readout (1) folds as a new gated exploratory row in
  `results_master.tsv` or stands as a separate small output.
