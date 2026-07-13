# Reconcile the contested peak/valley signature: rotation null vs FOV-pseudobulk

> **RESOLVED 2026-07-13 (Steps 0-2 executed on beast). Verdict: the 37-gene signature is an
> artifact; no dose-driven peak/valley signal survives.** The rotation null was never the culprit.
> Its own statistic, rerun, gives 38 survivors on the retired ImmGen `a` bucket (reproducing the
> "~37"), 7 after swapping to current atlas tumor cells + current QC (the IFN module collapses;
> MHC-I splits, only H2-D1 holds), and 0 once aggregated to the FOV replicate unit (paired
> `~ FOV + zone` DESeq2 + continuous phase regression, flat p-histogram, frac p<0.05 = 0.046).
> Two-stage death: cell-label/QC swap kills the inflammatory program, pseudoreplication kills the
> residual. A phase null was geometrically futile (tumor spans ~4.9 beam periods, p-floor
> ~0.2-0.47) and skipped. Lone near-miss Clu is a single-gene lead, not the claimed program.
> Figure `dev/if_channel_check/pv_reconcile_figure.png`; scripts `pv_step{0,1,2}_*.R`. Full
> per-step results recorded inline below. Retained as the executed record, not a live to-do.
>
> **Spec, 2026-07-13.** Two defensible methods give opposite answers on the 4h MBRT
> Block_21 tumor: the within-block rotation null reports a 37-gene MHC-I/IFN/damage peak
> signature, while paired FOV-pseudobulk group-mean DE on current QC'd cells is null (0 genes,
> flat p-histogram, positive control fires). This plan isolates *why* and decides whether any
> dose-driven peak/valley signal survives. Companion to the Peak/valley idea-node
> (`spatial-rads.org`, Project Map) and the committed null result (`** Results`, Post-IR
> dynamics). Spec at root per convention.

## The disagreement is not one difference but four (this is the crux)

The rotation null (`dev/if_channel_check/pv_rotation_null.R`) and the FOV-pseudobulk test differ
on four axes simultaneously, so a naive rerun confounds them. Each must be isolated:

1. **Cell population.** Rotation null uses the retired ImmGen `a` bucket (`main_type == "a"`),
   which is only ~32% epithelial (rest endothelial/fibroblast/SMC per
   `project_tumor_model`). FOV-pseudobulk used the current merged-scale atlas *tumor*
   compartment. Different cells entirely.
2. **Counts / QC.** Rotation null reads pre-re-base, pre-current-QC parquets
   (`/tmp/mutter01_*.parquet`, now deleted; regenerable from the input metadata + counts
   parquets). FOV-pseudobulk reads `norm.rds` (current 4-criteria QC).
3. **Zone labels.** Rotation null uses a wide peak band (`within_sd + 0.15`, ~0.27 mm
   half-width, plus a transition class). FOV-pseudobulk used tight core bands (dist<0.10 peak /
   >0.40 valley), sharpening the contrast.
4. **Statistic (and normalization).** Rotation null is a **pooled per-cell** log-CPM
   `rowMeans(peak) - rowMeans(valley)`, every cell weighted equally, no aggregation. FOV-pseudobulk
   sums to FOV×zone units and blocks on FOV, so the replicate unit is the FOV, not the cell. Note
   this bundles two separable things: the **aggregation unit** (cell vs FOV, the pseudoreplication
   axis) and the **normalization** (per-cell log-CPM vs pseudobulk library-size), which matter
   independently for sparse counts. Step 2 should hold normalization fixed while swapping only the
   aggregation unit, so a difference is attributed to the right cause.

**The signal to explain, stated carefully.** B2m and H2-K1 are peak-*up* in the rotation null but
carry small negative coefficients in the FOV test (-0.24, -0.23). Do **not** over-read this as a
biological sign flip. The FOV test is null with a flat p-histogram, so those two coefficients are
unstable draws from a null and "valley-leaning" is noise rather than direction. The thing that
actually needs explaining is why one method calls significance and the other does not, a
discrepancy in calibration or power, not a reversed effect. The prime suspects are #1 (cell
population) and #4 (pooled-mean weighting); a genuine direction reversal is the least likely
reading and is not assumed.

**Matched thresholding is a precondition, not a finding.** Before any confound is blamed, both arms
must apply identical multiple-testing correction. The "37 survivors" rotation cutoff and the FOV
BH-FDR are not obviously matched, so part of the disagreement could be a pure thresholding artifact
(lenient per-gene rotation quantile vs pooled BH). Every comparison below reports both arms under a
single agreed correction. The original 37-gene cutoff is recoverable from `pv_rotation_null.R` and
must be read out and stated before Step 1, so the port is faithful.

## Primary readout: effect-size stability, not a p-value verdict

This is n=1: one block, one animal, descriptive by project convention (report effect sizes, not
p-values). The primary readout is therefore **whether the 37-gene effect sizes are stable or
collapse toward zero as each confound is toggled**, not whether a gene "passes" a null. For each of
the 37 genes, report the peak-valley log2FC with a bootstrap CI (resampling cells/FOVs) under every
confound toggle: retired-vs-atlas labels, old-vs-current QC, per-cell-vs-FOV aggregation,
wide-vs-core bands. A signature that is an artifact of one confound collapses toward zero when that
confound is corrected; a real dose effect holds its magnitude and direction across toggles. This
framing sidesteps the rotation-vs-phase-null philosophical argument for the main verdict and matches
admissible inference for a single animal. The null models (rotation, phase) are retained as a
**secondary** calibration applied only to genes that survive the effect-size-stability screen, never
as the primary gate, and their p-values are read as effect-size-consistency checks, not confirmatory
significance.

## The deeper method question: rotation null vs phase null

The rotation null changes the stripe **angle**, which confounds beam direction with tumor
architecture. If the beam axis aligns with a core-to-edge expression gradient, the real
peak-valley delta is inflated by architecture; the +/-30 deg exclusion of near-angles then
*removes* the rotations that would carry that same architecture, deflating the null and
potentially manufacturing significance. The dose-isolating control instead holds the angle fixed
and shifts only the **phase** (which cells are labelled peak vs valley): a **phase null** preserves
both tumor architecture and the beam axis, varying only dose registration.

Two caveats keep this from being circular. First, the phase null only *out*-performs the rotation
null in the specific regime where tumor architecture aligns with the beam axis; if architecture is
isotropic or orthogonal to the beam, the two nulls agree and the preference is moot. That regime is
an empirical question about *this* block, not an a priori fact, so it is tested (the beam-axis
architecture diagnostic below), not assumed. Second, "the rotation null manufactures significance"
is a hypothesis, not an established fact: it must be *demonstrated* by a synthetic
architecture-only, zero-dose simulation (impose a real core-to-edge gradient along the beam axis
with no periodic dose term, then run the rotation null and measure its false-positive rate on genes
that carry only the gradient). Absent that run, the artifact verdict is asserted, not shown. The
FOV-paired test *approximates* a local phase null (it varies registration while holding the section
fixed), which is a reason to weight it, but it is spatial dilution, not a calibrated phase-null
test, so it does not by itself settle the question and the plan stays agnostic until Step 1 runs.

## Step 0 result (executed 2026-07-13): phase null is futile on this block

`dev/if_channel_check/pv_step0_powerfloor.R` (geometry-only, run on beast) reproduced the contested
stripe fit (theta=-25 deg, spacing 1.039 mm, peak_half 0.269 mm) and settled the gate:

- Tumor spans **4.91 beam periods** along the beam-perpendicular axis, so fully-independent phase
  offsets number ~5 at most.
- Phase-null p-floor is **~0.20 coarse** (1/periods) and **~0.47 optimistic** (peak-set
  decorrelation lag 0.245 mm, i.e. the wide 0.269 mm band barely changes membership as phase
  slides). Both are far above 0.05, so **no gene can clear a phase null regardless of true effect.**
- **Decision: skip the phase null (former Steps 3-3.5 phase-null arms).** The FOV-pseudobulk arm and
  the effect-size-stability screen carry the verdict; the continuous phase regression stays as the
  powerful periodic estimator but is judged on effect size + CI, not a phase-null p.
- FOV pseudobulk units on the real grid: **43 peak / 37 valley = 80 units**, ample for the FOV arm
  and the phase regression.
- Registration: only the 0.119 mm peak-FOV-centroid scatter proxy exists; the true
  H2AX-overlay-to-centroid uncertainty is still unmeasured and must gate any band width used.

The steps below retain the phase-null description for provenance, but the phase-null *execution* is
cancelled by this gate; read Steps 3-3.5 as superseded.

## Step 1 result (executed 2026-07-13): the signature is largely a label/QC artifact

`dev/if_channel_check/pv_step1_effectsize_swap.R` ran the original pooled-per-cell statistic + wide
real bands + rotation null on two cell sets on Block_21, holding geometry fixed.

- **Reproduction (Step 0b, folded in).** On the retired `a` bucket (102,421 cells) the original
  method regenerates **38 survivors** (the reported "~37"; the one-gene difference is filter/seed
  provenance, not decision-critical). The MHC-I/IFN/damage members are all present: B2m, H2-D1,
  H2-K1, Ifitm2, Ifitm3, Ifi27, Oas1a/g, Cdkn1a, Bax.
- **Swap to current atlas-tumor + current QC (57,039 cells).** Survivor count drops **38 -> 7**.
  Only **5 of the 38** original survivors are retained by the new-set rotation null; median survivor
  effect roughly halves (|delta| 0.074 -> 0.042, ratio 0.57).
- **The interferon module specifically collapses.** Ifitm2, Ifitm3, Ifi27, Oas1a/g all fall to
  bootstrap CIs spanning zero. MHC-I is split, not cleanly killed: H2-D1 survives (0.092, CI
  [0.04,0.14]), B2m shrinks but stays peak-up (0.079, CI [0.02,0.13]), H2-K1 collapses (0.042, CI
  [-0.01,0.09]).
- **Not pure noise.** Genome-wide peak-valley delta correlates old-vs-new at r=0.61 (Spearman 0.63),
  and a handful of genes are effect-size-stable across the swap: Clu (0.286 -> 0.281, essentially
  unchanged), Tpt1, Itm2b, H2-D1, Cdkn1a. These are the only real candidates for a residual
  dose-associated signal, and they are stress/secretory/housekeeping-leaning, not the claimed
  IFN/MHC-I inflammatory program.

**Verdict from Steps 0+1:** the 37-gene MHC-I/IFN/damage peak signature as a coherent program is
**substantially a retired-`a`-bucket + old-QC artifact** (the IFN component does not survive; MHC-I
is half-killed), and it cannot be rescued by a phase null (futile per Step 0). A weak, directionally
consistent residual (~5-7 genes, effect halved) remains but must still clear the pseudoreplication
test before any claim: the pooled statistic weights ~57k cells equally, and the FOV-paired test on
this block was null, so these residual survivors are prime pseudoreplication suspects. That is
Step 2 (FOV-pseudobulk on the residual survivors) and Step 3.5 (continuous phase regression at the
FOV unit, judged on effect size + CI, no phase-null p).

## Step 2 result (executed 2026-07-13): pseudoreplication is the killer; signal is null at the FOV unit

`dev/if_channel_check/pv_step2_fov_pseudobulk.R`, same atlas-tumor cells + current QC as Step 1,
same geometry, swapping only the aggregation unit to the FOV.

- **Design A (paired within-FOV, `~ FOV + zone` DESeq2).** 134 FOV-zone pseudobulk units across **67
  paired FOVs** (67 peak / 67 valley). **0 genome-wide survivors at padj<0.05.** Of the 7 Step-1
  residual survivors, **0 are significant**; of 11 IFN/MHC members, **0 have a peak-up CI excluding
  zero** (B2m 0.013, H2-K1 -0.003, H2-D1 -0.016, all CIs spanning 0).
- **Design B (continuous phase regression, 80 FOV units, limma-voom on mean distance-to-peak).**
  **0 genome-wide survivors.** The more powerful periodic estimator finds nothing either.
- **The lone near-miss is Clu**, not the inflammatory program: within-FOV peak-up log2FC 0.228 (CI
  [0.11,0.34]) but **padj 0.084 (n.s.)** and a null phase slope (padj 0.98). It is the one gene that
  was effect-stable across Step 1's label swap and it still does not clear FOV-level correction, so
  it is at most a single-gene lead, not a signature.
- Proliferation genes trend peak-*down* (Mki67 -0.187, Ube2c -0.089; both n.s.), the opposite of a
  damage/arrest read, consistent with the whole thing being architecture + pseudoreplication rather
  than dose biology.

**Confound attribution.** The 37-gene signature dies in two stages: the interferon module and half
of MHC-I are killed by the **cell-label + QC swap** (Step 1: 38 -> 7 survivors), and the entire
residual is killed by **pseudoreplication** (Step 2: 7 -> 0 at the FOV unit). Confounds #1/#2 and #4
each carry part of the artifact; the null model (rotation vs phase) never had to be adjudicated,
because the signature does not survive even the rotation null's own statistic once the cells and the
replicate unit are corrected, and the phase null was independently shown futile (Step 0).

## Plan

**Step 0 (power-floor gate, run first, five minutes).** Before any confound work, compute the
resolving power the phase null can achieve on this block. The block has ~4 stripe periods, so the
number of *independent* phase offsets is a handful; count the achievable distinct offsets and the
FOV/pseudobulk-unit total, and derive the smallest attainable permutation p-value (with only k
independent offsets the floor sits near 1/k, i.e. ~0.25 for four). If that floor exceeds any
biologically plausible per-gene effect, a "survives the phase null" verdict is structurally
unreachable regardless of the truth, and Steps 3-3.5 are archaeology: in that case skip the phase
null entirely and let the effect-size-stability screen (above) plus the FOV-pseudobulk result carry
the verdict. This gate is cheap, decisive for sequencing, and its inputs (offset count, unit count)
are knowable before a single test is run.

**Step 0b (rotation-null reproduction, only if cheap).** Regenerate the old `/tmp` intermediates
from the input parquets and rerun `pv_rotation_null.R` unchanged; confirm it reproduces ~37
surviving genes. If regeneration is more than a quick script (the SBRT-null sidecar
`tumor_pv_mbrt_vs_sbrt_null.tsv` must also be rebuilt), skip: the *ported* version (Step 1) is the
baseline, and exact reproduction is provenance, not decision-critical.

**Step 1: port the rotation null to current inputs, method unchanged.** Same rotation-null
statistic and wide bands, but on current atlas tumor cells + `norm.rds`. Isolates confounds #1+#2
(population and QC). Does the 37-gene signature survive the label/data swap? If it collapses here,
the signature was an artifact of the retired `a` bucket, and the method question is moot. **This is
the cheapest step and may close the question on its own:** collapse here delivers the only
decision-critical deliverable (stop citing the 37-gene signature), and by the plan's own logic
(the FOV-paired test already *is* a local phase null, and it is null) an architecture/label
artifact is the expected outcome. Run Step 1 first and gate the rest on its result rather than
committing to all five steps up front.

**Step 2: swap the statistic, hold the null.** On the ported inputs, replace the pooled-per-cell
mean with FOV-pseudobulk aggregation, still using the rotation null for calibration. Isolates
confound #4 (pseudoreplication). This is the single most likely explanation for the disagreement.
Pin the normalization so the swap is clean: aggregate raw counts to FOV×zone units and apply DESeq2
size factors at the pseudobulk unit (standard pseudobulk re-normalization), rather than carrying the
per-cell log-CPM forward. "Holding normalization fixed" here means both arms of the Step-2
comparison use this same pseudobulk normalization and differ *only* in the aggregation unit; the
per-cell-log-CPM-vs-pseudobulk-size-factor contrast is a separate axis reported alongside, not
conflated into the pseudoreplication test.

**Step 3: swap the null, hold the statistic.** Replace the rotation null with a **phase null**
(fixed real angle, random phase offset within the beam period; recompute the per-gene delta
distribution). Run it under both the pooled-mean and the FOV-pseudobulk statistic. The
pooled-mean-under-phase-null arm (3a) is the single load-bearing novel test in this plan: it is the
one thing the existing FOV-paired test has not already answered, because it isolates whether the
pooled per-cell statistic manufactures signal even under a legitimate spatial null. The dose-isolating
question is: does any gene's real peak-valley delta exceed a null that varies only dose registration,
holding architecture and beam axis fixed?

Two calibration requirements gate interpretation of a null phase-null result, without which "no
survivors" is uninterpretable (could be true absence or a dead/underpowered null):

- **Positive control, anchored to the observed effect.** Spike a synthetic periodic signal onto one
  or more genes at the real beam period, at a magnitude matched to the effect the rotation null
  claims (the observed B2m/H2-K1 peak-valley log2FC, not an arbitrarily large spike), and confirm the
  phase null recovers it with the correct sign. Report a small recovery curve across a magnitude
  range bracketing that observed effect, with an explicit success criterion (correct-direction
  recovery at the claimed magnitude). A spike large enough to pass trivially proves nothing about
  sensitivity at the biologically realistic magnitude, so a downstream null is only an "interpretable
  null" if the control fires *at the observed effect size*. A gene with a known non-beam gradient
  serves as the negative-direction control.
- **Achievable resolution.** State the number of *independent* phase offsets the block geometry
  admits (~4 stripe periods bounds this to a handful), since that integer, not band-width, floors
  the smallest attainable p-value. Report it up front so a null is read as "underpowered to N" not
  "absent."

**Step 4: architecture diagnostic (run early, it is logically prior).** The whole rotation-vs-phase
argument hinges on whether tumor architecture aligns with the beam axis, so this diagnostic gates
the phase-null preference and should run near the front, not last. Two parts:

- **Characterize the architecture.** Project the top tumor expression PCs onto the beam-parallel and
  beam-perpendicular axes and quantify any monotone core-to-edge gradient along the beam axis. This
  says whether the confounding regime (architecture aligned with beam) actually holds for this block.
- **Test periodicity, not axis loading.** A genuine periodic dose response is *by construction* also
  aligned with the beam-perpendicular axis, so a survivor loading on a beam-axis PC does **not**
  discriminate dose from architecture. The discriminating test is whether a survivor's expression is
  **periodic at the beam frequency** (power at the 1.02 mm period in a phase-folded profile or its
  spatial spectrum) versus a **monotone** gradient. Monotone-along-beam implicates architecture;
  periodic-at-beam-frequency implicates dose. Report the periodic-vs-monotone decomposition per
  survivor, not the PC-axis alignment alone.
- **Architecture-only null simulation** (the demonstration the "rotation null manufactures
  significance" claim needs). Simulate expression carrying only the measured core-to-edge gradient
  along the beam axis, with no periodic dose term, and run the rotation null on it. The false-positive
  rate on these gradient-only genes is the direct measure of how much significance the rotation null
  manufactures from architecture; a high rate demonstrates the artifact rather than asserting it.

**Step 3.5 (preferred primary, if Steps 1-3 leave the question open): continuous phase regression.**
Rather than carrying binary peak/valley bands plus a band-width sensitivity sweep through every
step, regress expression on **phase** (signed distance to nearest peak, modulo the beam period) as a
continuous covariate, at the FOV-pseudobulk unit. This dissolves the band-width nuisance axis
entirely, is strictly more powerful for a periodic estimand than a two-band contrast, and the
usual periodicity objection to spatially-variable-gene tests does not apply because the regressor is
phase, not raw distance. If adopted, it becomes the primary resolving analysis and the binary bands
demote to a robustness check. This is a lighter-weight relative of the deferred continuous-SVG
reframing, kept in scope precisely because regressing on phase (not distance) sidesteps that
method's periodicity problem.

**Step 5: verdict + doc.** State which confound (or the phase-null itself) explains why the two
methods disagree, judged primarily on whether the 37-gene effect sizes stay stable or collapse under
the confound toggles, and whether any dose-driven peak/valley signal remains under the phase-null
secondary calibration (if the power floor allowed it to run). Resolve the CONTESTED idea-node to a
finding (survives / artifact / underpowered), update `** Results`, CLAUDE.md, and
`project_mbrt_mechanism_status` (which still asserts the 37-gene signature as real). Band-width and
core-vs-wide are a sensitivity axis throughout, not a separate step.

## Interpretive guards

- Tumor-cell definition must be the current atlas compartment, not the `a` bucket, everywhere
  except the Step 0 reproduction (whose purpose is to match the old result).
- The phase null half-width must be >= the stripe-geometry registration uncertainty; below it,
  the null selects noise (same caveat as the M01 signature plan's core bands). The H2AX-to-cell
  registration uncertainty is not yet measured; measure it (H2AX-stripe-to-centroid offset spread)
  before setting the band, so the guard is enforced with a number rather than assumed satisfied.
- n=1 block throughout: even a phase-null survivor is descriptive for one animal, replication-gated
  on M02 H2AX (the M02 companion plan).
- No dose/sparing language (MBRT mean dose unrecorded).

## What this does NOT do

- Does not touch M02 (no H2AX overlay yet; companion plan).
- Does not revive the full continuous-SVG reframing (`plan-svg-v2.md`) with spectral/Moran's-test
  machinery; if both the phase null and FOV-pseudobulk are null, that SVG escalation is the next
  step but is out of scope here. Note the narrower continuous **phase** regression (Step 3.5) *is*
  in scope: it captures most of the SVG power gain for a periodic estimand without the periodicity
  pitfall, because it regresses on phase rather than raw `dist_to_peak`.
