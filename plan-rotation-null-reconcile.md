# Reconcile the contested peak/valley signature: rotation null vs FOV-pseudobulk

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
4. **Statistic.** Rotation null is a **pooled per-cell** log-CPM `rowMeans(peak) - rowMeans(valley)`,
   every cell weighted equally, no aggregation. FOV-pseudobulk sums to FOV×zone units and blocks
   on FOV, so the replicate unit is the FOV, not the cell.

**The signal to explain:** B2m and H2-K1 are peak-*up* in the rotation null but valley-*leaning*
in the FOV test (-0.24, -0.23). That is a sign flip, not a power difference, so at least one of
the four confounds inverts direction, not just significance. The prime suspects are #1 (cell
population) and #4 (pooled-mean weighting).

## The deeper method question: rotation null vs phase null

The rotation null changes the stripe **angle**, which confounds beam direction with tumor
architecture. If the beam axis aligns with a core-to-edge expression gradient, the real
peak-valley delta is inflated by architecture; the +/-30 deg exclusion of near-angles then
*removes* the rotations that would carry that same architecture, deflating the null and
manufacturing significance. The dose-isolating control instead holds the angle fixed and shifts
only the **phase** (which cells are labelled peak vs valley): a **phase null** preserves both
tumor architecture and the beam axis, varying only dose registration. The FOV-paired test is
already a local phase null, which is the a priori reason to trust it over the rotation null.

## Plan

**Step 0 (baseline, only if cheap).** Regenerate the old `/tmp` intermediates from the input
parquets and rerun `pv_rotation_null.R` unchanged; confirm it reproduces ~37 surviving genes. If
regeneration is more than a quick script (the SBRT-null sidecar `tumor_pv_mbrt_vs_sbrt_null.tsv`
must also be rebuilt), skip: the *ported* version (Step 1) is the baseline, and exact
reproduction is provenance, not decision-critical.

**Step 1: port the rotation null to current inputs, method unchanged.** Same rotation-null
statistic and wide bands, but on current atlas tumor cells + `norm.rds`. Isolates confounds #1+#2
(population and QC). Does the 37-gene signature survive the label/data swap? If it collapses here,
the signature was an artifact of the retired `a` bucket, and the method question is moot.

**Step 2: swap the statistic, hold the null.** On the ported inputs, replace the pooled-per-cell
mean with FOV-pseudobulk aggregation, still using the rotation null for calibration. Isolates
confound #4 (pseudoreplication). This is the single most likely explanation for the disagreement.

**Step 3: swap the null, hold the statistic.** Replace the rotation null with a **phase null**
(fixed real angle, random phase offset within the beam period; recompute the per-gene delta
distribution). Run it under both the pooled-mean and the FOV-pseudobulk statistic. This is the
decisive dose-isolating test: does any gene's real peak-valley delta exceed a null that varies
only dose registration, holding architecture and beam axis fixed?

**Step 4: architecture diagnostic.** Project the top tumor expression PCs onto the
beam-perpendicular axis. If the rotation-null "survivors" load on a PC that runs along the beam
axis, the signal is architecture-along-beam, not dose. Directly tests the confound the phase null
is designed to remove.

**Step 5: verdict + doc.** State which of the four confounds explains the sign flip, and whether
any dose-driven peak/valley signal survives the phase null. Resolve the CONTESTED idea-node to a
finding (survives / artifact / underpowered), update `** Results`, CLAUDE.md, and
`project_mbrt_mechanism_status` (which still asserts the 37-gene signature as real). Band-width and
core-vs-wide are a sensitivity axis throughout, not a separate step.

## Interpretive guards

- Tumor-cell definition must be the current atlas compartment, not the `a` bucket, everywhere
  except the Step 0 reproduction (whose purpose is to match the old result).
- The phase null half-width must be >= the stripe-geometry registration uncertainty; below it,
  the null selects noise (same caveat as the M01 signature plan's core bands).
- n=1 block throughout: even a phase-null survivor is descriptive for one animal, replication-gated
  on M02 H2AX (the M02 companion plan).
- No dose/sparing language (MBRT mean dose unrecorded).

## What this does NOT do

- Does not touch M02 (no H2AX overlay yet; companion plan).
- Does not revive the continuous-SVG reframing (`plan-svg-v2.md`), the principled third option
  (regress on `dist_to_peak`, spectral/Moran's test); if both the phase null and FOV-pseudobulk
  are null, SVG is the natural next escalation, but it is out of scope here.
