# Differential-layer robustness pass — design spec

**Date:** 2026-06-23
**Status:** design approved (brainstorming); implementation plan pending (`-impl.md`)
**Goal:** Close the confounding gaps a blind devils-advocate panel surfaced in the MBRT/SBRT/Control differential layer, so each reported effect is defensibly *regulation* (or *composition*, correctly labelled) and not an artifact of cell-state mixing, unlabelled cells, detection sparsity, or unrecorded dose. One consolidated rerun of the differential layer at the end.

**Scope:** the differential layer only — now `workflows/aggregate_differential.smk` (split out of the former `aggregate_differential.smk` on 2026-07-12; typing is locked and reproducible in `aggregate_typing.smk`, not touched). Inference cohort is Mutter_02 day-2, n=4/arm, balanced 3-conditions/slide. Mutter_01 stays descriptive (n=1).

---

## Background — the confounding modes

A four-lens blind panel (2026-06-22, logged) returned 3 solid-with-caveats / 1 needs-rework on the differential design. Six findings survived triage as valid-needs-action. The unifying defect: **within-coarse-label pseudobulk DE conflates a shift in cell-state mixture with transcriptional regulation**, and several secondary confounders (unlabelled cells, detection sparsity, unrecorded MBRT dose, a non-independent gene set) are unmodelled. The p21/`Cdkn1a` case made it concrete — an apparent ~2-fold within-fibroblast downregulation in SBRT was a 17%→7% drop in the *fraction* of expressing cells with per-expresser level flat, i.e. resting→myofibroblast state turnover, not regulation.

Each fix below is grounded in field-standard method precedent (citations in the final section), per project rule `feedback_verify_method_precedent`.

---

## Design

### Fix 1 — Unlabelled stroma (33%, 441k cells): model it, test it, sensitivity-check it

**Decision.** Treat `unassigned` as an explicit category, never silently dropped or redistributed. Composition is a closed (sum-to-one) system, so a treatment-associated shift in the unassigned fraction mechanically distorts every labelled type's proportion. Redistribution and further force-recovery are rejected at this sparsity (marker rescue already converged with 33% residual).

**Implementation** (`scripts/aggregate/composition.R`):
1. Run propeller with `unassigned` retained as a real category in the proportion model.
2. Add a primary test of the **unassigned fraction itself** across the three contrasts. If it does not differ by arm, the confound is largely defused; if it does, it becomes an interpretive backdrop for every labelled-type call.
3. Emit composition results **with and without** the unassigned category; flag any labelled-type call whose direction/significance flips between the two as confound-driven.
4. **Verified state of play** (not an assertion to add): each cell carries exactly one `cell_subtype`, so `unassigned` cells form their *own* pseudobulk group and do **not** contaminate any labelled type's pseudobulk — the per-type DE is clean. But `unassigned` is a string label (441k cells) and `deg_pseudobulk.R:45` loops over all `cell_type` values, so it is **already tested in DESeq2** and its (uninterpretable, grab-bag) DE rows are in `results_master`. Fix: tag every `cell_type == "unassigned"` DE row non-interpretable; keep `unassigned` only as the composition category from steps 1–3, never as a reported DE finding.

**Acceptance:** a `composition_unassigned_sensitivity.tsv` showing each labelled-type effect under both inclusion regimes + the unassigned-fraction test.

### Fix 2 — Detection vs level: pseudobulk stays primary, add a per-sample detection test

**Decision.** Keep pseudobulk DESeq2 as the headline DE (the correct multi-sample method). Add, at the n=4 sample level: (a) the detection pair — per-arm % expressing and mean-level-among-expressers — on every DE row; (b) a **limma empirical-Bayes moderated test on arcsin-sqrt-transformed per-sample detection fractions** (the propeller/Phipson variance-stabilized small-n machinery, which borrows strength across genes — this is what makes it calibrated at n=4), run in parallel with the same moderated test on the per-expresser level. A **beta-binomial GLM is rejected**: n=4/arm (12 obs, 3-level factor) leaves too few residual df for a stable dispersion estimate; naïve cell-as-replicate MAST is also rejected (FDR inflation).

**Classifier (sparsity-aware).** At median 58 genes/cell a genuine induction often presents as a 0→1 detection event, so a flat per-expresser level does **not** by itself disprove regulation. Therefore: default class is **`ambiguous`**; assign **`fraction_shift`** only when the per-expresser level is *reliably measurable* (mean ≥ ~3 counts among expressers) **and** flat **and** the detection fraction moved; assign **`regulation`** when the per-expresser level itself moves. p21/`Cdkn1a` qualifies as `fraction_shift` (per-expresser ~4.6–4.9 counts, well above single-count noise, flat; fraction 17→7%); a sparse 0→1 gene lands in `ambiguous`, not `fraction_shift`.

**Implementation:**
1. `pseudobulk_build.R:39` — emit a `pct_expr` matrix alongside the count sum via the existing group indicator: `pos <- as(cnt > 0, "dgCMatrix") %*% G; pct <- sweep(as.matrix(pos), 2, n_cells, "/")`. Store as a second SE assay; also retain per-group mean-among-expressers.
2. `deg_pseudobulk.R` — carry per-arm `pct_expr` and `mean_among_expr` for each gene×cell_type×contrast into the DE output.
3. New `scripts/aggregate/detection_test.R` — two limma-eBayes moderated fits (arcsin detection fraction; per-expresser level) per gene×cell_type×contrast + the sparsity-aware classifier; emit `detection_test_m02day2.tsv`.
4. `assemble_results.R` — add `pct_expr_*`, `mean_among_expr_*`, `detection_padj`, `level_padj`, and `call_class ∈ {regulation, fraction_shift, ambiguous}` to `results_master.tsv`.

**Acceptance:** every DE row carries its detection pair + both moderated p-values; the p21/`Cdkn1a` SBRT row classifies `fraction_shift`; a known sparse gene classifies `ambiguous` rather than `fraction_shift`.

### Fix 3 — Cell-state mixing: pre-registered fibroblast sub-states + propeller (primary); Milo optional

**Decision.** Resolve the fibroblast state-mix confound with a pre-registered sub-state split adjudicated by propeller — **not** Milo as the primary arbiter. All major states already carry tier-2 labels, so cluster-free DA is not needed to *find* states; the one unresolved confound is the resting↔myofibroblast continuum inside the Fibroblast label.

(a) Subcluster fibroblast into resting↔activated/myofibroblast, anchored on `Acta2`, `Tagln`, `Myh11`, core collagens **plus at least one fibroblast-specific activation marker absent from the SmoothMuscle panel — `Postn`, `Fap`, or `Col5a1` — as a required specificity anchor**. Without it, `Acta2`/`Tagln` (which also define SmoothMuscle and sit in the Fibrosis set) score SMC-contaminated cells as "activated"; this is the same specificity-anchor rule already used in the stroma rescue (a lineage-exclusive marker must individually clear the gate). Marker panels **pre-registered in config and committed before any DE is inspected** (guards post-selection / double-dipping); split gated by the detectability test. Myeloid M1/M2 already done. Nothing finer than resting/activated is attempted at this sparsity.

(b) **propeller on the resting/activated sub-states** across the three contrasts is the regulation-vs-composition arbiter: if the activated-fibroblast proportion shifts with SBRT and the collagen/`Acta2` DE tracks that proportion, the program is (partly) composition, not within-state regulation.

(c) **Milo is optional and off the critical path** — a label-free cross-check on a **~200k stratified subsample** (proportional by cell type), design **`~ condition` only** (slide is already regressed out of the scVI latent; adding `slide_id` double-corrects and costs df at n=4, and the balanced design makes slide orthogonal to condition anyway). Its first job is to confirm SpatialFDR is calibrated at our scale (Dann 2022 validated FDR at ~100k cells, not 3.3M) before any larger run. Run only if the propeller adjudication leaves residual ambiguity.

**Implementation:**
1. `config/substate_markers.yaml` (new) — frozen fibroblast resting/activated panels incl. the non-SMC specificity anchor; committed before the rerun.
2. New `scripts/aggregate/substate_split.R` — score + assign fibroblast sub-state, detectability-gated + specificity-anchor-gated; emit a sub-state label column joinable to `full_labels.parquet`.
3. `composition.R` — re-run propeller at sub-state resolution; cross-tabulate the activated-fibroblast proportion shift against the collagen/`Acta2` DE.
4. *(optional, deferred)* `scripts/aggregate/milo_da.py` (env `spatial-rads-scvi`, `milopy`) — ~200k stratified subsample, `~ condition`, SpatialFDR-calibration check + neighbourhood DA.

**Acceptance:** fibroblast resting/activated split passes both gates (detectability + a non-SMC specificity anchor); a propeller sub-state composition test exists; the SBRT collagen/`Acta2` signal is labelled regulation-vs-composition. Milo only if invoked, with its FDR-calibration check reported.

### Fix 4 — Claim scoping (no code, or annotation only)

- **MBRT dose unrecorded.** Every MBRT_vs_SBRT statement is descriptive and labelled dose-confounded; never "sparing." SBRT>MBRT stromal hits are explicitly noted as equally consistent with MBRT delivering lower mean dose. Add an open dosimetry ask to Mutter/Fazzari. (Annotation in `results_master.tsv` contrast metadata + the reanalysis narrative.)
- **M01↔M02.** Stop framing cross-cohort concordance as validation; the cohorts *anti*-correlate on the SBRT stromal signal. Present the actual M01-vs-M02 effect-size scatter, label the disagreement a falsification/caveat to explain (timepoint, cohort construction), not corroboration.

### Fix 5 — STING gene-set independence

**Decision.** STING (3 on-panel genes, 2 shared with the IFN sets) is correlated with the IFN scores by construction. Add a **gene-set overlap diagnostic** (Jaccard / shared-gene matrix across the curated sets) to the pathway outputs; stop counting STING as independent immune-activation evidence in H1. (`pathway_summary.R` or the gene-set reader; emit `geneset_overlap.tsv`.)

### Plumbing — make the rerun cheap and safe

1. `pathway_summary.R:190` empty-facet bug — **already fixed** (`tier == "primary"`, was retired `pathway_source == "project"`).
2. **Split UCell-score compute from plotting** — cache scored cells to a tsv so a plot-layer failure never forces the ~5h UCell recompute again. (`pathway_summary.R` → a `pathway_scores` compute rule + a separate plot rule.)
3. One consolidated `aggregate_differential.smk` rerun of the differential layer after all of the above → refreshed `results_master.tsv` + `detectability_summary.tsv`.

---

## Success criteria

- Composition + DE conclusions reported with explicit unassigned-sensitivity; the unassigned fraction's own differential is tested.
- Every DE row carries its detection pair and a `regulation | fraction_shift | ambiguous` class; p21/`Cdkn1a` lands as `fraction_shift`.
- Fibroblast resting/activated split passes detectability + specificity-anchor gates; a propeller sub-state composition test labels the SBRT collagen/`Acta2` signal regulation-vs-composition (Milo only if propeller leaves residual ambiguity).
- Gene-set overlap diagnostic emitted; STING demoted from independent evidence.
- MBRT_vs_SBRT and M01 framing corrected in outputs/narrative.
- A single rerun regenerates `results_master.tsv`; UCell compute is cached against plot-layer failure.

## Out of scope

- Within-slide / peak-valley spatial MBRT analysis (gated on absent M02 peak/valley staining; the arm-level design cannot resolve a peak-restricted MBRT signal — bounded, not fixed here).
- Re-running the locked typing chain.
- Full-scale Milo / DCATS / CNA. Propeller on the pre-registered sub-states is the primary arbiter; Milo is an optional subsample cross-check only (Fix 3c). Full-cohort cluster-free DA is deferred unless the subsample SpatialFDR check passes and propeller leaves ambiguity.

## Method precedent (evidence base)

- **Composition / unassigned:** Büttner et al. 2021 *Nat Commun* (scCODA, closed-composition coupling); Phipson et al. 2022 *Bioinformatics* (propeller, variance-stabilized small-n); Mangiola et al. 2023 *PNAS* (sccomp); Xenium annotation benchmark 2025 *BMC Bioinformatics* (sparse-panel unassigned + segmentation spillover). Sensitivity-with/without is the procedural standard; no published numeric unassigned cutoff.
- **Detection vs level:** Squair et al. 2021 *Nat Commun* (pseudobulk > cell-level); Soneson & Robinson 2018 *Nat Methods* (detection-rate artifacts); Finak et al. 2015 *Genome Biol* (MAST two-part hurdle); Zimmerman et al. 2021 *Nat Commun* (mixed hurdle = proportion + conditional magnitude); Junttila et al. 2022 *Brief Bioinform* (pseudobulk FDR calibration); Cable et al. 2022 *Nat Methods* (C-SIDE, cell-type-conditioned spatial DE).
- **State granularity / cluster-free DA:** Dann et al. 2022 *Nat Biotechnol* (Milo); Heumos et al. 2023 *Nat Rev Genet* (fixed labels inappropriate on continua → kNN DA; compositionality-aware tests); Luecken & Theis 2019 *Mol Syst Biol* (subclustering valid only when annotation-validated); Zhang et al. 2019 *Cell Systems* (post-clustering inference / double-dipping); Reshef et al. 2022 *Nat Biotechnol* (CNA, scalable alternative); Yu et al. 2023 *Genome Biol* (DA-method benchmark: Milo/CNA robust).

---

## Results — robustness pass (2026-06-23)

*Detection/sub-state/sensitivity outputs are final on disk; the assembled master is being refreshed by the consolidated rerun. DE-derived findings are stable (counts assay unchanged).*

**Detection-vs-level decomposition** (per cell type, n=4, limma-eBayes on arcsin per-sample fractions; BH within cell_type×contrast; symmetric level floor = 2 counts):
- SBRT_vs_Ctrl: **16 fraction_shift, 6 regulation**; MBRT_vs_Ctrl: **0 significant** (all ambiguous); MBRT_vs_SBRT: 4 regulation, 3 fraction_shift.
- **The SBRT collagen/contractile program is `fraction_shift`** — Col1a1, Col1a2, Col4a1, Col5a1, Acta2 across fibroblast/adipocyte/endothelial show significant *detection* increases with non-significant per-expresser level. The strong DESeq2 collagen p-values (≤1e-7) are the pseudobulk shadow of a cell-state shift, **not per-cell upregulation**.
- **SBRT tumor = proliferative skew:** Ccnd1 (+0.53 level), Hmgb2 (+0.44), Mki67 (det+0.13, padj 0.037) up as `fraction_shift`; Cdkn1a/p21 marginally down (det/level −0.13, mae 1.4 → `ambiguous`). Coherent with cycling regrowth of p53-mutant 4T1 survivors at 48h, not arrest (the p53→p21 axis is non-functional, so p21 induction was never expected).
- **p21/Cdkn1a in SBRT stroma reclassified `ambiguous`** — detection-significant but per-expresser level too sparse (~1.4 counts) to separate fraction from regulation. The earlier "definite fraction shift" rested on pooled log-normalized values that overstated it.

**Fibroblast sub-state** (pre-registered split; gate PASSED — 313,088 activated / 419,011 resting; specificity anchor Col5a1 enriched 2.6× in activated, confirming fibroblast-specific myofibroblasts, not Acta2-driven SMC bleed):
- Activated-fibroblast proportion **trends up in SBRT** (+0.51 logit, ~1.4×) with resting flat (−0.09), but **not FDR-significant at n=4** (padj 0.59). The fibrosis signal is a directionally-coherent myofibroblast-expansion trend, underpowered to confirm in this cohort.

**Confounder checks:**
- 33% unassigned stroma does **not** confound composition: its own fraction is arm-stable (SBRT −0.83 logit, padj 0.68) and no labelled-type call flips meaningfully with vs without it (only Endothelial MBRT_vs_SBRT, both null).
- STING demoted from independent H1 evidence: 2 of its 3 on-panel genes (Isg15, Cxcl10) are shared with the Type-I/II IFN sets.

**Net reframe:** the day-2 signal is **SBRT-driven fibroblast activation + tumor proliferative regrowth, both read as cell-state composition shifts**; MBRT is null at whole-compartment scale (the peak-restricted signature is spatially diluted, unrecoverable in this cohort). The strong original DESeq2/pathway hits are **recontextualised as composition, not overturned**.
