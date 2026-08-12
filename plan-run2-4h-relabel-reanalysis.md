# Run 2 timepoint relabel and 4h re-analysis — implementation plan

## Context

`build_metadata_xlsx.R` hardcodes all M02 samples as 48h harvest. The yH2AX layout doc shows slides 3-4 are 4h. Fix the labels and re-run the 4h analysis with all three MBRT animals (1 from Run 1, 2 from Run 2) integrated.

Rob's R01 hypothesis (from the Aug 9 email): MBRT creates spatially distinct injury/survival states at 4h — valley cells retain proliferative activity with less checkpoint/arrest signaling relative to peak cells, which show more direct damage/stress. He asks for:

1. **Cross-experiment integration** (not n=1 per condition)
2. **Peak-vs-valley pseudobulk expression matrix** for the DESeq2 table
3. **Three prespecified module scores** — damage/stress, checkpoint/arrest, proliferation — compared between peak and valley tumor cells. Test whether valleys have higher proliferation with lower checkpoint/arrest.
4. **Macrophage biology** (CCL8, antigen-presenting monocyte-derived populations) — relevant to the SPORE, worth flagging in the 4h results

## What the plan does

Fix the mislabeled metadata, then produce two deliverables for Rob:

**Deliverable A (fast):** Arm-level 4h tables — MBRT vs SBRT vs Control across all three integrated animals per treated arm. Both descriptive effect sizes (as in the tables already sent) and DESeq2 pseudobulk with real replicates.

**Deliverable B (gated on H2AX overlay):** Peak-vs-valley module-score analysis on the 4h MBRT sections, testing Rob's three modules (damage/stress, checkpoint/arrest, proliferation) at the FOV-pseudobulk level with animal as the replicate unit.

## Decisions already made

- **Engine cohort split:** MBRT_vs_SBRT at 4h goes into its own 2-condition cohort (no controls) so the DESeq2 script's sample-count floor doesn't block the 3-vs-3 contrast. Vs-Ctrl contrasts run separately as descriptive (n=2 control) plus inferential with caveats.
- **Day-2 at n=2:** Run both descriptive and inferential. Answers whether SBRT fibrosis signal holds after removing the mislabeled 4h slides.
- **Control sections:** Treat as independent. Add a note to verify block-to-mouse mapping with Fazzari.
- **Peak/valley framing:** Exploratory. The M01 single-section analysis found no significant peak-vs-valley differences at the gene level. With 3 animals integrated, power increases. Rob's module-based approach (damage/stress, checkpoint/arrest, proliferation) is the primary test, not single-gene DE.

---

## Phase 0: Verify the true labels (GATE)

Unchanged from draft. Read the layout doc, check RDS metadata for any timepoint field, get Rob's one-line confirmation. Output: `data/sources/mutter02_block_layout.csv`. Note to check block-to-mouse mapping with Fazzari (control provenance), but don't gate the plan on it — treat controls as independent with that caveat.

## Phase 1: Fix the metadata build

**Files:** `scripts/build_metadata_xlsx.R`, `data/metadata.xlsx`, `workflows/data_model.smk`.

1. Replace the inline M02 tribble with a read of the Phase 0 CSV.
2. Remove the hardcoded `timepoint_h = 48`; set from per-block value.
3. Fix condition mapping: `MBRT_4h`/`SBRT_4h` for slides 3-4, keep `_day2` for slides 1-2.
4. Guard: assert M02 timepoints not all-identical; assert every CSV block appears in built samples.
5. Rebuild data model on beast.
6. **Verify:** `samples.tsv` shows s3/s4 at `timepoint_h=4`, s1/s2 at `48`.

## Phase 2: Re-scope comparisons and fix engine policy

**Files:** `config/comparisons.yaml`, `scripts/build_comparisons.R`, `config/engine_params.yaml`.

1. `mutter02_day2` cohort → slides 1-2 only (n=2/arm).
2. New `combined_4h_mbrt_vs_sbrt` cohort: MBRT and SBRT only (no controls), M01+M02, n=3 vs n=3. Passes the engine floor.
3. New `combined_4h_vs_ctrl` entries: descriptive tier + inferential at n=2 (lower `min_samples` to 2 for this cohort, caveat in output).
4. Add registry/engine cross-check: warn when `inference_capable=TRUE` but realized n < engine's `min_samples`.
5. **Verify:** `comparisons.tsv` shows correct n and tiers for all new entries.

## Phase 3: Confirm cell typing invariant

Grep typing scripts for condition/timepoint keys. Expectation: invariant, no retype.

## Phase 3.5: Make producer scripts cohort-aware

**Files:** `scripts/aggregate/composition.R`, `niches.R`, `spatial_mixing.R`, `myeloid_polarization.R`, `substate_composition.R`, `pathway_scores.R`.

These scripts hardcode M02-specific condition levels (`"Control","MBRT_day2","SBRT_day2"`) and `n_samples_per_group = 4L`. Replace with reads from the comparison registry or the pseudobulk colData so they work for a 4h cohort.

Scope: replace hardcoded literals with data-driven values. Not a full refactor.

## Phase 4: Arm-level 4h analysis (Deliverable A)

1. Run MBRT_vs_SBRT through DESeq2 (3v3, own cohort).
2. Run vs-Ctrl descriptive effect sizes + inferential at n=2 with caveats.
3. Run producers: composition/propeller, GSEA, UCell (damage/stress, checkpoint/arrest, proliferation modules per Rob's ask), myeloid M1/M2 + CCL8/macrophage chemokine flagging.
4. Cross-run concordance: M01-4h vs M02-4h (same-timepoint replication, answering Rob's Q1).
5. **Verify:** 4h results_master partition; MBRT_vs_SBRT has DE rows; send tables to Rob.

## Phase 5: Day-2 re-scope

1. Re-run day-2 on slides 1-2 only (n=2/arm). Both descriptive and inferential.
2. Diff against prior n=4: which SBRT fibrosis hits survive vs which were inflated by the mislabeled 4h slides.
3. Caveat all prior "day-2 n=4" claims.
4. **Verify:** day-2 results_master partition at n=2 with written diff.

## Phase 6: Peak/valley on M02 4h slides (Deliverable B, gated on H2AX)

This is exploratory. The M01 single-section analysis found no gene-level peak-vs-valley differences. With 3 animals the test has more power, and Rob's module-based approach is the primary analysis.

1. Confirm M02 H2AX covers slides 3-4. Register beam tracks per slide (don't assume M01 geometry).
2. Call peak/valley zones; QC against H2AX.
3. **Rob's three prespecified modules:** Define damage/stress, checkpoint/arrest, and proliferation gene sets from the CosMx panel. Score per cell (UCell or AddModuleScore). Test peak vs valley in tumor cells at FOV-pseudobulk level, animal as replicate unit (3 MBRT animals).
4. **Rob's secondary tests:** (a) Whether valleys have a greater proportion of tumor cells with high proliferation but low checkpoint/arrest. (b) Whether proliferation relative to checkpoint engagement is higher in valley cells. (c) Whether valley cells retain proliferation at comparable damage/stress levels.
5. Send peak-vs-valley pseudobulk expression matrix to Rob (his Q2).
6. **Verify:** module scores + secondary tests with animal-level replication; comparison to M01 single-section result stated.

## Phase 7: Docs, figures, memory

1. Correct all "n=4" / "no 4h anchor" / "day-2 only" references in CLAUDE.md and plan docs.
2. Regenerate experimental-design figure (M02 now 2×2), sample-quality suite, arm figures.
3. Memory note: timepoint-hardcode bug and the guard now in place.
4. **Verify:** grep for stale M02 claims returns only corrected hits.

## Verification summary

- Phase 1: `samples.tsv` corrected
- Phase 2: `comparisons.tsv` correct n and tiers; engine cross-check passes
- Phase 4: DESeq2 produces non-empty 4h DE table; tables sent to Rob
- Phase 5: day-2 n=2 diff against old n=4 written
- Phase 6: module scores at animal-level replication
- Phase 7: no stale M02 references
