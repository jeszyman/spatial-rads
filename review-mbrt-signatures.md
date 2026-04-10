# Review Checklist: MBRT Signature Persistence and Replication

Execute this checklist after the autonomous run completes. Each section has specific files to check and criteria for pass/fail.

## 1. Orchestrator health

- [ ] Read `results/signatures/orchestrator.log`
  - All three agents ran? (executor, reviewer, reviser)
  - Exit codes all 0?
  - Total wall-clock time reasonable (<10h)?
- [ ] Read `results/signatures/watchdog.log`
  - Disk never exceeded 90%?
  - No sustained CPU idle (agent hung)?

## 2. Executor output validation

- [ ] Read `results/signatures/executor.log`
  - Any errors, tracebacks, or "Error in" lines?
  - Any "skipping" or "cannot find" warnings?
  - Did all 8 tasks complete?

- [ ] `results/signatures/data_inventory.tsv` exists and makes sense
  - All 4 Mutter_02 slides represented?
  - Cell counts plausible (>10K per slide)?
  - Gene panel overlap >90%?
  - Spatial coordinates present?

- [ ] File counts:
  - `find results/signatures/data -name "*.tsv" | wc -l` — expect >= 10
  - `find results/signatures/plots -name "*.png" | wc -l` — expect >= 12
  - `results/signatures/SUMMARY.md` exists and is >500 words?

## 3. Scientific intention assessment

- [ ] Read `results/signatures/SUMMARY.md`
  - Does it directly answer: "Do peak/valley signatures persist at later timepoints?"
  - Does it directly answer: "Do signatures replicate in Mutter_02?"
  - Does it directly answer: "Are signatures absent in SBRT?"
  - Is the interpretive boundary maintained (scoring, not classification, at >4h)?

## 4. Data quality checks

- [ ] `results/signatures/plots/mutter02_umap_landscape.png`
  - Clusters look reasonable? No single dominant cluster?
  - Cell types labeled?
- [ ] `results/signatures/plots/batch_umap_merged.png`
  - Is dataset (Mutter_01 vs 02) confounded with biology?
  - Was batch correction applied? Was it justified?
- [ ] `results/signatures/mutter02_qc_summary.tsv`
  - QC retention rates plausible (40-90%)?
  - No condition with <1000 cells post-QC?

## 5. Signature analysis checks

- [ ] `results/signatures/data/peak_signature_genes.tsv`
  - No circularity: Cdkn1a/p21 excluded? DDR classification genes excluded?
  - Genes are biologically interpretable?
- [ ] `results/signatures/plots/signature_kinetics_peak.png`
  - Kinetic curves show a trend (decay, persistence, or rebound)?
  - Error bars present? (SD, not SEM)
  - Control baseline flat?
- [ ] `results/signatures/plots/spatial_peak_sig_all_timepoints.png`
  - Same color scale across timepoints?
  - Any visible spatial structure at >4h? (report honestly either way)
- [ ] `results/signatures/plots/sbrt_vs_peak_valley_signatures.png`
  - SBRT distribution distinguishable from peak and/or valley?

## 6. Replication checks

- [ ] `results/signatures/plots/mutter02_signature_violins.png`
  - MBRT vs Control difference visible?
  - Effect direction consistent with Mutter_01?
- [ ] `results/signatures/data/concordance_m01_m02.tsv`
  - Direction concordance >60%?
  - Effect size correlation positive?
- [ ] `results/signatures/plots/concordance_scatter.png`
  - Points cluster along diagonal? Outliers explained?
- [ ] Pseudobulk results (if available)
  - Independent units truly independent?
  - FDR correction applied?
  - Top hits biologically plausible?

## 7. Review quality

- [ ] Read `results/signatures/REVIEW.md`
  - Numbered findings (R1, R2, ...)?
  - Critiques specific and actionable?
  - Did reviewer catch statistical issues?
  - Did reviewer assess the interpretive boundary?
- [ ] Read `results/signatures/reviewer.log`
  - Did reviewer actually read the figures?
  - Any errors in reviewer's own reasoning?

## 8. Reviser quality

- [ ] Read `results/signatures/SUMMARY.md` "Revision Notes" section
  - Every R# addressed?
  - Fixes verified (re-ran analysis, updated figures)?
  - Limitations acknowledged honestly?
- [ ] `git log --oneline` on feature branch
  - Clean commit history (executor commit, reviewer commit, reviser commit)?
  - No junk commits?

## 9. Git diff review

- [ ] `git diff main..HEAD --stat`
  - Only expected files changed?
  - No modifications to existing `dev/peak_valley_analysis/` files?
  - No secrets or large binaries committed?
- [ ] `git diff main..HEAD -- results/signatures/SUMMARY.md`
  - Final SUMMARY.md is coherent and complete?

## 10. Decision

- [ ] **Accept and merge**: Results answer the scientific intention, replication looks solid, limitations acknowledged.
- [ ] **Iterate**: Specific gaps identified → write a focused follow-up plan targeting those gaps, re-run Steps 2-3.
- [ ] **Expand**: Results suggest new questions worth pursuing (e.g., spatial neighborhood analysis, cell-type-specific signatures) → new plan.
