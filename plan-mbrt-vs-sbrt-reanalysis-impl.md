# MBRT vs SBRT reanalysis — implementation plan

> **For agentic workers:** implement task-by-task; each task ends with an independently
> verifiable deliverable. Steps use checkbox (`- [ ]`) syntax. Spec:
> `plan-mbrt-vs-sbrt-reanalysis.md`. This is a bioinformatics analysis (R + snakemake), so
> "verify" = dry-run / output-schema / value-check / figure-inspection, not pytest.

**Goal:** Refresh the arm-level MBRT/SBRT/Control differential on the reproducible labels and add
the devils-advocate fixes — symmetric panel-coverage, a magnitude-floored detectability view, and
an honest effect-size figure set across all three contrasts.

**Architecture:** The differential layer already exists in `workflows/aggregate.smk` and emits
`results/aggregate/results_master.tsv` via `scripts/aggregate/assemble_results.R`. We (a) verify
prereqs + regenerate it on the rebuilt labels, (b) extend `assemble_results.R` with a coverage
join + magnitude-floor `trend_call` + a `detectability_summary.tsv`, (c) add two figure scripts.

**Tech stack:** R (`data.table`, `ggplot2`), Snakemake (driver in `basecamp` env; R rules via
`conda run -n spatial-rads Rscript`). Heavy intermediates under
`/mnt/data/projects/spatial-rads/`; small outputs + plots under `results/aggregate/`.

## Global Constraints (verbatim from the spec; every task inherits these)

- All three contrasts (`MBRT_vs_Ctrl`, `MBRT_vs_SBRT`, `SBRT_vs_Ctrl`) are first-class;
  significance is shown but **never gates** what is reported.
- **Effect-size-forward**: every readout reports effect + 95% CI + n=4 MDE.
- **Symmetric panel coverage**: every gene-set's coverage fraction is shown (no set silently
  protected); coverage is from `results/data_model/gene_set_panel_coverage.tsv`.
- **Magnitude floor**: a direction is a "trend" only when `|effect| >= MDE_FLOOR * mde`
  (`MDE_FLOOR = 0.5`); below that it is `below-floor`, not a signal.
- **Dose language gated**: never the word "sparing"; MBRT-vs-SBRT is labeled dose-confounded.
- **M01 is descriptive only** (n=1); never a falsification/coherence arm.
- **No volcano plots** (sparse ~950-gene panel).
- Labels come from `results/aggregate/full_labels.parquet` (reproducible, bit-identical rebuild);
  gene-sets from `results/data_model/pathway_sets.tsv` + `gene_set_panel_coverage.tsv`.
- Snakemake: `TMPDIR=/mnt/data/projects/spatial-rads/tmp`, BLAS pinned, `threads: 1` on R rules.

## File Structure

- Modify `scripts/aggregate/assemble_results.R` — add coverage join (Task 3) + magnitude-floor
  `trend_call` and `detectability_summary.tsv` (Task 4).
- Modify `workflows/aggregate.smk` — `assemble_results` rule: add the coverage input + the new
  output (Task 4); add the two figures to `rule all` (Task 7).
- Create `scripts/aggregate/fig_program_panel.R` — a priori 8-program effect-size forest (Task 5).
- Create `scripts/aggregate/fig_contrast_effects.R` — 3-contrast effect view + SBRT fibrosis (Task 6).

---

### Task 1: Verify Step-6 prerequisites + dry-run the differential targets

**Files:** none modified (verification only).

- [ ] **Step 1: Confirm the three old prereqs are closed.**
  Run:
  ```bash
  cd /home/jeszyman/repos/spatial-rads
  grep -n "full_labels\|cell_type_atlas" scripts/aggregate/pseudobulk_build.R   # expect full_labels, NOT cell_type_atlas
  grep -n "x_slide_mm\|coords" scripts/aggregate/coords_necrosis.R              # coords producer exists
  ```
  Expected: `pseudobulk_build.R` reads `full_labels.parquet` (already confirmed); `coords_necrosis`
  produces `coords_necrosis.parquet`.
- [ ] **Step 2: Dry-run the differential chain to `results_master.tsv`.**
  Run:
  ```bash
  TMPDIR=/mnt/data/projects/spatial-rads/tmp conda run -n basecamp \
    snakemake -s workflows/aggregate.smk results/aggregate/results_master.tsv \
    --dry-run --rerun-incomplete 2>&1 | tail -30
  ```
  Expected: clean DAG, ends "This was a dry-run", no MissingInputException; the differential rules
  (composition, pseudobulk_build, deg_pseudobulk, gsea, pathway_summary, build_gene_sets,
  panel_coverage, power_mde, coords_necrosis, celltype_qc, niches, spatial_mixing,
  myeloid_polarization, concordance_m01_m02, assemble_results) appear. **If any arg desync errors,
  fix the offending rule's shell args to match the current script signature before proceeding.**

### Task 2: Regenerate `results_master.tsv` on the reproducible labels (Step 6 run)

**Files:** none modified (pipeline run).

- [ ] **Step 1: Run the differential chain (background; GPU not needed).**
  ```bash
  TMPDIR=/mnt/data/projects/spatial-rads/tmp conda run -n basecamp --no-capture-output \
    snakemake -s workflows/aggregate.smk results/aggregate/results_master.tsv \
    --cores 8 --rerun-incomplete 2>&1 | tee logs/aggregate/step6_refresh.log
  ```
- [ ] **Step 2: Verify the refreshed table reproduces the known result.**
  ```bash
  conda run -n spatial-rads-scvi python -c "
  import pandas as pd; d=pd.read_csv('results/aggregate/results_master.tsv',sep='\t',low_memory=False)
  print('rows', len(d))
  c=d[d.tier=='confirmatory']; c['hit']=c['padj_confirmatory']<0.05
  print(c.groupby('contrast')['hit'].sum())"
  ```
  Expected (direction must hold): `SBRT_vs_Ctrl` ≫ `MBRT_vs_SBRT` > `MBRT_vs_Ctrl` confirmatory
  hits (~27 / ~20 / ~2); SBRT-driven H3 stromal hits present. Commit nothing yet.

### Task 3: Symmetric panel-coverage join in `assemble_results.R`

**Files:** Modify `scripts/aggregate/assemble_results.R`; Modify `workflows/aggregate.smk`
(`assemble_results` rule input).

**Interfaces — Produces:** `results_master.tsv` gains 3 columns on every `pathway`/`gsea` row:
`n_panel` (int), `n_set_total` (int), `panel_cov_frac` (numeric, `n_panel/n_set_total`).

- [ ] **Step 1: Add the coverage arg + read.** After line 28 (`...sets_p <- a[9]; out_p <- a[10]`),
  shift `out_p` to `a[11]` and add `cov_p <- a[10]`. After line 39 (`gs_dt <- fread(sets_p)`), add:
  ```r
  cov_dt <- fread(cov_p)                       # set, tier, source, n_total, n_panel, usable, thin
  setnames(cov_dt, c("n_total","n_panel"), c("n_set_total","n_panel"), skip_absent = TRUE)
  ```
- [ ] **Step 2: Join coverage onto the master before write.** Immediately before
  `setcolorder(...)` (line 134), add:
  ```r
  master[, `:=`(n_panel = NA_integer_, n_set_total = NA_integer_, panel_cov_frac = NA_real_)]
  iset <- master$readout_class %in% c("pathway","gsea")
  mm <- cov_dt[match(master$feature[iset], cov_dt$set)]
  master[iset, `:=`(n_panel = mm$n_panel, n_set_total = mm$n_set_total,
                    panel_cov_frac = round(mm$n_panel / mm$n_set_total, 3))]
  ```
  and add `"n_panel","n_set_total","panel_cov_frac"` to the `setcolorder` tail.
- [ ] **Step 3: Wire the rule input + arg.** In `workflows/aggregate.smk` `assemble_results` rule,
  add input `cov = "results/data_model/gene_set_panel_coverage.tsv"` and insert `{input.cov}` in the
  shell **before** `{output.master}` (matching the new positional `cov_p <- a[10]`, `out_p <- a[11]`).
- [ ] **Step 4: Verify.** Re-run `assemble_results` only
  (`snakemake ... results/aggregate/results_master.tsv --force --dry-run` then run), then:
  ```bash
  conda run -n spatial-rads-scvi python -c "
  import pandas as pd; d=pd.read_csv('results/aggregate/results_master.tsv',sep='\t',low_memory=False)
  p=d[d.readout_class=='pathway']
  print(p[['feature','panel_cov_frac']].drop_duplicates().sort_values('panel_cov_frac').head(6))"
  ```
  Expected: TypeI_interferon ~0.37, STING ~0.67, Hypoxia ~0.38 all shown — **every** primary set
  carries a coverage fraction (none NA among primary sets).

### Task 4: Magnitude-floor `trend_call` + `detectability_summary.tsv`

**Files:** Modify `scripts/aggregate/assemble_results.R`; Modify `workflows/aggregate.smk`.

**Interfaces — Produces:** `results_master.tsv` gains `clears_mde` (logical) + `trend_call`
(`up`/`down`/`below-floor`/`na`); new file `results/aggregate/detectability_summary.tsv`
(readout_class × contrast × hypothesis → n, n_clears_mde, n_trend_up, n_trend_down).

- [ ] **Step 1: Add the columns** after the `abs_effect_lt_mde` line (132):
  ```r
  MDE_FLOOR <- 0.5
  master[, clears_mde := is.finite(mde) & abs(effect) >= mde]
  master[, trend_call := fifelse(!is.finite(mde) | abs(effect) < MDE_FLOOR * mde, "below-floor",
                          fifelse(effect > 0, "up", "down"))]
  master[is.na(effect), trend_call := "na"]
  ```
  Add `"clears_mde","trend_call"` to `setcolorder`.
- [ ] **Step 2: Emit the detectability summary** before the final `cat(...)` inspect block:
  ```r
  det <- master[, .(n = .N, n_clears_mde = sum(clears_mde, na.rm = TRUE),
                    n_trend_up = sum(trend_call == "up"), n_trend_down = sum(trend_call == "down")),
                by = .(readout_class, contrast, hypothesis, tier)]
  fwrite(det, sub("results_master", "detectability_summary", out_p), sep = "\t")
  ```
- [ ] **Step 3: Declare the new output** in the `assemble_results` rule
  (`detect = "results/aggregate/detectability_summary.tsv"`); the script derives its path from
  `out_p`, so no extra arg.
- [ ] **Step 4: Verify.** Re-run; then:
  ```bash
  column -t -s$'\t' results/aggregate/detectability_summary.tsv | grep -E "confirmatory|contrast" | head
  ```
  Expected: rows for all 3 contrasts; `SBRT_vs_Ctrl` H3 has `n_clears_mde > 0`; `MBRT_vs_Ctrl`
  confirmatory `n_clears_mde` small (≈0–2) — the diluted-but-shown MBRT result.

### Task 5: A priori program panel figure — `fig_program_panel.R`

**Files:** Create `scripts/aggregate/fig_program_panel.R`; Test output
`results/aggregate/plots/program_panel.png`.

**Interfaces — Consumes:** `results_master.tsv` (`readout_class=='pathway'`, primary sets);
**Produces:** one PNG.

- [ ] **Step 1: Write the script.** Effect (limma estimate) + 95% CI per program × cell type,
  faceted by the 8 primary programs, one colored point/interval per contrast; x-label annotated
  with `panel_cov_frac` per program; a vertical 0 line; points whose `trend_call=="below-floor"`
  drawn hollow/greyed. Use `ggplot2` `geom_pointrange(position=position_dodge)`, `facet_wrap(~feature)`,
  `coord_flip()`. Title carries no "sparing" language. Args: `<results_master.tsv> <out.png>`.
- [ ] **Step 2: Run.**
  ```bash
  conda run -n spatial-rads Rscript scripts/aggregate/fig_program_panel.R \
    results/aggregate/results_master.tsv results/aggregate/plots/program_panel.png
  ```
- [ ] **Step 3: Inspect the PNG (required).** `Read` `results/aggregate/plots/program_panel.png`;
  confirm: 8 program facets, 3 contrasts distinguishable, CIs visible, coverage fractions shown,
  below-floor points de-emphasized, no clipping/overplotting. Fix + re-render if any fail.

### Task 6: Contrast effect-size + SBRT fibrosis figure — `fig_contrast_effects.R`

**Files:** Create `scripts/aggregate/fig_contrast_effects.R`; Test outputs
`results/aggregate/plots/contrast_effects.png`, `results/aggregate/plots/sbrt_fibrosis.png`.

**Interfaces — Consumes:** `results_master.tsv`; **Produces:** two PNGs.

- [ ] **Step 1: Write the script.** (a) `contrast_effects.png`: per compartment, the distribution of
  DE effect sizes for the three contrasts (e.g. ridgeline / boxplot), showing SBRT_vs_Ctrl shifted
  vs the near-zero MBRT_vs_Ctrl — the honest "where the signal is" view. (b) `sbrt_fibrosis.png`:
  effect + CI of the Fibrosis_remodeling / Stromal_stress genes in stroma across the 3 contrasts.
  No volcano. Args: `<results_master.tsv> <out_contrast.png> <out_fibrosis.png>`.
- [ ] **Step 2: Run** (analogous `Rscript` call).
- [ ] **Step 3: Inspect both PNGs (required).** `Read` each; confirm the SBRT > MBRT fibrosis
  pattern is visible and the MBRT-vs-Ctrl near-zero shift is honestly shown; axes fit data; no
  clipping. Fix + re-render if needed.

### Task 7: Wire figures into the workflow + close out

**Files:** Modify `workflows/aggregate.smk`.

- [ ] **Step 1: Add `program_panel` + `fig_contrast_effects` rules** (R rules, `threads: 1`,
  consuming `results_master.tsv`), and add the 3 new PNGs + `detectability_summary.tsv` to
  `rule all`.
- [ ] **Step 2: Dry-run.**
  ```bash
  TMPDIR=/mnt/data/projects/spatial-rads/tmp conda run -n basecamp \
    snakemake -s workflows/aggregate.smk --dry-run 2>&1 | tail -20
  ```
  Expected: clean; the new figure rules appear.
- [ ] **Step 3: Commit.**
  ```bash
  repo-commit "feat(mbrt-reanalysis): symmetric coverage + detectability + honest contrast figures" \
    scripts/aggregate/assemble_results.R scripts/aggregate/fig_program_panel.R \
    scripts/aggregate/fig_contrast_effects.R workflows/aggregate.smk \
    results/aggregate/detectability_summary.tsv
  ```

---

## Out of this plan (follow-ons)

- **Prose writeup** of the honest narrative (SBRT fibrosis + MBRT spatial-dilution reframe + the two
  gating inputs) — a `sci-write` task, not implementation.
- **M02 peak/valley spatial** — blocked on staining (gating input); separate plan when it arrives.
- **MBRT dose physics** — lab input from Fazzari/Mutter; unblocks pattern-vs-magnitude.
- The `umap` integration-QC rule (refactor task #8) — wire alongside Task 7 if doing the refactor
  closeout in the same pass.

## Self-review

- Spec coverage: 3-contrast effect-size reporting (Tasks 2,5,6) ✓; symmetric coverage (Task 3) ✓;
  detectability table + magnitude floor (Task 4) ✓; SBRT-positive figure (Task 6) ✓; M01 descriptive
  (no task elevates it) ✓; reproducible-labels substrate (Tasks 1–2) ✓; gating inputs documented
  (out-of-plan) ✓. The decision-gate is realized as the detectability_summary + trend_call, not a
  binary label — matches the revised spec.
- No pytest TDD: intentional — this is an analysis pipeline; verification is dry-run/value/figure.
