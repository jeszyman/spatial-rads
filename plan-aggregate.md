# `aggregate.smk` — design, decisions, status

Cross-sample / cross-dataset analysis workflow consuming the 23 per-sample
scored RDS produced by `processing.smk`. Reads the flank cohort (20 of 23
samples), produces an integrated cell-type embedding plus a Banksy spatial-
niche embedding, and runs three parallel analysis tracks (composition,
cell-type-resolved expression, spatial structure) under the discipline that
the pipeline emits richly-annotated tables rather than pre-filtered results.

This file supersedes the earlier handoff "plan for a plan." Companion to
`plan-processing-pipeline.md` (the per-sample upstream that aggregate
consumes). Brainstormed and devils-advocate-reviewed 2026-05-31.

## v2.0 — Atlas re-typing precedes this workflow (2026-05-31)

Every `cell_type`-keyed stage here (Stage 1a labels, composition, pseudobulk DE,
GSEA, pathway, per-cell DE, M01/M02 concordance) consumes per-sample labels that
are being **replaced** by the v2.0 external-atlas re-typing — all cells typed
directly against the **NanoString CellProfileLibrary mouse mammary atlas** +
ImmGen, no M01 -> M02 transfer (see `plan-processing-pipeline.md` Atlas stem).
**Consequence:** any cell-type output already produced here is on old labels and
must re-run after re-typing. Stage 1a's "labels remain per-sample (Yi-ImmGen /
TransferData)" note is v1 and will change. The **Option C** integrate-then-transfer
typing path is superseded by the external-atlas approach.

---

## Goal

Produce reproducible cross-condition (MBRT vs SBRT vs Control) analyses for
the flank cohort, with formal statistical inference on the n=4 replicated
M02 day-2 stratum and honest descriptive analysis on the n=1 M01 timecourse
stratum, joined by an explicit cross-dataset concordance lens. Best-practice
CosMx aggregation: integrated embedding for visualization and clustering QC,
Banksy niche embedding for spatial structure, all R-native.

Peak/valley analysis is explicitly **out of scope** (stays in
`dev/peak_valley_analysis/`); the M01 tongue cohort is **out of scope**
(separate biology, n=1 throughout).

---

## Cohort

20 flank samples from `processing/scored/`:

| Group | n | Notes |
|---|---|---|
| M01 flank timecourse | 8 | 1h, 4h, day2, day6 × {MBRT, SBRT}; single Control sam0001 at t=0 |
| M02 flank day2 | 12 | 4 slides × {MBRT, SBRT, Control}, n=4 per condition |

Total ≈ 3.3M cells × 950 genes (common panel). M02 day2 is the only
inference-bearing stratum (replicates); M01 timecourse is descriptive only.

---

## Workflows

```
metadata.xlsx
  -> data_model.smk        -> data/data_model.rda + samples.tsv
  -> processing.smk        -> 23 × scored/{sample}.scored.rds        [done]
  -> aggregate.smk         -> merged/integrated object + 3 analysis tracks  [this plan]
```

Heavy intermediates: `/mnt/data/projects/spatial-rads/aggregate/`. Small
tabular outputs and plots: `results/aggregate/`. Snakemake driver in
`basecamp`; R rules via `conda run -n spatial-rads Rscript`. BLAS pinned
at workflow level (`OMP/OPENBLAS/MKL_NUM_THREADS=1` in `shell.prefix`);
every R rule declares `threads: 1`.

---

## `aggregate.smk` stages

### Stage 0 — Merge

**`merge_pilot.R`** (new — added per devils-advocate critique #2)

Memory budget check before the full merge. Loads 5 M01 + 5 M02 scored RDS
(small-to-medium sized samples), runs Seurat v5 `merge()`, records peak RSS
via `Rprof(memory.profiling=TRUE)` or `gc()` snapshots. Emits
`results/aggregate/merge_pilot_memory.tsv` with: n_samples, n_cells,
peak_rss_gb, final_object_gb. Decision rule: if peak RSS scaled to 20
samples projects under 100 GB, use single-pass merge; otherwise fall back to
per-dataset-then-cross.

**`merge.R`**

Single-pass Seurat v5 `merge()` on all 20 flank scored.rds, sparse-only,
no densification. Carries all per-cell metadata + spatial coordinates +
existing pathway score columns. Validates FK against `samples.tsv`.

- Heavy output: `/mnt/data/projects/spatial-rads/aggregate/merged.rds`
  (~3.3M cells × 950 genes, sparse).
- Test artifact: `results/aggregate/merge_summary.tsv` —
  sample_id × {n_cells, n_genes_detected, mean_counts, max_counts,
  frac_negative_probes, ram_peak_gb}; final row aggregates totals.

Fallback (plan B if pilot blows budget): `merge_m01.R` + `merge_m02.R` +
`merge_cross.R`. Documented in code comments; not pre-tangled.

### Stage 1a — `embed_celltype.R`

Standard cell-type integration. Skips HVG selection (curated 950-gene panel
is already biologically selected).

- `ScaleData(features=all_950)` → `RunPCA(npcs=30)`
- `RunHarmony(group.by.vars = c("dataset", "slide_id"), theta=c(2, 2))`
  [updated per devils-advocate critique #1: include `dataset` covariate,
  not slide_id alone]
- `RunUMAP(dims=1:30, n.neighbors=30, min.dist=0.3, reduction="harmony")`
- `FindNeighbors(reduction="harmony", dims=1:30) → FindClusters` Leiden
  resolution 0.5 (production); sweep {0.3, 0.5, 0.8, 1.2} during prototype.
- **Cell type labels remain the per-sample assignments** (Yi-ImmGen for M01,
  Seurat anchor TransferData for M02). Leiden clusters are QC only —
  marker-concordance scored per cluster (Krt8/Epcam → tumor; Cd3e/Cd8a →
  T; Cd19/Ms4a1 → B; etc.).
- Heavy output: `merged_celltype.rds`.
- Test artifact: `results/aggregate/celltype_embed_qc.tsv` — n_pcs,
  harmony_iters, silhouette(dataset)/silhouette(slide_id) pre/post,
  LISI(dataset)/LISI(slide_id) pre/post, n_clusters,
  per-cluster {n_cells, modal_cell_type, marker_concordance_score}.

**Formal LISI(dataset) checkpoint:** downstream tracks must not run if
LISI(dataset) < 1.8 post-Harmony (indicates dataset axis not adequately
mixed; would invalidate the cross-dataset concordance lens). Failure mode
surfaces as an explicit rule failure with message and remediation hint
(re-run with `theta=c(4,2)` or revisit batch covariates).

### Stage 1b — `embed_niche.R`

Banksy spatial niche embedding.

- Convert `merged.rds` → `SpatialExperiment` (one-line via
  `as.SingleCellExperiment` then add `spatialCoords`).
- `runBanksy(lambda=0.8, k_geom=c(6,18), M=1, group="sample_id",
  use_agf=TRUE)` — per-sample neighbor graphs (no cross-sample neighbors).
- `runPCA(npcs=20)` on the BANKSY matrix.
- `RunHarmony(group.by.vars = c("dataset", "slide_id"))` on the BANKSY PCs
  [applies devils-advocate critique #1 here as well].
- Leiden on Harmony-corrected Banksy PCs, resolution 0.4 (production);
  sweep {0.2, 0.4, 0.6, 1.0} during prototype.
- Niche names manually assigned after inspecting
  `niche_celltype_composition.tsv` (a small lookup table in
  `config/niche_names.tsv` keyed by numeric niche_id).
- Heavy output: `merged_niche.rds`.
- Test artifacts:
  - `results/aggregate/niche_embed_qc.tsv` — n_niches, harmony iters,
    silhouette/LISI(dataset, slide_id) pre/post, per-niche {n_cells,
    modal_cell_type, composition_entropy, mean_pathway scores}.
  - `results/aggregate/niche_spatial_qc.tsv` — per-sample bbox of
    (x_slide_mm, y_slide_mm), bbox overlap matrix for samples on the same
    slide [added per devils-advocate critique #3 — confirms per-sample
    subset before k-NN is doing its job on shared-slide M01 samples].

Memory pre-check: expected peak 50–80 GB; run with `--cores 1` (no
concurrent jobs during this rule).

### Track 1 — `composition.R`

Cell-type composition per condition.

- M02 day2 inference: `propeller` (speckle BioC) with design
  `~ 0 + condition + slide_id`. Three contrasts: MBRT_vs_Ctrl,
  SBRT_vs_Ctrl, MBRT_vs_SBRT. Logit transform. BH padj across
  (cell_type × contrast).
- M01 timecourse descriptive: per-sample proportions only, no test.
- No silent drops except on method failure; emit
  `composition_dropped_celltypes.tsv` listing skip reason and the quality
  metrics that drove it.

Outputs:

| File | Schema |
|---|---|
| `composition_by_sample.tsv` | sample_id, cell_type, n_cells, fraction, condition, timepoint_h, dataset, slide_id |
| `composition_test_m02day2.tsv` | cell_type, contrast, log2FC_logit, t_stat, pvalue, padj, method="propeller", n_samples_per_group, mean_n_cells, dataset="M02" |
| `composition_dropped_celltypes.tsv` | cell_type, reason, mean_n_cells |
| `plots/composition_m02day2_bars.png` | stacked bars per sample, faceted by condition |
| `plots/composition_m02day2_forest.png` | per-celltype forest of the three contrasts with 95% CI |
| `plots/composition_m01_timecourse.png` | line plot, fraction vs timepoint, one panel per cell type |

### Track 2 inference — `pseudobulk_build.R` + `deg_pseudobulk.R`

**`pseudobulk_build.R`** (M02 day2 only):

- For each (sample × cell_type): sum raw counts of those cells → one column.
- Build `SummarizedExperiment`: assay = counts; colData = sample_id,
  cell_type, condition, slide_id, timepoint_h=48, n_cells, total_counts,
  mean_libsize; rowData = gene_symbol.
- No filtering at build — quality columns attached for downstream choice.
- Heavy output: `pseudobulk_se.rds`.
- Test artifact: `pseudobulk_qc.tsv` — sample × cell_type × {n_cells,
  total_counts, mean_libsize, frac_genes_with_count_gt_5}.

**`deg_pseudobulk.R`**:

- Per cell type: subset SE, gene-level filter (≥10 counts in ≥4 samples,
  emit drop log), `DESeqDataSetFromMatrix(design = ~ slide_id + condition)`,
  `DESeq()`, three contrasts via `results()` + `lfcShrink(type="apeglm")`.
- BH padj per (cell_type × contrast). Cell-type-level skip only when DESeq2
  cannot fit (zero cells of that type in any condition); emit reason.
- Attach quality columns to every row: n_samples_used,
  min_counts_per_sample, mean_n_cells_per_sample.

Outputs:

| File | Schema |
|---|---|
| `degs_pseudobulk_m02day2.tsv` | cell_type, contrast, gene, log2FC, lfcSE, stat, pvalue, padj, baseMean, n_samples_used, min_counts_per_sample, mean_n_cells_per_sample, dataset="M02" |
| `deg_summary_m02day2.tsv` | cell_type × contrast → n_genes_tested, n_padj_05, n_padj_05_lfc_1, top5_up, top5_down |
| `degs_pseudobulk_skipped.tsv` | cell_type × contrast → reason, n_cells_in_failed_group |
| `plots/volcano_{cell_type}_{contrast}.png` | volcano per (cell_type × contrast) |

### Track 2 descriptive — `deg_percell.R`

M01 flank timecourse only.

- Subset `merged.rds` to M01 flank cells.
- Per (cell_type × contrast × timepoint): `FindMarkers(test="wilcox",
  logfc.threshold=0, min.pct=0, max.cells.per.ident=20000)`
  [downsampling cap added per devils-advocate critique #4 — caps per-group
  scan at 20k cells while preserving statistical resolution].
- 10 contrasts per cell type (MBRT_*h vs Ctrl_0; SBRT_*h vs Ctrl_0;
  MBRT_*h vs SBRT_*h matched at 4h/day2/day6).
- Skip (cell_type × contrast) pairs where either group has <20 cells; log
  to `degs_percell_skipped.tsv`.
- Per-call wall-clock timing recorded in summary artifact for runtime
  visibility.

Outputs:

| File | Schema |
|---|---|
| `degs_percell_m01.tsv` | cell_type, contrast, timepoint_h, gene, log2FC, pct.1, pct.2, pct.diff, cohens_d, n_cells_1, n_cells_2, pvalue, padj_bh, ref_is_baseline_t0, dataset="M01" |
| `deg_percell_summary_m01.tsv` | cell_type × contrast × timepoint_h → n_genes_tested, n_log2FC_gt_1, n_log2FC_lt_neg1, top5_up, top5_down, runtime_seconds |
| `degs_percell_skipped.tsv` | cell_type × contrast → reason, n_cells_1, n_cells_2 |

**p-value caveat documented in `SCHEMA.md`:** Wilcoxon p-values across
millions of cells from n=1 samples are not valid biological inference
(Squair 2021). pvalue/padj_bh kept for ranking and completeness, not
interpretation. Inferential effect-size columns are log2FC, pct.diff, and
cohens_d.

### Track 2 cross-dataset lens — `concordance_m01_m02.R`

Inner-joins M01 day2 results (timepoint_h=48) with M02 day2 pseudobulk
results on `(cell_type, contrast_key, gene)`. Maps:

| M01 contrast | M02 contrast | contrast_key | ref_asymmetric |
|---|---|---|---|
| MBRT_day2 vs Ctrl_0 | MBRT_vs_Ctrl | MBRT_vs_Ctrl | TRUE |
| SBRT_day2 vs Ctrl_0 | SBRT_vs_Ctrl | SBRT_vs_Ctrl | TRUE |
| MBRT_day2 vs SBRT_day2 | MBRT_vs_SBRT | MBRT_vs_SBRT | FALSE |

Outputs:

| File | Schema |
|---|---|
| `concordance_m01_m02_day2.tsv` | cell_type, contrast_key, gene, m01_log2FC, m01_pct_diff, m01_cohens_d, m02_log2FC, m02_padj, sign_agree, m02_sig (padj<.05 bool), ref_asymmetric |
| `concordance_summary_day2.tsv` | cell_type × contrast_key → primary_rho (populated iff contrast_key=MBRT_vs_SBRT), approximate_rho_asymmetric (populated iff ref_asymmetric=TRUE), frac_sign_agree_overall, frac_sign_agree_in_m02_sig, n_m02_sig, n_m02_sig_concordant, top10_M02_sig_with_M01_rank, ref_asymmetric |
| `plots/concordance_{cell_type}.png` | M01 vs M02 log2FC scatter per (cell_type × contrast_key), points colored by `m02_padj`, asymmetric-ref panels labeled |

**Two-tier concordance metric** [added per devils-advocate critique #5]:
`primary_rho` is computed only on `MBRT_vs_SBRT` rows (symmetric refs in
both cohorts — the clean cross-dataset reproducibility metric).
`approximate_rho_asymmetric` is computed only on `ref_asymmetric=TRUE` rows
(time-confounded in M01 — descriptive only). The two coexist as separate
columns; consumers headline `primary_rho`.

### Pathway tracks — `pathway_summary.R` + `gsea.R`

**`pathway_summary.R`:**

- Per-cell UCell + AddModuleScore on:
  - 4 project-priority pathways (TypeI-IFN, TypeII-IFN, DDR, STING from
    `config/pathway_gene_lists.yaml`; tier="primary")
  - MSigDB Hallmark (~50 sets via `msigdbr` package; tier="exploratory")
- Summarize per (sample × cell_type × pathway × score_type): mean, sd,
  median, n_cells.
- M02 day2 inference: `limma` on per-sample means, design
  `~ slide_id + condition`. BH padj across (cell_type × pathway ×
  score_type × contrast).
- M01 timecourse descriptive: per-sample mean only, no test.
- UCell vs AMS concordance: Pearson r per (cell_type × pathway × dataset).
- **No coverage filtering** — quality columns (n_set_genes, n_panel_genes,
  panel_coverage_frac) attached per row; downstream applies thresholds.

Outputs:

| File | Schema |
|---|---|
| `pathway_scores_summary.tsv` | sample_id, cell_type, pathway_name, pathway_source, tier, score_type, mean, sd, median, n_cells, condition, timepoint_h, dataset, slide_id, n_set_genes, n_panel_genes, panel_coverage_frac |
| `pathway_test_m02day2.tsv` | cell_type, pathway_name, pathway_source, tier, score_type, contrast, estimate, se, t_stat, pvalue, padj_bh, n_samples_per_group, n_panel_genes, panel_coverage_frac, dataset="M02" |
| `pathway_ucell_ams_concordance.tsv` | cell_type, pathway_name, dataset, pearson_r, n_samples |
| `plots/pathway_heatmap_m02day2.png`, `plots/pathway_timecourse_m01.png`, `plots/pathway_ucell_vs_ams_scatter.png` | per stratum |

**`gsea.R`:**

- `fgsea` on DESeq2 `stat`-ranked gene lists per (cell_type × contrast),
  against the same 4+50=54 gene sets.
- Output: NES, pvalue, padj_bh, leading-edge genes, leading-edge size,
  n_set_genes, n_panel_genes per (cell_type × contrast × pathway).

Outputs:

| File | Schema |
|---|---|
| `gsea_pseudobulk_m02day2.tsv` | cell_type, contrast, pathway_name, pathway_source, tier, NES, pvalue, padj_bh, leading_edge_genes (comma-sep), leading_edge_size, n_set_genes, n_panel_genes, dataset="M02" |

### Track 3 spatial — `niche_composition.R` + `colocalization.R`

**`niche_composition.R`:** mirrors `composition.R` on Banksy niches.

- M02 day2: propeller on niche fractions, design
  `~ 0 + condition + slide_id`, three contrasts.
- M01 timecourse: descriptive niche-fraction-by-sample only.
- Plus niche-by-celltype composition (defines niche identity) and
  niche-by-pathway mean scores.

Outputs:

| File | Schema |
|---|---|
| `niche_composition_by_sample.tsv` | sample_id, niche_id, niche_name, n_cells, fraction, condition, timepoint_h, dataset, slide_id |
| `niche_celltype_composition.tsv` | niche_id, niche_name, cell_type, mean_fraction, sd_fraction, entropy |
| `niche_test_m02day2.tsv` | niche_id, niche_name, contrast, log2FC_logit, t_stat, pvalue, padj_bh, n_samples_per_group, method="propeller", dataset="M02" |
| `plots/niche_composition_m02day2_bars.png`, `plots/niche_composition_m02day2_forest.png`, `plots/niche_composition_m01_timecourse.png` | per stratum |

**`colocalization.R`:** pairwise cell-cell colocalization via `spicyR`.

- Per-sample spatial graphs from (x_slide_mm, y_slide_mm) (per-sample
  subset enforced upstream of k-NN; spicyR L-function or weighted-pair
  correlation metric per (cell_type_A × cell_type_B × sample)).
- M02 day2 inference: spicyR hierarchical model (slide as random effect).
- M01 timecourse: per-sample values only, no test.

Outputs:

| File | Schema |
|---|---|
| `colocalization_by_sample.tsv` | sample_id, cell_type_A, cell_type_B, coloc_score, n_A, n_B, n_pairs_within_radius, dataset, condition, timepoint_h, slide_id |
| `colocalization_test_m02day2.tsv` | cell_type_A, cell_type_B, contrast, estimate, se, t_stat, pvalue, padj_bh, n_samples_per_group, method="spicyR", dataset="M02" |
| `colocalization_skipped.tsv` | cell_type_A × cell_type_B → reason |
| `plots/colocalization_m02day2_heatmap.png`, `plots/colocalization_m01_timecourse.png` | per stratum |

---

## Key decisions and rationale

1. **Cohort = flank-only (20 samples).** Tongue is different biology, n=1
   throughout, M01-only; would muddy the cell-type landscape figure if
   integrated. Tongue gets its own future aggregate run if needed.

2. **Two-track DE.** M02 day2 pseudobulk DESeq2 is the only honest
   inference path (slide as replicate, n=4); M01 per-cell Wilcoxon is
   descriptive only. **Rejected:** pooled cross-dataset pseudobulk
   (confounds CosMx run batch with replication); uniform per-cell Wilcoxon
   (n=1 inflated-FDR fallacy at scale, per Squair 2021).

3. **Hybrid M01-vs-M02 framing.** M02 runs unbiased on full panel; separate
   `concordance_m01_m02.R` provides the discovery-validation lens.
   **Rejected:** pre-restricting M02 to M01 nominations (loses anything M01
   missed).

4. **Cell type labels canonical from per-sample assignments.** Stage 1a
   clusters are QC-only via marker concordance. **Rejected:** cluster-
   majority relabeling (washes out rare populations); fresh anchor transfer
   on integrated space (extra compute for likely-same answer).

5. **Banksy λ=0.8 only.** Cell-type-aware λ=0.2 mode would duplicate Stage
   1a; skip it.

6. **All-R stack.** Banksy-R + spicyR for spatial methods, no Python /
   AnnData / Squidpy sidecar — avoids dev friction during prototype-fold
   cycles. .h5ad export deferrable to a one-shot if ever needed.

7. **Pipeline emits richly-annotated tables; silent drops only on method
   failure.** Quality columns (panel_coverage_frac, n_samples_used,
   n_cells_per_group, mean_n_cells_per_sample) attached so downstream
   analysis applies its own thresholds. No upstream coverage cutoffs, no
   upstream rare-cell-type cutoffs.

8. **MSigDB Hallmark only for broad pathway scoring** (no KEGG, no
   Reactome, no GO BP). Keeps the exploratory tier focused; can expand
   later if reviewers demand.

### Devils-advocate adjustments (2026-05-31, sonnet review)

The five action items returned by blind review are folded in above:

- **#1 Harmony covariates.** `c("dataset", "slide_id")` in both 1a and 1b
  (uniform; the slide_id-only choice under-corrects the dominant M01/M02
  axis confirmed by InSituType failure in `processing.smk`). LISI(dataset)
  gate ≥1.8 before downstream tracks run.
- **#2 Merge memory pilot.** `merge_pilot.R` precedes `merge.R`; the
  pilot's `merge_pilot_memory.tsv` informs the single-pass vs per-dataset
  decision before committing to the 20-sample run. 15 GB swap already used
  → effective free RAM is ~105 GB, not ~120 GB; Seurat v5 sparse join can
  spike 2–3× final object size.
- **#3 Banksy multi-region slides.** Explicit per-sample subset enforced
  before k-NN; bbox-overlap check on shared-slide M01 samples in
  `niche_spatial_qc.tsv`.
- **#4 deg_percell.R runtime.** `max.cells.per.ident=20000` downsampling +
  per-call wall-clock in the summary artifact.
- **#5 Concordance two-tier metric.** `primary_rho` on symmetric
  MBRT_vs_SBRT only; `approximate_rho_asymmetric` on time-confounded
  vs-Ctrl rows. Headline figure uses `primary_rho`.

---

## Development process

Mirrors `processing.smk`'s two-phase pattern (per `feedback_dev_workflow.md`).

**Per-stage cycle:**

1. Sketch rule directly in `workflows/aggregate.smk` (not org).
2. Prototype `scripts/aggregate/<stage>.R` as standalone with hardcoded
   args; iterate via `conda run -n spatial-rads Rscript ...` on a 1–2
   sample subset.
3. Generalize to `commandArgs(trailingOnly=TRUE)` matching the rule's
   positional shell.
4. Dry-run `snakemake -s workflows/aggregate.smk --dry-run`.
5. Execute the rule; produces canonical object + small test TSV.
6. Verify: `Read` the test TSV, inspect key metrics, view any plots via
   `Read`. Fix and loop back to (2) if wrong.
7. **Fold to org only after the rule runs end-to-end and verifies.** Move
   the rule into a `*** aggregate.smk` block under a new `** Aggregate
   analysis` heading in `spatial-rads.org`; move the R script body into a
   `*** <stage>.R` block with `:tangle scripts/aggregate/<stage>.R`.
   Tangle. `diff` against the working files to confirm byte-identical.
8. Add the stage's outputs to `rule all:` once folded.

**Snakemake hygiene baked from line 1** (per `feedback_smk_thread_hygiene`):

- `shell.prefix("export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1; ")`
- `threads: 1` on every R rule
- Every R script declared as an `input:` to its rule (so edits trigger
  reruns — per `feedback_snakemake_script_tracking`)
- Run command:
  `TMPDIR=/mnt/data/projects/spatial-rads/tmp conda run -n basecamp snakemake -s workflows/aggregate.smk --cores N`

**Output discipline:**

- Heavy intermediates → `/mnt/data/projects/spatial-rads/aggregate/`.
- Small tabular outputs → `results/aggregate/`.
- Plots → `results/aggregate/plots/` (one PNG per panel; not bundled into
  HTML).
- Long-format tables with explicit FKs: `sample_id` joins
  `results/data_model/samples.tsv`; `gene` joins
  `results/processing/common_genes.tsv`.
- A `results/aggregate/SCHEMA.md` declares column meanings, units, FK
  references, and interpretation caveats (e.g., p-value non-inference in
  `degs_percell_m01.tsv`, ref-asymmetry semantics in concordance).

---

## Build order

Stages designed so each can be tested in isolation before downstream
depends on it.

1. `merge_pilot.R` — memory characterization
2. `merge.R` — verify cell count, gene panel, FK integrity
3. `embed_celltype.R` — standard embedding; marker-concordance QC;
   LISI(dataset) gate
4. `composition.R` — cheap, immediate science output, validates merge
   metadata
5. `pseudobulk_build.R` — build SE only; verify counts per
   (cell_type × sample)
6. `deg_pseudobulk.R` — DESeq2 on SE
7. `deg_percell.R` — independent of pseudobulk path
8. `concordance_m01_m02.R` — depends on 6 + 7
9. `pathway_summary.R` — independent of DE
10. `gsea.R` — depends on 6
11. `embed_niche.R` — Banksy + Harmony
12. `niche_composition.R` — analogous to (4), on niches
13. `colocalization.R` — spicyR

---

## Status

Designed and devils-advocate-reviewed 2026-05-31. No code written yet.
Next action: writing-plans skill → implementation plan keyed to this spec.

---

## Key files

| What | Where |
|---|---|
| This spec | `plan-aggregate.md` (this file) |
| Companion upstream spec | `plan-processing-pipeline.md` |
| Per-sample inputs | `/mnt/data/projects/spatial-rads/processing/scored/*.scored.rds` (23) |
| Sample sheet | `results/data_model/samples.tsv` |
| Common gene panel | `results/processing/common_genes.tsv` (950) |
| Project-priority pathway lists | `config/pathway_gene_lists.yaml` |
| Workflow (to be written) | `workflows/aggregate.smk` |
| R scripts (to be written) | `scripts/aggregate/*.R` |
| Heavy intermediates | `/mnt/data/projects/spatial-rads/aggregate/` |
| Tabular results | `results/aggregate/` |
| Plots | `results/aggregate/plots/` |
| Output schema declarations | `results/aggregate/SCHEMA.md` |
| Literate notebook (canonical after fold) | `** Aggregate analysis` heading in `spatial-rads.org` |
