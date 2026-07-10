# plan-preprocessing.md — invariant per-sample preprocessing workflow

**Status: as-built (split executed 2026-07-10).** The per-sample pipeline was split into this invariant preprocessing workflow; per-sample cell typing + pathway scoring were **deleted** (canonical typing is cohort-scale in `aggregate.smk`). Supersedes `plan-processing-pipeline.md` and folds the adapter/control facts from `plan-mutter01-controls.md` (both deleted). Companion to `plan-aggregate.md` / `plan-aggregate-refactor.md` (the cohort-scale consumer).

## What it is

`workflows/preprocessing.smk`, tangled from `spatial-rads.org` `**** preprocessing.smk`. A deterministic function of `(raw input, pinned config)`: dataset adapter → 4-criteria QC filter → LogNormalize, plus report-only QC. **Per-sample terminus = `norm.rds`.** ~30-min run at `--cores 8` (23 samples).

## Invariant contract

The per-sample output is a pure function of raw input + pinned config — no cohort-scale or method-laden step. Onboarding a new dataset = write **one adapter** emitting the common-panel Seurat schema; nothing downstream changes. This is why typing/pathway (method-laden, cohort-superseded) were removed from the per-sample side.

## Rules (10)

- **Adapters (dataset-specific):** `adapt_mutter01` (M01 per-slide raw RDS + Yi/Condition parquet join, split by `Condition`), `adapt_mutter02` (M02 per-slide raw RDS, negprobe/falsecode dropped, deterministic k-means y-band region assignment).
- **Shared:** `qc_filter` (4-criteria: min_counts/min_features/max_prop_negative/area_nmads), `normalize` (LogNormalize).
- **Report-only QC (never drops cells/genes):** `samplesheet`, `probe_qc`, `control_qc` (negprobe+falsecode characterization + per-cell `cell_controls.parquet` sidecar + flag-only per-FOV falsecode QC), `contamination_qc` (SpatialQM MECR, sample + FOV grain), `sample_metrics` (SpatialQM per-sample technical QC → `sample_tech_metrics.tsv`: TPC/features/sparsity/SNR/specificityFDR; see `plan-qc-metrics.md`), `qc_report`, `qc_plots`.

## Key decisions (folded from plan-processing-pipeline.md)

- **Probe-QC is a diagnostic, not a filter** — NanoString guidance; dropping low-expression genes discards rare-cell markers (an earlier drop-based version wrongly removed Pecam1/Ms4a1/Ncr1/Cd4/Cd163).
- **QC harmonization** — `propNegative` recomputed; `qcFlagsCell` not used as a filter.
- **Contamination QC report-only** — MECR over `config/lineage_markers.yaml`; per-FOV triage cols separate segmentation contamination from necrosis vs mixed biology. Cohort clean (MECR medians M01 0.030 / M02 0.021; 5/2469 FOVs flagged). FastReseg infeasible (no transcript coords).
- **Thread hygiene** — `shell.prefix("set -euo pipefail; export OMP/OPENBLAS/MKL_NUM_THREADS=1")` (strict mode + BLAS pin) + per-rule `threads: 1`; `--cores N` gives N concurrent single-threaded R jobs.
- **Style-guide compliant** — `D_`/`R_` path constants, `preproc_` rule prefix, `message:` directives, scripts as `input:` deps, standard directive order; `conda run` invocation retained (not `--use-conda`).

## M01 raw-RDS ingestion (folded from plan-mutter01-controls.md)

`adapt_mutter01` ingests the 4 per-slide raw RDS (counts byte-identical to the legacy parquet; invariance verified 2026-06-22). Control sidecar + flag-only per-FOV falsecode QC implemented 2026-06-15. **Deferred, owned by `plan-aggregate-refactor.md` Step 3:** negprobe retirement (blocked — 0% barcode overlap between the raw-RDS negprobes and the scored cells). **Optional future items here:** retain negprobes/falsecode assays on per-sample objects; promote the falsecode-FOV-exclusion flag to a decision.

## Handoff to the aggregate workflow

The aggregate merge reads `norm.rds` (re-pointed from `scored.rds` 2026-07-10 across `merge.R`/`merge_pilot.R`/`coords_necrosis.R`/`recover_negprobes.R` + the `aggregate.smk` expands) and re-types all cells at merged scale (canonical `full_labels.parquet`). The vestigial per-sample `cell_type` requirement in the merge/export was removed and the integration-QC benchmark label nulled (`full_export.R`) — **re-base that benchmark on the merged coarse labels at the next aggregate rebuild.**

## Files

- Literate source: `spatial-rads.org` `*** Initial processing and QC` → `**** preprocessing.smk` + per-script `**** <name>.R`.
- `config/config.yaml` (qc, normalize); `config/lineage_markers.yaml` + `config/pathway_gene_lists.yaml` (shared with `contamination_qc` / `data_model` gene-sets — retained).
- Outputs: `{datadir}/processing/{raw,qc,norm}/*.rds` + `results/processing/*.tsv`, `plots/`.
