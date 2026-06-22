# Aggregate typing — provenance manifest (Step 1 of plan-aggregate-refactor.md)

> Records exactly what produced the **locked** `results/aggregate/full_labels.parquet`, so the
> A2 re-wire (Steps 3–4) reproduces it and the Step 5 re-validation has a fixed baseline.
> Acceptance for the rebuild is **re-validation of the four Q4 gates + marker recall against
> the baseline below, NOT byte-identity** (scVI-GPU + Leiden-igraph are not byte-reproducible).
> Captured 2026-06-22 at repo HEAD `f4cf34e` (audit before any wiring change).

## 1. Hand-run chain that produced the locked labels

Order grounded against `scripts/aggregate/` + `{AGG}/full/` artifact mtimes (`{AGG}` =
`/mnt/data/projects/spatial-rads/aggregate`):

| # | Script | Env | Input → Output |
|---|---|---|---|
| 1 | `merge.R` | spatial-rads | 20 flank `scored.rds` → `{AGG}/merged.rds` (May 31; **stale vs the Jun-15 re-based scored.rds, but count-invariant** per the verified M01 re-base) |
| 2 | `recover_negprobes.R` (parquet-sourced) | spatial-rads | M01 metadata parquet + M02 raw RDS → `{AGG}/cell_neg.tsv` |
| 3 | `full_export.R` | spatial-rads | `merged.rds`,`cell_neg.tsv` → `{AGG}/full/mtx/{counts.mtx,features,barcodes}`, `{AGG}/full/obs.parquet` |
| 4 | **scVI training — UNTRACKED** | spatial-rads-scvi | `full/mtx/*`,`obs.parquet` → `{AGG}/full/scvi_latent.parquet`, `{AGG}/full/scvi_model/model.pt` (Jun 1; **no committed script / log** — only `pilot_scvi.py` exists. Step 3 authors `full_scvi.py` from §3 below.) |
| 5 | `full_cluster.py` | spatial-rads-scvi | latent + `lineage_markers.yaml` → `cluster_checkpoint.h5ad`, per-res labels/gates/recall, `full_sweep_summary.tsv` |
| 6 | `finalize_tier1.py` | spatial-rads-scvi | sweep (res=3.0 adopted) → `results/aggregate/full_coarse_labels.parquet`, `full_coarse_summary.tsv`, `full_gates.json` |
| 7a | `tier2_immune_subcluster.py` → `tier2_singler.R` → `tier2_immune_rescue.py` | scvi / R / scvi | coarse immune → `{AGG}/full/immune_subtypes_rescued.parquet` |
| 7b | `tier2_stroma_subcluster.py` → `tier2_stroma_ucell.R` → `tier2_stroma_rescue.py` | scvi / R / scvi | coarse stroma → `{AGG}/full/stroma_subtypes_rescued.parquet` |
| 8 | `unify_labels.py` | spatial-rads-scvi | all three tiers → `results/aggregate/full_labels.parquet` (3,277,090 rows; `compartment` + `cell_subtype` + `subtype_source` + `rescued`) |

Helpers (cited, not re-run as DAG nodes): `tier2_detectability.py`, `spatialqm_metrics.R`,
`tier2_marker_check.py`, `bcell_diagnostic.py`.

## 2. Step 1.1 live/dead audit — `merged_typed.rds` consumers (re-point to `merged.rds`)

`merged_typed.rds` carries the **dead** per-cell InSituType `cell_type`. Verified from each
script's `commandArgs`: **all three consumers ignore it** (labels come from
`full_labels.parquet`); they use the object only for raw counts + base per-cell metadata
(`dataset`, `sample_id`, `condition`, `slide_id`, `timepoint_h`) — all present on `merged.rds`.

| Consumer | Uses merged_typed for | Reads InSituType `cell_type`? | Re-point safe? |
|---|---|---|---|
| `pseudobulk_build.R` | counts + meta; overwrites `cell_type` from labels (L25–26) | No | **Yes** |
| `pathway_summary.R` | counts + meta; labels joined from parquet; header already says `<merged.rds>` | No | **Yes** (intended) |
| `panel_coverage.R` | counts + `md$dataset` only; `compartment`/`subtype` from labels | No | **Yes** |

⇒ Deleting `embed_celltype` / `prepare_reference` / `typing_insitutype` + the `merged_typed.rds`
output (Step 3) drops no live information. `recover_negprobes`→`cell_neg.tsv` is **live** (feeds
the scVI export) and is retired-not-deleted (re-based onto the M01 RDS in Step 3).

## 3. scVI hyperparameters (recovered from the locked `scvi_model/model.pt`)

Model init (`attr_dict.init_params_`): `n_hidden=128, n_latent=30, n_layers=2,
dropout_rate=0.1, dispersion="gene", gene_likelihood="nb", use_observed_lib_size=True,
latent_distribution="normal"`. `setup_anndata`: `batch_key="slide_id",
continuous_covariate_keys=["neg"]` (no labels / categorical / size-factor keys). Training:
**50 epochs** (per plan v2.5; not stored in `model.pt`). `full_scvi.py` (Step 3) must encode
exactly these.

## 4. Seeds

`full_cluster.py` `SEED=0`; iLISI subsample `N_SUB=50000`; Leiden `flavor="leidenalg",
n_iterations=2`. **scVI training seed was not journaled** for the locked run (reproducibility
gap this manifest closes): `full_scvi.py` pins `scvi.settings.seed = 0` going forward.

## 5. Environment versions (captured 2026-06-22)

- **spatial-rads-scvi** (tier-1): python 3.11.15 · scvi-tools 1.4.2 · torch 2.10.0 (CUDA 12.9,
  GPU available) · scanpy 1.11.5 · anndata 0.12.16 · numpy 2.4.6 · leidenalg 0.12.0 / igraph 1.0.0.
- **spatial-rads** (tier-2): Seurat 5.5.0 · SingleR 2.12.0 · celldex 1.20.0 · UCell 2.14.0 ·
  arrow 22.0.0 · Matrix 1.7.5 · data.table 1.17.8.
- A conda lock for both envs is committed as the frozen-env record at the Step 5 re-validation.

## 6. Re-validation baseline (locked result the rebuild must reproduce within tolerance)

- **Coarse (tier-1):** tumor 1,577,685 (48.1%) / stroma 1,331,601 (40.6%) / immune 367,800
  (11.2%) / 4 unassigned. Scaled iLISI(dataset) **0.67** (floor 0.5, ceiling 1.70 for the 29/71
  split). All four Q4 gates ALL_PASS at res=3.0.
- **Immune (tier-2, rescued):** Macrophage 300,975 / ILC 34,340 / T 16,046 / DC 5,315 / NK 4,991
  / Plasma 3,573 / Mast 1,121 / Neutrophils 259 / epithelial-contam 1,180. B cells near-absent.
- **Stroma (tier-2, rescued):** Fibroblast 732,099 / unassigned 441,227 / SmoothMuscle 72,124 /
  Adipocyte 49,214 / Endothelial 36,937. Pericyte 0 (unresolvable de novo).
- **Differential headline (Step 6 must preserve direction):** day-2 = SBRT-driven stromal
  fibrosis (H3 = 20 hits SBRT_vs_Ctrl, `padj_confirmatory` to ~2e-11); MBRT ~null at
  whole-compartment scale (MBRT_vs_Ctrl 2/207 confirmatory hits). `results_master.tsv` 32,961
  rows (621 confirmatory / 32,340 exploratory).

Re-validation tolerance is a stated band on these (compartment fractions ±~2 pp; gates must
re-pass; headline direction unchanged), reflecting scVI/Leiden non-determinism — not byte-identity.
