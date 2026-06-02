# spatial-rads

Spatial transcriptomics analysis of microbeam radiation therapy (MBRT) using CosMx single-cell data from the Mutter Lab (Mayo Clinic).

## Scientific Context

First-in-field spatial transcriptomics study of MBRT peak/valley biology. MBRT delivers dose through ~1mm-spaced parallel beams creating alternating high-dose peaks and low-dose valleys. Key question: do peaks show direct damage signatures while valleys show immune priming / bystander effects, and how does this compare to uniform SBRT?

## Datasets

- **Mutter_01** (current): 11 conditions (MBRT/SBRT/Control × timepoints), 971K cells post-QC, CosMx ~1000-gene panel. Input parquets at `/mnt/data/projects/spatial-rads/inputs/mutter01/`
- **Mutter_02**: 4 slides, 12 samples (3 per slide: Control/MBRT/SBRT at 2d), 2.35M cells post-QC. Raw Seurat RDS at `/mnt/data/projects/spatial-rads/inputs/mutter02/`. Fully processed via `processing.smk` (2026-05-30): sample assignment, QC, LogNormalize, cell typing via Seurat anchor TransferData against a pooled M01 cell reference (replaces InSituType, which failed on cross-dataset batch shift in this per-sample application — see `plan-processing-pipeline.md`; note the `aggregate.smk` layer re-types all cells with InSituType at merged scale, a different application), and UCell/AddModuleScore pathway scoring.

## Analysis Architecture

**Reproducible pipeline (per-sample / pre-aggregate)** — `workflows/processing.smk`, completed 2026-05-30. Adapter (raw → common 950-gene panel Seurat) → 4-criteria QC filter → LogNormalize → probe-QC diagnostic → cell typing (Yi labels + tumor_epithelial relabel for M01; Seurat anchor TransferData for M02) → UCell + AddModuleScore pathway scoring. 23 per-sample scored.rds at `/mnt/data/projects/spatial-rads/processing/scored/`. See `plan-processing-pipeline.md`.

**Typing rewrite staged (v2.0, uncommitted, NOT yet run — 2026-05-31).** `scripts/celltype.R` + `config/lineage_markers.yaml` replace the two-branch v1 typing above with *de-anchored* scoring: both datasets typed **identically** by UCell against shared coarse lineage markers (~15 lineages; Epithelial→`tumor_epithelial`), per-cell argmax→`cell_type` with a min-score/margin→`unassigned` cutoff, **no M01→M02 transfer** (cross-dataset labels become structural, not inherited). The 23 scored.rds on disk (2026-05-30) still carry **v1** labels — M01 = Yi ImmGen (`yi_celltype`→`cell_type`), M02 = Seurat TransferData (`cell_type` + `celltype_prob`), M02 sharing M01's vocabulary because it inherited it. Rerun `processing.smk` to adopt v2.0. See the v2.0 banner in `plan-processing-pipeline.md`.

**Cross-sample workflow (`aggregate.smk`)** — in active build (2026-06-01); current design/status in `plan-aggregate.md`. **STATUS 2026-06-01 — the per-cell InSituType typing described below was RUN and FAILED validation** (84.6% of cells mislabeled `tumor_epithelial`, immune gold-marker recall <5%); root cause is the per-cell *unit of inference* at this panel sparsity, NOT the `EPI_GATE`/tumor anchor — do NOT re-tune those. Coarse typing has pivoted to **cluster-then-annotate** (`plan-aggregate.md` v2.3; memory `project_insitutype_aggregate_failure.md`). **No cell-type labels — per-sample TransferData OR aggregate InSituType — are currently trustworthy; do not run composition / DEG / pathway / neighborhood on them until v2.3 lands.** The InSituType narrative that follows is retained as the record of the failed detour: Merge → **InSituType cell typing run once on merged M01+M02 raw counts** against an external NanoString mammary + ImmGen atlas (semi-supervised, `update_reference_profiles=TRUE`, per-cell negprobe modeling; de-novo clusters absorb the 4T1 tumor, Epcam/Krt8 → `tumor_epithelial`) — this **re-types all cells and supersedes the per-sample TransferData labels**, making cross-dataset typing structural with no M01→M02 transfer → Harmony embedding + **Leiden clustering via `igraph::cluster_leiden`** (standard tools; `lisi` package for the batch-mixing gate) as UMAP / clustering-QC substrate → Banksy niche → hierarchical immune second pass (SingleR/celldex ImmGen on the immune subset) → composition + pseudobulk cross-dataset DEGs. The earlier Option C integrate-then-transfer typing path is superseded. *(As-built 2026-06-01: `embed_celltype.R` groups Harmony on `slide_id` alone with interim LISI gate 1.5 — Harmony 2.0.2 LAPACK-errors on ≥2 covariates, `project_harmony_lapack`; `slide_id` nests `dataset` so the dataset axis is still corrected; the planned `c("dataset","slide_id")` + gate 1.8 awaits a LAPACK-enabled Harmony. The expensive embedding is checkpointed to `aggregate/embed_checkpoint.rds` after FindNeighbors so re-clustering costs minutes; the gate runs before the viz-only UMAP. Clusters are QC-only — InSituType types from raw counts.)*

**Scientific analysis layers** — exploratory work in `dev/`:

1. **QC & Landscape** (complete) — `dev/peak_valley_analysis/00-03_*.R`: load, filter, normalize, cluster, cell type validation.
2. **MBRT vs SBRT Kinetics** (complete) — `dev/02_deg_analysis.R` etc.: DEGs, pathway kinetics, composition, spatial neighborhoods.
3. **Peak/Valley at 4h** (complete) — `dev/peak_valley_analysis/05-12_*.R`: stripe detection, H2AX validation, classification, peak-vs-valley DEGs/pathways/composition/neighborhoods, SBRT comparison.
4. **Signature Projection** (planned, runs on jeff-frag-test VM) — project peak/valley signatures across all timepoints; see `plan-mbrt-signatures.md`.

## Code Organization

- `spatial-rads.org` — main literate notebook (~3900 lines), all analysis + documentation. Top `* Plans` heading indexes active plan docs.
- `workflows/*.smk` — Snakemake workflows. `data_model.smk`, `processing.smk` (tangled from org); `aggregate.smk` in active build (2026-06-01).
- `scripts/*.R` — R steps invoked by snakemake rules. The `processing.smk` steps are tangled from org; the `scripts/aggregate/*.R` steps are authored directly as standalone scripts (not in the org notebook).
- `plan-processing-pipeline.md` — full design/decisions/status of the per-sample reproducible pipeline.
- `plan-aggregate.md` — handoff doc for the next cross-sample / cross-dataset workflow build.
- `plan-mbrt-signatures.md` — Layer-4 signature projection plan (autonomous VM runs).
- `results/processing/` — small TSVs and reports from `processing.smk` (qc_summary, celltype_summary, probe_qc_report, plots).
- `dev/peak_valley_analysis/` — scientific Layer-3 work, exploratory.
- `config/spatial-rads-conda-env.yaml` — conda environment definition.
- `config/config.yaml` — workflow config (datadir, samplesheet, qc thresholds, normalize scale_factor).
- `config/pathway_gene_lists.yaml` — canonical pathway gene sets (IFN-I/II, DDR, STING).
- `data/metadata.xlsx` — experimental design metadata (single source of truth). Relational sheets: `mice`/`datasets`/`slides`/`samples` plus `if_channels` (one row per dataset×morphology-IF stain). The `datasets` sheet carries dataset-level RNA-panel design (`panel_base`, `panel_custom`, `assay_type`=RNA, `panel_n`, `n_negprobe`, `n_falsecode`); `if_channels` captures the 5 CosMx morphology stains (DAPI/CD298.B2M/PanCK/CD45/G) with segmentation/lineage `role`. IF channels are morphology stains for segmentation+lineage, NOT a protein-expression assay. Built from a legacy flat source by `scripts/build_metadata_xlsx.R`.
- `data/metadata_schema.yaml` — Frictionless-style schema (PK/enum/FK/required) consumed by `make_data_model.R` for validation. One resource per xlsx sheet.
- `data/data_model.rda` — relational model (list of tibbles: mice/datasets/slides/samples/if_channels) emitted by `data_model.smk`; `results/data_model/samples.tsv` is the sample-grain denormalized join (dataset-grain dims like `if_channels` are intentionally not joined in).
- `data/sources/` — write-protected vendor documents (NanoString panel files, quotes).

## Tech Stack

- **R** (primary): Seurat v5, tidyverse, data.table, arrow, ComplexHeatmap, RANN, patchwork
- **Conda env**: `spatial-rads` (R ≥4.4, conda-forge + bioconda)
- **Compute**: jeff-beast (48 cores, 124 GB RAM). Memory, not CPU, is the practical bound on Seurat operations over the merged ~3.32M-cell cohort.
- **Data storage**: `/mnt/data/projects/spatial-rads/` (local disk; mirrored from GCS). Heavy intermediates (raw/qc/norm/typed/scored RDS, reference) live here; small TSVs + plots in `results/`.
- **Snakemake**: workflows under `workflows/*.smk`, driver in `basecamp` env, R steps via `conda run -n spatial-rads Rscript`. **Always set `TMPDIR=/mnt/data/projects/spatial-rads/tmp`** (not root `/tmp`). The processing.smk workflow bakes `shell.prefix("export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1; ")` + per-rule `threads: 1` to prevent R BLAS thread oversubscription (see memory `feedback_smk_thread_hygiene`); aggregate.smk should do the same.

## Key Conventions

- **Snakemake invocation**: `TMPDIR=/mnt/data/projects/spatial-rads/tmp conda run -n basecamp snakemake -s workflows/<name>.smk --cores N`. Dry-run before execute (`--dry-run`).
- **Per-sample data layout** (post-`processing.smk`): `processing/raw/{s}.raw.rds` → `processing/qc/{s}.qc.rds` → `processing/norm/{s}.norm.rds` → `processing/typed/{s}.typed.rds` → `processing/scored/{s}.scored.rds`, all under `/mnt/data/projects/spatial-rads/`.
- **Cell typing**: M01 keeps Yi ImmGen labels with confident-epithelial `a`/`b`/NA cells relabeled to `tumor_epithelial`; M02 uses Seurat anchor TransferData against a 20k natural-abundance pooled-M01 cell reference. **For "tumor compartment" cuts, group `a` + `tumor_epithelial`** — M02 tumor cells land mostly in `a` faithfully to M01. *(These are the per-sample processing labels; `aggregate.smk` re-types all cells with InSituType against the external atlas — see Analysis Architecture — emitting `tumor_epithelial` directly, so the `a`-bucket grouping applies only to per-sample / legacy outputs.)*
- **Exploratory R scripts** (legacy): `conda run -n spatial-rads Rscript dev/peak_valley_analysis/<script>.R`
- **Stats**: n=1 per condition in Mutter_01 → report effect sizes (log2FC, % expressed), NOT p-values. Mutter_02 brings replicates for formal inference at 2d.
- Peak/valley ground truth from H2AX IHC at 4h only.
- Stripe model: 15-deg tilt, 1.02 mm spacing, 4 peaks.
- Tongue and flank models analyzed separately.

## Key Contacts

- Yi Liu — CosMx data provider, cell type annotations, pathway scores
- Rob Mutter — Lab PI
- Jenn Fazzari — slide layout coordination

## Autonomous VM Runs

Layer 4 (signature projection) runs on **jeff-frag-test** GCP VM via the virtual-scientist skill. See memory `reference_vm_frag_test.md` for gcloud coords (project `aif-usr-p-chaudhuri-lab-83f0`, zone `us-west4-b`, user `ext_szymanski_jeffrey_mayo_edu`).

- **Orchestrator**: `run-mbrt-signatures.sh` — three sequential `claude -p` phases (executor → reviewer → reviser) with validation gates
- **Scratch data on VM**: `/mnt/data/spatial-rads/` (staged with `cp`, NOT symlink from `/mnt/gcs/`)
- **Mutter_02 slides on VM**: `/mnt/data/spatial-rads/mutter02/seuratObject_0[1-4]_Mutter_02_CosMmR.RDS`
- **Results**: written to `results/signatures/` on the VM's feature branch, retrieved via `git format-patch` + `scp` (Mayo VMs cannot push)
- **Max plan overage**: `claude -p` must include `--max-budget-usd 50` or the headless API stream closes with `UND_ERR_SOCKET` when the overage wall is hit
- **Observability**: `executor.log` is unreliable (Node buffers stdout). Ground truth is the session jsonl at `~/.claude/projects/<escaped-cwd>/*.jsonl`; the launcher captures the active jsonl pointer to `results/signatures/<phase>.jsonl.path`. Tail filter is documented in the virtual-scientist skill.

## Statistical Caveats

- All Mutter_01 analyses are descriptive (n=1). Formal testing awaits Mutter_02 replicates.
- 1000-gene panel has blind spots; negative results need cautious interpretation.
- "Signature persistence" at later timepoints may reflect cell turnover, not intrinsic memory.
