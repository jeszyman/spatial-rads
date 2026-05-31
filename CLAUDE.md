# spatial-rads

Spatial transcriptomics analysis of microbeam radiation therapy (MBRT) using CosMx single-cell data from the Mutter Lab (Mayo Clinic).

## Scientific Context

First-in-field spatial transcriptomics study of MBRT peak/valley biology. MBRT delivers dose through ~1mm-spaced parallel beams creating alternating high-dose peaks and low-dose valleys. Key question: do peaks show direct damage signatures while valleys show immune priming / bystander effects, and how does this compare to uniform SBRT?

## Datasets

- **Mutter_01** (current): 11 conditions (MBRT/SBRT/Control × timepoints), 971K cells post-QC, CosMx ~1000-gene panel. Input parquets at `/mnt/data/projects/spatial-rads/inputs/mutter01/`
- **Mutter_02**: 4 slides, 12 samples (3 per slide: Control/MBRT/SBRT at 2d), 2.35M cells post-QC. Raw Seurat RDS at `/mnt/data/projects/spatial-rads/inputs/mutter02/`. Fully processed via `processing.smk` (2026-05-30): sample assignment, QC, LogNormalize, cell typing via Seurat anchor TransferData against a pooled M01 cell reference (replaces InSituType, which failed on cross-dataset batch shift — see `plan-processing-pipeline.md`), and UCell/AddModuleScore pathway scoring.

## Analysis Architecture

**Reproducible pipeline (per-sample / pre-aggregate)** — `workflows/processing.smk`, completed 2026-05-30. Adapter (raw → common 950-gene panel Seurat) → 4-criteria QC filter → LogNormalize → probe-QC diagnostic → cell typing (Yi labels + tumor_epithelial relabel for M01; Seurat anchor TransferData for M02) → UCell + AddModuleScore pathway scoring. 23 per-sample scored.rds at `/mnt/data/projects/spatial-rads/processing/scored/`. See `plan-processing-pipeline.md`.

**Cross-sample workflow (`aggregate.smk`)** — not yet built. Bootstrap design in `plan-aggregate.md`; covers merge → Harmony integration → cluster → UMAP → cell-type landscape → cross-dataset DEGs. Option C (integrate-then-transfer for typing) deferred here.

**Scientific analysis layers** — exploratory work in `dev/`:

1. **QC & Landscape** (complete) — `dev/peak_valley_analysis/00-03_*.R`: load, filter, normalize, cluster, cell type validation.
2. **MBRT vs SBRT Kinetics** (complete) — `dev/02_deg_analysis.R` etc.: DEGs, pathway kinetics, composition, spatial neighborhoods.
3. **Peak/Valley at 4h** (complete) — `dev/peak_valley_analysis/05-12_*.R`: stripe detection, H2AX validation, classification, peak-vs-valley DEGs/pathways/composition/neighborhoods, SBRT comparison.
4. **Signature Projection** (planned, runs on jeff-frag-test VM) — project peak/valley signatures across all timepoints; see `plan-mbrt-signatures.md`.

## Code Organization

- `spatial-rads.org` — main literate notebook (~3900 lines), all analysis + documentation. Top `* Plans` heading indexes active plan docs.
- `workflows/*.smk` — Snakemake workflows (tangled from org). Currently: `data_model.smk`, `processing.smk`. Future: `aggregate.smk`.
- `scripts/*.R` — R steps invoked by snakemake rules (tangled from org).
- `plan-processing-pipeline.md` — full design/decisions/status of the per-sample reproducible pipeline.
- `plan-aggregate.md` — handoff doc for the next cross-sample / cross-dataset workflow build.
- `plan-mbrt-signatures.md` — Layer-4 signature projection plan (autonomous VM runs).
- `results/processing/` — small TSVs and reports from `processing.smk` (qc_summary, celltype_summary, probe_qc_report, plots).
- `dev/peak_valley_analysis/` — scientific Layer-3 work, exploratory.
- `config/spatial-rads-conda-env.yaml` — conda environment definition.
- `config/config.yaml` — workflow config (datadir, samplesheet, qc thresholds, normalize scale_factor).
- `config/pathway_gene_lists.yaml` — canonical pathway gene sets (IFN-I/II, DDR, STING).
- `data/metadata.xlsx` — experimental design metadata (single source of truth).
- `data/data_model.rda` — relational sample/slide/dataset model emitted by `data_model.smk`.
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
- **Cell typing**: M01 keeps Yi ImmGen labels with confident-epithelial `a`/`b`/NA cells relabeled to `tumor_epithelial`; M02 uses Seurat anchor TransferData against a 20k natural-abundance pooled-M01 cell reference. **For "tumor compartment" cuts, group `a` + `tumor_epithelial`** — M02 tumor cells land mostly in `a` faithfully to M01.
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
