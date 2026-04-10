# spatial-rads

Spatial transcriptomics analysis of microbeam radiation therapy (MBRT) using CosMx single-cell data from the Mutter Lab (Mayo Clinic).

## Scientific Context

First-in-field spatial transcriptomics study of MBRT peak/valley biology. MBRT delivers dose through ~1mm-spaced parallel beams creating alternating high-dose peaks and low-dose valleys. Key question: do peaks show direct damage signatures while valleys show immune priming / bystander effects, and how does this compare to uniform SBRT?

## Datasets

- **Mutter_01** (current): 11 conditions (MBRT/SBRT/Control × timepoints), 971K cells post-QC, CosMx ~1000-gene panel. Input parquets on GCS: `/mnt/gcs/jeszyman/projects/spatial-rads/inputs/`
- **Mutter_02** (incoming): 4 slides, biological replicates + new tumor model. Raw Seurat objects available from Yi Liu. Not yet processed.

## Analysis Architecture

4 progressive layers, each building on the last:

1. **QC & Landscape** (complete) — scripts 00-03: load, filter, normalize, cluster, cell type validation
2. **MBRT vs SBRT Kinetics** (complete) — script 04: DEGs, pathway kinetics, composition, spatial neighborhoods
3. **Peak/Valley at 4h** (complete) — scripts 05-12: stripe detection, H2AX validation, classification, peak vs valley DEGs/pathways/composition/neighborhoods, SBRT comparison
4. **Signature Projection** (planned) — project peak/valley signatures across all timepoints

## Code Organization

- `spatial-rads.org` — main literate programming file (~2700 lines), all analysis + documentation
- `dev/peak_valley_analysis/` — tangled R scripts (00-12), `run_all.R` master runner
- `dev/peak_valley_analysis/data/` — output TSV files (14 tables)
- `dev/peak_valley_analysis/plots/` — output figures (25 PNGs)
- `config/spatial-rads-conda-env.yaml` — conda environment definition
- `docs/superpowers/` — design specs and implementation plans
- `metadata.xlsx` — experimental design metadata

## Tech Stack

- **R** (primary): Seurat v5, tidyverse, data.table, arrow, ComplexHeatmap, RANN, patchwork
- **Conda env**: `spatial-rads` (R ≥4.4, conda-forge + bioconda)
- **Compute**: jeff-beast (48 cores), future multicore parallelization
- **Data storage**: `/mnt/gcs/jeszyman/projects/spatial-rads/` (GCS mount, read-only unless explicit permission)
- **No Snakemake yet** — deferred until analyses stabilize + replicates arrive

## Key Conventions

- Run R scripts via: `conda run -n spatial-rads Rscript dev/peak_valley_analysis/<script>.R`
- Cached RDS objects avoid recomputation between scripts
- n=1 per condition: report effect sizes (log2FC, % expressed), NOT p-values
- Peak/valley ground truth from gamma-H2AX IHC at 4h only
- Stripe model: 15-deg tilt, 1.02mm spacing, 4 peaks
- Tongue and flank models analyzed separately

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
