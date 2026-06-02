# spatial-rads

Spatial transcriptomics analysis of microbeam radiation therapy (MBRT) using CosMx single-cell data from the Mutter Lab (Mayo Clinic).

## Scientific Context

First-in-field spatial transcriptomics study of MBRT peak/valley biology. MBRT delivers dose through ~1mm-spaced parallel beams creating alternating high-dose peaks and low-dose valleys. Key question: do peaks show direct damage signatures while valleys show immune priming / bystander effects, and how does this compare to uniform SBRT?

## Datasets

- **Mutter_01** (current): 11 conditions (MBRT/SBRT/Control × timepoints), 971K cells post-QC, CosMx ~1000-gene panel. Input parquets at `/mnt/data/projects/spatial-rads/inputs/mutter01/`
- **Mutter_02**: 4 slides, 12 samples (3 per slide: Control/MBRT/SBRT at 2d), 2.35M cells post-QC. Raw Seurat RDS at `/mnt/data/projects/spatial-rads/inputs/mutter02/`. Fully processed via `processing.smk` (2026-05-30): sample assignment, QC, LogNormalize, cell typing via Seurat anchor TransferData against a pooled M01 cell reference (replaces InSituType, which failed on cross-dataset batch shift in this per-sample application — see `plan-processing-pipeline.md`; note the `aggregate.smk` layer re-types all cells at merged scale via scVI integration + cluster-then-annotate, a different application that supersedes these per-sample labels for cross-dataset work — see Analysis Architecture), and UCell/AddModuleScore pathway scoring.

## Analysis Architecture

**Reproducible pipeline (per-sample / pre-aggregate)** — `workflows/processing.smk`, completed 2026-05-30. Adapter (raw → common 950-gene panel Seurat) → 4-criteria QC filter → LogNormalize → probe-QC diagnostic → cell typing (Yi labels + tumor_epithelial relabel for M01; Seurat anchor TransferData for M02) → UCell + AddModuleScore pathway scoring. 23 per-sample scored.rds at `/mnt/data/projects/spatial-rads/processing/scored/`. See `plan-processing-pipeline.md`.

**Typing rewrite staged (v2.0, uncommitted, NOT yet run — 2026-05-31).** `scripts/celltype.R` + `config/lineage_markers.yaml` replace the two-branch v1 typing above with *de-anchored* scoring: both datasets typed **identically** by UCell against shared coarse lineage markers (~15 lineages; Epithelial→`tumor_epithelial`), per-cell argmax→`cell_type` with a min-score/margin→`unassigned` cutoff, **no M01→M02 transfer** (cross-dataset labels become structural, not inherited). The 23 scored.rds on disk (2026-05-30) still carry **v1** labels — M01 = Yi ImmGen (`yi_celltype`→`cell_type`), M02 = Seurat TransferData (`cell_type` + `celltype_prob`), M02 sharing M01's vocabulary because it inherited it. Rerun `processing.smk` to adopt v2.0. See the v2.0 banner in `plan-processing-pipeline.md`.

**Cross-sample workflow (`aggregate.smk`)** — cell typing **executed and locked 2026-06-02**; current design/status in `plan-aggregate.md` (v2.5 banner). Cross-dataset typing assigns identity **once on the merged 3,277,090-cell, 950-gene object** (no M01→M02 transfer; typing is structural). Coarse typing is **cluster-then-annotate, not per-cell**: an initial per-cell InSituType run failed validation (84.6% mislabeled `tumor_epithelial`, immune gold-marker recall <5%) because the per-cell *unit of inference* is unreachable at this panel sparsity (median 58/950 genes per cell) — that was the root cause, NOT the tumor anchor (memory `project_insitutype_aggregate_failure.md`). **Tier-1 (coarse):** scVI integration [Lopez 2018] (30 latent dims, 2 layers, NB likelihood, `batch_key=slide_id`, per-cell negprobe as continuous covariate, 50 epochs, GPU-direct on the RTX A4000 — scVI replaced Harmony, which 2.0.2-LAPACK-errors on ≥2 covariates) → 15-NN graph → **Leiden** [Traag 2019] (`flavor="leidenalg"`, `n_iterations=2`; igraph over-partitions, run-to-convergence intractable at 3.27M nodes) with a swept resolution (0.5 failed premise; 1.0/2.0/3.0 swept; **res=3.0 adopted**, finest stable partition) → cluster-level marker annotation (≥100 cells + signal-to-background ≥2). All four pre-registered Q4 gates PASSED (premise; batch-mixing by imbalance-scaled iLISI 0.67 vs 0.5 floor, ceiling 1.70 for the 29/71 split, per-cluster minority demoted to a QC flag; marker-recall; negprobe QC). Result: **tumor 48.1% (1,577,685) / stroma 40.6% (1,331,601) / immune 11.2% (367,800), 4 unassigned.** **Tier-2 (immune):** cluster-level **SingleR** [Aran 2019] vs celldex **ImmGen** (label.main) on the 367,800 immune cells (subtyped at res=3.0; rare lineages re-examined at res=6.0). SingleR assigns no B/Plasma/Neutrophil class at this sparsity, so a **marker-rescue step** (`tier2_immune_rescue.py` — the Mast-rule "marker evidence overrides reference confidence", precedent Cheng et al. BMC Bioinformatics 2025 Xenium benchmark) overrides SingleR where a res=6.0 subcluster's **identity markers** (lineage-defining, non-secreted, non-shared subset of each panel) win by argmax AND ≥2 corroborate. Recovers **Plasma 3,573** (the SingleR-"DC" subcluster 21: Mzb1 17.8×/Xbp1 19.7× bg) and **Neutrophils 259** (subcluster 272: Elane/Prtn3/Mpo, S100a8 260×); only 3,832 cells (1.0% of immune) reassigned. A plain signal-to-background argmax was rejected for *assignment* (it stays the detectability *gate*): ambient immunoglobulin (Jchain/Igkc) and non-specific genes (Cd37 pan-leukocyte, Xbp1 secretory) corrupt set-means — naive argmax falsely flipped an 8k T-cell subcluster (Cd3e 44×) to Plasma. **Final immune: macrophages 300,975 / ILC 34,340 / T 16,046 / DC 5,315 / NK 4,991 / Plasma 3,573 / Mast 1,121 / Neutrophils 259 (+1,180 epithelial-contamination).** **B cells confirmed near-absent**: 0/298 res=6.0 subclusters have B as the argmax identity (B's signal always rank ≥2); B markers on-panel but diffuse (2.4–3.3% of cells carry any B count). (An in-session "B RETAIN" reading was an artifact of per-lineage-best detectability output — B's panel clears background at subcluster 3, but that subcluster is dominantly Plasma/T — corrected here.) Final labels at `aggregate/full/immune_subtypes_rescued.parquet` (pre-rescue SingleR-only in `immune_subtypes.parquet`). These labels **supersede the per-sample TransferData labels** for all cross-dataset analysis. **Tier-2 (stroma):** UCell cluster-level argmax+margin (`tier2_stroma_ucell.R`) against the `config/lineage_markers.yaml` stromal sets (NOT SingleR/MouseRNAseqData, whose bulk reference lacks pericyte/smooth-muscle classes) → Fibroblast 796,933 / SmoothMuscle 73,570 / unassigned 461,098 (pre-rescue UCell). A **marker-rescue step** (`tier2_stroma_rescue.py`, the same identity-marker override built for immune, examining r6 subclusters) then recovered the minor lineages UCell's fibroblast-washed rank score missed but the detectability gate confirmed present: **Endothelial 36,937** (its best subcluster 8 had been 99.5% mislabeled Fibroblast by UCell; anchored by Cdh5/Pecam1) and **Adipocyte 49,214** (anchored by lipid-droplet Cidea; present in BOTH cohorts so not a single-dataset artifact). A **specificity anchor** (≥1 lineage-exclusive marker must individually clear the gate) blocks shared-gene false flips — Fabp4 is adipocyte+endothelial+macrophage (s2b 254 on the adipocyte cluster but also 39 on the endothelial one), so a Fabp4-driven call with no Cidea is rejected; this caught 2 false subclusters (Cidea 0.83/1.96). 86,151 cells (6.5% of stroma) reassigned, only out of Fibroblast/SMC/unassigned. **Pericyte stays 0** (best subcluster 8 shared with and dominated by Endothelial; only 3 panel markers, Pdgfrb fibroblast-shared — present but unresolvable de novo, folded into the perivascular endothelial call). **Final stroma: Fibroblast 732,099 / unassigned 441,227 / SmoothMuscle 72,124 / Adipocyte 49,214 / Endothelial 36,937.** Final labels at `aggregate/full/stroma_subtypes_rescued.parquet` (pre-rescue UCell-only in `stroma_subtypes.parquet`). **Unified per-cell label table** joining all three tiers (`unify_labels.py` → `results/aggregate/full_labels.parquet`, 3,277,090 rows; `compartment` + final `cell_subtype`) is the single canonical input for downstream analysis. Still parked: BANKSY niche, composition + pseudobulk cross-dataset DEGs, tumor-state pathway scoring. Conda: `spatial-rads-scvi` (Python: scvi-tools/scanpy/leidenalg/GPU torch) for tier-1, `spatial-rads` (R: SingleR/celldex, needs `scrapper`) for tier-2. Tier scripts: `scripts/aggregate/{full_export.R,full_cluster.py,finalize_tier1.py,tier2_immune_subcluster.py,tier2_singler.R,tier2_immune_rescue.py,tier2_stroma_subcluster.py,tier2_stroma_ucell.R,tier2_stroma_rescue.py,tier2_detectability.py,bcell_diagnostic.py,unify_labels.py}`. Published-QC helpers vendored (pinned to SpatialQM commit 36e1a59d4ca3): `spatialqm_metrics.R` (`sqm_mecr` sample-level segmentation-contamination QC; `sqm_pseudobulk_cor` immune-vs-ImmGen cross-check) — `tier2_detectability.py` is the per-lineage signal-to-background gate (= SpatialQM `getMeanSignalRatio`, cited not re-implemented).

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
- **Compute**: jeff-beast (48 cores, 124 GB RAM, **NVIDIA RTX A4000 16 GB GPU**). Memory, not CPU, is the practical bound on Seurat operations over the merged ~3.32M-cell cohort. The GPU drives scVI integration (env `spatial-rads-scvi`), which trains GPU-direct — no CPU sketch-train workaround needed.
- **Data storage**: `/mnt/data/projects/spatial-rads/` (local disk; mirrored from GCS). Heavy intermediates (raw/qc/norm/typed/scored RDS, reference) live here; small TSVs + plots in `results/`.
- **Snakemake**: workflows under `workflows/*.smk`, driver in `basecamp` env, R steps via `conda run -n spatial-rads Rscript`. **Always set `TMPDIR=/mnt/data/projects/spatial-rads/tmp`** (not root `/tmp`). The processing.smk workflow bakes `shell.prefix("export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1; ")` + per-rule `threads: 1` to prevent R BLAS thread oversubscription (see memory `feedback_smk_thread_hygiene`); aggregate.smk should do the same.

## Key Conventions

- **Snakemake invocation**: `TMPDIR=/mnt/data/projects/spatial-rads/tmp conda run -n basecamp snakemake -s workflows/<name>.smk --cores N`. Dry-run before execute (`--dry-run`).
- **Per-sample data layout** (post-`processing.smk`): `processing/raw/{s}.raw.rds` → `processing/qc/{s}.qc.rds` → `processing/norm/{s}.norm.rds` → `processing/typed/{s}.typed.rds` → `processing/scored/{s}.scored.rds`, all under `/mnt/data/projects/spatial-rads/`.
- **Cell typing**: M01 keeps Yi ImmGen labels with confident-epithelial `a`/`b`/NA cells relabeled to `tumor_epithelial`; M02 uses Seurat anchor TransferData against a 20k natural-abundance pooled-M01 cell reference. **For "tumor compartment" cuts, group `a` + `tumor_epithelial`** — M02 tumor cells land mostly in `a` faithfully to M01. *(These are the per-sample processing labels; `aggregate.smk` re-types all cells at merged scale via scVI cluster-then-annotate — see Analysis Architecture — emitting coarse `tumor`/`stroma`/`immune` (plus tier-2 immune subtypes) directly, so the `a`-bucket grouping applies only to per-sample / legacy outputs. For cross-dataset work, use the **unified per-cell label table** at `results/aggregate/full_labels.parquet` (one row per cell: `compartment` = tumor/stroma/immune + final `cell_subtype` joining all three tiers, plus `subtype_source` and `rescued`; built by `scripts/aggregate/unify_labels.py`). Its per-tier components: aggregate coarse labels at `results/aggregate/full_coarse_labels.parquet`, immune subtypes at `/mnt/data/projects/spatial-rads/aggregate/full/immune_subtypes_rescued.parquet` (marker-rescued final; `immune_subtypes.parquet` is the pre-rescue SingleR-only intermediate), and stroma subtypes at `aggregate/full/stroma_subtypes_rescued.parquet` (marker-rescued final; `stroma_subtypes.parquet` is the pre-rescue UCell-only intermediate).)*
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
