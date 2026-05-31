# `aggregate.smk` — session handoff (plan for a plan)

Bootstrap material for the next session, which will fully design and build the
**aggregate** workflow (per-sample → cross-sample / cross-dataset). Not the final
spec — read this first, then brainstorm the actual design with `devils-advocate`,
then write the implementation plan.

Companion to `plan-processing-pipeline.md` (the per-sample upstream that aggregate
consumes).

---

## What to read before starting

1. **`CLAUDE.md`** — project conventions, tech stack, contacts, the n=1 statistical caveats.
2. **`plan-processing-pipeline.md`** — what's upstream, key decisions (especially the
   Seurat anchor TransferData substitution and the Option-C note for cross-dataset typing).
3. **Memory files** at `~/.claude/projects/-home-jeszyman-repos-spatial-rads/memory/`:
   - `feedback_smk_thread_hygiene.md` — **mandatory**: bake `shell.prefix` BLAS pinning +
     `threads: 1` into `aggregate.smk` from line 1, not as a retrofit.
   - `feedback_dev_workflow.md` — two-phase: standalone R for iteration → fold into
     literate org once stable.
   - `feedback_local_data.md` — `gsutil cp` GCS inputs to local disk before R reads.
   - `feedback_no_local_tmp.md` — `TMPDIR=/mnt/data/projects/spatial-rads/tmp` always.
   - `feedback_no_polling.md` — background jobs + notification, never poll.
   - `project_tumor_cell_annotation.md`, `project_tumor_model.md` — biology context.
4. **`spatial-rads.org`** — `** Initial processing and QC` for the working literate
   pattern (smk org block + R-script headings tangling byte-identical), and the
   `* Project Map` heading for current scientific framing.
5. **`dev/` exploratory code** — prior analytical patterns to lift, NOT re-invent:
   - `dev/02_deg_analysis.R` — DEG patterns for cross-condition contrasts.
   - `dev/explore_composition.R` — composition matrices, plots.
   - `dev/explore_landscape.R` — cell-type landscape.
   - `dev/explore_spatial_nn.R` — spatial neighborhoods.
   - `dev/peak_valley_analysis/` — Layer 3 work; aggregate should leave room for
     peak/valley integration (probably its own future smk; do not absorb).

Also worth checking: `git log --oneline -10`, `git status`, `ls /mnt/data/projects/spatial-rads/processing/scored/ | wc -l` (expect 23).

---

## Project context (terse)

CosMx spatial transcriptomics, two datasets (Mutter_01, Mutter_02), MBRT-vs-SBRT
radiotherapy of 4T1 flank tumors in Balb/c. Goal: cross-dataset reproducible
MBRT-vs-SBRT findings, where "reproducible" = same pipeline through both datasets.
Mutter_01 is n=1 per condition (descriptive only); Mutter_02 brings replicates for
formal inference. ~3.32M cells post-QC across the 23 samples.

---

## What `processing.smk` hands to `aggregate.smk`

23 per-sample scored Seurat RDS at
`/mnt/data/projects/spatial-rads/processing/scored/<sample_id>.scored.rds`:

- **Assay**: RNA counts + data (LogNormalize, sf=1e4) on the common **950-gene** panel
  (Mutter_01 1000 ∩ Mutter_02 972).
- **Metadata columns**: `sample_id`, `dataset` (Mutter_01/Mutter_02), `treatment`
  (NT/MBRT/SBRT), `condition`, `timepoint_h`, `model` (flank/tongue), `cell_type`
  (Yi ImmGen + `tumor_epithelial` relabel for M01; Seurat anchor TransferData against
  the M01 reference for M02), `celltype_prob` (M02 only — `prediction.score.max`),
  pathway columns `TypeI_interferon_UCell` / `_AMS`, `TypeII_interferon_*`,
  `DNA_Damage_Repair_*`, `STING_*`, and the spatial coordinates `x_slide_mm`,
  `y_slide_mm`, plus morphology / IF channel means.
- **Sample-level master sheet**: `results/data_model/samples.tsv` (23 rows; FK-validated).
- **Reference**: `/mnt/data/projects/spatial-rads/processing/celltype_reference.rds`
  (20k pooled M01 cells, natural-abundance; available for re-validation, not strictly
  needed at aggregate stage).

Composition reality after the Seurat-transfer fix: `a` bucket dominates the M02 tumor
compartment (~55%) faithfully to M01 labeling; group `a`+`tumor_epithelial` as
"tumor compartment" for analysis.

---

## What `aggregate.smk` is *for* and what it is *not*

**Is for**: anything that needs all cells in one place — merge, batch assessment,
integration (Harmony or similar), joint embedding/clustering, cross-sample
composition, cross-dataset DEGs by celltype × condition × timepoint, the
cell-type landscape that frames everything downstream.

**Is not**: per-sample work (lives in `processing.smk`), peak/valley analysis
(stays in `dev/peak_valley_analysis/` for now, possibly its own future smk),
formal statistical inference on Mutter_01 alone (n=1 per condition; descriptive
stats only — defer p-values to Mutter_02 replicate strata).

---

## Sketched stages — these are starting points, not the design

Use brainstorming to firm these up. Probably:

1. **Merge** — pool all 23 scored.rds into a single Seurat (or AnnData) object.
   ~3.32M cells × 950 genes. Sparse-only; never densify the full matrix at once.
   Watch memory; this might need disk-backed (HDF5/anndata) or per-dataset
   intermediate merges.
2. **Pre-integration batch assessment** — PCA on merged; quantify dataset/slide
   batch effect (silhouette / kBET / LISI) to motivate integration choice.
3. **Integration** — Harmony is the default expectation (handles batch covariates
   cleanly, scales well). Batch covariate likely `dataset` (M01 vs M02), possibly
   `slide_id`. Output `harmony` reduction.
4. **Embedding** — UMAP on harmony dims. Maybe per-dataset UMAPs too for sanity.
5. **Clustering** — Leiden at one or a small set of resolutions. Annotate by
   majority `cell_type` for sanity (validates the transfer).
6. **Cell-type landscape** — composition matrices stratified by
   `dataset × treatment × timepoint × model`, plus plots (bar/heatmap).
7. **Cross-dataset DEGs** — per `cell_type`, MBRT-vs-Control, SBRT-vs-Control,
   MBRT-vs-SBRT at each timepoint. Pseudobulk where sample replicates exist
   (Mutter_02); per-cell with descriptive language on Mutter_01.
8. **Spatial neighborhoods** — optional at aggregate stage; the per-sample
   coordinates are preserved, so per-sample neighborhood graphs can be built here
   or kept separate. Cross-sample neighborhood comparisons probably belong here.
9. **Reports + figures** — landscape composition, UMAPs by stratum, DEG volcanoes,
   pathway score distributions by stratum.

**Option C (deferred from processing.smk):** integrate M01+M02 first (Harmony),
then transfer labels in the integrated space. Could *replace* per-sample anchor
transfer or *coexist* as a sanity check (compare cluster majority labels vs the
per-sample transfer labels). Decide during brainstorm.

---

## Open design questions to settle in brainstorming

- Single merged object vs per-dataset-then-merge? Memory budget will force this.
- Harmony batch covariates: `dataset` only? `dataset + slide_id`?
  `dataset + slide_id + timepoint_h`? (Don't regress out treatment — that's signal.)
- Clustering: single resolution or sweep? Annotate with the existing `cell_type`
  via majority vote vs re-cluster cell types?
- DEG method: pseudobulk (DESeq2/edgeR) where replicates exist vs per-cell
  (Wilcoxon/MAST)? n=1 in Mutter_01 forces descriptive-only there.
- Do we run Option C as primary, fallback, or sanity-check-only for typing?
- What gets exported to the report — composition only, or composition + UMAP +
  DEG tables in a self-contained HTML?
- Where does Layer 3 (peak/valley) wire in? Aggregate should leave a clean
  output (typed + scored + clustered + integrated object) that Layer 3 can
  consume sample-by-sample.
- Spatial neighborhoods: per-sample, then aggregate? Or skip at this stage?

---

## Development conventions

- **Literate org**: rules and R scripts under a new `** Aggregate analysis`
  heading in `spatial-rads.org`, mirroring the `** Initial processing and QC`
  structure. Each rule + each script gets its own `*** <name>.R` (or `.smk`)
  with PROPERTIES-level `:tangle` and a single `#+begin_src` block.
- **Two-phase dev**: prototype in `scripts/` first (fast iteration; per-script
  conda-run-Rscript), fold into the org once stable. Watch for the
  script-not-tracked landmine — `shell:`-invoked R scripts are NOT inputs to
  their snakemake rules and don't trigger reruns on edit. Either declare them as
  `input:` deps in `aggregate.smk` from the start, or use the `script:` directive
  (per `~/.claude/rules/snakemake.md`).
- **Snakemake**: workflow at `workflows/aggregate.smk`, driver in `basecamp` env,
  R steps via `conda run -n spatial-rads Rscript`. **Bake from line 1**:
  ```
  shell.prefix("export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1; ")
  ```
  and `threads: 1` on every R rule. Run with
  `TMPDIR=/mnt/data/projects/spatial-rads/tmp conda run -n basecamp snakemake -s workflows/aggregate.smk --cores N`.
  Always dry-run before execute.
- **Big files**: heavy intermediates (merged.rds, harmonized.rds, the multi-million
  cell embedding) go to `/mnt/data/projects/spatial-rads/aggregate/` (or similar),
  not `results/`. Small TSVs and plots go in `results/aggregate/`.
- **Verification**: per the global verification rule, every "done"/"fixed" claim
  needs same-turn evidence (run output, file listing, marker-concordance check).
  Run rendered figures through `Read` before claiming they look right.
- **Process discipline**: brainstorming → devils-advocate critique → spec →
  writing-plans → executing-plans. Use `show-then-execute` for small code chunks.

---

## Hard-won lessons from this session — apply to aggregate

1. **BLAS oversubscription is the silent killer**: `--cores 4` × default
   BLAS threading = swamp. Pin BLAS at the workflow level via `shell.prefix`
   from day one. With BLAS pinned, `--cores 4–8` is safe on jeff-beast.
2. **Memory bounds the cores ceiling, not CPU**: this session, Seurat anchor
   transfer on big M02 queries hit ~68 GB resident per job. For aggregate's
   merge / integrate / cluster on ~3.32M cells, expect even larger spikes. Plan
   for sparse-only and chunked processing; pilot per-stage memory budget before
   choosing `--cores N`.
3. **OOM-killed R leaves a 0-byte log + snakemake error "Killed"**: when a
   concurrent job dies and snakemake SIGTERMs the others, all the per-job
   logs end up empty. Diagnose by reproducing in isolation (`--cores 1`).
4. **Cross-dataset is hard**: InSituType failed on M02 because count-likelihood
   matching can't absorb per-gene platform/efficiency differences between two
   CosMx runs. Use rank/correlation/anchor methods (Seurat anchors at
   per-sample, Harmony at aggregate). Validate post-integration with
   marker-gene concordance (Krt8+/Epcam+ → tumor-compartment cluster), not just
   label distribution.
5. **Validate methods with markers, not labels**: ground-truth-free, robust to
   typing artifacts. Worked for the typing-bug investigation; will work for
   integration QC ("do M01 epithelial + M02 epithelial co-cluster?").
6. **Snakemake doesn't track R scripts referenced via `shell:`** — code edits
   don't trigger reruns. Either declare scripts as `input:` deps from the start
   or use `script:` directive. Costly in iteration if missed.
7. **`pkill -f` matches its own command line** — self-kill bug burned this
   session. Use explicit PID kills or `pgrep | grep -v $$ | xargs kill`.
8. **Background, don't poll**: long jobs → `run_in_background=true`, single
   completion notification, never poll. Per `feedback_no_polling`.
9. **`a`-bucket interpretability** is a known wart of the Yi labeling that
   aggregate inherits. Group `a` + `tumor_epithelial` for any "tumor" cuts;
   note in figures. Refining the taxonomy (split `a` into tumor vs
   stromal-microenv) is a possible aggregate-stage task if needed.

---

## Bootstrap procedure for next session (in order)

1. Sanity-check state: `git log --oneline -5`, `git status`, confirm 23 scored
   RDS exist, confirm `feedback_smk_thread_hygiene` memory still applies.
2. Read this file, then `CLAUDE.md`, then `plan-processing-pipeline.md`.
3. Read the relevant exploratory code in `dev/` (composition, landscape, DEG,
   spatial) to inventory what patterns already exist.
4. **Brainstorm the design** with the brainstorming skill — work through the
   Open Design Questions above. Don't anchor on the sketched stages; restate
   the job-to-be-done with implementation choices stripped out (per the
   `de-anchor before divergent work` rule in `CLAUDE.md`).
5. **Devils-advocate** the proposed design before writing the spec.
6. Write the spec (`docs/superpowers/specs/YYYY-MM-DD-aggregate-design.md` or
   equivalent), then the implementation plan, then execute.
7. Add `aggregate.smk` and its R scripts directly under
   `** Aggregate analysis` in `spatial-rads.org` from the start (literate from
   day one, no standalone phase). Tangle + dry-run after every edit.
8. First end-to-end milestone: a merged, Harmonized, UMAP'd, clustered object
   + a cell-type landscape composition TSV. Everything else (DEGs, spatial,
   reports) layered on after that's stable.

---

## Open items to flag at the start of next session

- `probe_qc_report.tsv` flagged candidates (11 probes at background in both
  datasets, including Cd4) — review and decide smooth / flag-only / drop. Not
  blocking aggregate but should be settled before publication-grade analysis.
- The Option C decision — semi-supervised typing in the integrated space —
  should happen during aggregate brainstorm.
- Whether to refine the `a`-bucket label (split tumor vs stromal-TME) — affects
  every cross-condition tumor stratification downstream.

---

## Key paths cheat-sheet

| What | Where |
|---|---|
| Per-sample scored RDS | `/mnt/data/projects/spatial-rads/processing/scored/<sample_id>.scored.rds` |
| Cell-typing reference | `/mnt/data/projects/spatial-rads/processing/celltype_reference.rds` |
| Master sample sheet | `results/data_model/samples.tsv` |
| Common gene panel | `results/processing/common_genes.tsv` |
| Pathway gene lists | `config/pathway_gene_lists.yaml` |
| Processing workflow | `workflows/processing.smk` + `** Initial processing and QC` in `spatial-rads.org` |
| Literate notebook | `spatial-rads.org` |
| Aggregate scratch | `/mnt/data/projects/spatial-rads/aggregate/` (create) |
| Aggregate results | `results/aggregate/` (create) |
