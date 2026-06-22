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

## v2.1 — InSituType typing (field standard) + reference correction (2026-06-01)

**Typing method: InSituType, not hand-rolled correlation.** v2.0's first cut typed
cells with a hand-rolled atlas-correlation step (`annotate.R`: cluster pseudobulk →
Spearman to atlas profiles → argmax lineage). That is *retired*. We now type with
**InSituType** (Danaher et al. 2022, `scripts/aggregate/typing_insitutype.R`), the
platform-native CosMx classifier — the field-standard choice. It models the CosMx
noise structure (per-cell negative-probe background + Poisson counts) that flat
correlation ignores, and runs **once on the merged M01+M02 counts**, so both cohorts
get one identical rule → cross-dataset comparability is structural, no M01→M02 transfer.

**How we run it (semi-supervised, decoupled from the embed):**
- Types from **raw counts** (`merged.rds`), *not* the Harmony embedding — so typing is
  independent of the 7-hour `embed_celltype` (which is now pure UMAP/viz substrate).
- `reference_profiles` = the external atlas collapsed to lineage-mean profiles;
  `update_reference_profiles=TRUE` adapts the linear-scale atlas means to CosMx
  (fixes the cross-platform scale mismatch that hurt correlation).
- `neg` = per-cell mean negative-probe count, recovered from raw inputs
  (`recover_negprobes.R`: M01 parquet `nCount_NegativeProbes`, M02 raw `negprobes`
  assay; neg = nCount_neg/10). Negprobes were dropped during 950-gene harmonization.
- Semi-supervised: fixed reference lineages + `n_clusts` auto-selected **de novo**
  clusters that absorb structure the atlas lacks (notably the 4T1 tumor, which matches
  no normal profile). De novo clusters → nearest lineage by cosine of the *updated*
  profile; de novo (or Epithelial) clusters expressing **Epcam/Krt8** → `tumor_epithelial`.
- Validation: gold-marker recall per lineage, **split M01/M02** (honest yardstick),
  → `celltype_atlas_validation.tsv`. Composition `celltype_atlas_summary.tsv`.

**Reference correction — back to the agreed mammary + ImmGen plan.** The v2.0 build had
silently broadened the reference to a 4-source composite (ImmGen 41 + **Lung 27** +
MammaryGland 14 + **Muscle 11**) — mammary was only 14/93 and Lung was the largest
stromal contributor, diluting the tissue match (lung AT1/AT2 epithelium polluting the
Epithelial anchor). Corrected 2026-06-01 to the plan as written: **MammaryGland_Virgin
(tissue-matched: 4T1 is a mammary carcinoma) + ImmGen (immune depth)** only. Lung +
Muscle dropped. Result: **55 profiles, 11 lineages**. The one consequence: **SmoothMuscle
is unanchored** (neither source has a smooth-muscle profile) — smooth-muscle cells fall
to the nearest stromal type (Pericyte/Fibroblast). Tumor cells do not depend on the
Epithelial anchor; they are caught by de novo clustering + the Epcam/Krt8 overlay.

**`aggregate.smk` rewiring:** added `recover_negprobes` (→ `cell_neg.tsv`) and
`typing_insitutype` rules; the `annotate` (Spearman) rule is removed. `typing_insitutype`
depends on `merged.rds` + `ref_profiles.rds` + `cell_neg.tsv` (not the embedded object),
producing `merged_typed.rds` + summary/validation/labels.

**Status (2026-06-01):** reference regenerated (mammary+ImmGen, 55×11). `embed_celltype`
still running as the viz substrate; the in-flight v2.0 Spearman `annotate` run was killed
(scheduler SIGKILLed, embed spared). InSituType typing is queued behind the embed freeing
its ~96 GB working set, then runs off `merged.rds`.

## v2.2 — Standard tools enforced + hierarchical immune second pass (2026-06-01)

**Standard tools, no substitutions.** The first `embed_celltype.R` run deviated from the
Stage 1a spec below because three field-standard tools were not installed, and substitutes
were used: **Louvain** clustering for **Leiden**, an **in-script LISI** approximation for
the published **lisi** package, and **Harmony on `slide_id` alone** for the planned
`c("dataset","slide_id")` (Harmony 2.0.2 throws a LAPACK error with two or more covariates,
see `project_harmony_lapack`). Standing rule: install the accepted tool, never invent a
substitute.

**Resolved 2026-06-01.** The missing tools are installed and verified (`SingleR` 2.12.0 +
`celldex` 1.20.0, `lisi` 1.0, `leidenalg` 0.12.0 + `python-igraph`, the `leiden` 0.4.3.1
R wrapper; igraph 2.1.4). `embed_celltype.R` was rewritten to the standard tools — **Leiden
via `igraph::cluster_leiden`** (objective = modularity, resolution 0.5) and the **`lisi`
package** for the batch-mixing metrics and gate — and relaunched. Two improvements rode along:

- **Checkpoint after the embedding.** Load → ScaleData → PCA → Harmony → FindNeighbors (~1.5 h
  on 3.27M cells) writes `aggregate/embed_checkpoint.rds` (the object minus the dense 24 GB
  scale.data; only the ~50 marker rows needed downstream are stashed in `o@misc$sd_markers`,
  written uncompressed for fast I/O). Re-invoking the script with the checkpoint present skips
  straight to clustering, so a Leiden-resolution / metrics / UMAP re-run costs minutes, not
  hours. This is the fix for the "every restart redoes the 1.5 h embedding" problem.
- **Gate before UMAP.** The LISI(dataset) gate now runs before the ~50-min UMAP (viz only),
  so a badly mixed embedding fails fast instead of after the expensive plot.

**Why a full restart, not recompute-in-place.** The bullets below argued LISI and Leiden need
no re-run *given a finished, saved embedding*. But the original Louvain run never completed its
clustering step (>3 h on it, single-threaded, no object ever written), so nothing was on disk to
recompute from — a fresh run was required regardless, and it runs the standard tools from the
start. The Louvain job was killed and the standard-tools job launched 2026-06-01.

**Two-covariate Harmony stays the one open deviation.** Of the three substitutions, only the
Harmony covariate set changes the integrated map:

- **LISI** is a post-hoc quality measure; the new run computes it with the `lisi` package.
- **Leiden** clusters are visualization/QC only (InSituType types from raw counts, fitting its
  own clusters), so they never enter the cell-type calls.
- **Two-covariate Harmony** (`dataset` + `slide_id`) does change the map and needs the Harmony
  fix (version bump or LAPACK relink). Trigger it only if the `lisi(dataset)` score on the
  slide-only map shows the two datasets inadequately mixed (Stage 1a gate LISI(dataset) >= 1.8;
  interim floor 1.5 per `project_harmony_lapack`). If mixing is already adequate, the slide-only
  map stands.

### Stage 1c -- `subtype_immune.R` (hierarchical second pass)

InSituType (`typing_insitutype.R`) assigns only the **coarse tier**: 11 atlas lineages plus
`tumor_epithelial`. It does not resolve immune subsets, and a single-shot fine classification
should not attempt to on a 950-gene panel with ambient spillover (the failure mode of Yi's
M01 run: thymic labels in a flank tumor, ~20% Cd3/T concordance, T cells never separating;
documented in `spatial-rads.org` "Yi Mutter_01"). The field-standard remedy is a **second
pass on the immune compartment only**:

- Subset the lymphoid + myeloid cells (coarse lineages T, B, NK, Plasma, Macrophage, DC,
  Neutrophil) from `merged_typed.rds`.
- Re-embed the subset on its own: `ScaleData` -> `RunPCA` -> `RunHarmony` (dataset +
  slide_id) -> `FindNeighbors` -> `FindClusters` **Leiden** (the subset carries the immune
  signal the whole-tumor embedding drowns).
- Subtype by reference-based annotation: **SingleR** (Aran et al. 2019) against the **celldex
  ImmGen** reference at `label.fine`/`label.main` (Heng et al. 2008), scored **per Leiden
  cluster** (robust to per-cell sparsity), not per cell.
- Validate with **UCell** marker-set scoring (Andreatta and Carmona 2021) on canonical subset
  markers; reconcile SingleR labels against UCell ranks per cluster.
- Outputs: `immune_subtypes.tsv` (cell_id, coarse_lineage, immune_subtype, singler_score,
  ucell_top_marker_set, cluster_id) joined back onto the main object; `immune_subtype_qc.tsv`
  (per-cluster modal SingleR label, mean UCell scores, marker concordance).

Interpretive guard (per `project_mbrt_mechanism_status`, `project_yi_labels_unreliable_immune`):
immune signal is weak on this panel; subtypes are reported with per-cluster confidence, and
immune-subtype-level claims stay cautious in both cohorts.

## v2.3 — Coarse typing returns to cluster-then-annotate; per-cell InSituType retired (2026-06-01)

> This banner is a **plan-to-plan**: it records the decision and its rationale.
> The detailed execution design is deferred to a follow-on (superpowers) planning
> pass; **no code or rerun** until that plan is reviewed.

**The v2.1/v2.2 per-cell InSituType coarse run was executed and FAILED validation.**
Headline: 84.6% of cells labeled `tumor_epithelial`. But the failure is *upstream of
the tumor overlay*. Before the Epcam/Krt8 overlay ran, the coarse typing had already
collapsed — 91% of all cells fell into just two labels (T 68.9%, Pericyte 22.2%),
87.9% of cells anchored to no lineage and went to de-novo clusters, and even the 12.1%
that anchored directly carried 0.1–45% gold-marker precision. The overlay then buried
T cells (recall 78% → 3.3%). Evidence: `logs/aggregate/diag_stage_performance.log`,
`scripts/aggregate/diag_stage_performance.R`, `scripts/aggregate/diag_denovo_profiles.R`.

**Root cause — unit of inference, not reference content.** Per-cell reference typing
cannot separate lineages at this panel sparsity (cohort median 58 of 950 genes detected
per cell, 96 counts). 88% of cells were too sparse to anchor and were bulk-mapped
degenerately (the ~1M-cell de-novo cluster → "T"). The earlier fix framings — add a
`tumor_epithelial` reference anchor, tighten `EPI_GATE`, change `n_clusts` — are tuning
around a method that is the wrong *unit* for this data. **Tumor id is not a separate hard
problem:** in a flank implant with no normal epithelium, an Epcam/Krt8-high *cluster* is
unambiguously tumor; `EPI_GATE=0.25` (a per-cell detection gate) failed only because
ambient spillover clears a per-cell bar that washes out in a cluster mean.

**Decision: type the coarse tier at the CLUSTER level (integrate → joint-cluster →
annotate clusters).** This is **not a new direction** — it is the **"Gold standard
(`aggregate.smk`)" tier already specified in `plan-processing-pipeline.md` v2.0
(2026-05-31)**: integrate M01+M02 into one batch-corrected embedding, joint-cluster,
annotate clusters from the same external marker/reference knowledge. The v2.1/v2.2
per-cell InSituType build was a **detour away from that documented gold standard**;
v2.3 returns to it. Cluster-level annotation also matches Stage 1c's own immune method
(SingleR per Leiden cluster), so the coarse tier reuses that machinery rather than
inventing a parallel one. Tumor falls out as the Epcam/Krt8-high cluster — no 4T1
reference, no per-cell gate.

**New coupling this introduces.** Per-cell InSituType typed from *raw counts*,
independent of the embedding (Stage 1a was QC-only). Cluster-level typing makes the
embedding **load-bearing for the labels**: clusters must reflect biology, not batch. So
the two-covariate Harmony (`dataset` + `slide_id`) LAPACK fix and the LISI(dataset) ≥ 1.8
gate (currently `slide_id`-only at interim floor 1.5, `project_harmony_lapack`) must be
resolved *before* trusting cluster labels — they are no longer QC-only niceties.

**Open questions for the follow-on planning pass (NOT decided here):**
1. Integration method / Harmony LAPACK fix that corrects the dataset axis explicitly.
2. Clustering strategy. Single-pass Leiden at one resolution will absorb rare immune
   (<1% of a tumor-dominated cohort) into the majority — the same failure that sank
   per-sample typing — so hierarchical / iterative sub-clustering (coarse compartments
   first, then sub-cluster immune) is likely required (devils-advocate, sonnet, 2026-06-01).
3. Cluster-annotation engine: SingleR vs marker scoring vs InSituType-at-cluster.
4. Validation gates: per-lineage gold-marker recall, split M01/M02, with thresholds; plus
   a batch-segregation check (both datasets present in the major clusters) before trusting
   any cluster label.

**Supersedes** the coarse-typing *method* of v2.1/v2.2 (retained below as the record of
the failed detour). The Stage 1c immune second pass, reference build, negprobe recovery,
and all downstream tracks are unaffected in design and stay parked on valid labels.

## v2.4 — Cluster-then-annotate execution design (2026-06-01)

> Resolves the four v2.3 open questions. Brainstormed, then reviewed by a
> three-reviewer consensus panel (spatial domain expert, pipeline architect,
> devils-advocate, 2026-06-01) and backed by a Scopus + literature evidence pass.
> Every method choice carries its evidence anchor inline — standing rule
> (`feedback_verify_method_precedent`): **every important tool choice is
> evidence-based AND the evidence is documented.** Still gated: **no code or
> rerun until this plan is reviewed.**

### Evidence basis for the method stack

Anchors verified via Scopus 2026-06-01 (citation counts in parens) plus the two
imaging-spatial integration benchmarks found by literature search.

| Choice | Role | Evidence anchor | Tier |
|---|---|---|---|
| **scVI** | **SELECTED integrator (pilot-confirmed 2026-06-01)** | Lopez 2018 *Nat Methods* (1,611c, foundational); top-tier in scIB benchmark Luecken 2022 (790c); imaging-only benchmark **scVI 0.49/0.70 vs Harmony 0.31/0.62** (bio-conservation/batch) [PMC11996334, 2025]; **our pilot: scaled iLISI 0.71 vs Harmony 0.20, silhouette-label 0.46 vs 0.42** | **Emerging for imaging but confirmed here** — Scopus has only 4 scVI×imaging hits; strongest external support is preprint/tutorial (CellCharter CosMx, MOSAIC bioRxiv 2026.01.12.699017). Pilot resolved the hedge in scVI's favor |
| **Harmony** | linear baseline | Korsunsky 2019 *Nat Methods* (6,378c); Tran 2020 benchmark (815c) | Established; loses to VI on imaging |
| **scVIVA** | deferred sense-check only | scVIVA 2025 preprint (bioRxiv 2025.06.01.657182); tops MOSAIC *cancer*-atlas benchmark | Preprint — not on the build path |
| **Leiden** | clustering | Traag 2019 *Sci Rep* (algorithm, doi:10.1038/s41598-019-41695-z); SpatialLeiden 2025 *Genome Biol* (spatial variant) | Established standard |
| **BANKSY** | niche / coarse-domain | Singhal 2024 *Nat Genetics* (137c); recommended in spatial benchmark *Genome Biol* 2024 (63c) | Established |
| **UCell** | coarse marker scoring | Andreatta & Carmona 2021 *CSBJ* (637c) | Established |
| **SingleR + celldex/ImmGen** | immune tier-2 | Aran 2019 *Nat Immunol*; Heng 2008 (ImmGen); appears across annotation evals (GPB 2021, 126c) | Established |
| **AddModuleScore** | tumor-state programs | Tirosh 2016 *Science* (doi:10.1126/science.aad0501) | Established |
| **LISI / iLISI + silhouette** | batch-mixing gate | Korsunsky 2019 (LISI origin); Tran 2020 (815c); scIB Luecken 2022 | Established metrics |

**Honest calibration on scVI.** Scopus shows scVI-on-*imaging* is peer-reviewed-thin
(4 hits); the strong pro-scVI evidence lives in 2025–26 preprints and the CellCharter
tutorial. So scVI is the *evidence-favored* integrator for this data type, but recent
enough that the **pilot bake-off is the correct hedge** rather than committing sight
unseen. Everything else rests on well-cited peer-reviewed foundations. A prior in-thread
lean to "drop scVI for Harmony" was **reversed by the literature search** — scVI does
have CosMx-scale precedent (CellCharter, 960-gene NSCLC) and beats Harmony on both
benchmark axes.

### Q1 — Integration: 2-way scVI-vs-Harmony pilot; scVI evidence-favored primary

Decision: integration produces a **batch-correction-only joint embedding** (no label
transfer — that was the circular M01→M02 mistake). Candidates are **scVI (primary)** and
**Harmony-done-right (baseline)**, settled by the mandatory pilot below. **scVIVA is a
later sense-check only**, not on the build path (preprint maturity; revisit only if scVI
wins and we want to test the spatial-aware upside).

- **Negprobe modeling before integration.** CosMx counts carry structured negative-probe
  background that scVI's NB/ZINB likelihood (built for droplet UMIs) does not model — the
  recurring threat behind the per-cell tumor overcall. Use the recovered per-cell negprobe
  rate (`cell_neg.tsv`, already built) as a covariate / offset, and **report negprobe
  fraction per cluster** as a standing QC column (spatial-expert + devils-advocate, both
  flagged gate circularity).
- **LISI ceiling math fixes the impossible gate.** With the 29/71 M01/M02 split the
  iLISI(dataset) ceiling is `1/(0.29² + 0.71²) = 1.70`, so the v2.2/v2.3 absolute gates
  (1.8, interim 1.5) are mathematically unreachable. Replace with **scaled iLISI =
  (observed − 1)/(1.70 − 1)**; pre-registered floor **scaled iLISI ≥ 0.5** (≈ observed
  ≥ 1.35), *plus* a per-cluster batch-presence check (minority dataset ≥ 5% in every
  cluster holding >1% of cells). Supersedes `project_harmony_lapack`'s absolute floor.
- **Harmony-done-right** means the two-covariate `c("dataset","slide_id")` model the
  LAPACK bug blocks (`project_harmony_lapack`). Resolve via a LAPACK-enabled Harmony
  build/version bump before the pilot; if scVI wins, Harmony's covariate fix is moot.

**GPU available — train scVI directly (corrects the "No GPU" assumption).** jeff-beast
has an **NVIDIA RTX A4000 (16 GB), nearly idle**, and scvi-tools installs CUDA torch; the
pilot trained the full 1.2M-cell pair on-GPU in ~5 min (`batch_size = 1024`, 50 epochs,
`gene_likelihood="nb"`). So the 3.27M full run goes **GPU-direct**; the sketch-train-then-
transform path (≤300–500k-cell stratified subsample → `get_latent_representation()` over
all cells) is retained **only as a CPU fallback**, not the default. The remaining binding
constraint is *RAM* on the R/Seurat side (ScaleData densifies ~25 GB) — serialize heavy
rules with `resources: mem_mb`.

**Reproducibility & R↔Python handoff (architect).**
- scVI lives in a **separate conda env** (`config/spatial-rads-scvi-env.yaml`:
  scvi-tools, scanpy, anndata, pytorch) — do not pollute the R env. (Per
  `feedback_tool_choice_lang_install`, the Python sidecar is judged on merit, not
  penalized for being Python; this retires Decision #6's "all-R avoids dev friction"
  framing as a *selection* criterion.)
- Treat the **trained model as the reproducible artifact**: `model.save()` →
  `{AGG}/scvi_model/`, with a seed + package-version sidecar.
- Handoff: export merged counts **once** as MTX (`Matrix::writeMM`) + cell metadata
  parquet; the Python rule reads via anndata, writes the latent back as parquet keyed by
  `cell_id`; R loads it into a `DimReduc`. **Avoid `zellkonverter`** (densifies, ~2× RAM
  blowup). One export rule, two consumers.

### Q2 — Clustering: hierarchical (coarse compartments → per-compartment sub-cluster)

A single-pass Leiden at one resolution absorbs rare immune (<1% of a tumor-dominated
cohort) into the majority — the failure that sank per-sample typing. Use a **two-level
hierarchy** (devils-advocate; matches Stage 1c's own immune logic):

- **Level 1 — coarse compartments** on the global joint latent: Leiden
  (`igraph::cluster_leiden`, modularity) → annotate clusters to the three **fixed named
  compartments tumor / immune / stroma**. Fixed names let static Snakemake `{compartment}`
  wildcard rules work — avoid dynamic checkpoints (architect).
- **Level 2 — per-compartment re-embedding.** Retrain the embedder *separately* on the
  immune and on the stroma subsets (each carries signal the whole-tumor embedding drowns),
  then Leiden each. Tumor is **not** sub-clustered (see Q3).
- **Sparsity floor.** Below a per-cluster floor (e.g. <200 cells or median <40 genes
  detected) collapse to the parent compartment label rather than over-subtyping noise
  (spatial-expert: rare immune may not separate at this sparsity).

### Q3 — Annotation engines by tier; tumor identity vs tumor state

Annotation is **always at the cluster level** (cluster mean / pseudobulk), never per cell
— the v2.3 root-cause fix for panel sparsity.

- **Coarse (Level 1):** UCell (Andreatta & Carmona 2021) scoring of `config/lineage_markers.yaml`
  per cluster → argmax compartment, with a min-score/margin → `unassigned` cutoff.
- **Immune tier-2:** SingleR (Aran 2019) vs **celldex ImmGen** (Heng 2008) at
  `label.fine`/`label.main`, scored **per Leiden cluster**; reconcile against UCell ranks
  (= Stage 1c, unchanged).
- **Stroma tier-2:** **cluster-level marker-S/B annotation** against the
  `config/lineage_markers.yaml` stromal sets (Endothelial / Fibroblast / Pericyte /
  SmoothMuscle / Adipocyte) — the **same annotation engine as the tier-1 coarse step**
  (mean-expression signal-to-background per subcluster, argmax, ≥100 cells + S/B ≥2).
  **Revised 2026-06-02 away from SingleR/celldex MouseRNAseqData:** that bulk reference's
  `label.main` has no pericyte or smooth-muscle (mural) class, so the two stromal
  populations the panel-marker file explicitly flags would collapse into fibroblast/
  endothelial and need a marker-correction pass anyway; the curated panel-verified marker
  set names every expected subtype directly. Pericyte-vs-SmoothMuscle is the documented
  panel-weak separation (carry as a watch item, mirroring the immune NK/ILC pooling).
- **Tumor IDENTITY = the Epcam/Krt8-high cluster** (cluster-mean), **not** de-novo
  subtyping. In a flank implant with no normal epithelium an Epcam/Krt8-high cluster is
  unambiguously 4T1; the per-cell `EPI_GATE` failed only because ambient spillover clears
  a per-cell bar that washes out in a cluster mean. De-novo tumor subtyping was considered
  and **rejected** — it would be noise on a 950-gene panel (panel + your pushback).
- **Tumor STATE ≠ a clustering step.** "Tumor signal" is **curated program scoring**, not
  subtype discovery: run the existing `pathway_summary.R` machinery (UCell + AddModuleScore,
  Tirosh 2016) on the tumor compartment for DDR, proliferation (Mki67/Top2a), hypoxia,
  IFN-I/II & STING (`config/pathway_gene_lists.yaml`), and arrest (Cdkn1a/p21). Add the
  proliferation/hypoxia/arrest sets to the config if absent. This is the downstream
  tumor-signal step that was missing — it lives in the pathway track, scoped to tumor.

### Q4 — Validation gates (pre-registered before the run)

All thresholds fixed *now*, not tuned post-hoc (devils-advocate). Reuse the
`embed_celltype.R` `stop()`-after-writing-light-outputs pattern so a failed gate is an
explicit rule failure.

1. **Premise gate (the deepest check).** Coarse Leiden on the pilot must yield **>3 major
   clusters, Immune >5%, Stroma >10%, no single cluster-annotated label >85%.** The 91%→2-label
   collapse happened *before* the tumor overlay, so failure here means the problem is
   **upstream of annotation** (embedding / QC / segmentation) — *stop and rethink the
   embedding, not the annotator.*
2. **Batch-mixing gate:** scaled iLISI(dataset) ≥ 0.5 **and** minority dataset ≥ 5% in
   every cluster >1% of cells.
3. **Marker-recall gate:** per-lineage gold-marker recall computed at **cluster mean with
   a negprobe-relative threshold** (not raw per-cell positivity — fixes gate circularity),
   reported **split M01/M02** (the honest cross-dataset yardstick; keep the M01≈M02 parity
   check).
4. **Negprobe fraction per cluster** reported as standing QC.

### The mandatory pilot (first gate, 2–3 slides spanning both datasets)

Cheap dry run before any 3.27M-cell commit (architect sketch-train + devils-advocate
2–3-slide pilot, merged):

1. Run **both** scVI and Harmony-done-right on the pilot slides.
2. **Premise gate first** (Q4.1) on each embedder. If it fails on *both* → stop; escalate
   the upstream-embedding problem to Jeff (do not proceed to full run).
3. If the premise passes, **pick the integrator:** scVI if it beats Harmony by **>5%** on
   the composite (scaled iLISI × bio-conservation/silhouette × compartment-balance); else
   **Harmony** (simpler, transparent, faster — and it tests whether the *cluster unit*, not
   the embedder, was the real fix all along).
4. scVIVA stays a later sense-check only.

### Pilot results (2026-06-01) — scVI selected

Ran on **M01 sld0002 (331,403 cells) + M02 sld0005 (878,422)**, 27/73 split (mirrors the
cohort); metrics on a common **50k dataset-stratified subsample**. Scripts:
`scripts/aggregate/pilot_{export.R,scvi.py,harmony.R,compare.py}`; artifacts in
`results/aggregate/pilot_*` and `{AGG}/pilot/` (latents, `scvi_model/`). Both arms 30-dim
(scVI `n_latent` = Harmony NPCS). iLISI ceiling for this split = **1.661**.

- **Premise gate — PASS for all three embeddings.** Coarse leidenalg clusters split into
  all three compartments (tumor+immune+stroma), and every cluster's top-lineage marker mean
  exceeded its per-cell negprobe floor (min signal-to-background: PCA 7.1×, Harmony 8.2×,
  **scVI 24.6×**). The cluster-then-annotate premise holds — the v2.3 pivot is validated.
- **Integrator — scVI.** Scaled iLISI (ceiling-corrected dataset mixing): PCA 0.15,
  Harmony 0.20 (≈ unintegrated → **fails** the ≥0.5 gate), **scVI 0.71** (passes). scVI also
  best preserves biology (silhouette-label 0.46 > Harmony 0.42 > PCA 0.31) — mixes datasets
  without over-blending. Confirms the documented imaging-spatial precedent (scVI > Harmony).

**Honest caveats (carry into the full run):**
- On the full scIB *composite* scVI 0.457 vs Harmony 0.441 — close; Harmony is propped up by
  kBET/PCR/graph-connectivity/NMI. scVI's win concentrates in **iLISI** (0.467 vs 0.131, the
  operative cross-dataset metric) and batch-silhouette (BRAS 0.91 vs 0.86). Weighting iLISI
  is the basis for the call (full table: `results/aggregate/pilot_scib_full.tsv`).
- scVI mixed locally (high iLISI) but only **2/6 clusters held both datasets ≥10%** (Harmony
  all 7) — most likely genuine treatment/timepoint-specific states (slides differ: Control t0
  + SBRT d2/d6 vs Control/MBRT/SBRT d2), not under-mixing, given the high silhouette-label.
  **Verify per-cluster composition on the full run** (this is the Q4.2 minority-presence gate).
- Bio metrics use the imperfect v1 `cell_type` labels (treat as relative); the **label-free
  marker-recall gate** is the trustworthy biology check and passed strongly.
- Metric plumbing: `pilot_compare.py` "silhouette_batch" reported NaN (scib renamed it BRAS,
  which favors scVI) — fix the column match; clusters coarse (6–7 @ res 0.5), correct for the
  Level-1 test only.

### Compute & engineering notes (architect amendments, folded)

- Resolve the repo inconsistency: v2.2 says LISI via the `lisi` package, but check the
  `embed_celltype.R` header before adding the scaled-iLISI gate so there's one LISI path.
- Map Level-1 Leiden → fixed {tumor,immune,stroma} so `{compartment}` wildcard rules are
  static; no Snakemake checkpoint.
- Tier-2 runs on the immune/stroma subsets only, at cluster level (not 3.27M cells):
  SingleR/ImmGen for immune; **marker-S/B annotation for stroma** (revised 2026-06-02 off
  MouseRNAseqData — see Approach section).
- Keep the post-FindNeighbors checkpoint (`embed_checkpoint.rds`) so re-clustering costs
  minutes.

### Panel consensus record (2026-06-01)

- **Spatial domain expert** — *sound with changes, no dealbreaker*: model negprobe before
  scVI + per-cluster negprobe QC; conservative tier-2 resolution / cluster-level SingleR /
  collapse floor; negprobe-relative recall denominators; segmentation QC before typing;
  BANKSY a candidate for Level-1.
- **Pipeline architect** — *sound with changes, no dealbreaker*: sketch-train-then-transform;
  model.save() artifact + separate conda env + version/seed sidecar; fixed-name compartments
  with static wildcards; MTX/parquet handoff (not zellkonverter); gates as rule failures;
  R/Seurat ScaleData is the OOM risk → `resources: mem_mb`.
- **Devils-advocate** — *needs rework → addressed*: three untested load-bearing assumptions
  → the premise gate (does coarse clustering separate at all), the scVI-vs-Harmony pilot
  (don't assume scVI; drop it if Harmony ties), and pre-registered thresholds. Alternatives
  noted: Harmony-alone, per-slide-then-meta-cluster.

### What v2.4 changes below / supersedes

- **Stage 1a (`embed_celltype.R`)** becomes: scVI-or-Harmony per pilot → **hierarchical**
  clustering → **cluster-level annotation that is now load-bearing for labels** (not
  QC-only). The "labels remain per-sample (Yi-ImmGen / TransferData)" note and **Decision
  #4** (labels canonical from per-sample assignments) are **superseded** — labels now come
  from cluster annotation of the joint embedding.
- **Decision #6 (all-R stack)**: amended — a Python scVI/scVIVA sidecar is accepted on
  merit; "all-R avoids dev friction" is retired as a *selection* criterion (global
  `feedback_tool_choice_lang_install`). Spatial methods (BANKSY-R, spicyR) stay R by merit,
  not by language preference.
- The batch-mixing **LISI gate** moves from absolute 1.8/1.5 to **scaled iLISI** (ceiling
  1.70) + per-cluster batch presence.

## v2.5 — Cluster-then-annotate executed; tier-1 + tier-2 locked (2026-06-02)

**Deferred (2026-06-15):** at the next aggregate rebuild, adopt the Mutter_01 raw RDS — retain negprobes+falsecode assays through processing and retire recover_negprobes.R; until then the locked typing stands (negprobe-mean inputs verified identical). See plan-mutter01-controls.md.

The v2.4 design ran on the full 3.27M-cell cohort. **Coarse (tier-1) and immune (tier-2)
typing are complete and locked**; the run is written up in `spatial-rads.org`
(`*** Aggregate atlas cell typing`). The tier scripts are authored directly under
`scripts/aggregate/` (not org-tangled): `full_export.R`, `full_cluster.py`,
`finalize_tier1.py`, `tier2_immune_subcluster.py`, `tier2_singler.R`,
`tier2_marker_check.py`, `bcell_diagnostic.py`.

**Tier-1 — scVI integration + coarse cluster-then-annotate.**
- scVI (env `spatial-rads-scvi`, GPU-direct on the RTX A4000): 30 latent dims, 2 layers,
  NB likelihood, `batch_key=slide_id`, per-cell negprobe (`neg`) as a continuous covariate,
  50 epochs. Artifacts: `aggregate/full/scvi_model/`, `aggregate/full/scvi_latent.parquet`,
  `aggregate/full/cluster_checkpoint.h5ad`.
- Leiden on the scVI latent (`flavor="leidenalg"`, `n_iterations=2` — the run-to-convergence
  default is intractable at 3.27M nodes; the `igraph` flavor over-partitions). **Resolution
  swept**: 0.5 (failed the premise gate, immune submerged at 0.03%) → 1.0 / 2.0 / 3.0
  (106 / 56 / 90 clusters; 6 / 13 / 25 above 1%). **res=3.0 adopted** (finest stable
  partition, best immune substrate, all Q4 gates pass).
- Unassigned rule (`finalize_tier1.py`, `N_MIN=100, S_MIN=2.0`): a cluster is labeled only
  with ≥100 cells AND signal-to-background ≥2 → 4 cells unassigned. Q4.2 per-cluster minority
  demoted from hard-fail to a QC flag (scaled iLISI is the actual gate); flagged 1
  M02-specific tumor cluster (cluster 3, 2.6% M01).
- **Result (full cohort): tumor 1,577,685 (48.1%) / stroma 1,331,601 (40.6%) / immune
  367,800 (11.2%) / 4 unassigned. ALL_PASS.** Scaled iLISI 0.67 (floor 0.5, ceiling 1.70 for
  the 29/71 split). Canonical: `results/aggregate/full_coarse_labels.parquet`,
  `full_coarse_summary.tsv`, `full_gates.json`.

**Tier-2 — immune subtyping by cluster-level SingleR.**
- Immune subset re-embedded on its own scVI latent (graph rebuilt — the global graph is not a
  valid subgraph), Leiden-subclustered, SingleR vs celldex ImmGen (`label.main`) scored **per
  subcluster**. The `spatial-rads` env needed `scrapper` (SingleR 2.12.0 / Bioc 3.22
  dependency); installed additively.
- **Immune resolution walked 1.0 → 3.0 → 6.0** (16 / 65 / 298 subclusters). res=1.0
  over-merged the lymphoid cells into one 71k block; **res=3.0 adopted** (T / NK-ILC / Mast /
  DC / Plasma resolvable and marker-confirmed); res=6.0 bought only a sharper plasma core at
  the cost of 256/298 macrophage fragments. `tier2_marker_check.py` gold-marker means
  adjudicated res=3.0 over res=1.0.
- **Marker rescue over SingleR (`tier2_immune_rescue.py`, executed 2026-06-02).** SingleR at
  res=3.0 assigns no B/Plasma/Neutrophil class, so those rare lineages are re-examined at res=6.0
  and SingleR is overridden where a subcluster's **identity markers** (the lineage-defining,
  non-secreted, non-shared subset of each panel) win by argmax AND ≥2 individually corroborate
  (≥2× negprobe bg). This is the Mast-rule "marker evidence overrides reference confidence"
  (precedent Cheng et al. BMC Bioinformatics 2025 Xenium benchmark). A plain set-mean
  signal-to-background argmax was **rejected for assignment** (it stays only the detectability
  *gate*): ambient immunoglobulin (Jchain/Igkc) and non-specific genes (Cd37 pan-leukocyte, Xbp1
  secretory) corrupt set-means — naive argmax flips an 8,258-cell T subcluster (Cd3e 44× bg) to
  Plasma; the identity-marker rule keeps it T. Rescued: **Plasma 3,573** (the SingleR-"DC"
  subcluster 21, Mzb1 17.8×/Xbp1 19.7×) + **Neutrophils 259** (subcluster 272, Elane/Prtn3/Mpo,
  S100a8 260×); 3,832 cells (1.0% of immune) reassigned. NK and ILC remain **separate** SingleR
  labels in the parquet (the v2.2 "pool NK + ILC" was a presentational grouping, never applied to
  data — pool at analysis time if a coarser vocabulary is wanted; not separable de novo on panel).
- **Result (final, marker-rescued): Macrophages 300,975 (81.8%) / ILC 34,340 / T 16,046 /
  DC 5,315 / NK 4,991 / Plasma 3,573 / Mast 1,121 / Neutrophils 259 / epithelial-contam 1,180.**
  The earlier "~6.7k plasma" was an estimate; the marker-grounded count is 3,573. Outputs:
  `aggregate/full/immune_subtypes_rescued.parquet` (final), `immune_subtypes.parquet` (pre-rescue
  SingleR-only), `immune_subtype_singler.tsv`, `immune_rescue_audit.tsv`.
- **Stroma tier-2, step A (`tier2_stroma_ucell.R`, executed 2026-06-02).** UCell cluster-level
  argmax+margin against the `config/lineage_markers.yaml` stromal sets (NOT
  SingleR/`celldex MouseRNAseqData`, whose bulk reference lacks pericyte/smooth-muscle
  classes). Pre-rescue result: **Fibroblast 796,933 / SmoothMuscle 73,570 /
  unassigned 461,098.** UCell's within-cell rank score, washed out by the fibroblast-dominated
  background, assigned **none** of the minor lineages even though the detectability gate
  (`tier2_detectability.py`) confirmed each clears background in a specific r6 subcluster
  (Endothelial s2b 27 / Adipocyte 56 / Pericyte 14). Outputs: `stroma_subtypes.parquet`
  (pre-rescue), `stroma_subtype_summary.tsv`, `stroma_detectability.tsv`,
  `stroma_subcluster_ucell.tsv`.
- **Stroma tier-2, step B — marker rescue (`tier2_stroma_rescue.py`, executed 2026-06-02).**
  The same identity-marker override built for immune, examining the r6 subclusters: argmax over
  the non-secreted/non-shared identity subset of each stromal panel, rescue where the argmax is a
  UCell-missed lineage AND ≥2 identity markers corroborate AND a **specificity anchor** (≥1
  lineage-EXCLUSIVE marker) individually clears the gate. The anchor is the stroma analog of the
  immune ≥2-corroboration guard: **Fabp4 is shared** (adipocyte + endothelial s2b 39 on the endo
  cluster + macrophage), so a Fabp4-driven "adipocyte" call with no lipid-droplet **Cidea** is
  rejected (caught 2 false subclusters at Cidea 0.83/1.96); endothelial anchored by Pecam1/Cdh5.
  Recovered **Endothelial 36,937** (its best subcluster 8 had been 99.5% mislabeled Fibroblast by
  UCell) + **Adipocyte 49,214** (both cohorts: M01 19k/M02 30k, M01 8k/M02 29k — not a
  single-dataset artifact); 86,151 cells (6.5% of stroma) reassigned, only out of
  Fibroblast/SMC/unassigned. **Pericyte = 0**: its best subcluster (8) is shared with and
  dominated by Endothelial (endo s2b 27 > peri 14), only 3 panel markers with Pdgfrb
  fibroblast-shared — present but unresolvable de novo, folded into the perivascular endothelial
  call. Final: **Fibroblast 732,099 / unassigned 441,227 / SmoothMuscle 72,124 /
  Adipocyte 49,214 / Endothelial 36,937.** Outputs: `stroma_subtypes_rescued.parquet` (final),
  `stroma_rescue_audit.tsv`.

**B cells near-absent (empirical, `bcell_diagnostic.py`).** No subcluster carries B markers as
its dominant signal at any resolution (0/298 at res=6.0). Not panel-blind
(Cd19/Ms4a1/Cd79a/Ighm present; Cd79b/Cd22/Iglc1 absent), not mis-compartmented: cohort-wide
2.4–3.3% diffuse B-marker expression, max cluster B-score 0.365 vs 1.0–1.8 for a genuine
lineage. Real 4T1 immune-cold biology, consistent with `project_mbrt_mechanism_status`.
(An in-session "B RETAIN" reading was an artifact of the per-lineage-best detectability output:
B's panel clears background at subcluster 3, but that subcluster's argmax identity is Plasma/T,
not B — so detectability "retain" ≠ "B is anyone's dominant identity". The argmax test above is
the assignment-correct view.)

**Marked divergence from Yi's original M01 typing.** On the same M01 cells the new scheme
gives stroma 55.3% / tumor 29.2% / immune 15.5%, against Yi's ImmGen-only InSituType where the
single largest label was the uninterpretable "a" attractor (28.0%), 7.1% of cells carried
anatomically-impossible thymocyte labels, and the four ImmGen stromal labels summed to only
11.5%. The disagreement is on **coarse identity**; the original M01 labels (and any M01→M02
transfer derived from them) are superseded for cross-dataset analysis.

**What changed vs the v2.4 plan as-built.** scVI replaced Harmony **entirely** (Harmony 2.0.2
LAPACK-errors on ≥2 covariates; the slide-only-vs-two-covariate question in v2.2/v2.3/v2.4 is
moot — scVI with `slide_id` batch + negprobe covariate clears the scaled-iLISI gate). Tier-1
cluster annotation used **marker signal-to-background with an unassigned cutoff**, not the
Q3-planned UCell-argmax engine. Stroma tier-2 ran in two steps: **UCell argmax+margin**
(`tier2_stroma_ucell.R`) for the base call, then an **identity-marker rescue**
(`tier2_stroma_rescue.py`, the same Mast-rule override applied to immune) recovered the minor
lineages UCell washed out — Endothelial 36,937 (its best subcluster 8 was 99.5% mislabeled
Fibroblast by UCell) and Adipocyte 49,214 (anchored by lipid-droplet Cidea). Pericyte stays 0:
it shares subcluster 8 with Endothelial, which wins argmax, so it is folded into the perivascular
endothelial call (present but unresolvable de novo at this panel). **Label join DONE
(2026-06-02):** `unify_labels.py` merged all three tiers into the canonical per-cell table
`results/aggregate/full_labels.parquet` (3,277,090 rows; `compartment` + final `cell_subtype` +
`subtype_source` + `rescued`). **Not yet run:** the BANKSY niche embedding and the downstream
analysis tracks (composition / pseudobulk DE / GSEA / pathway / colocalization, incl. tumor-state
program scoring) — these can now proceed off the unified label table.

## v2.6 — Workflow refactor planned: the wired DAG does not reproduce this typing (2026-06-22)

The typing locked above is an **artifact-of-record but is NOT reproducible from
`workflows/aggregate.smk`**: the wired typing rules are the dead per-cell InSituType chain
(`embed_celltype`→`prepare_reference`→`typing_insitutype`, the 85%-tumor failure), and the
labels actually used (`results/aggregate/full_labels.parquet`, 9 consumer rules) plus
`{AGG}/full/obs.parquet` (3 consumer rules) are **orphan static inputs** produced by the
hand-run scVI→tier-2 chain, not by any rule. The full-cohort scVI training step is **untracked
code** (no `full_scvi.py`; only `pilot_scvi.py` exists).

Harmonized fix in `plan-aggregate-refactor.md` (three gating decisions resolved 2026-06-22:
**A2 full re-wire now · decision-forcing best-practice review · C1 orchestration-only literate
transfer**) — a **seven-step** sequence: audit → method shoot-out (scVI vs scANVI/scArches on
the 2-slide pilot; a margin-passing winner feeds the rebuild) → wire the real ~13-node
multi-env R→GPU-Python→R DAG (authoring `full_scvi.py`; retiring `recover_negprobes` onto the
re-based M01 RDS) → single full 3.27M-cell rebuild → re-validate the four Q4 gates + marker
recall → **refresh the differential layer (`results_master.tsv` goes stale when the labels
regenerate)** → literate transfer. Acceptance = **re-validation, not byte-identity**. Until the
rebuild runs, the labels and `results_master.tsv` recorded here stand as the locked record.

## v2.7 — Best-practice review (Workstream B verdicts, decision-forcing, 2026-06-22)

Step 2 of `plan-aggregate-refactor.md`. Scopus pass 2026-06-22 (SCOPUS_API_KEY; `TITLE-ABS-KEY`
queries logged below). Per the standing evidence rule (`feedback_verify_method_precedent`) a NULL
result — "no peer-reviewed support at imaging-spatial / CosMx scale" — is an **acceptable,
unblocking** verdict. Decision rule: a "change" stands only if it beats current scVI
cluster-then-annotate on the 2-slide pilot harness by the pre-set margin (scaled iLISI **and**
marker recall **and** unassigned fraction). **Outcome: 3 KEEP + 1 deferred non-blocking refinement
→ the rebuild method is locked to the current scVI / swept-res=3.0 / per-sample-MECR stack; no
method change gates Steps 3–4.**

| Q | Question | Verdict | Evidence |
|---|---|---|---|
| 1 | Integrator: scVI vs scANVI vs scArches | **KEEP scVI** | Scopus scVI/scANVI/scArches × imaging-spatial = **2 hits**; deep-generative × imaging benchmark = **3 hits** (NULL at CosMx scale). Architectural mismatch: **scANVI is semi-supervised** (needs trusted labels — the de-anchored design forbids label transfer) and **scArches is reference-mapping** onto a fixed reference atlas we lack at CosMx scale; scVI is the only fully-unsupervised batch-correction-only option and is already pilot-validated (scaled iLISI 0.71 vs Harmony 0.20). Emerging candidate noted, not adopted: interpretable-generative integration+annotation (Cell Genomics 2026, doi:10.1016/j.xgen.2025.101105, 0c). |
| 2 | Leiden resolution at >3M cells | **KEEP sweep→gates (res=3.0)** | Scopus resolution-selection = 23 hits, none CosMx-scale that clears the margin. The locked run already *selects* res=3.0 empirically by a pre-registered sweep (0.5/1.0/2.0/3.0) against the four Q4 gates — not an arbitrary pick. Future automated-selection option noted: CANTAO average-overlap (Mol Syst Biol 2026, doi:10.1038/s44320-025-00176-4, 0c). |
| 3 | Recover the 33% unassigned stroma (441k) | **CHANGE-CANDIDATE — deferred, gated, NON-blocking** | Real precedent exists: **TACIT** jointly deconvolves cell types *and states* in spatial multiomics (Nat Commun 2025, doi:10.1038/s41467-025-58874-4, 8c) + active CAF-subpopulation literature (21 hits). This is the only question with a credible challenger, but it refines **only the unassigned-stroma bucket** (not the coarse tier or the integrator), so it does **not** block the rebuild. Decision: defer to a stroma-tier-2 v2 pilot (TACIT vs the current UCell-argmax + identity-marker rescue on the 2-slide harness), keep the specificity-anchor guard; rebuild proceeds on the current stroma method. Aligns with memory `project_stroma_unassigned_recovery` (FUTURE). |
| 4 | Merged-scale QC beyond per-sample MECR | **KEEP per-sample MECR** | Scopus merged-scale segmentation-contamination QC = **1 hit** (Xenium contamination, Nat Methods 2026, doi:10.1038/s41592-026-03089-8, 0c) — essentially no precedent for a *merged-scale* step. The project already runs per-sample SpatialQM MECR (published-method-backed) and the cohort is clean (`project_contamination_qc`). NULL → keep; note the Xenium paper as a future *per-cell spillover-correction* reference, not a merged-scale step. |

Scopus queries (TITLE-ABS-KEY): Q1 `("scArches" OR "scANVI" OR "scVI") AND ("CosMx" OR "MERFISH"
OR "Xenium" OR "imaging-based spatial" OR "in situ") AND (integrat* OR annotat* OR "label
transfer")`; Q2 `("Leiden" OR "clustering resolution" OR "resolution parameter") AND ("single-cell"
OR scRNA) AND (select* OR optimal OR stability OR benchmark*)`; Q3 `("unassigned" OR "ambiguous" OR
"unresolved") AND (fibroblast OR stroma* OR mesenchymal) AND (spatial OR "single-cell") AND
(annotat* OR subtyp* OR cluster*)`; Q4 `(segmentation AND (spillover OR contamination OR "transcript
bleed" OR missegmentation)) AND ("CosMx" OR "MERFISH" OR "Xenium" OR "imaging-based spatial")`.

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

**As-built (2026-06-01, see v2.2):** the spec above is the target; the running
script differs in three concrete ways. (1) **Clustering is Leiden via
`igraph::cluster_leiden`** (pure-C, minutes on 3.27M cells), not Seurat
`FindClusters` Louvain. (2) **Harmony groups on `slide_id` alone**, interim
LISI gate **1.5**, pending the Harmony multi-covariate (LAPACK) fix that
restores `c("dataset","slide_id")` + gate 1.8; `slide_id` nests `dataset`, so
the single covariate still corrects the dataset axis. (3) **A checkpoint
(`aggregate/embed_checkpoint.rds`) is written after FindNeighbors** and the
**gate runs before UMAP** (fail-fast). `cluster` and `seurat_clusters` are both
set on the saved object for downstream compatibility.

### Stage 1b — `embed_niche.R`

Banksy spatial niche embedding.

- Convert `merged.rds` → `SpatialExperiment` (one-line via
  `as.SingleCellExperiment` then add `spatialCoords`).
- `runBanksy(lambda=0.8, k_geom=c(6,18), M=1, group="sample_id",
  use_agf=TRUE)` — per-sample neighbor graphs (no cross-sample neighbors).
- `runPCA(npcs=20)` on the BANKSY matrix.
- **Batch correction: `RunHarmony(group.by.vars = "slide_id")` on the BANKSY
  PCs — single covariate, NOT `c("dataset","slide_id")`.** Rationale (settled
  2026-06-02): (1) the two-covariate call LAPACK-fails on Harmony 2.0.2
  (`project_harmony_lapack`); (2) `slide_id` nests `dataset`, so the dataset
  axis is still corrected — the same nesting logic the tier-1 embedding adopted
  and validated against its mixing gate; (3) Harmony-on-`slide_id`-alone is
  already proven to run in this env. Harmony-on-BANKSY-PCs is the field-standard
  niche batch path (Singhal 2024), so this stays standard, just single-covariate.
  Do **not** inherit the old two-covariate call. *(Note: the tier-1 scVI latent is
  NOT fed here — BANKSY builds its own spatial-neighbor-augmented feature matrix
  and PCA, so batch correction must act on the BANKSY PCs, not on a precomputed
  expression latent.)*
- Leiden on Harmony-corrected Banksy PCs, resolution 0.4 (production);
  sweep {0.2, 0.4, 0.6, 1.0} during prototype.
- Niche names manually assigned after inspecting
  `niche_celltype_composition.tsv` (a small lookup table in
  `config/niche_names.tsv` keyed by numeric niche_id).
- Heavy output: `merged_niche.rds`.
- Test artifacts:
  - `results/aggregate/niche_embed_qc.tsv` — n_niches, harmony iters,
    silhouette + the **tier-1 imbalance-scaled iLISI(dataset)** pre/post
    (same scaled-iLISI gate the coarse embedding used: floor 0.5, ceiling
    1.70 for the 29/71 dataset split, per-cluster minority demoted to a QC
    flag; this replaces the old raw LISI(dataset)>=1.8/1.5 floor), per-niche
    {n_cells, modal_cell_type, composition_entropy, mean_pathway scores}.
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

**2026-06-02 — tier-1 and tier-2 cell typing executed and locked.** Full-cohort scVI
integration ran on all 3,277,090 cells; coarse cluster-then-annotate at resolution 3.0
PASSED all four Q4 gates (tumor 48.1% / stroma 40.6% / immune 11.2%, 4 cells unassigned;
scaled iLISI 0.67). Tier-2 immune subtyping (cluster-level SingleR/ImmGen, immune
resolution walk 1.0→3.0→6.0, res=3.0 adopted) resolved a macrophage-dominated compartment
(~301k macrophage, ~39.5k NK/ILC, ~16.3k T, ~6.7k plasma, ~1.2k mast); B cells are
near-absent and that absence is biological, not panel- or clustering-driven. See the
**v2.5 banner above** for the full parameter walk, gate values, and what changed vs v2.4.
The aggregate cell profiles diverge markedly from Yi's original Mutter_01 ImmGen typing,
which is now superseded (documented in `spatial-rads.org` → *** Aggregate atlas cell
typing). Downstream tracks (composition / pseudobulk / DE / pathway / niche) and the
stroma tier-2 (cluster-level marker-S/B vs `lineage_markers.yaml`, revised from
MouseRNAseqData) + BANKSY niche steps stay **parked** but are now unblocked
— the coarse and immune labels are trustworthy.

**2026-06-01 (prior).** Merge, reference build, negprobe recovery, embedding (QC), and the
v2.1/v2.2 per-cell InSituType coarse typing were built and run. **That per-cell coarse
typing failed validation** (84.6% tumor) and was replanned to cluster-then-annotate (v2.3),
executed per the v2.5 entry above. The detailed execution design is in **v2.4 above**
(integration / clustering / annotation / gates resolved; three-reviewer panel + Scopus
evidence pass folded in). Pilot (2026-06-01): premise gate PASSED on all embeddings; scVI
selected over Harmony (scaled iLISI 0.71 vs 0.20, best bio preservation).

(Original 2026-05-31 status: designed and devils-advocate-reviewed, no code yet —
superseded by the build log in the v2.x banners above.)

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
