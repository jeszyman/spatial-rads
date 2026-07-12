# Aggregate workflow — harmonized refactor plan (decisions resolved 2026-06-22)

> **Status: Steps 1–5 EXECUTED and committed 2026-06-22.** The real scVI cluster-then-annotate
> typing chain is wired (`full_scvi.py` authored, dead InSituType rules deleted); the full
> 3.27M-cell rebuild reproduced the locked labels bit-identically (all four Q4 gates pass),
> exceeding the re-validation acceptance bar. **Remaining: Step 6** (refresh `results_master.tsv`
> on the rebuilt labels), **Step 7** (literate transfer), and the `recover_negprobes` retirement
> (still blocked — 0% barcode overlap between the raw-RDS negprobes and the scored cells). A devils-advocate panel
> (2026-06-22, run `wf_c66fbfa4-c5f`) split the original *investigate → update → transfer*
> bundle into three decoupled workstreams (A honesty / B best-practice / C literate). The
> user resolved the three gating decisions 2026-06-22:
> **A2 (full re-wire now) · B decision-forcing · C1 (orchestration-only + document).**
> This file now harmonizes those choices into one sequenced execution plan. Companion to
> `plan-aggregate.md` (the as-built typing record).

## Decisions resolved (2026-06-22)

| # | Decision | Chosen | Consequence |
|---|---|---|---|
| 1 | Workflow-honesty target | **A2 — full re-wire now** | The locked typing is *rebuilt* as a real DAG, not quarantined. Folds in the deferred M01 raw-RDS re-base. Acceptance = re-validation, not byte-identity. |
| 2 | Investigation stance | **Decision-forcing** | B's four questions get keep/change verdicts; a "change" that beats current scVI on the pilot harness by the pre-set margin **feeds the rebuild** (so we rebuild once, with the winning method). |
| 3 | Literate transfer | **C1 — orchestration-only + document** | `aggregate.smk` + thin R wrappers go to org; heavy GPU-Python / tier-2 scripts are *documented* (pointer to `scripts/aggregate/`), not tangled. Keeps the mature-code exemption. |

**Why these cohere.** A2 + decision-forcing B means B's method shoot-out must complete *before*
the integration method is committed — otherwise we rebuild twice. The harmonized sequence
therefore runs B's cheap 2-slide shoot-out first, locks the integrator/clustering choice, then
does the single full 3.27M-cell rebuild. The literate transfer runs last (it depends on the rules
existing and verifying).

---

## Context — the workflow-reality gap (grounded 2026-06-22, verified against live `aggregate.smk`)

`workflows/aggregate.smk` does **not** reproduce the locked results it appears to:

- **Dead chain wired.** The typing rules in the DAG are the **superseded, failed per-cell
  path**: `embed_celltype` (Harmony) → `prepare_reference` → `recover_negprobes` →
  `typing_insitutype` (InSituType) → `merged_typed.rds`. That path validated at **85%
  mistyped tumor** and was abandoned. A clean `snakemake -s workflows/aggregate.smk` would
  silently re-type with it, no error.
- **Orphan static inputs (consumed, produced by NO rule) — there are TWO, not one:**
  - `results/aggregate/full_labels.parquet` — input to **9 rules** (`composition`,
    `pseudobulk_build`, `pathway_summary`, `panel_coverage`, `celltype_qc`, `niches`,
    `spatial_mixing`, `myeloid_polarization`, `concordance_m01_m02`). Producer:
    `unify_labels.py` (hand-run). **Confirmed on disk** (25 MB, 2026-06-02).
  - `{AGG}/full/obs.parquet` — input to `composition`, `niches`, `spatial_mixing`.
    Producer: `full_export.R` (hand-run). **Confirmed on disk** (21 MB, 2026-06-01).
- **`scvi_latent.parquet` has no committed full-run producer (new finding).** `full_cluster.py`
  *reads* it; `full_export.R` does not write it; only `pilot_scvi.py` exists under
  `scripts/aggregate/`. The full-cohort scVI training was run off-DAG, so the integrator step —
  the most expensive and method-defining node — is currently **untracked code**. The rebuild must
  author it (`full_scvi.py`) as a first-class rule.
- **`merged_typed.rds` is read for counts only** by `pseudobulk_build`, `pathway_summary`,
  `panel_coverage` (they label from `full_labels.parquet`); its counts equal `merged.rds`'s.
  Re-point them to `merged.rds`.

The hand-run chain that actually produced the locked labels (grounded against the scripts +
`{AGG}/full/`):

```
merge.R ─► merged.rds
recover_negprobes(.R, parquet-sourced) ─► cell_neg.tsv
full_export.R(merged.rds, cell_neg.tsv) ─► full/mtx/{counts.mtx,features,barcodes}, full/obs.parquet
[UNTRACKED scVI train] ─► full/scvi_latent.parquet, full/scvi_model/
full_cluster.py(full/, lineage_markers.yaml) ─► cluster_checkpoint.h5ad, per-res labels/gates/recall, full_sweep_summary
finalize_tier1.py ─► full_coarse_labels.parquet, full_coarse_summary.tsv, full_gates.json
  ├─ immune: tier2_immune_subcluster.py ─► immune_checkpoint.h5ad, immune_mtx, immune_subclusters.parquet
  │          tier2_singler.R ─► immune_subtype_singler.tsv, immune_subtypes.parquet
  │          tier2_immune_rescue.py ─► immune_subtypes_rescued.parquet
  └─ stroma: tier2_stroma_subcluster.py ─► stroma_checkpoint.h5ad, stroma_mtx, stroma_r*_subclusters.parquet
             tier2_stroma_ucell.R ─► stroma_subtypes.parquet (+ ucell tsvs)
             tier2_stroma_rescue.py ─► stroma_subtypes_rescued.parquet
unify_labels.py ─► results/aggregate/full_labels.parquet   (the canonical per-cell table)
```

Helpers (cited, not re-implemented): `tier2_detectability.py` (= SpatialQM signal-to-background
gate), `spatialqm_metrics.R` (MECR / pseudobulk-cor), `tier2_marker_check.py` (resolution
adjudication), `bcell_diagnostic.py` (B-cell absence check).

## Panel-forced corrections (still binding)

1. **"Locked" ≠ bit-identical.** scVI GPU training + Leiden `igraph` are not byte-reproducible
   across library/driver versions. Reproducibility = the committed `full_labels.parquet` as
   artifact-of-record + a frozen env/seed record + **re-validation of the four pre-registered
   Q4 gates + marker recall** on any re-run. **Not** byte-identity. (This is the rebuild's
   acceptance test.)
2. **Literate transfer of the heavy chain conflicts with the project's own guideline** ("mature/
   complex code may not benefit from a literate approach"). C1 honors it — orchestration to org,
   heavy scripts documented-not-tangled.

---

## Harmonized execution sequence

Seven steps. Envs: `basecamp` (snakemake driver), `spatial-rads` (R: Seurat/SingleR/celldex/
UCell — needs `scrapper`), `spatial-rads-scvi` (GPU Python: scvi-tools/scanpy/leidenalg/torch).
BLAS pinned at workflow level; every R rule `threads: 1`; GPU rules `resources: gpu=1`.
`TMPDIR=/mnt/data/projects/spatial-rads/tmp` always.

```
1 audit ─► 2 method shoot-out ─► 3 wire DAG ─► 4 full rebuild ─► 5 re-validate ─► 6 refresh results ─► 7 literate
           (decision-forcing,     (method-agnostic   (adopts the    (4 Q4 gates    (results_master.tsv   (org orchestration
            sets the method)       skeleton + the     shoot-out      + recall)      on rebuilt labels)    + Methods pointer)
                                   method commit)     winner)
```

### Step 1 — Audit + provenance capture (cheap; do first; blocks nothing downstream)
Establishes the exact rule I/O contract the wiring must reproduce, and the manifest of record.
- **1.1** For each `merged_typed.rds` consumer, confirm from its `commandArgs` that it ignores
  the `cell_type` column (labels come from `full_labels.parquet`) → re-point to `merged.rds` is
  safe. *Acceptance:* 1-page live/dead table; every consumer shown count-only.
- **1.2** Record the exact hand-run order + commit SHA + seeds + conda env versions that produced
  the current locked `full_labels.parquet` → `typing_provenance.md`. This is the **baseline** the
  rebuild's re-validation compares against (composition, gate values).
- **1.3** Pin down the **untracked scVI step**: recover the exact parameters used for the locked
  run (30 latent dims, 2 layers, NB likelihood, `batch_key=slide_id`, neg covariate, 50 epochs,
  seed) from logs / `scvi_model/`, so `full_scvi.py` reproduces them.

### Step 2 — Method shoot-out (Workstream B; decision-forcing; sets the method before any full run)
Bounded research pass — same shape as this session's processing-QC research-swarm. **Four fixed
questions (do not expand):**
1. scVI vs **scANVI** vs reference-mapping/**scArches** for cross-dataset CosMx typing.
2. Leiden resolution selection at >3M cells (current: res=3.0 adopted by sweep).
3. Recovering the **33% unassigned stroma** (441k cells) — precedent-first (memory
   `project_stroma_unassigned_recovery`); keep the specificity-anchor guard.
4. Whether a **merged-scale contamination/QC** step adds value beyond the per-sample MECR now in
   `processing.smk` (memory `project_contamination_qc`).
- **Evidence rule:** every finding cited (Scopus/lit per `feedback_verify_method_precedent`). A
  NULL result ("no peer-reviewed support at CosMx scale") is an **acceptable, unblocking** verdict.
- **Decision rule (pre-committed):** a "change" verdict stands **only if** the alternative beats
  current scVI cluster-then-annotate on the existing 2-slide pilot harness (`pilot_*`) by a pre-set
  margin on the named metric set — **imbalance-scaled iLISI AND marker recall AND unassigned
  fraction**. Winner → its config is what Steps 3–4 wire and run. No winner → Step 4 reproduces the
  current scVI stack (still a valid reproducibility win).
- **Output:** a `plan-aggregate.md` "best-practice review" banner — one keep/change row per
  question, with citations and (if change) the margin-passing pilot result.
- *Acceptance (Step 2):* banner exists with all four verdicts; the integrator + clustering config
  for Step 4 is locked.

### Step 3 — Wire the real chain as a DAG (method-agnostic skeleton; commit Step 2's method)
Replace the dead InSituType rules with the real ~13-node multi-env chain. Skeleton can be authored
in parallel with Step 2; only the integrator rule's method/params wait on Step 2.

| Rule | Script | Env | Key input → output | GPU |
|---|---|---|---|---|
| `neg_recover` | re-based neg step (**retires** parquet path) | spatial-rads | M01 **re-based RDS** negprobes assay + M02 raw RDS → `cell_neg.tsv` | — |
| `full_export` | `full_export.R` | spatial-rads | `merged.rds`,`cell_neg.tsv` → `full/mtx/*`,`full/obs.parquet` | — |
| `full_scvi` | **`full_scvi.py` (NEW — author it)** | spatial-rads-scvi | `full/mtx/*`,`full/obs.parquet` → `scvi_latent.parquet`,`scvi_model/` | ✔ |
| `full_cluster` | `full_cluster.py` | spatial-rads-scvi | latent + markers → labels/gates/recall, `cluster_checkpoint.h5ad` | — |
| `finalize_tier1` | `finalize_tier1.py` | spatial-rads-scvi | sweep → `full_coarse_labels.parquet`,`full_gates.json` | — |
| `tier2_immune_*` | subcluster.py → singler.R → rescue.py | scvi / R / scvi | coarse immune → `immune_subtypes_rescued.parquet` | — |
| `tier2_stroma_*` | subcluster.py → ucell.R → rescue.py | scvi / R / scvi | coarse stroma → `stroma_subtypes_rescued.parquet` | — |
| `unify_labels` | `unify_labels.py` | spatial-rads-scvi | all three tiers → `results/aggregate/full_labels.parquet` | — |

- **Delete** rules `embed_celltype`, `prepare_reference`, `typing_insitutype` and the
  `merged_typed.rds` output. **Re-point** `pseudobulk_build` / `pathway_summary` / `panel_coverage`
  to `merged.rds`.
- **Retire `recover_negprobes.R`'s parquet path:** `neg_recover` sources M01 from the re-based
  per-slide RDS negprobes assay (the deferred `plan-mutter01-controls.md` item); M02 already uses
  raw RDS. Its output (`cell_neg.tsv`) interface is unchanged, so `full_export.R` is untouched.
- **Handoff contract = MTX + parquet** (the existing hand-run interfaces). Export counts once
  (`Matrix::writeMM`); Python reads via anndata, writes latent back keyed by `cell_id`. **Avoid
  `zellkonverter`** (densifies, ~2× RAM).
- Snakemake hygiene from line 1: every script an `input:` dep (`feedback_snakemake_script_tracking`);
  `threads: 1` + BLAS pin; serialize heavy rules with `resources: mem_mb` (R ScaleData/merge is the
  OOM risk).
- *Acceptance (Step 3):* `snakemake -n` shows the full typing chain producing `full_labels.parquet`
  (no orphan inputs, no dead rules); the differential rules resolve against it.

### Step 4 — Single full 3.27M-cell rebuild (GPU-direct on the RTX A4000)
Run Step 3's DAG end-to-end with Step 2's committed method, frozen seed/env. Reuses the
`cluster_checkpoint.h5ad` pattern so re-tries cost minutes.

### Step 5 — Re-validation (acceptance gate — NOT byte-identity)
- Re-run the **four pre-registered Q4 gates + per-lineage marker recall** (`full_gates.json`).
- Sanity-compare new coarse composition against the locked baseline (tumor 48.1 / stroma 40.6 /
  immune 11.2; scaled iLISI 0.67). Expect **close, not identical**.
- Commit the reproducibility artifact: conda lock + scVI/torch/leidenalg/CUDA versions + seeds.
  `full_labels.parquet` is the committed artifact-of-record.
- *Acceptance (Step 5):* ALL_PASS on the Q4 gates; composition within a stated tolerance of
  baseline; frozen-env record committed.

### Step 6 — Refresh the differential layer (in scope BECAUSE A2 re-types)
The 13 downstream rules (`composition` … `assemble_results` → `results_master.tsv`) are **already
built and working** (2026-06-03), but they consumed the *hand-run* `full_labels.parquet`. The Step 4
rebuild regenerates that file, so the committed `results_master.tsv` is now stale.
- Re-run the differential layer on the rebuilt labels — **refresh, not redesign** (no science
  change; only the label input is refreshed and now reproducible).
- *Acceptance (Step 6):* `results_master.tsv` regenerated end-to-end from the wired DAG; spot-check
  the headline result (day-2 SBRT-driven stromal fibrosis; MBRT ~null) is unchanged in direction.

### Step 7 — Literate transfer (C1; depends on Step 3 rules existing + Steps 4–5 verifying)
- **Verify-first:** inventory the org's block types (`:tangle` target vs `jupyter-R` vs
  `:tangle no`) and confirm a plain `org-babel-tangle` over the mixed set is safe **before** moving
  any block (tangle hazard; the org file is 4284 lines).
- Author `aggregate.smk` + thin R **wrapper** rules as org-babel blocks under a new
  `** Aggregate analysis` heading; fold each block **only after its rule runs end-to-end and
  verifies** (the pre-req gate from `plan-aggregate.md`).
- **Document** the heavy scripts (GPU-Python tier-1, tier-2 R/Py) in an org Methods entry that
  *points at* `scripts/aggregate/` rather than tangling them. Keeps the guideline; literate index
  without doubling the file.

**2026-07-12 addendum — executed as a two-file split, not one `aggregate.smk`.** The literate
transfer described above (one `** Aggregate analysis` heading tangling one `aggregate.smk`) was
carried out instead as **two** org-tangled workflows, split at the label-handoff seam:
`aggregate_typing.smk` (15 rules, terminus `full_labels.parquet` + `merged.rds`) and
`aggregate_differential.smk` (21 rules, consumes those three leaf inputs, terminus
`results_master.tsv`). Rationale: the typing/differential boundary is the biological seam —
typing is a one-time structural identity assignment on the merged cohort, while the differential
layer (and the future peak/valley spatial analysis) only ever *consume* labels, never recompute
them — and the split isolates the GPU/multi-env typing DAG (scVI, Leiden, SingleR, tier-2 rescue
scripts) from the fast-iterating R differential work, so a differential-side edit no longer
forces a re-think of the typing DAG's environment/GPU dependencies. A rule-conservation diff
confirmed no rule was lost or duplicated across the split (36 named rules old == 15 + 21 new),
resolved output-path values were confirmed unchanged (`{AGG}`/`{FULL}` renamed to
`{D_AGG}`/`{D_FULL}` in both new files with identical resolved strings), and both new workflows
dry-run clean (only the expected `MissingInputException` for not-yet-materialized upstream
intermediates). The monolithic `workflows/aggregate.smk` was removed after these checks passed.

---

## Cross-cutting — reproducibility artifact

The Step 4 rebuild and any B-driven switch run through: (a) `full_labels.parquet` = committed
artifact-of-record; (b) frozen env (conda lock + recorded scVI/torch/leidenalg/CUDA versions);
(c) recorded seeds; (d) acceptance = **re-validation of the four Q4 gates + marker recall**, not
byte-identity. The rebuild **is gated on / includes** the deferred M01 raw-RDS re-base (Step 3's
`neg_recover` retires the parquet path onto it).

## Out of scope

- The differential **science design** (complete 2026-06-03). Step 6 *refreshes its outputs* on the
  rebuilt labels but does not change the analysis.
- M01 deferred items in `plan-mutter01-controls.md` **except** the `recover_negprobes` retirement,
  which Step 3 absorbs.
- The tongue cohort (n=1, separate biology).
- The BANKSY niche embedding (`embed_niche.R`) — the wired niches are the k-means/kNN-composition
  niches in the differential layer, not BANKSY; BANKSY stays a future item per `plan-aggregate.md`.

## Provenance

Restructured from a bundled investigate→update→transfer proposition by a blinded devils-advocate
panel (2 needs-rework / 2 solid-with-caveats; all four lenses confirmed the workflow-reality gap;
shared blind spot = wire-vs-quarantine), run `wf_c66fbfa4-c5f`, 2026-06-22. Three open decisions
resolved by the user 2026-06-22 (A2 / decision-forcing / C1) and harmonized into the seven-step
sequence above; ground-truthed against live `aggregate.smk` + `scripts/aggregate/` + `{AGG}/full/`
the same day, which added the `obs.parquet` second-orphan and the untracked-`scvi_latent.parquet`
findings.
