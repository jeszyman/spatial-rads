# Aggregate workflow — reproducibility-wiring, best-practice review, literate transfer (PLAN ONLY)

> **PLAN ONLY — for later execution.** A devils-advocate panel (2 needs-rework / 2 solid-with-caveats, 2026-06-22) restructured this from a single bundled *investigate → update → transfer* plan into **three decoupled, independently-executable workstreams**. No workstream gates another except C depends on A. Each can be approved and run on its own. Supersedes nothing committed; companion to `plan-aggregate.md` (the as-built record).

## Context — the workflow-reality gap

`workflows/aggregate.smk` does **not** reproduce the locked results it appears to. Grounded 2026-06-22:

- The **wired** typing rules are the **superseded, failed** per-cell path: `embed_celltype` (Harmony) → `prepare_reference` → `recover_negprobes` → `typing_insitutype` (InSituType), which validated at **85% mistyped tumor** and was abandoned.
- The labels the project **actually uses**, `results/aggregate/full_labels.parquet`, are consumed as a **static input by 9 rules but produced by NO rule**. They came from ~10 **hand-run** scripts absent from the DAG: `full_export.R` → `pilot_scvi.py` (scVI, GPU) → `full_cluster.py` → `finalize_tier1.py` → `tier2_immune_subcluster.py` → `tier2_singler.R` → `tier2_immune_rescue.py` → `tier2_stroma_{subcluster.py,ucell.R,rescue.py}` → `unify_labels.py`.
- Consequence: a clean `snakemake -s workflows/aggregate.smk` **silently re-types with the wrong (failed) chain**, no error. The locked typing is not reproducible from the workflow.

Live/dead classification (grounded — do NOT assume):
- `merged_typed.rds` carries **dead InSituType labels** but is read by `pseudobulk_build`, `pathway_summary`, `panel_coverage` **for counts only** (they label from `full_labels.parquet`). Its counts equal `merged.rds`'s.
- `recover_negprobes` → `cell_neg.tsv` feeds the **live** hand-run `full_export.R` (per-cell negprobe covariate for scVI), so it is **not** purely dead.
- So "retire one script" understates the job (panel finding 1): it is a careful re-point + quarantine, or a full multi-environment (R → GPU-Python → R) DAG build.

This session also re-based M01 onto per-slide RDS (`scored.rds` now invariant; aggregate "armed" to rebuild) and added per-sample contamination QC to `processing.smk` — both inform Workstream B.

## Panel-forced corrections (baked into this plan)

1. **Decoupled** the three workstreams; A is method-agnostic and runs first.
2. **"Locked" ≠ bit-identical.** scVI GPU training + Leiden `igraph` are not byte-reproducible across library/driver versions. Reproducibility here = the committed `full_labels.parquet` as artifact-of-record + a frozen-env/seed record + **re-validation of the four pre-registered Q4 gates + marker recall** on any re-run. Not byte-identity.
3. **Literate transfer of the heavy chain conflicts with the project's own guideline** ("mature/complex code may not benefit from a literate approach"). Descoped + flagged as a decision, not assumed.

---

## Workstream A — make `aggregate.smk` honest (reproducibility). INDEPENDENT · do first · method-agnostic

Two targets; **A1 recommended now**, A2 deferred to a deliberate rebuild.

### A1 (recommended) — quarantine + lock ("honest by subtraction")
- **Task A1.0 — audit.** For each rule, confirm whether it uses `merged_typed.rds` for counts (re-pointable to `merged.rds`) or for InSituType labels (none expected — verify against each script's `commandArgs`). Output: a 1-page live/dead table. *Acceptance:* every `merged_typed.rds` consumer shown to ignore its `cell_type` column.
- **Task A1.1 — re-point counts.** Change `pseudobulk_build`, `pathway_summary`, `panel_coverage` to read `merged.rds` (counts identical, typing dropped). *Acceptance:* scripts read counts + join `full_labels.parquet`; no reference to `merged_typed.rds`'s typing.
- **Task A1.2 — delete the dead InSituType chain.** Remove rules `embed_celltype`, `prepare_reference`, `typing_insitutype` and the `merged_typed.rds` output. Keep `recover_negprobes`/`cell_neg.tsv` (feeds the locked typing); mark it locked (below).
- **Task A1.3 — lock the typing artifacts.** Declare `merged.rds`, `cell_neg.tsv`, and `full_labels.parquet` as `protected()`/`ancient()` (or a sentinel) so a clean run will not regenerate them, and commit a **provenance manifest** (`typing_provenance.md`): the exact hand-run script order + commit SHA + seeds + conda env versions that produced `full_labels.parquet`.
- *Acceptance (A1):* `snakemake -n` shows **no typing rule** wanting to run; differential rules resolve against `merged.rds` + the locked `full_labels.parquet`; dead rules gone; manifest committed. Locked results untouched.

### A2 (deferred — at the deliberate seed-frozen rebuild) — full from-scratch wiring
- Wire the real chain as rules with per-rule conda envs + a GPU resource: `full_export.R` (env `spatial-rads`) → scVI (env `spatial-rads-scvi`, `resources: gpu=1`) → `full_cluster.py` → `finalize_tier1.py` → `tier2_*` → `unify_labels.py`. Rule I/O contract = MTX/parquet (the existing hand-run interfaces). Retire `recover_negprobes.R` here (source `neg` from the re-based M01 RDS assay — the deferred item from `plan-mutter01-controls.md`).
- *Gate:* run only as a deliberate rebuild that **re-validates** the four Q4 gates + marker recall (results will not be byte-identical — accept by re-validation, not diff). Triggered by either the M01 re-base rebuild or a Workstream-B method switch.

## Workstream B — best-practice investigation. INDEPENDENT · parallel · decision-forcing BY RULE

Bounded research pass (same shape as the processing-QC research-swarm this session). **Questions (fixed, do not expand):** (1) scVI vs scANVI vs reference-mapping/scArches for cross-dataset CosMx typing; (2) Leiden resolution selection at >3M cells; (3) recovering the 33% unassigned stroma (precedent-first — see memory `project_stroma_unassigned_recovery`); (4) whether a merged-scale contamination/QC step adds value beyond the per-sample MECR now in `processing.smk`.

- **Evidence rule:** every method finding cited (Scopus/lit per memory `feedback_verify_method_precedent`); a **NULL result ("no peer-reviewed support at CosMx scale") is an acceptable, unblocking outcome.**
- **Output form:** a `plan-aggregate.md` "best-practice review" banner — one decision row per question (keep / change + the evidence).
- **Pre-committed decision rule (resolves the locked-results tension):** adopt an alternative method **only if** it beats the current scVI cluster-then-annotate on the **existing 2-slide bake-off harness** (`pilot_*`) by a pre-set margin on a named metric set (imbalance-scaled iLISI batch-mixing **and** marker recall **and** unassigned fraction). If it wins → triggers an **A2** seed-frozen rebuild with re-validation. If not → documented, **no change to locked results.**
- *Acceptance (B):* the banner exists with a keep/change verdict + citations for all four questions; if any verdict is "change," a margin-passing bake-off result is attached.

## Workstream C — literate transfer. INDEPENDENT · SEPARATE · depends on A · **DECISION REQUIRED**

**Active policy conflict:** the project guideline exempts mature/complex code from the literate approach, and `plan-aggregate.md` chose standalone authoring *deliberately*. The aggregate chain is 12 GPU-Python + complex tier-2 R — exactly that class. The org notebook is **4284 lines**; a full transfer ~doubles it and mixes tangle-target R with non-tangleable `jupyter-R`/`:tangle no` blocks (tangle hazard).

- **C1 (recommended) — transfer the orchestration only.** Author `aggregate.smk` + thin R *wrapper* rules as org-babel blocks; **document** the heavy scripts (GPU-Python tier-1, tier-2 R/Py) in an org Methods entry that *points at* `scripts/aggregate/` rather than tangling them. Keeps the guideline intact; gives a literate index without the doubling.
- **C2 — full transfer of all ~40 scripts** (the original ask). Only viable after splitting the org file into sub-notebooks; reverses the guideline; tangle-feasibility must be proven first on a mixed-block subset.
- **Pre-req gate (from `plan-aggregate.md`):** fold a script to org **only after its rule runs end-to-end and verifies** — so C depends on **A** (rules existing), not on B.
- **Verify-first task:** inventory the org's block types (`:tangle` target vs `jupyter-R` vs `:tangle no`) and confirm a plain `org-babel-tangle` over the mixed set is safe before moving any block.

## Cross-cutting — reproducibility artifact (resolves the locked-results tension)

Any A2 rebuild or B-driven switch runs through: (a) `full_labels.parquet` is the committed artifact-of-record; (b) a frozen env (conda lock + recorded scVI/torch/leidenalg/CUDA versions); (c) recorded seeds; (d) **acceptance = re-validation of the four Q4 gates + marker recall**, NOT byte-identity. State explicitly whether the rebuild is gated on (or excluded from) the deferred M01 raw-RDS re-base.

## Open decisions (for the user, before execution)

1. **Workflow-honesty target:** A1 quarantine+lock now (recommended) vs A2 full-wire now.
2. **Literate transfer:** C1 orchestration-only + document (recommended) vs C2 full transfer (needs org split + tangle proof).
3. **Investigation stance:** decision-forcing (pre-set switch margin, as written) vs exploratory-only (informs, never auto-triggers a rebuild).

## Out of scope

- The differential layer (complete 2026-06-03; `results_master.tsv`).
- M01 deferred items in `plan-mutter01-controls.md` except where A2 absorbs them (retire `recover_negprobes.R`).
- The tongue cohort (n=1, separate biology).

## Provenance

Restructured from a bundled investigate→update→transfer proposition after a blinded devils-advocate panel (2 needs-rework / 2 solid-with-caveats; all four lenses confirmed the workflow-reality gap; shared blind spot = wire-vs-quarantine, resolved here toward A1). Panel run `wf_c66fbfa4-c5f`, 2026-06-22.
