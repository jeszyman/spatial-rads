# Pathway gene-sets + panel generation — design/spec

Single-source-of-truth refactor for the common gene panel and the tiered pathway gene
sets, plus collapse of the redundant double pathway scoring. Companion to
`plan-processing-pipeline.md` and `plan-aggregate-refactor.md`.

## Context / problem

Four coupled defects in how gene sets and the panel are handled today:

1. **Broken panel script.** `scripts/common_gene_panel.R` reads Mutter_01 genes from the
   counts **parquet** (`counts_path`), but the 2026-06-22 rds re-base set M01
   `format=rds` and `counts_path=NA` for both datasets. It also filters on a `dataset`
   column the sample sheet doesn't have (it is `name`). The rule would now fail.
2. **Duplicated gene-set definitions.** The 2-tier structure (curated confirmatory +
   ~50 MSigDB Hallmark exploratory) already exists, but is **re-derived independently**
   in ~8 scripts (`gsea.R`, `pathway_summary.R`, `build_gene_sets.R`, `panel_coverage.R`,
   `assemble_results.R`, plus yaml-reading `contamination_qc.R`, `qc_plots.R`,
   `tier2_stroma_ucell.R`). Each `gsea.R`/`pathway_summary.R` calls `msigdbr` **live and
   unpinned** — a reproducibility hole (a version bump silently shifts the set list).
3. **Redundant double scoring.** Per-cell program scores are computed twice: per-sample
   (`processing.smk/pathway_score.R`, curated 8) and again merged
   (`aggregate.smk/pathway_summary.R`, curated + Hallmark, "recomputes"). UCell is
   rank-based **within a cell** (cohort-independent), so per-sample UCell == merged UCell
   exactly; the per-sample values are recomputed and discarded. No inference consumer
   reads the per-sample `_UCell`/`_AMS` columns (`assemble_results.R`, `power_mde.R` read
   only `pathway_summary.R`'s merged output).
4. **No tier-1 provenance / freeze.** The curated sets lack recorded provenance and have
   no guard that their membership stays frozen.

## Goal

One generated, reproducible, panel-aware, **tier-structured** gene-set artifact plus the
common-gene panel, both produced in `data_model.smk` (Stage A); every consumer reads
them instead of re-deriving; pathway scoring happens **once**, at merged scale.

## Design

### Artifacts (produced by `data_model.smk`)

- **`results/data_model/common_genes.tsv`** — the 950-gene panel (Mutter_01 ∩ Mutter_02),
  built by reading RNA-assay rownames from one per-slide **RDS** of each dataset (both are
  `format=rds`), via `raw_input_path` (uniform M01/M02 read; drops the parquet/`dataset`
  paths). The file is git-tracked. **Invariance guard:** the freshly-built list must
  `setequal` the currently-committed version of this same file before it is overwritten,
  else hard-error (identity, not cardinality). The committed file is the frozen reference
  that protects the locked aggregate typed on this panel. (The existing
  `results/processing/common_genes.tsv` is moved here, carrying its git history as the
  initial snapshot.)
- **`results/data_model/pathway_sets.tsv`** — long/tidy: columns `set, tier, source, gene`.
  - `tier ∈ {primary, exploratory}` (matches the labels the scoring scripts already use;
    `primary` = the curated confirmatory sets).
  - tier=primary from `config/pathway_gene_lists.yaml` (the human-maintained source, kept
    git-tracked and authoritative), each with a `source:` provenance string.
  - tier=exploratory = MSigDB mouse Hallmark, pulled **once** via the existing pattern
    `msigdbr(species = "Mus musculus", collection = "H")`. The **msigdbr package version**
    is recorded (header comment + a `_provenance` sidecar row/file).
  - Hallmark set names are namespaced (e.g. `HALLMARK_*`) so they can never collide with a
    primary set name.
- **`results/data_model/gene_set_panel_coverage.tsv`** — per set: `set, tier, source,
  n_total, n_panel, usable`. **Every** set is listed, including excluded ones (nothing
  silently dropped). `usable = n_panel >= min_panel_genes` (config `pathway.min_panel_genes`,
  default 5, documented). Primary sets are always `usable=TRUE` even if thin, flagged
  `thin = n_panel < min_panel_genes`.

### Guards (the three the devils-advocate panel required)

1. **Panel invariance** — `setequal` against the committed snapshot, hard-error (above).
2. **Tier-1 freeze** — generated tier=primary membership must `setequal` the committed
   `config/pathway_gene_lists.yaml`; mismatch hard-errors (pre-registration integrity).
3. **Completeness** — post-write assertion that every primary set appears in
   `pathway_sets.tsv`; missing → error (honors "nothing silently dropped").

### Rewire consumers (read the artifact; none call `msigdbr`)

| Script | Change |
|---|---|
| `gsea.R` | Replace inline `read_yaml + msigdbr` (L28-36) with a read of `pathway_sets.tsv`; filter by `tier`. Net **removal** of duplicated logic. |
| `pathway_summary.R` | Same: read `pathway_sets.tsv` instead of yaml+live msigdbr (L52-58). Becomes the **single** scoring path. |
| `assemble_results.R` | Read tier=primary gene membership for the H1/H2/H3 families from the artifact (was `read_yaml`). Confirmatory family hypothesis lookup is unchanged. |
| `build_gene_sets.R` | Promoted from `scripts/aggregate/` to `scripts/build_gene_sets.R` and extended from audit-only to the **generator**: emits `pathway_sets.tsv` + `gene_set_panel_coverage.tsv` from the curated yaml + pinned msigdbr, against the panel. The `aggregate.smk` audit-only `build_gene_sets` rule is retired. |
| `panel_coverage.R` | Read `pathway_sets.tsv` (was yaml). |
| `contamination_qc.R`, `qc_plots.R`, `tier2_stroma_ucell.R` | Confirm whether they need pathway sets or `lineage_markers.yaml`; repoint pathway-set reads only. (lineage markers are a separate file, untouched.) |

### Collapse the double scoring

- **Delete** `processing.smk/pathway_score.R` and the per-sample `scored` stage. Per-sample
  pipeline now terminates at `typed.rds` (the `scored` stage was `typed` + redundant
  pathway columns).
- `pathway_summary.R` (merged scale, sourcing the artifact) is the single scoring path.
- **Aggregate merge input** `*.scored.rds` → `*.typed.rds`: update `aggregate.smk`
  (`SCORED`/`PILOT`/`FLANK` expands) and the `.scored.rds`-stripping scripts (`merge.R`,
  `merge_pilot.R`, `coords_necrosis.R`, `recover_negprobes.R`, `cell_assignment_map.R`).

### Path-move rewiring (panel)

`results/processing/common_genes.tsv` → `results/data_model/common_genes.tsv` in:
`processing.smk` (`GENES`), `aggregate.smk` (`PANEL`), `adapt_mutter01.R`, `adapt_mutter02.R`,
`probe_qc.R`, `build_celltype_reference.R`, `prepare_reference.R`, and the
`lineage_markers.yaml` comment. The committed snapshot stays at the old path as the frozen
invariance reference (or is moved + re-committed; decide at implementation).

### `data_model.smk` rule graph

`make_data_model` → `common_gene_panel` → `build_gene_sets`; `rule all` requires
`common_genes.tsv` + `pathway_sets.tsv` + `gene_set_panel_coverage.tsv`. Note: Stage A now
loads two per-slide RDS (one per dataset) for rownames — heavier than its current
metadata-only profile, but bounded (two files, not the full cohort).

## Sequencing / locked-state safety

The scored→typed rename changes an aggregate **merge input**, so it cannot be applied to the
locked aggregate piecemeal. It is sequenced **with the pending aggregate rebuild**
(`merged.rds` is already missing; dead InSituType chain — see `plan-aggregate-refactor.md`).
Two phases:

- **Phase 1 (no aggregate re-run):** add the two artifacts + guards in `data_model.smk`;
  rewire the gene-set readers and the panel path; pin msigdbr. Run targeted
  (`snakemake results/data_model/...`) + `--touch` any otherwise-flagged locked outputs so
  nothing recomputes gratuitously.
- **Phase 2 (with the aggregate rebuild):** delete per-sample `pathway_score.R` / `scored`
  stage and switch the merge to `typed.rds`. Folded into `plan-aggregate-refactor.md`'s
  rebuild so the merge re-runs intentionally, once.

## Out of scope

`config/lineage_markers.yaml` (typing markers, separate concern); the aggregate typing
method itself; adding new pathways beyond the existing curated + Hallmark tiers.

## Verification

1. `data_model.smk` dry-run, then targeted run: `pathway_sets.tsv` has both tiers,
   `_provenance` records the msigdbr version; `gene_set_panel_coverage.tsv` lists all sets
   incl. excluded with `usable`/`thin`.
2. Panel invariance: regenerated `common_genes.tsv` `setequal` committed snapshot (assert in
   the rule; 950 genes).
3. Tier-1 freeze: corrupt one curated gene in a scratch copy → build hard-errors.
4. `grep -rn "msigdbr" scripts/` → only the generator calls it; `grep -rn
   "pathway_gene_lists.yaml" scripts/` → only the generator (+ the committed source).
5. `grep -rn "common_genes.tsv\|\.scored\.rds" workflows/ scripts/` → no stale refs after
   each phase.
6. Phase-1 reproducibility: `gsea.R` / `pathway_summary.R` outputs unchanged vs current when
   run on the same inputs (same gene sets, just centrally sourced) — diff the result tables.
