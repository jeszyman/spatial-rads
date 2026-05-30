# Reproducible processing pipeline — design, decisions, status

Session record for the consolidation of spatial-rads into reproducible Snakemake
workflows. Captures architecture and the (heavily deliberated) decisions so they
aren't lost. Companion to `plan-mbrt-signatures.md`.

## Goal

Consolidate the ad hoc analyses into reproducible workflows so **Mutter_01 and
Mutter_02 yield comparable SBRT-vs-MBRT findings** — a true cross-dataset
reproducibility check, not a comparison of two different pipelines. Both datasets
are rebuilt **from raw through identical code**; Yi Liu's pre-computed Mutter_01
vendor layers are held only as a **concordance check**, never as the substrate.

## Workflows

- **`data_model.smk`** — relational `data/metadata.xlsx` (sheets: mice, datasets,
  slides, samples) + Frictionless `data/metadata_schema.yaml`, validated by
  `make_data_model.R` → `data/data_model.rda` + master `results/data_model/samples.tsv`
  (23 samples). *Committed.*
- **`processing.smk`** — the **per-sample / pre-aggregate** stage (below).
- **`aggregate.smk`** — **NOT YET BUILT.** Merge → batch/integration (Harmony) →
  clustering → UMAP → cell-type landscape/composition → cross-dataset DEGs. The
  clustering/UMAP/cell-type-validation/batch blocks deleted from the old Methods
  section belong here (preserved in git history as reference).

## `processing.smk` stages (per-sample)

The dividing line: *anything you can do to one sample by itself* lives here;
anything that needs all cells merged is `aggregate.smk`.

1. **adapter** (`adapt_mutter01`/`adapt_mutter02`) — raw → common-format sparse
   Seurat on the common **~950-gene** panel (M01 1000 ∩ M02 972). M01 from the
   counts parquet (transpose + join on `cell_id`, vendor cols sequestered to
   `yi_reference.tsv`); M02 from per-slide RDS with **seeded** y-band k-means
   sample assignment.
2. **probe-QC diagnostic** (`probe_qc`) — **report-only**, flags candidate failed
   probes (`probe_qc_report.tsv`); **does not drop genes** (see Decisions).
3. **QC filter** (`qc_filter`) — four-criteria gate (`nCount_RNA>20`,
   `nFeature_RNA>10`, `propNegative<0.5`, cell-area MAD outlier), identical on both.
   `propNegative` recomputed for M02 from its negprobes; `qcFlagsCell` is M01-only
   and used only for the Yi concordance.
4. **normalize** (`normalize`) — `LogNormalize`(1e4) + variable features.
5. **cell typing** (`build_celltype_reference` + `celltype`) — reference =
   pooled M01 cells (uniform random subsample to 20k preserving natural
   abundance) labeled with `cell_type` = Yi ImmGen label, confident-epithelial
   `a`/`b`/NA → `tumor_epithelial`. M01 keeps Yi's labels (+ relabel); M02 =
   **Seurat anchor label transfer** (`FindTransferAnchors` reduction=CCA +
   `TransferData` weight.reduction=CCA) against the reference. Replaces
   InSituType, which collapsed M02's tumor compartment into a Lymphatic.endo
   sink under the cross-dataset platform shift (see Decisions).
6. **pathway scoring** (`pathway_score`) — `UCell` primary + `AddModuleScore`
   (seed=42) secondary, from `config/pathway_gene_lists.yaml` (IFN-I/II, DDR, STING).

Reports: `qc_report`, `qc_plots`, `yi_concordance`, `celltype_report`.

## Key decisions and rationale

- **Rebuild M01 from raw** (vs reuse Yi's layers): otherwise the cross-dataset
  comparison is two pipelines. Verified M01 counts parquet is raw integers; M02 is
  raw RDS. Yi's layers → `yi_reference.tsv`, concordance only.
- **QC harmonization**: M02 was delivered raw and **lacks `qcFlagsCell`/`propNegative`**
  (never requested — confirmed in Gmail). So: recompute `propNegative` for M02 from
  negprobes; drop `qcFlagsCell` as a *filter* (keep as M01 check); identical
  four-criteria gate on both. (Earlier run: retention M01 78.8% / M02 83.2%, the
  latter matching the documented ~2.35M.)
- **Typing method: Seurat anchor TransferData, not InSituType (Option B).**
  The first build used InSituType against an M01-derived 44-profile mean-count
  reference. On M02 it confidently mis-typed: Lymphatic.endothelial **38%**
  sink, tumor compartment **0.15%**, prob median 0.99. **Decisive control**:
  same code + reference typed an M01 sample correctly (tumor compartment 29%
  recovered, lymph 4%), proving the failure was cross-dataset batch shift —
  InSituType's count-likelihood matching can't absorb the per-gene
  platform/efficiency differences between two CosMx runs. Fixed by switching
  to Seurat `FindTransferAnchors` (reduction=CCA) + `TransferData`
  (weight.reduction=CCA): rank/correlation-based, batch-robust, the standard
  tool for cross-dataset label transfer. Validated by marker concordance —
  Krt8+/Epcam+ M02 cells → tumor compartment **64–71% across all 12 M02
  samples** (vs 0.1% under InSituType), Lymphatic.endo **<1%** (vs ~37%).
  The InSituType `NEG_N=50`-vs-actual-10 background defect surfaced during
  this investigation becomes moot — InSituType is no longer used. **Option C
  (integrate M01+M02 with Harmony then transfer)** noted as the future ceiling
  and deferred to `aggregate.smk`.
- **Reference = pooled M01 cells, natural abundance, not profiles.** Built by
  `build_celltype_reference.R`: pool the 11 M01 norm objects' common-panel
  counts, attach `cell_type` = Yi's ImmGen label with confident-epithelial
  `a`/`b`/NA → `tumor_epithelial` (matches `celltype.R`'s M01 branch logic so
  cross-dataset labels are consistent). **Uniform random subsample to 20,000
  cells preserving natural abundance** — not per-label balancing, which
  flattens priors and over-represents rare ImmGen subsets (thymocytes, spleen
  subtypes), causing epithelial scatter (smoke test: only 8% of Krt8+ cells
  reached `tumor_epithelial` under per-label balancing). 20k also bounds
  per-job memory — at 40k cells the Seurat CCA OOM-killed even `--cores 1` on
  sam0023 (375k cells). Background: Yi's `ImmuneAtlas_ImmGen` is the public
  ImmGen ontology, but the panel-projected matrix isn't in NanoString's
  CellProfileLibrary (mouse = brain only); we reverse-engineer it from Yi's
  labels. `tumor_epithelial` is added because ImmGen has no tumor type (4T1
  cells got dumped into Yi's de-novo `a` bucket; only ~32% of `a` is actually
  epithelial). M02 tumor cells end up mostly in `a` (~55%, faithful to M01) +
  some `tumor_epithelial` (~2-3%); for analysis we group `a`+`tumor_epithelial`
  as the "tumor compartment".
- **Pathway scoring**: `UCell` primary (rank-based, no random background, dataset-
  independent) over `AddModuleScore` (thin ~950-gene background pool is biologically
  contaminated). Resolved the conflicting IFN-I gene lists to the documented
  NanoString-module set, committed to `config/pathway_gene_lists.yaml`.
- **probe-QC is a diagnostic, not a filter.** NanoString (Danaher, CosMx Analysis
  Scratch Space) advises **against dropping low-expression genes** — it discards
  rare-cell-type markers — and InSituType models background internally. An earlier
  drop-based version wrongly removed lineage markers (Pecam1, Ms4a1, Ncr1, Cd4,
  Cd163). Failed probes (Ozirmak Lermi 2025, 8–31.9% of a human panel at neg-control
  levels) are a narrow flag-and-review case. Documented in the org "Background and
  probe QC" section with citations.
- **per-sample vs aggregate split**: clustering, UMAP, batch assessment, landscape
  all require the merged object → `aggregate.smk`, deliberately deferred.

## Upstream + repo changes

- **science.org** `Data layout` standard: added the **workflow-linked sample sheet**
  convention; **harmonized `data_model.rda` → `data/` root**; added a consolidation
  TODO.
- **Literate consolidation**: workflows + scripts tangle byte-for-byte from
  `spatial-rads.org`; old superseded QC code blocks deleted (in git history).
- Five commits on `main` (data model, processing, org consolidation, gitignore, QC
  cleanup); the cell-typing/pathway/probe-QC scripts + config + org doc edits are
  currently **uncommitted** (pending verified run).

## Status / next

- **`processing.smk` full re-run complete 2026-05-30** (49 jobs, `--cores 1`
  + BLAS threads pinned). All 23 per-sample scored.rds produced. Marker-
  concordance validation across the previously-failing big M02 samples
  (sam0012/13/14/16/22/23): tumor compartment 64–71%, Lymphatic.endo
  0.4–0.6%, prob median 0.40–0.42. UCell + AMS pathway columns present.
- **Memory/thread lesson**: Seurat-transfer per-job memory on M02 queries
  (~300k+ cells) is heavy enough that `--cores 2` OOM-killed an R process
  at the original 40k reference; even `--cores 1` OOM-killed sam0023 at 40k,
  fixed by halving the reference to 20k cells (natural-abundance subsample).
  Two safe-defaults to bake into the next org fold:
  `shell.prefix("export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1; ")`
  and `threads: 1` per R rule, so future runs can use `--cores 4` without
  thinking about env vars. (Captured in memory `feedback_smk_thread_hygiene`.)
- **Next**:
  1. **Fold the cell-typing/pathway/probe-QC layer into the literate
     `spatial-rads.org`** — currently lives only in `scripts/` + `workflows/`;
     the org's `processing.smk` block is the **stale 10-rule version** and
     would *regress* the workflow back to no-celltype/no-pathway if tangled
     (landmine, removed by the fold).
  2. **Build `aggregate.smk`** — merge → Harmony → cluster → UMAP → cell-type
     landscape → cross-dataset DEGs. Also where **Option C** for cross-dataset
     typing lives if we ever want to escalate from per-sample anchor transfer.
  3. **Review `probe_qc_report.tsv`** flagged candidates (smooth, not drop).

## Key files

`workflows/{data_model,processing}.smk` · `scripts/*.R` ·
`config/{config,common,pathway_gene_lists}.yaml` ·
`data/{metadata.xlsx,metadata_schema.yaml,data_model.rda}` · `results/processing/*`
