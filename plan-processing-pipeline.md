# Reproducible processing pipeline — design, decisions, status

Session record for the consolidation of spatial-rads into reproducible Snakemake
workflows. Captures architecture and the (heavily deliberated) decisions so they
aren't lost. Companion to `plan-mbrt-signatures.md`.

## v2.0 supersession — Atlas re-typing (2026-05-31)

Cell typing (stage 5) is being **replaced**, not extended. v1 typed Mutter_01
with Yi's ImmGen labels and transferred them M01 -> M02 by Seurat anchors; it
leaves **48.5% of the 3.27M flank cells unclassified** because ImmGen is
immune-only. v2.0 types **all cells (M01 + M02) directly against one fixed
external reference** (profile-native, supervised) — no M01 -> M02 transfer, so
batch-robust by construction. Decided reference (**option A**): ImmGen (immune;
panel-space CSV already on disk) + GeoMx MCA `MammaryGland` projected onto our
950 panel (non-immune) + marker-based malignant 4T1, with **marker-based lineage
scoring (option B) as the panel-native fallback**.

Load-bearing finding (verified 2026-05-31): **there is no panel-native mouse
mammary profile.** NanoString ships two easily-conflated repos — (1)
**CosMx-Cell-Profiles** (panel-native, InSituType-ready) whose mouse tier is
**Brain only**, on a *different* 1k mouse panel (**177/950 gene overlap**,
missing Epcam/Krt8/Cdh5), so unusable here; (2) **CellProfileLibrary** (GeoMx-WTA)
which has the MCA mammary profiles but in whole-transcriptome space and
spatially-naive. The mammary reference must therefore come from repo (2) and be
projected onto our panel.

**Downstream impact:** re-typing invalidates every cell-type-keyed output — the
processing typing stage AND all of `aggregate.smk`'s cell-type tracks
(composition, pseudobulk DE, GSEA, pathway, per-cell DE, concordance) must re-run
on the new labels. Companion plans bumped to v2.0: `plan-aggregate.md`,
`plan-outcomes-de.md`, `plan-mbrt-signatures.md`. The build is scoped in the
**"Atlas re-typing — work plan (v2.0 stem)"** section at the bottom of this file.

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
   sink under the cross-dataset platform shift (see Decisions). **v2.0: this stage
   is being replaced by external-atlas profile-native typing — see the Atlas stem.**
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
  ImmGen ontology (Yoshida et al. Cell 2019, PMID 30686579); NanoString's
  CellProfileLibrary ships it, but in GeoMx-WTA gene space rather than projected
  onto the CosMx 950 panel, so we reverse-engineer the panel-space profile from
  Yi's labels to keep M01/M02 labels consistent. (Note two easily-conflated
  NanoString repos: the on-disk **CosMx-Cell-Profiles** is mouse-**brain-only**;
  the separate **CellProfileLibrary** GeoMx library has mammary — that is the one
  v2.0 uses. See the v2.0 banner + Atlas stem.)
  `tumor_epithelial` is added because ImmGen has no tumor type (4T1
  cells got dumped into Yi's de-novo `a` bucket; only ~32% of `a` is actually
  epithelial). M02 tumor cells end up mostly in `a` (~55%, faithful to M01) +
  some `tumor_epithelial` (~2-3%); for analysis we group `a`+`tumor_epithelial`
  as the "tumor compartment".
- **Atlas expansion for the non-immune lineages (committed 2026-05-31, not yet
  built).** ImmGen is immune-only, so 48.5% of the 3.27M-cell flank cohort
  (1.59M cells) sits in the InSituType de novo `a` bucket and the M02 tumor
  compartment is almost entirely unassigned (1.28M of 1.59M); well-matched
  reference mapping leaves single-digit-to-15% unassigned, so this is reference
  mismatch, not noise. **Decision: type against the NanoString CellProfileLibrary
  mouse mammary atlas** (markers are only a break-glass contingency if projection
  fails, not a co-equal choice). Load-bearing caveat (verified 2026-05-31): the
  panel-native NanoString library already on disk (**CosMx-Cell-Profiles**) is
  **mouse-brain-only** on a different 1k panel (**177/950 gene overlap**, missing
  Epcam/Krt8/Cdh5), so it is unusable; the mammary profiles live in the *separate*
  GeoMx-WTA **CellProfileLibrary** (github.com/Nanostring-Biostats/CellProfileLibrary,
  `Mouse/Adult`; MCA = Han et al., "Mapping the Mouse Cell Atlas by Microwell-Seq,"
  Cell 172:1091 (2018), PMID 29474909), which is whole-transcriptome and must be
  projected onto our 950 panel (expected to retain most genes, unlike the brain
  library's 177). ImmGen handles the immune compartment (panel-space CSV already
  on disk at `analysis/objects/mbrt_vs_sbrt/ImmuneAtlas_ImmGen_derived_from_M01.csv`);
  malignant 4T1 has no turn-key atlas and stays marker-called. InSituType
  (supervised mode, consumes profile matrices) is installed; SingleR/celldex are
  not. Full build steps in the **v2.0 Atlas stem** at the bottom of this file.
  Plan: type the non-immune compartment with a **profile-native typer**
  (SingleR correlation, or InSituType's *supervised* mode — both consume reference
  profile matrices directly) against **ImmGen (immune) + MCA `MammaryGland`
  non-immune profiles**, plus **marker-based assignment for the malignant 4T1
  identity**. This is deliberately close to "use the NanoString library as-is":
  the library already ships ImmGen, so the combined reference is just more profile
  columns in one matrix, not a custom-engineered object. (Tabula Muris Senis,
  Nature 583:590 (2020), is the broader-but-generic alternative to MCA.) The
  malignant compartment has no off-the-shelf reference (MCA mammary profiles are
  normal luminal/basal/stromal, not carcinoma); a **4T1-specific scRNA reference
  is deferred to future work**, so near-term malignant cells stay marker-called
  (`tumor_epithelial`).

  **What "assembling" the reference actually requires** (two mandatory, one
  optional — there is far less custom "harmonization" than first sketched):
  1. *Gene-space subset (mandatory).* The library profiles are GeoMx-WTA (~20k
     genes); intersect to our 950 CosMx panel before scoring. Not optional — a
     format requirement for any use of these profiles.
  2. *Concatenate per-tissue RData (mandatory, mechanical).* The library ships
     separate per-tissue objects, so even "full NanoString" means `load()`-ing the
     chosen tissues + ImmGen and `cbind`-ing their profile columns into one labeled
     matrix.
  3. *Dedupe immune (optional polish).* ImmGen gives fine labels (`T.8.naive`); an
     MCA tissue gives coarse `T cell`/`macrophage`. Keeping both splits immune
     cells across two vocabularies. Sourcing immune only from ImmGen (drop the
     atlas's coarse immune columns) keeps one clean ontology — label hygiene, not
     correctness.

  **One scope decision, settled on the pilot: all ~20 MCA tissues vs the
  mammary-relevant subset.** Extra tissues (testis/placenta/etc.) rarely
  mis-attract under correlation (an irrelevant type only wins a cell if it
  genuinely correlates best, which itself flags a no-fit), but more candidate
  types = more spurious near-ties on marginal cells. Pilot both, compare the
  unassigned/confidence distributions.

  **Method change vs current pipeline.** M01 currently keeps Yi's ImmGen labels
  and M02 inherits via Seurat anchor `TransferData`. Profile-native typing scores
  **each cell independently against the fixed external profiles**, which sidesteps
  the cross-dataset batch shift that broke InSituType de-novo in the first place —
  so we can type **M01 and M02 both directly** against the same reference (no
  M01→M02 transfer needed), or keep the M01-then-inherit path; settle on the
  pilot. Panel adequacy for the target lineages is confirmed against
  `results/processing/common_genes.tsv` — epithelial 6/8, endothelial 7/9,
  fibroblast 8/8, smooth-muscle 4/5 markers present, **pericyte the thin spot at
  3/8** (Rgs5 / Pdgfrb / Notch3 present; Cspg4 / Mcam / Des / Kcnj8 / Abcc9
  absent), so pericyte-vs-smooth-muscle is the weakest call. The `aggregate.smk`
  Stage 1a `embed_celltype` de-novo-cluster + marker-concordance step is the
  reference-free cross-check, not the fix.
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
- **Next** (priority re-ordered 2026-05-31 — **atlas re-typing is now #1**: all
  cell-type-resolved aggregate outputs already built, i.e. composition, pseudobulk
  DE, GSEA, pathway, embed_celltype, are computed on the current ImmGen+`a` labels
  and must be re-run once the types change, so lock typing before stacking more
  analysis on labels that will move):
  1. **Atlas-expand the typing reference** (see Key decisions + the v2.0 Atlas
     stem) — type the non-immune compartment with a **profile-native typer**
     (supervised InSituType) against **ImmGen + MCA `MammaryGland` profiles**, markers for
     malignant 4T1 (4T1-specific reference = future work). Mandatory prep: subset
     library profiles to the 950 CosMx panel and `cbind` the per-tissue RData +
     ImmGen into one labeled matrix; optionally dedupe immune (ImmGen only). Pilot
     resolves two forks: (a) all ~20 MCA tissues vs mammary-relevant subset, and
     (b) type M01+M02 both direct vs M01-then-inherit. Pericyte resolution is the
     known weak spot (3/8 markers). **Benchmark target**: cut unassigned from
     48.5% toward the well-matched <~15% norm. Pilot on one tumor-heavy M02 sample
     before the full 23-sample re-run; report the stromal/vascular rescue
     separately from how many malignant cells anchor to MCA epithelial (the
     malignant majority is the uncertain part).
  2. **Fold the cell-typing/pathway/probe-QC layer into the literate
     `spatial-rads.org`** — currently lives only in `scripts/` + `workflows/`;
     the org's `processing.smk` block is the **stale 10-rule version** and
     would *regress* the workflow back to no-celltype/no-pathway if tangled
     (landmine, removed by the fold). Fold the atlas re-typing in at the same time.
  3. **Build `aggregate.smk`** — merge → Harmony → cluster → UMAP → cell-type
     landscape → cross-dataset DEGs. Also where **Option C** for cross-dataset
     typing lives if we ever want to escalate from per-sample anchor transfer.
  4. **Review `probe_qc_report.tsv`** flagged candidates (smooth, not drop).

## Key files

`workflows/{data_model,processing}.smk` · `scripts/*.R` ·
`config/{config,common,pathway_gene_lists}.yaml` ·
`data/{metadata.xlsx,metadata_schema.yaml,data_model.rda}` · `results/processing/*`

## Atlas re-typing — work plan (v2.0 stem)

**Status: stem.** The direction is decided; expand into a committed sub-plan on
resume. Supersedes v1 stage 5 (Yi-ImmGen + M01 -> M02 anchor transfer).

### Job to be done
Assign every flank cell (M01 + M02, 3.27M) to a lineage by scoring it directly
against one fixed external reference — cutting the 48.5% ImmGen-only unassigned
rate toward the well-matched <~15% norm — so all downstream cell-type-resolved
analyses run on honest labels.

### Decided approach (one path)
- **Type against the NanoString CellProfileLibrary mouse mammary atlas.** Profile-
  native, supervised typing applied to **all cells directly** (no M01 -> M02
  transfer); batch-robust because the labels are defined by a fixed external
  reference, not by either dataset.
- **Reference = ImmGen (immune) + GeoMx MCA `MammaryGland` projected to our 950
  panel (non-immune) + marker-called malignant 4T1.**
- **Markers are a break-glass contingency only** — used for the non-immune
  compartment *only if* projecting the GeoMx atlas onto the 950 panel leaves too
  few genes to separate lineages (not expected: whole-transcriptome should retain
  most of our 950, unlike the brain library's 177). Not a co-equal option.
- Typer: supervised InSituType (installed, consumes profile matrices). SingleR
  optional (needs a Bioconductor install).

### Inputs / what's on disk
- ImmGen panel-space profile (immune, ready): `analysis/objects/mbrt_vs_sbrt/ImmuneAtlas_ImmGen_derived_from_M01.csv` (1000 genes × Yi types incl `a`,`b`).
- CosMx-Cell-Profiles repo — **NOT usable** (mouse = brain, 177/950 overlap): `/mnt/data/projects/spatial-rads/CosMx-Cell-Profiles`.
- GeoMx MCA `MammaryGland` — **NOT yet downloaded**: github.com/Nanostring-Biostats/CellProfileLibrary `Mouse/Adult`.
- Our 950 panel: `results/processing/common_genes.tsv`.
- 23 per-sample objects: `/mnt/data/projects/spatial-rads/processing/{norm,typed,scored}/`.

### Step skeleton
1. **Download the good atlas (FIRST ACTION).** The mammary profile is NOT on disk
   — clone the GeoMx CellProfileLibrary and pull `Mouse/Adult/MammaryGland`:
   ```
   git clone https://github.com/Nanostring-Biostats/CellProfileLibrary.git \
     /mnt/data/projects/spatial-rads/refs/CellProfileLibrary
   # mouse profiles under CellProfileLibrary/Mouse/Adult/
   #   MammaryGland_*.RData  (optionally all ~20 Mouse/Adult tissues for the pilot)
   ```
2. Project the profile onto the 950 panel; QC gene coverage + lineage separability. **Gate:** if it can't separate epithelial/endo/fibro/pericyte, drop to the marker contingency.
3. Assemble one labeled profile matrix: ImmGen immune columns + MCA non-immune columns (malignant handled by markers, separately).
4. **Pilot** on one tumor-heavy M02 sample: supervised type, measure unassigned drop vs 48.5%, report stromal/vascular rescue separately from malignant epithelial.
5. If the pilot passes, revise `build_celltype_reference.R` / `celltype.R` to the external-atlas path; re-run processing typing for all 23 samples.
6. Re-run `aggregate.smk` cell-type tracks on the new labels.

### Open decisions to settle on the pilot
1. All ~20 MCA Mouse/Adult tissues vs MammaryGland-only (extra tissues broaden coverage but add spurious near-ties).
2. Type M01+M02 both-direct vs keep M01-then-inherit (both-direct is cleaner/batch-robust).
3. Malignant 4T1 marker set + threshold (generalize the existing `tumor_epithelial` relabel). 4T1-specific scRNA reference = future work.

### Panel adequacy (already verified, `common_genes.tsv`)
Epithelial 6/8, endothelial 7/9, fibroblast 8/8, smooth-muscle 4/5, **pericyte
3/8** (Rgs5/Pdgfrb/Notch3 present; Cspg4/Mcam/Des/Kcnj8/Abcc9 absent — weakest
call). 48% unassigned is a *reference* gap, not a *panel* gap.

### Success criteria
Unassigned < ~15% cohort-wide; non-immune lineages populated at plausible
fractions; tumor compartment coherent; M01/M02 label distributions concordant
where biology matches.

### Downstream ordering
Re-typing moves labels, so it must precede `aggregate.smk` cell-type tracks,
`plan-outcomes-de.md` findings, and `plan-mbrt-signatures.md` per-cell-type
signatures. Lock typing before stacking analysis on labels that will move.
