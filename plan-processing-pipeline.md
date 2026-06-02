# Reproducible processing pipeline — design, decisions, status

Session record for the consolidation of spatial-rads into reproducible Snakemake
workflows. Captures architecture and the (heavily deliberated) decisions so they
aren't lost. Companion to `plan-mbrt-signatures.md`.

## v2.0 — De-anchored re-typing (2026-05-31)

**The defect is M01-anchoring, not just the unassigned rate.** v1 types Mutter_01
with Yi's ImmGen labels and makes Mutter_02 *inherit* them via Seurat CCA transfer
— so M02's labels are **derived from M01**. For a cross-dataset SBRT-vs-MBRT
comparison that is circular: regression-to-M01 artifacts (proportions, fine immune
subtypes) can masquerade as biological concordance. The 48.5%-unassigned blob is
the *symptom* in the non-immune compartment; the disease is that the labels are
anchored to one of the two datasets being compared. Measured 2026-05-31: pooled
M02 transfer `prediction.score.max` median **0.41** (32% >0.5, 8% >0.7), and 55%
of M02 cells sit in M01's de-novo `a` bucket — the labeling is soft *and* M01-shaped.

**The fix: anchor every cell to external reference knowledge, not to M01.** Type
M01 and M02 *identically* against the same fixed external knowledge, **no M01→M02
transfer** — comparability then falls out, because neither dataset defines the
other's labels.

**"Markers vs atlas" is a false binary — drop it.** Canonical markers *are*
external reference knowledge (ImmGen *is* the immune atlas; Epcam/Pecam1/Col1a1/
Rgs5/Acta2/Adipoq are atlas+literature distilled to a few diagnostic genes).
Markers and a full atlas are the low- and high-resolution ends of the *same* thing
— using external priors. The real decision is two orthogonal axes:
1. **Anchor** — M01 (circular) vs external (de-anchored). *This* is the fix.
2. **Resolution / cost** — a few curated marker genes (panel-native, dataset-
   neutral, *no projection*, coarse) ↔ a full external profile through a classifier
   (finer in principle, but needs the profile *in panel space* — a GeoMx→CosMx
   projection; on 950 genes the extra resolution barely materializes).

Caveat that pins the cheap end: the on-disk `ImmuneAtlas_ImmGen_derived_from_M01.csv`
is **reverse-engineered from M01** (M01's own mean expression under ImmGen labels),
*not* a clean external reference — which is exactly why direct InSituType against it
failed on M02 (see Key decisions). A genuinely-external panel-space ImmGen needs the
same projection as MCA. So the only *free*, fully-de-anchored, panel-native option is
**curated markers**.

**Three tiers (decided direction):**
- **Near-term (decided, fits `processing.smk`):** curated **marker lineage scoring
  applied identically to M01 + M02** — coarse immune lineages (T/B/NK/myeloid/DC) +
  non-immune lineages (epithelial, endothelial, fibroblast, pericyte, smooth-muscle,
  adipocyte) + marker-called 4T1 (`tumor_epithelial`) + honest `unassigned`. Dataset-
  neutral, panel-native, no transfer, no projection. De-anchors all three compartments.
  **Trades fine-subtype resolution for honesty** — fixes the bias, not the panel-limited
  softness, so trust *coarse lineage*, not fine subsets (which the transfer already
  typed weakest).
- **Gold standard (`aggregate.smk`):** integrate M01+M02 into one batch-corrected
  embedding → joint-cluster → annotate clusters from the *same* external marker/
  reference knowledge. Identity is defined **once on the joint object**, so cross-
  dataset comparability is structural. Heavier (3.3M-cell memory; r-harmony LAPACK
  constraint, memory `project_harmony_lapack`).
  **(2026-06-01: this is the authoritative cross-dataset typing direction. `aggregate.smk`
  first detoured to per-cell InSituType, which failed validation; aggregate typing is
  returning to this joint-cluster-then-annotate approach. Single source of truth for the
  method is now `plan-aggregate.md` v2.3.)**
- **Optional upgrade (prepublication):** project genuinely-external profiles onto the
  950 panel for *named, finer* types — GeoMx MCA `MammaryGland` (non-immune subtypes)
  and/or external ImmGen (fine immune) — **gated on a projection-quality check**. MCA is
  normal mammary (no carcinoma), so 4T1 stays marker-called regardless and the gain is
  marginal on this panel. Pair with a **4T1-specific subtype reference/signature** (the
  malignant compartment is one marker-called bucket today). This is the "full atlas +
  tumor-subtype" manuscript scope, not the near-term fix.

**Reference-availability facts (verified 2026-05-31, still true):** there is no panel-
native mouse mammary profile. NanoString's on-disk **CosMx-Cell-Profiles** mouse tier
is **Brain only** (different 1k panel, **177/950 overlap**, missing Epcam/Krt8/Cdh5) —
unusable; the MCA mammary profiles live in the separate **CellProfileLibrary** (GeoMx-
WTA), so any atlas use (tier 3) requires the GeoMx→CosMx projection.

**Downstream impact:** re-typing moves every cell-type-keyed output — the processing
typing stage AND all of `aggregate.smk`'s cell-type tracks (composition, pseudobulk
DE, GSEA, pathway, per-cell DE, concordance) must re-run on the new labels. Companion
plans to bump: `plan-aggregate.md`, `plan-outcomes-de.md`, `plan-mbrt-signatures.md`.
Build scoped in the **v2.0 stem** at the bottom of this file.

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
   sink under the cross-dataset platform shift (see Decisions). **v2.0: this
   transfer is M01-anchored (M02 labels are derived from M01 → circular for a
   cross-dataset comparison); being replaced by de-anchored marker lineage scoring
   applied identically to M01 + M02, no transfer. See the v2.0 banner + stem.**
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
- **Atlas expansion for the non-immune lineages — SUPERSEDED 2026-05-31** by the
  de-anchored framing (see v2.0 banner + stem): "markers vs atlas" was a false
  binary (both are external knowledge), so the MCA atlas is demoted to the optional
  prepublication tier and dataset-neutral markers are the near-term path. The
  *verified facts* below remain valid. ImmGen is immune-only, so 48.5% of the 3.27M-cell flank cohort
  (1.59M cells) sits in the InSituType de novo `a` bucket and the M02 tumor
  compartment is almost entirely unassigned (1.28M of 1.59M); well-matched
  reference mapping leaves single-digit-to-15% unassigned, so this is reference
  mismatch, not noise. **Decision (original, now superseded — see v2.0 banner +
  stem): type against the NanoString CellProfileLibrary mouse mammary atlas**
  (markers cast as only a break-glass contingency). **Corrected:** markers and
  atlas are one external-knowledge spectrum, not a binary — near-term uses dataset-
  neutral curated markers applied identically to M01 + M02; the MCA atlas is the
  optional prepublication upgrade. Load-bearing caveat (verified 2026-05-31): the
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
- **Next** (priority re-ordered 2026-05-31 — **de-anchored re-typing is now #1**: all
  cell-type-resolved aggregate outputs already built, i.e. composition, pseudobulk
  DE, GSEA, pathway, embed_celltype, are computed on the current ImmGen+`a` labels
  and must be re-run once the types change, so lock typing before stacking more
  analysis on labels that will move):
  1. **De-anchor the typing** (see v2.0 banner + stem) — replace the M01→M02 anchor
     transfer with **dataset-neutral marker lineage scoring** (`config/lineage_markers.yaml`,
     to author) applied **identically to M01 + M02**: coarse immune lineages +
     non-immune lineages (epithelial, endothelial, fibroblast, pericyte, smooth-
     muscle, adipocyte) + marker-called 4T1 (`tumor_epithelial`) + honest
     `unassigned`. No transfer, no projection. NB the on-disk ImmGen CSV is
     **M01-derived/tainted**, so a genuinely-external panel-space ImmGen (and MCA)
     is the tier-3 upgrade, not the near-term fix. Pericyte is the known weak spot
     (3/8 markers). **Benchmark**: report coarse-lineage fractions M01 vs M02 (they
     must no longer be M01-derived); cut unassigned toward the well-matched <~15%
     norm. Pilot on one tumor-heavy M02 sample before the full 23-sample re-run.
  2. **Fold the cell-typing/pathway/probe-QC layer into the literate
     `spatial-rads.org`** — currently lives only in `scripts/` + `workflows/`;
     the org's `processing.smk` block is the **stale 10-rule version** and
     would *regress* the workflow back to no-celltype/no-pathway if tangled
     (landmine, removed by the fold). Fold the de-anchored re-typing in at the same time.
  3. **Build `aggregate.smk`** — merge → Harmony → cluster → UMAP → cell-type
     landscape → cross-dataset DEGs. This is also the **gold-standard typing tier**
     (integrate → joint-cluster → annotate from the same external marker/reference
     knowledge; v2.0 stem tier 2) — identity defined once on the joint object.
  4. **Review `probe_qc_report.tsv`** flagged candidates (smooth, not drop).

## Key files

`workflows/{data_model,processing}.smk` · `scripts/*.R` ·
`config/{config,common,pathway_gene_lists}.yaml` ·
`data/{metadata.xlsx,metadata_schema.yaml,data_model.rda}` · `results/processing/*`

## De-anchored re-typing — work plan (v2.0 stem)

**Status: stem.** Direction decided; expand into a committed sub-plan on resume.
Supersedes v1 stage 5 (Yi-ImmGen on M01 + M01→M02 anchor transfer) and the earlier
"Atlas re-typing" stem (which framed markers-vs-atlas as a binary — it isn't).

### Job to be done
Assign every flank cell (M01 + M02, 3.27M) to a lineage by scoring it against
**fixed external reference knowledge**, applied **identically to M01 and M02 with no
M01→M02 transfer** — so cross-dataset comparability falls out structurally (neither
dataset defines the other's labels) instead of being assumed.

### The two axes (don't conflate them)
1. **Anchor** — M01 (current, circular) vs external (de-anchored). *This is the fix.*
2. **Resolution / cost** — curated marker genes (panel-native, dataset-neutral, no
   projection, coarse) ↔ a full external profile via classifier (finer in principle,
   needs the profile *in panel space*; on 950 genes the extra resolution is marginal).

"Markers vs atlas" is a false binary: both are external prior knowledge — markers are
the atlas/literature distilled to a few diagnostic genes. The only *free*, fully-de-
anchored, panel-native option is **curated markers**; anything richer needs a
GeoMx→CosMx projection (incl. a genuinely-external ImmGen — the on-disk CSV is
M01-derived, see Inputs).

### Compartment-specific anchoring (what de-anchoring buys, per compartment)
- **Immune** — currently anchored to M01's ImmGen labels, yet already the weakest-
  resolved (fine subtypes typed at lowest confidence under transfer). Markers give
  honest coarse lineages; fine immune subtypes need external ImmGen (tier 3).
- **Non-immune** — inherited M01's failure-to-resolve (the 48.5% `a` blob). Markers
  de-anchor *and populate* these lineages directly — the biggest near-term gain.
- **Tumor (4T1)** — marker-called regardless (no off-the-shelf carcinoma reference);
  de-anchoring doesn't change the mechanism, only removes the M01-bucket dependence.

### Tier 1 — decided near-term path (fits `processing.smk`)
- **Curated marker lineage scoring** (UCell/AddModuleScore over
  `config/lineage_markers.yaml`, to author), applied identically to M01 + M02:
  coarse immune (T/B/NK/myeloid/DC) + non-immune (epithelial, endothelial,
  fibroblast, pericyte, smooth-muscle, adipocyte) + marker-called `tumor_epithelial`
  + honest `unassigned`.
- No transfer, no projection, no M01 dependence. Trades fine-subtype resolution for
  honesty — trust *coarse lineage*, not fine subsets.
- Revise `celltype.R` (and largely retire `build_celltype_reference.R`'s pooled-M01
  reference) so scoring replaces the anchor transfer on both datasets.

### Tier 2 — gold standard (`aggregate.smk`)
Integrate M01 + M02 into one batch-corrected embedding → joint-cluster → annotate
clusters from the same external marker/reference knowledge. Identity defined **once on
the joint object** → comparability is structural. Heavier (3.3M-cell memory; r-harmony
LAPACK ≤1 group var, LISI gate ~1.5 — memory `project_harmony_lapack`). The deferred
"Option C".

### Tier 3 — optional prepublication upgrade (gated on projection QC)
Project genuinely-external profiles onto the 950 panel for *named, finer* types:
- **External ImmGen** (fine immune) — needs the same GeoMx→CosMx projection as MCA,
  because the on-disk CSV is M01-derived/tainted.
- **GeoMx MCA `MammaryGland`** (Han et al. Cell 2018, PMID 29474909; non-immune
  subtypes) — normal mammary, **no carcinoma**, so 4T1 stays marker-called and the
  gain is marginal on this panel.
- **4T1-specific subtype reference/signature** — the malignant compartment is one
  marker-called bucket today; a 4T1 scRNA reference would split it. Future work.

### Inputs / what's on disk
- ImmGen panel-space CSV — **M01-derived/tainted, NOT clean external**:
  `analysis/objects/mbrt_vs_sbrt/ImmuneAtlas_ImmGen_derived_from_M01.csv`.
- CosMx-Cell-Profiles repo — **NOT usable** (mouse = brain, 177/950 overlap):
  `/mnt/data/projects/spatial-rads/CosMx-Cell-Profiles`.
- GeoMx MCA `MammaryGland` — **NOT downloaded** (tier 3 only):
  github.com/Nanostring-Biostats/CellProfileLibrary `Mouse/Adult`.
- Marker sets — from literature/atlas knowledge → `config/lineage_markers.yaml` (to author).
- Our 950 panel: `results/processing/common_genes.tsv`.
- 23 per-sample objects: `/mnt/data/projects/spatial-rads/processing/{norm,typed,scored}/`.

### Measured evidence (2026-05-31)
- Current transfer is soft AND M01-shaped: pooled M02 `prediction.score.max` median
  **0.41** (32% >0.5, 8% >0.7); **55%** of M02 cells fall in M01's de-novo `a` bucket.
- InSituType-direct against the M01-derived profile failed on M02 (Lymphatic.endo 38%
  sink, tumor 0.15%, prob median 0.99) but typed an M01 control correctly — i.e. the
  failure was the M01-tainted reference under cross-dataset shift, not the method alone.

### Open decisions to settle on the pilot
1. Coarse-lineage marker sets + scoring thresholds (`lineage_markers.yaml`).
2. `unassigned` policy — score-margin / minimum-score cutoff.
3. Malignant 4T1 marker set + threshold (generalize the existing `tumor_epithelial` relabel).
4. Whether tier 3 (external ImmGen / MCA projection) is worth it for this manuscript.

### Panel adequacy (verified, `common_genes.tsv`)
Epithelial 6/8, endothelial 7/9, fibroblast 8/8, smooth-muscle 4/5, **pericyte 3/8**
(Rgs5/Pdgfrb/Notch3 present; Cspg4/Mcam/Des/Kcnj8/Abcc9 absent — weakest call,
pericyte-vs-smooth-muscle). The unassigned rate is a *reference/anchor* gap, not a
*panel* gap.

### Success criteria
**M02 labels are NOT derived from M01** (the load-bearing one). Coarse lineages
populated at plausible fractions on both datasets; unassigned toward the well-matched
<~15% norm; tumor compartment coherent; M01/M02 distributions concordant where biology
matches — now a real concordance check, not a transfer artifact.

### Downstream ordering
Re-typing moves labels, so it must precede `aggregate.smk` cell-type tracks,
`plan-outcomes-de.md` findings, and `plan-mbrt-signatures.md` per-cell-type
signatures. Lock typing before stacking analysis on labels that will move.
