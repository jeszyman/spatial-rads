# MBRT vs SBRT (vs Control) analysis — implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the cross-dataset MBRT-vs-SBRT-vs-Control differential analysis defined in `plan-mbrt-vs-sbrt.md`, consuming the locked unified per-cell label table, and emit one tier-tagged master results table plus figures.

**Architecture:** Standalone arg-driven R/Python scripts in `scripts/aggregate/`, runnable directly and wired into `workflows/aggregate.smk`. All cell-type identity comes from `results/aggregate/full_labels.parquet` (3,277,090 cells: `cell`, `compartment` = tumor/stroma/immune, `cell_subtype`, `subtype_source`, `rescued`). Counts come from `/mnt/data/projects/spatial-rads/aggregate/merged_typed.rds`. Inference is on Mutter_02 (M02) day-2 only, n=4/arm, randomized complete block (`~ slide_id + condition`); Mutter_01 (M01) is descriptive effect sizes only.

**Tech Stack:** R (Seurat v5, data.table, arrow, DESeq2, apeglm, speckle/propeller, limma, UCell, fgsea, RANN), conda env `spatial-rads`; Python (scanpy/scvi only where already used) env `spatial-rads-scvi`; snakemake driver in `basecamp`.

---

## Current state (what already exists — do not rebuild)

Verified on disk 2026-06-02:

- **`composition.R` — DONE and correct.** Reads `obs.parquet` + `full_labels.parquet` (uses `cell_subtype`), propeller `~0+condition+slide_id`, three contrasts (`MBRT_vs_Ctrl`, `SBRT_vs_Ctrl`, `MBRT_vs_SBRT`) with 95% CIs, global BH. Output `composition_test_m02day2.tsv` dated today (ran standalone). The M02 `condition` factor levels are literally `Control` / `MBRT_day2` / `SBRT_day2`.
- **`deg_pseudobulk.R` — correct, consumes the SummarizedExperiment only.** DESeq2 `~slide_id+condition`, three contrasts, apeglm shrinkage, abundance floor (≥10 cells in ≥3/4 samples per condition per cell type), gene filter (≥10 counts in ≥4 samples), per-contrast BH.
- **`gsea.R` — correct.** fgsea on DESeq2 Wald stat, the 4 YAML sets + 50 MSigDB Hallmark (mouse), BH within contrast. Reads `config/pathway_gene_lists.yaml`.
- **`pathway_summary.R` — running standalone now** on `merged.rds` + `full_labels.parquet`: UCell + AddModuleScore per cell → per-sample×cell-type means → limma (M02), plus M01 timecourse + concordance scatter.

**The one stale dependency:** `pseudobulk_se.rds` (and therefore `degs_pseudobulk_m02day2.tsv`, `gsea_*.tsv`) was built before typing was locked, from `pseudobulk_build.R:23` reading `cell_type_atlas` — the **rejected** per-cell InSituType labels (84.6% mistyped tumor). The DE chain must be rebuilt on `full_labels.parquet` (Task 1).

**`dev/mbrt_vs_sbrt/` is NOT reused as code.** Those scripts read an old `seurat_filtered.rds` and an old 9-category `cell_type_major` typing, hardcode paths, and define gene sets inline. They are cited only as **method provenance**; the spatial/niche/mixing/myeloid steps are re-implemented clean on the canonical inputs (Tasks 9–12), and their prior findings are re-tested, not imported.

---

## Phase 0 — Rebuild the confirmatory DE chain on correct labels; reconcile the workflow

### Task 1: Re-point pseudobulk to the unified labels

**Files:**
- Modify: `scripts/aggregate/pseudobulk_build.R` (args + lines 8, 16-19, 22-23)

- [ ] **Step 1: Add a labels arg and join `cell_subtype` instead of `cell_type_atlas`.**

Change the arg header comment (line 8) to:
```r
# Args: <merged_typed.rds> <full_labels.parquet> <out_se.rds> <out_qc.tsv>
```
Add `library(arrow)` to the `suppressPackageStartupMessages` block. Replace the arg parse (lines 16-19) with:
```r
args        <- commandArgs(trailingOnly = TRUE)
merged      <- args[1]
labels_path <- args[2]
out_se      <- args[3]
out_qc      <- args[4]
```
Replace lines 22-23 (`md <- o@meta.data` / `md$cell_type <- md$cell_type_atlas`) with:
```r
md  <- o@meta.data
lab <- as.data.table(read_parquet(labels_path))[, .(cell, cell_subtype)]
md$cell_type <- lab$cell_subtype[match(rownames(md), lab$cell)]
```

- [ ] **Step 2: Rebuild the pseudobulk SummarizedExperiment (standalone).**

Run:
```bash
TMPDIR=/mnt/data/projects/spatial-rads/tmp OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
conda run -n spatial-rads Rscript scripts/aggregate/pseudobulk_build.R \
  /mnt/data/projects/spatial-rads/aggregate/merged_typed.rds \
  results/aggregate/full_labels.parquet \
  /mnt/data/projects/spatial-rads/aggregate/pseudobulk_se.rds \
  results/aggregate/pseudobulk_qc.tsv
```
Expected stdout: `pseudobulk: <N> M02 cells -> <G> (sample x cell_type) columns x 950 genes | <K> cell types, 12 samples`. K should match the M02 cell_subtypes present (≈15–20), NOT a single tumor-dominated bucket.

- [ ] **Step 3: Sanity-check the new colData cell types.**

Run:
```bash
conda run -n spatial-rads Rscript -e 'se<-readRDS("/mnt/data/projects/spatial-rads/aggregate/pseudobulk_se.rds"); print(sort(table(SummarizedExperiment::colData(se)$cell_type)))'
```
Expected: multiple cell types (immune subtypes, fibroblast, endothelial, tumor, etc.), not ~100% one label.

- [ ] **Step 4: Commit.**
```bash
git add scripts/aggregate/pseudobulk_build.R
git commit -m "fix(aggregate): pseudobulk reads unified full_labels, not rejected atlas labels"
```

### Task 2: Re-run the DE + GSEA chain

**Files:**
- Output: `results/aggregate/degs_pseudobulk_m02day2.tsv`, `deg_summary_m02day2.tsv`, `degs_pseudobulk_skipped.tsv`, `gsea_pseudobulk_m02day2.tsv`

- [ ] **Step 1: Run DESeq2.** No volcano plots: the -log10 p axis + padj<0.05 cutoff reintroduce significance-first framing, which is invalid on this sparse targeted panel (no valid padj floor when most genes sit near the detection floor). Effect-size + CI + MDE carry the inference instead.
```bash
TMPDIR=/mnt/data/projects/spatial-rads/tmp OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
conda run -n spatial-rads Rscript scripts/aggregate/deg_pseudobulk.R \
  /mnt/data/projects/spatial-rads/aggregate/pseudobulk_se.rds \
  results/aggregate/degs_pseudobulk_m02day2.tsv \
  results/aggregate/deg_summary_m02day2.tsv \
  results/aggregate/degs_pseudobulk_skipped.tsv
```
Expected: `deg_summary_m02day2.tsv` lists per cell type × contrast counts of significant genes; `degs_pseudobulk_skipped.tsv` lists cell types failing the abundance floor.

- [ ] **Step 2: Run GSEA.**
```bash
conda run -n spatial-rads Rscript scripts/aggregate/gsea.R \
  results/aggregate/degs_pseudobulk_m02day2.tsv \
  config/pathway_gene_lists.yaml \
  results/aggregate/gsea_pseudobulk_m02day2.tsv
```
Expected: rows for each (cell_type × contrast × pathway) with NES + padj_bh.

- [ ] **Step 3: Commit results.**
```bash
git add results/aggregate/degs_pseudobulk_m02day2.tsv results/aggregate/deg_summary_m02day2.tsv results/aggregate/degs_pseudobulk_skipped.tsv results/aggregate/gsea_pseudobulk_m02day2.tsv
git commit -m "results(aggregate): rebuild M02 day2 pseudobulk DE + GSEA on unified labels"
```

### Task 3: Reconcile `aggregate.smk` rule args with the current script signatures

**Files:**
- Modify: `workflows/aggregate.smk` rules `composition` (lines 188-206), `pseudobulk_build` (209-220), `pathway_summary` (256-274)

- [ ] **Step 1: Fix the `composition` rule inputs/shell.** It currently passes `cell_metadata.tsv` (TSV) where the script reads a parquet, and `cell_atlas_labels.tsv` (rejected labels) where it expects the unified parquet. Set:
```python
    input:
        script = "scripts/aggregate/composition.R",
        obs    = f"{AGG}/full/obs.parquet",
        labels = "results/aggregate/full_labels.parquet",
    ...
    shell:
        "{RSCRIPT} {input.script} {input.obs} {input.labels} {output.by_sample} "
        "{output.test} {output.dropped} {output.bars} {output.forest} "
        "{output.timecourse} > {log} 2>&1"
```

- [ ] **Step 2: Fix the `pseudobulk_build` rule** to pass the labels parquet (matches Task 1):
```python
    input:
        script = "scripts/aggregate/pseudobulk_build.R",
        rds    = f"{AGG}/merged_typed.rds",
        labels = "results/aggregate/full_labels.parquet",
    ...
    shell:
        "{RSCRIPT} {input.script} {input.rds} {input.labels} {output.se} {output.qc} > {log} 2>&1"
```

- [ ] **Step 3: Fix the `pathway_summary` rule** to pass the labels parquet as arg 2 (the running standalone invocation does; the rule omits it):
```python
    input:
        script = "scripts/aggregate/pathway_summary.R",
        rds    = f"{AGG}/merged_typed.rds",
        labels = "results/aggregate/full_labels.parquet",
        yaml   = "config/pathway_gene_lists.yaml",
    ...
    shell:
        "{RSCRIPT} {input.script} {input.rds} {input.labels} {input.yaml} {output.summary} "
        "{output.test} {output.conc} {output.heatmap} {output.timecourse} "
        "{output.scatter} > {log} 2>&1"
```
Confirm the actual `pathway_summary.R` arg order by reading its `commandArgs` block and match exactly.

- [ ] **Step 4: Dry-run.**
```bash
TMPDIR=/mnt/data/projects/spatial-rads/tmp conda run -n basecamp snakemake -s workflows/aggregate.smk --dry-run composition pseudobulk_build deg_pseudobulk gsea pathway_summary
```
Expected: DAG resolves, no missing-input or wildcard errors.

- [ ] **Step 5: Commit.**
```bash
git add workflows/aggregate.smk
git commit -m "fix(aggregate.smk): reconcile composition/pseudobulk/pathway rule args with rewritten scripts"
```

---

## Phase 1 — Build the missing analysis inputs

### Task 4: Build the four missing gene sets, panel-filtered

**Files:**
- Modify: `config/pathway_gene_lists.yaml` (append 4 sets)
- Create: `scripts/aggregate/build_gene_sets.R`

- [ ] **Step 1: Append candidate sets to the YAML** (panel-filtered in Step 2; these are candidates):
```yaml
Angiogenesis: [Vegfa, Vegfc, Flt1, Kdr, Pecam1, Cdh5, Vwf, Angpt1, Angpt2, Tek, Pdgfb, Notch1, Dll4, Esm1, Apln, Aplnr, Nrp1, Cldn5, Tie1]
Hypoxia: [Hif1a, Vegfa, Slc2a1, Car9, Ldha, Pgk1, Bnip3, Ndrg1, Pdk1, Eno1, Aldoa, Hk2, Egln3]
Fibrosis_remodeling: [Col1a1, Col1a2, Col3a1, Col5a1, Acta2, Tagln, Fn1, Tgfb1, Tgfbr1, Ccn2, Postn, Sparc, Timp1, Mmp2, Mmp9, Lox, Fap, Vim]
Stromal_stress_senescence: [Cdkn1a, Cdkn2a, Trp53, Glb1, Serpine1, Il6, Cxcl1, Ccl2, Mmp3, Igfbp3, Gdf15, Bcl2l1]
```

- [ ] **Step 2: Write the coverage-audit script.** DEVIATION FROM ORIGINAL PLAN (2026-06-02): the script does **not** rewrite the YAML. Both scorers (`pathway_score.R`, `pathway_summary.R`) already drop off-panel genes at scoring time, and the YAML's shipped convention is full curated lists (the header documents drop-at-scoring). Permanently filtering the YAML would conflict with that tested pattern and stale the provenance comments (module sizes 19/17/18). So the script only emits the per-gene coverage audit; the YAML keeps full lists. It also reads the panel from `results/processing/common_genes.tsv` (the `PANEL` var in the smk) instead of loading the 1.3 GB `merged_typed.rds` for a name lookup.
```r
#!/usr/bin/env Rscript
# Audit each gene set in config/pathway_gene_lists.yaml against the common 950-gene
# panel; log each gene's panel membership. Does NOT rewrite the YAML (full curated
# lists kept; off-panel drops at scoring time).
# Args: <common_genes.tsv> <pathway_gene_lists.yaml> <out_coverage.tsv>
suppressPackageStartupMessages({library(yaml); library(data.table)})
a     <- commandArgs(trailingOnly = TRUE)
panel <- readLines(a[1])
sets  <- yaml::read_yaml(a[2])
cov <- rbindlist(lapply(names(sets), function(s) data.table(
  set = s, gene = sets[[s]], in_panel = sets[[s]] %in% panel)))
fwrite(cov, a[3], sep = "\t")
```

- [ ] **Step 3: Run it.**
```bash
conda run -n spatial-rads Rscript scripts/aggregate/build_gene_sets.R \
  results/processing/common_genes.tsv \
  config/pathway_gene_lists.yaml \
  results/aggregate/gene_set_panel_coverage.tsv
```
Realized coverage: Angiogenesis 14/19, Fibrosis_remodeling 13/18, Stromal_stress_senescence 8/12, **Hypoxia 5/13** (panel-thin). Hypoxia limitation recorded in `plan-mbrt-vs-sbrt.md` (Panel blind spots) and the YAML comment.

- [ ] **Step 4: Commit.**
```bash
git add config/pathway_gene_lists.yaml scripts/aggregate/build_gene_sets.R results/aggregate/gene_set_panel_coverage.tsv
git commit -m "feat(aggregate): add angiogenesis/hypoxia/fibrosis/stromal-stress gene sets, panel-filtered"
```

### Task 5: Per-cell-type detection report for every readout gene

**Files:**
- Create: `scripts/aggregate/panel_coverage.R`
- Output: `results/aggregate/readout_detection_m02.tsv`

- [ ] **Step 1: Write the script** — for every gene in every YAML set, report fraction of M02 cells with count > 0, overall and within the compartment the hypothesis targets (immune for H1, stroma/endothelial for H2/H3). Use `full_labels.parquet` for compartment, `merged_typed.rds` counts. Output one row per (gene, set, compartment, detect_frac, mean_count).

- [ ] **Step 2: Run + inspect.** Confirm each confirmatory readout has genes detected in >5% of its target compartment cells; genes below that are flagged "near-floor" so a null on them is not read as biological.

- [ ] **Step 3: Commit.**
```bash
git add scripts/aggregate/panel_coverage.R results/aggregate/readout_detection_m02.tsv
git commit -m "feat(aggregate): per-compartment detection report for all readout genes"
```

### Task 6: Power / minimum-detectable-effect table at n=4

**Files:**
- Create: `scripts/aggregate/power_mde.R`
- Output: `results/aggregate/power_mde.tsv`

- [ ] **Step 1: Write the script.** For the blocked design (n=4/arm, 6 residual df) compute, at 80% power, α=0.05: (a) composition — minimum detectable logit-proportion difference from the observed per-cell-type between-sample SD in `composition_by_sample.tsv`; (b) pseudobulk DE — minimum detectable log2FC at the median per-gene dispersion from the DESeq2 fit (read dispersions from `pseudobulk_se.rds` via `DESeq2::estimateDispersions` or reuse the fit) using a two-sample t approximation on log-CPM; (c) program scores — minimum detectable score delta from the observed per-sample×cell-type score SD. Use `pwr::pwr.t.test` (n=4 per group). One row per (readout_class, cell_type_or_gene_summary, MDE, observed_SD).

- [ ] **Step 2: Run + inspect.** Sanity: MDEs are finite and larger for sparse cell types. This table is the yardstick every null is reported against.

- [ ] **Step 3: Commit.**
```bash
git add scripts/aggregate/power_mde.R results/aggregate/power_mde.tsv
git commit -m "feat(aggregate): power / minimum-detectable-effect table for n=4 blocked design"
```

### Task 7: Per-cell coordinate + necrosis-flag table

**Files:**
- Create: `scripts/aggregate/coords_necrosis.R`
- Modify: `scripts/aggregate/merge.R:100-101` (add coords to the durable slim cache)
- Output: `/mnt/data/projects/spatial-rads/aggregate/coords_necrosis.parquet` (heavy per-cell parquet → data disk, not git-tracked `results/`)

- [ ] **Step 1: Write the extraction + necrosis script.** Iterate the 20 per-sample flank `scored.rds`, build the global cell id `paste0(sample_id, "_", colnames(counts))` (matching `merge.R:50`), pull `x_slide_mm`, `y_slide_mm`. For each `sample_id` (each a contiguous tissue region — the right unit; physical slides carry multiple non-overlapping M02 regions), compute the mean distance to the 20 nearest neighbors with `RANN::nn2(coords, k=21)` (drop self), and set `necrosis_zone = mean_knn_dist > quantile(., 0.90)` within that sample. This is the dev necrosis rule (`dev/mbrt_vs_sbrt/00_load_and_filter.R`) re-keyed to sample-level coordinate space. Write `cell, sample_id, x_slide_mm, y_slide_mm, mean_knn_dist, necrosis_zone` to parquet.

```r
#!/usr/bin/env Rscript
# Args: <samples.tsv> <out_coords_necrosis.parquet> <scored.rds...>
suppressPackageStartupMessages({library(Seurat); library(data.table); library(arrow); library(RANN)})
a <- commandArgs(trailingOnly = TRUE)
out <- a[2]; rds <- a[-(1:2)]
res <- rbindlist(lapply(rds, function(p){
  s <- sub("\\.scored\\.rds$", "", basename(p)); o <- readRDS(p); m <- o@meta.data
  d <- data.table(cell = paste0(s, "_", rownames(m)), sample_id = s,
                  x_slide_mm = m$x_slide_mm, y_slide_mm = m$y_slide_mm)
  xy <- as.matrix(d[, .(x_slide_mm, y_slide_mm)])
  nn <- RANN::nn2(xy, k = min(21, nrow(xy)))
  d[, mean_knn_dist := rowMeans(nn$nn.dists[, -1, drop = FALSE])]
  d[, necrosis_zone := mean_knn_dist > quantile(mean_knn_dist, 0.90)]
  d
}))
write_parquet(res, out)
cat(sprintf("coords: %d cells, %.1f%% necrosis-flagged\n", nrow(res), 100*mean(res$necrosis_zone)))
```

- [ ] **Step 2: Run it over the flank scored.rds.**
```bash
conda run -n spatial-rads Rscript scripts/aggregate/coords_necrosis.R \
  results/data_model/samples.tsv \
  results/aggregate/coords_necrosis.parquet \
  /mnt/data/projects/spatial-rads/processing/scored/sam000{1..8}.scored.rds \
  /mnt/data/projects/spatial-rads/processing/scored/sam00{12..23}.scored.rds
```
(Use the 20 flank sample ids; exclude the 3 tongue samples sam0009-0011.) Expected: row count ≈ 3.27M, ~10% necrosis-flagged per sample by construction.

- [ ] **Step 3: Verify join key matches the labels.**
```bash
conda run -n spatial-rads Rscript -e 'library(arrow); c<-read_parquet("results/aggregate/coords_necrosis.parquet"); l<-read_parquet("results/aggregate/full_labels.parquet"); cat("overlap:", mean(l$cell %in% c$cell), "\n")'
```
Expected: overlap ≈ 1.0 (every labeled cell has coordinates).

- [ ] **Step 4: Patch `merge.R` for durability** — add `"x_slide_mm", "y_slide_mm"` to `meta_cols` (line 100-101) so future merges carry coords in the slim cache. (No re-run of the 17 GB merge needed now; the parquet from Step 2 is the working source.)

- [ ] **Step 5: Commit.**
```bash
git add scripts/aggregate/coords_necrosis.R scripts/aggregate/merge.R results/aggregate/coords_necrosis.parquet
git commit -m "feat(aggregate): per-cell coords + per-sample necrosis flag for spatial tracks"
```

---

## Phase 2 — Re-implement the spatial / structural tracks on canonical inputs

Each of these is a clean re-implementation. Common input contract: join `results/aggregate/full_labels.parquet` (compartment, cell_subtype) + `/mnt/data/projects/spatial-rads/aggregate/coords_necrosis.parquet` (x/y, necrosis_zone) by `cell`; pull metadata (sample_id, dataset, condition, slide_id, timepoint_h) from `/mnt/data/projects/spatial-rads/aggregate/full/obs.parquet`. Method/parameters are ported from the cited dev script; the object and labels are new, so each prior finding is **re-tested, not imported**. Convention: heavy per-cell parquets land on the data disk (`/mnt/data/projects/spatial-rads/aggregate/`); only small TSVs + plots go to git-tracked `results/aggregate/`.

### Task 8: Cell-type-label QC (gate before anything rests on the labels)

**Files:**
- Create: `scripts/aggregate/celltype_qc.R`
- Output: `results/aggregate/plots/celltype_qc_dotplot.png`, `results/aggregate/celltype_qc_markers.tsv`

- [ ] **Step 1: Write the script.** Subsample merged to ~100k cells, set `Idents` to `cell_subtype` from `full_labels.parquet`, DotPlot of canonical lineage markers (e.g. Epcam/Krt8 tumor; Cd3e/Cd8a T; Cd19/Ms4a1 B; Mzb1/Jchain Plasma; Lyz2/Itgam macrophage; Pecam1/Cdh5 endothelial; Col1a1/Pdgfra fibroblast; Acta2/Myh11 smooth muscle; Cidea adipocyte). Provenance: `dev/mbrt_vs_sbrt/03_cell_type_validation.R` style.

- [ ] **Step 2: Run + VIEW the PNG.** Read the dotplot; confirm each marker enriches in its expected subtype. This is the go/no-go that the unified labels are sound for differential work.

- [ ] **Step 3: Commit.**
```bash
git add scripts/aggregate/celltype_qc.R results/aggregate/celltype_qc_markers.tsv results/aggregate/plots/celltype_qc_dotplot.png
git commit -m "feat(aggregate): cell-type-label QC dotplot on unified labels"
```

### Task 9: Cellular niches

**Files:**
- Create: `scripts/aggregate/niches.R`
- Output: `/mnt/data/projects/spatial-rads/aggregate/niche_per_cell.parquet` (heavy → data disk); `results/aggregate/niche_centroids.tsv`, `niche_frequency.tsv`, `niche_test_m02day2.tsv`, `plots/niche_centroids_heatmap.png`, `plots/niche_frequency_m02.png`

- [ ] **Step 1: Write the script (port method from `dev/mbrt_vs_sbrt/07_niche_clustering.R`).** Per sample, `RANN::nn2(xy, k=21)`; for each cell build the composition vector = fraction of each `compartment`-level type (tumor/stroma/immune; or a coarse cell_subtype grouping) among its 20 neighbors; stack all cells; `kmeans(comp, centers=6, nstart=10, iter.max=50, seed=42)` → niches N1–N6. Emit per-cell niche, per-niche centroid composition, and niche frequency per sample. **New vs dev:** uses unified labels, runs per sample over all 20 flank samples.

- [ ] **Step 2: Add the M02 niche-frequency test.** propeller on niche labels across the three arms (same `getTransformedProps` + `~0+condition+slide_id` + 3 contrasts pattern as `composition.R`), output `niche_test_m02day2.tsv` with CIs.

- [ ] **Step 3: Run + VIEW** the centroid heatmap (confirm 6 interpretable niches) and the M02 frequency plot.

- [ ] **Step 4: Commit.**
```bash
git add scripts/aggregate/niches.R results/aggregate/niche_*.tsv results/aggregate/niche_per_cell.parquet results/aggregate/plots/niche_*.png
git commit -m "feat(aggregate): spatial niches (k=20 NN composition, k-means K=6) + M02 arm test"
```

### Task 10: Immune-neighbor fraction + Keren tumor-immune mixing

**Files:**
- Create: `scripts/aggregate/spatial_mixing.R`
- Output: `results/aggregate/spatial_mixing_per_sample.tsv`, `spatial_mixing_test_m02day2.tsv`, `plots/mixing_m02.png`; `/mnt/data/projects/spatial-rads/aggregate/spatial_mixing_per_cell.parquet` (heavy → data disk)

- [ ] **Step 1: Write the script (port from `dev/mbrt_vs_sbrt/05_spatial_nn.R` + `12_mixing_score.R`).** Per sample, k=20 NN; per cell compute immune-neighbor fraction (fraction of 20 neighbors with `compartment=="immune"`). Keren 2018 mixing score per sample = (tumor↔immune neighbor edges) / (immune↔immune edges). **Necrosis exclusion:** drop `necrosis_zone==TRUE` cells before aggregating at day-2 (the dev convention at timepoint_h ∈ {48,144}). Aggregate to per-sample×condition.

- [ ] **Step 2: Add the M02 test.** limma on per-sample immune-neighbor fraction and mixing score, `~0+condition+slide_id`, three contrasts, CIs → `spatial_mixing_test_m02day2.tsv`.

- [ ] **Step 3: Run + VIEW** `mixing_m02.png` (mixing score and immune-neighbor fraction by arm).

- [ ] **Step 4: Commit.**
```bash
git add scripts/aggregate/spatial_mixing.R results/aggregate/spatial_mixing_*.tsv results/aggregate/spatial_mixing_per_cell.parquet results/aggregate/plots/mixing_m02.png
git commit -m "feat(aggregate): immune-neighbor fraction + Keren mixing score + M02 arm test"
```

### Task 11: Myeloid M1/M2 polarization

**Files:**
- Create: `scripts/aggregate/myeloid_polarization.R`
- Output: `results/aggregate/myeloid_m1m2_scores.tsv`, `myeloid_m1m2_test_m02day2.tsv`, `plots/myeloid_m1m2_m02.png`

- [ ] **Step 1: Write the script (port from `dev/mbrt_vs_sbrt/11_m1_m2_polarization.R`).** Subset to macrophage/myeloid `cell_subtype` cells; UCell score the M1 and M2 marker panels (the 16+16 gene lists in the dev script; panel-filter against the 950 panel and log the kept genes). Per-cell M1/M2 scores + M2/M1 ratio; aggregate to per-sample×condition means.

- [ ] **Step 2: Add the M02 test.** limma on per-sample M1, M2, and M2/M1 ratio, `~0+condition+slide_id`, three contrasts, CIs. This is where the dev "MBRT day-2 M2 skew" finding gets re-tested on correct typing.

- [ ] **Step 3: Run + VIEW** `myeloid_m1m2_m02.png`.

- [ ] **Step 4: Commit.**
```bash
git add scripts/aggregate/myeloid_polarization.R results/aggregate/myeloid_m1m2_*.tsv results/aggregate/plots/myeloid_m1m2_m02.png
git commit -m "feat(aggregate): myeloid M1/M2 polarization scoring + M02 arm test"
```

---

## Phase 3 — Cross-cutting, assembly, narrative

### Task 12: M01↔M02 day-2 concordance on unified labels

**Files:**
- Create: `scripts/aggregate/concordance_m01_m02.R`
- Output: `results/aggregate/concordance_m01_m02.tsv`, `plots/concordance_scatter.png`

- [ ] **Step 1: Write the script.** Compute per-cell-type day-2 effect sizes (log2FC, MBRT-vs-Control and SBRT-vs-Control) for M01 (descriptive, n=1 — simple mean-ratio of normalized expression, no p-values) and for M02 (from `degs_pseudobulk_m02day2.tsv`). Spearman correlation of the per-gene log2FC between cohorts, per cell type; scatter with the diagonal. Supersedes `dev/.../08_set2_validation.R` (which re-typed M02 independently — not needed now that both are jointly typed). This is a **qualitative concordance gate**, not validation.

- [ ] **Step 2: Run + VIEW** the scatter; print per-cell-type Spearman rho.

- [ ] **Step 3: Commit.**
```bash
git add scripts/aggregate/concordance_m01_m02.R results/aggregate/concordance_m01_m02.tsv results/aggregate/plots/concordance_scatter.png
git commit -m "feat(aggregate): M01-vs-M02 day2 effect-size concordance on unified labels"
```

### Task 13: Tier-tagged master results table (tiered FDR)

**Files:**
- Create: `scripts/aggregate/assemble_results.R`
- Output: `results/aggregate/results_master.tsv`

- [ ] **Step 1: Define the confirmatory family explicitly in the script** (a fixed lookup, not inferred):
  - **H1 (immune activation):** composition of the immune compartment + immune subtypes; pseudobulk DE within immune cell types restricted to the TypeI_interferon / TypeII_interferon / STING genes; program scores for those three sets in immune cells.
  - **H2 (vascular/oxygenation):** composition of endothelial cells; pseudobulk DE within endothelial restricted to Angiogenesis genes; program scores for Angiogenesis + Hypoxia in endothelial/tumor.
  - **H3 (stromal sparing):** composition of stromal subtypes; pseudobulk DE within fibroblast/stroma restricted to Fibrosis_remodeling + Stromal_stress_senescence genes; program scores for those in stroma.
  - Each across the three contrasts.

- [ ] **Step 2: Assemble.** Read every results table (composition, pseudobulk DE, GSEA, pathway, niche, mixing, myeloid). Tag each row `tier = "confirmatory"` if it matches the family lookup, else `"exploratory"`. Re-compute BH **within the confirmatory family** (across its rows only) → `padj_confirmatory`; leave each exploratory analysis's own within-analysis BH as `padj_exploratory`. Carry effect size + 95% CI + the matching MDE from `power_mde.tsv`. Emit `results_master.tsv` (one row per readout × cell type × contrast).

- [ ] **Step 3: Run + inspect.** Confirm every row has tier, effect size, CI, an FDR appropriate to its tier, and an MDE. Print counts of confirmatory hits at `padj_confirmatory < 0.05`.

- [ ] **Step 4: Commit.**
```bash
git add scripts/aggregate/assemble_results.R results/aggregate/results_master.tsv
git commit -m "feat(aggregate): tier-tagged master results table with confirmatory-family FDR"
```

### Task 14: Wire the new scripts into `aggregate.smk` and full dry-run

**Files:**
- Modify: `workflows/aggregate.smk` (add rules: `build_gene_sets`, `panel_coverage`, `power_mde`, `coords_necrosis`, `celltype_qc`, `niches`, `spatial_mixing`, `myeloid_polarization`, `concordance_m01_m02`, `assemble_results`; extend `rule all`)

- [ ] **Step 1: Add one rule per new script**, following the existing pattern (`script` as an `input:` dep per `feedback_snakemake_script_tracking`; `threads: 1`; BLAS-pinned via the existing `shell.prefix`; `threads: 4` only for the UCell scripts `pathway_summary`/`myeloid_polarization`/`niches` if they fork). Wire `results_master.tsv` as the terminal target in `rule all`.

- [ ] **Step 2: Full dry-run.**
```bash
TMPDIR=/mnt/data/projects/spatial-rads/tmp conda run -n basecamp snakemake -s workflows/aggregate.smk --dry-run
```
Expected: DAG builds to `results_master.tsv` with every rule's args matching its script's `commandArgs`.

- [ ] **Step 3: Commit.**
```bash
git add workflows/aggregate.smk
git commit -m "feat(aggregate.smk): wire gene-sets/coverage/power/coords/niches/mixing/myeloid/concordance/assembly rules"
```

### Task 15: Write limitations after looking; update the notebook

**Files:**
- Modify: `plan-mbrt-vs-sbrt.md` (Limitations section), `spatial-rads.org` (Methods/Results layer)

- [ ] **Step 1:** After the confirmatory run, fill the spatial-dilution limitation with the *realized* magnitude (e.g., how much the whole-compartment MBRT-vs-Control effect attenuates relative to per-niche or per-subtype), the dose-mismatch status (whether the MBRT mean dose arrived from Fazzari/Mutter), and which confirmatory readouts fell below their MDE (genuine nulls vs underpowered).

- [ ] **Step 2:** Add a durable Methods/Results entry to `spatial-rads.org` summarizing the pipeline and pointing at `results_master.tsv`. Use the `sci-write` skill for any manuscript-grade prose; keep the org notebook entry in the lab-notebook register.

- [ ] **Step 3: Commit.**
```bash
git add plan-mbrt-vs-sbrt.md spatial-rads.org
git commit -m "docs(mbrt-vs-sbrt): limitations written after first confirmatory run + notebook entry"
```

---

## Self-review (writing-plans checklist)

- **Spec coverage:** H1/H2/H3 confirmatory (composition T1/T9, DE T1-2, programs T4-5 + pathway_summary running) ✓; exploratory unbiased DE + 50 Hallmark (T2 gsea) ✓; DNA-damage exploratory (DDR already scored; in gsea/pathway) ✓; spatial niches/mixing/neighbor promoted early (T9-10) ✓; M2 polarization (T11) ✓; direct MBRT-vs-SBRT contrast (every test has it) ✓; M01 timecourse (composition/pathway) + concordance (T12) ✓; tiered FDR + power/MDE + effect sizes/CIs + necrosis exclusion + panel coverage + build gene sets + broaden stroma + pool endothelium (T4-6, T10, T13) ✓; pipeline fixes (T1, T3, T7) ✓; dev reuse-as-method + cell-type QC (T8, Phase-2 provenance notes) ✓; limitations look-first (T15) ✓.
- **Gaps to confirm at execution:** "pool endothelial cells to lift immature-vessel markers" is handled by per-compartment program scoring in endothelial (T5/pathway) — verify endothelial n is adequate when running T5; if not, pool endothelial across samples explicitly in the H2 program test.
- **Placeholders:** none — every code step has runnable content or a precise port spec naming the source script and the exact changes.
- **Naming consistency:** label column = `cell_subtype` / `compartment` throughout; join key = `cell`; contrasts = `MBRT_vs_Ctrl` / `SBRT_vs_Ctrl` / `MBRT_vs_SBRT`; condition factor levels = `Control` / `MBRT_day2` / `SBRT_day2` (matches `composition.R`).
