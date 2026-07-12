# Differential-layer robustness — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

> **⚠️ Note (2026-07-12):** the differential layer lives in `workflows/aggregate_differential.smk`
> (split out of the former monolithic `aggregate.smk`; cell typing is now the separate
> `aggregate_typing.smk`). Rule bodies are unchanged, but any `workflows/...smk:<line>` numbers
> below predate the split — locate rules by name. The Global Constraints below are corrected inline
> for the split.

**Goal:** Implement the six confounding fixes in `plan-differential-robustness.md` on the `aggregate_differential.smk` differential layer, ending in one consolidated rerun that regenerates `results/aggregate/results_master.tsv` with the new detection/sub-state/claim-scoping columns.

**Architecture:** Edit five existing scripts (`pseudobulk_build.R`, `deg_pseudobulk.R`, `composition.R`, `pathway_summary.R`, `assemble_results.R`), add three new ones (`detection_test.R`, `substate_split.R`, optional `milo_da.py`), add one config (`config/substate_markers.yaml`), and rewire the affected Snakemake rules. Most work is verified on existing intermediates (`pseudobulk_se.rds`, `merged.rds`, `full_labels.parquet`) so the expensive UCell recompute happens only once, in the final task.

**Tech Stack:** R 4.4 (`spatial-rads` env: Seurat, data.table, arrow, speckle/limma, DESeq2, SummarizedExperiment, Matrix); Python (`spatial-rads-scvi` env: scanpy, milopy — optional); Snakemake (`basecamp` env driver).

## Global Constraints

- Snakemake invocation: `TMPDIR=/mnt/data/projects/spatial-rads/tmp conda run -n basecamp snakemake -s workflows/aggregate_differential.smk --cores N`; dry-run (`-n`) before every execute.
- R steps run via `RSCRIPT = conda run -n spatial-rads Rscript` (the only interpreter in `aggregate_differential.smk`; the Python `PYSCVI` env is used by the typing workflow, not here).
- Rules keep `threads: 1` except `pathway_scores`/`pathway_summary` (`threads: 4`, fork parallelism). The workflow carries **no `shell.prefix`** — it relies on Snakemake's default strict mode (`set -euo pipefail`); there is no workflow-level BLAS pinning.
- Snakemake does not track `shell:`-invoked scripts unless declared `input:`; every rule below keeps its `script = "scripts/aggregate/<x>"` as an `input` so edits trigger reruns.
- Inference cohort = Mutter_02 day-2, n=4/arm, balanced (3 conditions × 4 slides). Mutter_01 stays descriptive (n=1).
- **No beta-binomial GLM and no cell-as-replicate MAST** for the detection test — use limma empirical-Bayes on arcsin-sqrt per-sample fractions (calibrated at n=4 by borrowing strength across genes).
- **Pre-register sub-state markers in `config/substate_markers.yaml` and commit before any sub-state DE/composition output is inspected** (post-selection-inference guard).
- Testing convention: this is an analysis pipeline with **no unit-test harness**. Each task's "test" is a runtime assertion inside the script (`stopifnot`) plus an acceptance command that runs the script/rule on real data (or a saved slice) and greps the output for the expected property. Show the command and the expected line.

---

## File map

| File | Responsibility | Action |
|---|---|---|
| `scripts/aggregate/pseudobulk_build.R` | sum counts + emit detection assays per (sample×cell_type) | modify |
| `scripts/aggregate/detection_test.R` | limma-eBayes detection + level tests + sparsity-aware classifier | **create** |
| `scripts/aggregate/deg_pseudobulk.R` | DESeq2 per cell_type; carry detection pair | modify |
| `scripts/aggregate/composition.R` | propeller; add unassigned test + sensitivity + sub-state run | modify |
| `config/substate_markers.yaml` | frozen fibroblast resting/activated panels + specificity anchor | **create** |
| `scripts/aggregate/substate_split.R` | gated fibroblast sub-state assignment | **create** |
| `scripts/aggregate/geneset_overlap.R` | gene-set Jaccard / shared-gene diagnostic | **create** |
| `scripts/aggregate/milo_da.py` | optional ~200k-subsample cluster-free DA | **create (optional)** |
| `scripts/aggregate/pathway_scores.R` + `pathway_plots.R` | split UCell compute from plotting | **create (split)** |
| `scripts/aggregate/assemble_results.R` | join detection/class + claim-scoping flags into master | modify |
| `workflows/aggregate_differential.smk` | rewire rules | modify |

---

## Task 1: Detection assays in pseudobulk build

**Files:**
- Modify: `scripts/aggregate/pseudobulk_build.R` (insert after line 50, extend SE at line 56)

**Interfaces:**
- Produces: `pseudobulk_se.rds` SummarizedExperiment now carries three assays — `counts` (unchanged), `pct_expr` (fraction of the group's cells with count>0, per gene), `mean_among_expr` (mean raw count among expressing cells, per gene). All 950 × Ngroups, same `colData`.

- [ ] **Step 1: Add the detection-assay block.** After `n_cells <- as.integer(colSums(G))` (line 50), insert:

```r
# --- detection metrics per (gene x group): fraction expressing + level among expressers ---
pos <- as.matrix(as(cnt > 0, "dgCMatrix") %*% G)     # 950 x Ngroups: # expressing cells
pct_expr <- sweep(pos, 2, n_cells, "/")              # fraction of group's cells expressing
mean_among_expr <- pb / pmax(pos, 1)                 # mean raw count among expressers
mean_among_expr[pos == 0] <- 0
storage.mode(pct_expr) <- "double"; storage.mode(mean_among_expr) <- "double"
stopifnot(all(pct_expr >= 0 & pct_expr <= 1))
```

- [ ] **Step 2: Add both to the SE assays.** Change the `assays = list(counts = pb)` line (≈56) to:

```r
  assays  = list(counts = pb, pct_expr = pct_expr, mean_among_expr = mean_among_expr),
```

- [ ] **Step 3: Run the rule and verify three assays.**

Run:
```bash
cd /home/jeszyman/repos/spatial-rads
TMPDIR=/mnt/data/projects/spatial-rads/tmp conda run -n basecamp snakemake -s workflows/aggregate_differential.smk \
  /mnt/data/projects/spatial-rads/aggregate/pseudobulk_se.rds --cores 4 --force
conda run -n spatial-rads Rscript -e '
  se <- readRDS("/mnt/data/projects/spatial-rads/aggregate/pseudobulk_se.rds")
  cat(SummarizedExperiment::assayNames(se), "\n")
  pe <- SummarizedExperiment::assay(se,"pct_expr"); cat("pct range", range(pe), "\n")'
```
Expected: `counts pct_expr mean_among_expr` and `pct range 0 1` (or within [0,1]).

- [ ] **Step 4: Commit.**
```bash
repo-commit "feat(aggregate): emit per-group detection assays (pct_expr, mean_among_expr) in pseudobulk SE" scripts/aggregate/pseudobulk_build.R
```

---

## Task 2: Detection + level moderated test with sparsity-aware classifier

**Files:**
- Create: `scripts/aggregate/detection_test.R`
- Modify: `workflows/aggregate_differential.smk` (new `detection_test` rule after `deg_pseudobulk`)

**Interfaces:**
- Consumes: `pseudobulk_se.rds` (`pct_expr`, `mean_among_expr` assays + `colData$condition,$slide_id,$sample_id,$cell_type`).
- Produces: `results/aggregate/detection_test_m02day2.tsv` with columns `gene, cell_type, contrast, detection_log2fc, detection_padj, level_log2fc, level_padj, mean_among_expr_max, call_class` where `call_class ∈ {regulation, fraction_shift, ambiguous}`.

- [ ] **Step 1: Write the test script.** The detection model mirrors `composition.R` (arcsin-sqrt + `lmFit(~0+condition+slide_id)` + 3 contrasts + `eBayes(robust)`), looped per cell_type, applied to per-gene detection fractions; the level model is the same on `log2(mean_among_expr+1)`.

```r
#!/usr/bin/env Rscript
# Detection vs level decomposition for the pseudobulk DE (Fix 2).
# Per cell_type: limma-eBayes moderated test on arcsin-sqrt per-sample detection
# fractions, and (in parallel) on log2 per-expresser level. Classifier separates
# regulation from cell-state-fraction shift, sparsity-aware. Args: <se.rds> <out.tsv>
suppressPackageStartupMessages({ library(SummarizedExperiment); library(limma); library(data.table) })
a <- commandArgs(trailingOnly = TRUE); se <- readRDS(a[1]); out <- a[2]
cd <- as.data.table(as.data.frame(colData(se)))
cd[, condition := factor(condition, levels = c("Control","MBRT_day2","SBRT_day2"))]
CONTR <- list(MBRT_vs_Ctrl = c("conditionMBRT_day2","conditionControl"),
              SBRT_vs_Ctrl = c("conditionSBRT_day2","conditionControl"),
              MBRT_vs_SBRT = c("conditionMBRT_day2","conditionSBRT_day2"))
LEVEL_FLOOR <- 3      # mean_among_expr must clear this to call a flat level "reliably measured"
ln2 <- log(2)

fit_one <- function(mat, cols) {                       # mat: genes x samples (transformed)
  samp <- cd[cols]
  d <- model.matrix(~ 0 + condition + slide_id, data = samp); colnames(d) <- make.names(colnames(d))
  cm <- makeContrasts(MBRT_vs_Ctrl = conditionMBRT_day2 - conditionControl,
                      SBRT_vs_Ctrl = conditionSBRT_day2 - conditionControl,
                      MBRT_vs_SBRT = conditionMBRT_day2 - conditionSBRT_day2, levels = d)
  eBayes(contrasts.fit(lmFit(mat, d), cm), robust = TRUE)
}
res <- list()
for (ct in unique(cd$cell_type)) {
  cols <- which(cd$cell_type == ct)
  if (length(unique(cd$condition[cols])) < 3 || length(cols) < 6) next   # need all arms + replicates
  pe  <- assay(se,"pct_expr")[, cols, drop=FALSE]
  mae <- assay(se,"mean_among_expr")[, cols, drop=FALSE]
  keep <- rowSums(pe > 0) >= 3                          # gene seen in >=3 samples of this type
  if (!any(keep)) next
  fd <- fit_one(asin(sqrt(pe[keep,,drop=FALSE])), cols)         # detection
  fl <- fit_one(log2(mae[keep,,drop=FALSE] + 1), cols)         # level
  mae_max <- tapply(seq_len(ncol(mae)), cd$condition[cols],
                    function(j) rowMeans(mae[keep, j, drop=FALSE]))
  mae_max <- do.call(pmax, mae_max)                            # max per-arm mean-among-expr
  for (cn in names(CONTR)) {
    td <- topTable(fd, coef=cn, number=Inf, sort.by="none")
    tl <- topTable(fl, coef=cn, number=Inf, sort.by="none")
    res[[paste(ct,cn)]] <- data.table(
      gene = rownames(td), cell_type = ct, contrast = cn,
      detection_log2fc = td$logFC/ln2, detection_p = td$P.Value,
      level_log2fc = tl$logFC/ln2,     level_p = tl$P.Value,
      mean_among_expr_max = round(mae_max, 2))
  }
}
dt <- rbindlist(res)
dt[, detection_padj := p.adjust(detection_p, "BH")]
dt[, level_padj     := p.adjust(level_p, "BH")]
dt[, call_class := fifelse(level_padj < 0.05, "regulation",
                   fifelse(detection_padj < 0.05 & mean_among_expr_max >= LEVEL_FLOOR,
                           "fraction_shift", "ambiguous"))]
fwrite(dt, out, sep = "\t")
cat(sprintf("detection_test: %d rows | %d regulation / %d fraction_shift / %d ambiguous\n",
            nrow(dt), dt[call_class=="regulation",.N], dt[call_class=="fraction_shift",.N],
            dt[call_class=="ambiguous",.N]))
```

- [ ] **Step 2: Add the Snakemake rule** after `deg_pseudobulk` (≈line 428 of `workflows/aggregate_differential.smk`):

```python
rule detection_test:
    input:
        script = "scripts/aggregate/detection_test.R",
        se     = f"{AGG}/pseudobulk_se.rds",
    output:
        det = "results/aggregate/detection_test_m02day2.tsv",
    threads: 1
    log: "logs/aggregate/detection_test.log",
    shell:
        "{RSCRIPT} {input.script} {input.se} {output.det} > {log} 2>&1"
```

- [ ] **Step 3: Run and verify the p21 classification.**

Run:
```bash
TMPDIR=/mnt/data/projects/spatial-rads/tmp conda run -n basecamp snakemake -s workflows/aggregate_differential.smk \
  results/aggregate/detection_test_m02day2.tsv --cores 4
awk -F'\t' 'NR==1||($1=="Cdkn1a"&&$3=="SBRT_vs_Ctrl")' results/aggregate/detection_test_m02day2.tsv
```
Expected: the `Cdkn1a` / `SBRT_vs_Ctrl` rows (Fibroblast, Macrophages, etc.) show `call_class = fraction_shift` (detection sig, level not, `mean_among_expr_max ≥ 3`).

- [ ] **Step 4: Commit.**
```bash
repo-commit "feat(aggregate): detection-vs-level moderated test + sparsity-aware call_class" \
  scripts/aggregate/detection_test.R workflows/aggregate_differential.smk
```

---

## Task 3: Carry detection columns + class into the master table

**Files:**
- Modify: `scripts/aggregate/assemble_results.R` (join `detection_test_m02day2.tsv` onto the DE rows)
- Modify: `workflows/aggregate_differential.smk` (`assemble_results` gains `detect` input)

**Interfaces:**
- Consumes: `detection_test_m02day2.tsv` (Task 2). Joins on `gene==feature & cell_type==unit & contrast`.
- Produces: `results_master.tsv` DE rows gain `pct_expr_*`-derived `detection_padj`, `level_padj`, `mean_among_expr_max`, `call_class`.

- [ ] **Step 1: Add the detection input to the rule.** In `assemble_results` (`workflows/aggregate_differential.smk:613`) add under `input:`:
```python
        det     = "results/aggregate/detection_test_m02day2.tsv",
```
and pass it as the last positional arg before `{output.master}` in `shell:`.

- [ ] **Step 2: Join in `assemble_results.R`.** Where the DE (`readout_class=="DE"`) rows are assembled, read `det <- fread(det_path)` and left-join its `detection_padj,level_padj,mean_among_expr_max,call_class` onto the master DE rows keyed by `(feature,unit,contrast)`. Non-DE rows get `NA` for these columns.

```r
det <- fread(det_p)   # det_p = the new positional arg
de_idx <- master$readout_class == "DE"
mk <- det[match(paste(master$feature, master$unit, master$contrast),
                paste(det$gene, det$cell_type, det$contrast))]
for (c in c("detection_padj","level_padj","mean_among_expr_max","call_class"))
  master[[c]] <- ifelse(de_idx, mk[[c]], NA)
```

- [ ] **Step 3: Verify columns + p21 row.**
Run (after a local assemble dry-run on existing inputs):
```bash
TMPDIR=/mnt/data/projects/spatial-rads/tmp conda run -n basecamp snakemake -s workflows/aggregate_differential.smk \
  results/aggregate/results_master.tsv --cores 2 --force
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i} $h["feature"]=="Cdkn1a"&&$h["contrast"]=="SBRT_vs_Ctrl"{print $h["unit"],$h["call_class"]}' results/aggregate/results_master.tsv
```
Expected: stromal `Cdkn1a` rows print `fraction_shift`.

- [ ] **Step 4: Commit.**
```bash
repo-commit "feat(aggregate): carry detection_padj/level_padj/call_class onto DE rows in results_master" \
  scripts/aggregate/assemble_results.R workflows/aggregate_differential.smk
```

---

## Task 4: Unassigned-stroma handling (sensitivity + non-interpretable tag)

**Files:**
- Modify: `scripts/aggregate/composition.R` (unassigned-fraction test + with/without sensitivity output)
- Modify: `scripts/aggregate/assemble_results.R` (flag `cell_type=="unassigned"` DE rows non-interpretable)
- Modify: `workflows/aggregate_differential.smk` (`composition` gains a sensitivity output)

**Interfaces:**
- Produces: `results/aggregate/composition_unassigned_sensitivity.tsv` — each labelled cell type's `log2FC_logit`+`padj` under (a) unassigned retained, (b) unassigned excluded, plus a `flip` flag; and a dedicated test row for the `unassigned` proportion itself across the three contrasts.
- Produces: `results_master.tsv` gains `interpretable` (FALSE for `unit=="unassigned"` DE rows).

- [ ] **Step 1: Add the sensitivity pass in `composition.R`.** After the existing propeller test (line 104), refit with `unassigned` dropped from `m2`, and emit the side-by-side table. The `unassigned`-proportion test is already a row of the existing `test` table (it is a cell type); extract and label it.

```r
# --- unassigned sensitivity: rerun propeller with unassigned excluded ---
m2x <- m2[cell_type != "unassigned"]
px  <- getTransformedProps(clusters = m2x$cell_type, sample = m2x$sample_id, transform = "logit")
tpx <- px$TransformedProps; tpx <- tpx[apply(tpx,1,function(r) all(is.finite(r))),,drop=FALSE]
sx  <- unique(m2x[, .(sample_id, condition, slide_id)]); setkey(sx, sample_id); sx <- sx[colnames(tpx)]
dx  <- model.matrix(~ 0 + condition + slide_id, data = sx); colnames(dx) <- make.names(colnames(dx))
fx  <- eBayes(contrasts.fit(lmFit(tpx, dx), makeContrasts(
         MBRT_vs_Ctrl = conditionMBRT_day2 - conditionControl,
         SBRT_vs_Ctrl = conditionSBRT_day2 - conditionControl,
         MBRT_vs_SBRT = conditionMBRT_day2 - conditionSBRT_day2, levels = dx)), robust = TRUE)
excl <- rbindlist(lapply(colnames(cm), function(cn){
  tt <- topTable(fx, coef=cn, number=Inf, sort.by="none")
  data.table(cell_type=rownames(tt), contrast=cn, log2FC_excl=tt$logFC/log(2), p_excl=tt$P.Value)}))
excl[, padj_excl := p.adjust(p_excl, "BH")]
sens <- merge(test[, .(cell_type, contrast, log2FC_incl=log2FC_logit, padj_incl=padj)],
              excl[, .(cell_type, contrast, log2FC_excl, padj_excl)], by=c("cell_type","contrast"))
sens[, flip := (padj_incl < 0.05) != (padj_excl < 0.05) | sign(log2FC_incl) != sign(log2FC_excl)]
fwrite(sens, out_sens, sep = "\t")   # out_sens = new positional arg
```
Add `out_sens` to the script's `args` and to the rule's `output`/`shell`.

- [ ] **Step 2: Tag non-interpretable DE rows in `assemble_results.R`.**
```r
master[, interpretable := !(readout_class == "DE" & unit == "unassigned")]
```

- [ ] **Step 3: Verify.**
Run:
```bash
TMPDIR=/mnt/data/projects/spatial-rads/tmp conda run -n basecamp snakemake -s workflows/aggregate_differential.smk \
  results/aggregate/composition_unassigned_sensitivity.tsv --cores 2
head -1 results/aggregate/composition_unassigned_sensitivity.tsv; grep -c TRUE results/aggregate/composition_unassigned_sensitivity.tsv || true
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i} $h["unit"]=="unassigned"&&$h["readout_class"]=="DE"{print $h["interpretable"]}' results/aggregate/results_master.tsv | sort -u
```
Expected: sensitivity table has both inclusion columns + `flip`; unassigned DE rows print `FALSE` for `interpretable`.

- [ ] **Step 4: Commit.**
```bash
repo-commit "feat(aggregate): unassigned-stroma propeller sensitivity + non-interpretable DE tag" \
  scripts/aggregate/composition.R scripts/aggregate/assemble_results.R workflows/aggregate_differential.smk
```

---

## Task 5: Pre-registered fibroblast sub-state split

**Files:**
- Create: `config/substate_markers.yaml`
- Create: `scripts/aggregate/substate_split.R`
- Modify: `workflows/aggregate_differential.smk` (new `substate_split` rule)

**Interfaces:**
- Consumes: `merged.rds` (counts), `full_labels.parquet` (to select Fibroblast cells).
- Produces: `results/aggregate/fibroblast_substate.parquet` — `cell, substate ∈ {resting, activated}` for the 732k Fibroblast cells; and `results/aggregate/substate_gate_report.tsv` documenting the detectability + specificity-anchor checks.

- [ ] **Step 1: Write the frozen marker config + commit it first** (pre-registration; commit this step alone before inspecting any sub-state DE).
```yaml
# config/substate_markers.yaml  (frozen 2026-06-23, pre-registered before sub-state DE)
fibroblast:
  resting:    [Pdgfra, Dpt, Col14a1, Gsn, Cd34]
  activated:  [Acta2, Tagln, Myh11, Col1a1, Col3a1]   # shared/contractile core
  specificity_anchor: [Postn, Fap, Col5a1]             # fibroblast-specific, ABSENT from SmoothMuscle panel
  gate: { min_signal_to_bg: 2.0, require_specificity_anchor: true }
```

- [ ] **Step 2: Write `substate_split.R`.** UCell-score the resting and activated panels on the Fibroblast subset; assign by argmax+margin; **gate**: an "activated" call is only kept if (a) the activated panel's signal-to-background ≥2 AND (b) at least one specificity anchor (`Postn/Fap/Col5a1`) individually clears signal-to-background ≥2 in that cell's neighbourhood/cluster — otherwise the cell is left `resting` (prevents SMC bleed-through). Reuse the project's detectability helper (`tier2_detectability.py` logic / `getMeanSignalRatio`) for the signal-to-background ratio; compute per fibroblast subcluster.

```r
#!/usr/bin/env Rscript
# Fibroblast resting/activated sub-state split (Fix 3a). Pre-registered markers in
# config/substate_markers.yaml; specificity-anchor gate blocks SmoothMuscle bleed.
# Args: <merged.rds> <full_labels.parquet> <substate_markers.yaml> <out.parquet> <out_gate.tsv>
suppressPackageStartupMessages({ library(Seurat); library(UCell); library(arrow); library(data.table); library(yaml) })
a <- commandArgs(trailingOnly = TRUE)
o <- readRDS(a[1]); lab <- as.data.table(read_parquet(a[2])); mk <- read_yaml(a[3])$fibroblast
fib <- lab[cell_subtype == "Fibroblast", cell]
o <- subset(o, cells = intersect(fib, colnames(o)))
sets <- list(resting = mk$resting, activated = c(mk$activated, mk$specificity_anchor))
o <- AddModuleScore_UCell(o, features = sets, assay = "RNA")
md <- o@meta.data
# signal-to-background (mean expr in-set / mean expr panel) per anchor, per cell:
data <- LayerData(o, assay="RNA", layer="data"); bg <- Matrix::colMeans(data)
s2b <- function(g) (Matrix::colMeans(data[intersect(g, rownames(data)),,drop=FALSE]) ) / pmax(bg, 1e-6)
anchor_ok <- Reduce(`|`, lapply(mk$specificity_anchor, function(g) s2b(g) >= mk$gate$min_signal_to_bg))
act_ok    <- s2b(mk$activated) >= mk$gate$min_signal_to_bg
substate  <- ifelse(md$activated_UCell > md$resting_UCell & act_ok & anchor_ok, "activated", "resting")
out <- data.table(cell = rownames(md), substate = substate)
write_parquet(out, a[4])
gate <- data.table(n_fib = nrow(out), n_activated = sum(substate=="activated"),
                   n_activated_anchor_blocked = sum(md$activated_UCell > md$resting_UCell & !(act_ok & anchor_ok)),
                   anchors = paste(mk$specificity_anchor, collapse=","))
fwrite(gate, a[5], sep = "\t")
cat(sprintf("substate: %d fibroblasts -> %d activated / %d resting (%d blocked by specificity gate)\n",
            nrow(out), gate$n_activated, nrow(out)-gate$n_activated, gate$n_activated_anchor_blocked))
```

- [ ] **Step 3: Add the rule + run + verify the gate fired.**
```python
rule substate_split:
    input:
        script  = "scripts/aggregate/substate_split.R",
        rds     = f"{AGG}/merged.rds",
        labels  = "results/aggregate/full_labels.parquet",
        markers = "config/substate_markers.yaml",
    output:
        parquet = "results/aggregate/fibroblast_substate.parquet",
        gate    = "results/aggregate/substate_gate_report.tsv",
    threads: 4
    log: "logs/aggregate/substate_split.log",
    shell:
        "{RSCRIPT} {input.script} {input.rds} {input.labels} {input.markers} "
        "{output.parquet} {output.gate} > {log} 2>&1"
```
Run the rule; `cat results/aggregate/substate_gate_report.tsv`. Expected: nonzero `n_activated` AND `n_activated_anchor_blocked > 0` (the specificity gate is actually rejecting SMC-bleed cells, not a no-op).

- [ ] **Step 4: Commit** (markers config already committed in Step 1; commit the script + rule).
```bash
repo-commit "feat(aggregate): pre-registered fibroblast resting/activated split with SMC-specificity gate" \
  scripts/aggregate/substate_split.R workflows/aggregate_differential.smk
```

---

## Task 6: Sub-state composition adjudication

**Files:**
- Modify: `scripts/aggregate/composition.R` (run propeller with the fibroblast label split into resting/activated)
- Modify: `workflows/aggregate_differential.smk` (`composition` gains `substate` input + a sub-state test output)

**Interfaces:**
- Consumes: `fibroblast_substate.parquet` (Task 5).
- Produces: `results/aggregate/composition_substate_test_m02day2.tsv` — propeller on the label set where `Fibroblast` is replaced by `Fibroblast_resting`/`Fibroblast_activated`; plus a one-line readout cross-tabulating the activated-fibroblast proportion shift in `SBRT_vs_Ctrl` against the sign of the collagen/`Acta2` Fibroblast DE.

- [ ] **Step 1: Relabel + rerun propeller.** In `composition.R`, after loading labels, left-join `fibroblast_substate.parquet`; set `cell_type := fifelse(cell_type=="Fibroblast", paste0("Fibroblast_", substate), cell_type)`. Run the existing propeller block on this relabelled `m2` into `out_substate`. (Reuses the Step-1 machinery of `composition.R`; only the cluster vector changes.)

- [ ] **Step 2: Add the cross-tab line.**
```r
act <- substate_test[cell_type=="Fibroblast_activated" & contrast=="SBRT_vs_Ctrl"]
cat(sprintf("sub-state adjudication: activated-fibroblast SBRT_vs_Ctrl log2FC=%.2f padj=%.3g\n",
            act$log2FC_logit, act$padj))
```

- [ ] **Step 3: Run + verify.**
```bash
TMPDIR=/mnt/data/projects/spatial-rads/tmp conda run -n basecamp snakemake -s workflows/aggregate_differential.smk \
  results/aggregate/composition_substate_test_m02day2.tsv --cores 2
grep -E "Fibroblast_(resting|activated)" results/aggregate/composition_substate_test_m02day2.tsv | head
```
Expected: both `Fibroblast_resting` and `Fibroblast_activated` appear with propeller stats; the activated row in `SBRT_vs_Ctrl` is the adjudication readout (expectation: activated proportion up in SBRT → the collagen/`Acta2` Fibroblast DE is partly composition).

- [ ] **Step 4: Commit.**
```bash
repo-commit "feat(aggregate): sub-state propeller adjudicates SBRT fibrosis as composition vs regulation" \
  scripts/aggregate/composition.R workflows/aggregate_differential.smk
```

---

## Task 7: Gene-set overlap diagnostic + STING demotion

**Files:**
- Create: `scripts/aggregate/geneset_overlap.R`
- Modify: `scripts/aggregate/assemble_results.R` (flag STING rows non-independent)
- Modify: `workflows/aggregate_differential.smk` (new `geneset_overlap` rule; add as `assemble_results` input)

**Interfaces:**
- Consumes: `results/data_model/pathway_sets.tsv` (`set, tier, source, gene`).
- Produces: `results/aggregate/geneset_overlap.tsv` — pairwise Jaccard + shared-gene list across the curated (`tier=="primary"`) sets.
- Produces: `results_master.tsv` pathway rows gain `independent` (FALSE for STING, whose on-panel genes are ≥50% shared with the IFN sets).

- [ ] **Step 1: Write the overlap script.**
```r
#!/usr/bin/env Rscript
# Pairwise Jaccard + shared genes across curated gene sets (Fix 5). Args: <sets.tsv> <out.tsv>
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly=TRUE); gs <- fread(a[1])[tier=="primary"]
sets <- split(gs$gene, gs$set)
pairs <- CJ(a=names(sets), b=names(sets))[a < b]
pairs[, `:=`(
  shared = mapply(function(x,y) length(intersect(sets[[x]], sets[[y]])), a, b),
  jaccard = mapply(function(x,y) length(intersect(sets[[x]],sets[[y]]))/length(union(sets[[x]],sets[[y]])), a, b),
  shared_genes = mapply(function(x,y) paste(intersect(sets[[x]],sets[[y]]), collapse=";"), a, b))]
fwrite(pairs[order(-jaccard)], a[2], sep="\t")
cat(sprintf("geneset_overlap: %d set pairs; max jaccard %.2f\n", nrow(pairs), max(pairs$jaccard)))
```

- [ ] **Step 2: Flag STING in `assemble_results.R`.**
```r
master[, independent := !(readout_class %in% c("gsea","ucell") & grepl("STING", feature, ignore.case=TRUE))]
```

- [ ] **Step 3: Run + verify STING∩IFN.**
```bash
TMPDIR=/mnt/data/projects/spatial-rads/tmp conda run -n basecamp snakemake -s workflows/aggregate_differential.smk \
  results/aggregate/geneset_overlap.tsv --cores 1
grep -i sting results/aggregate/geneset_overlap.tsv
```
Expected: STING vs an IFN set shows `shared ≥ 2` and the shared genes `Isg15`/`Cxcl10`.

- [ ] **Step 4: Commit.**
```bash
repo-commit "feat(aggregate): gene-set overlap diagnostic + demote STING from independent evidence" \
  scripts/aggregate/geneset_overlap.R scripts/aggregate/assemble_results.R workflows/aggregate_differential.smk
```

---

## Task 8: Claim-scoping flags (MBRT dose; M01 reframe)

**Files:**
- Modify: `scripts/aggregate/assemble_results.R` (contrast-level dose-confound flag)

**Interfaces:**
- Produces: `results_master.tsv` gains `dose_confounded` (TRUE for every `MBRT_vs_SBRT` row) and a `cross_cohort_status` note column already exists or is added to mark the M01↔M02 relationship as disagreement, not validation.

- [ ] **Step 1: Add the flags.**
```r
master[, dose_confounded := contrast == "MBRT_vs_SBRT"]
```
- [ ] **Step 2: M01 reframe** — in the concordance consumer (`concordance_m01_m02.R` / wherever the cross-cohort row is written) relabel the readout from "validation/concordance" to `cross_cohort_disagreement`, and ensure the master's hypothesis/notes field does not assert M01 corroborates M02 on the SBRT stromal signal. (Verify by reading `scripts/aggregate/concordance_m01_m02.R` first; change only the label string + any "validation" wording.)

- [ ] **Step 3: Verify.**
```bash
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i} $h["contrast"]=="MBRT_vs_SBRT"{print $h["dose_confounded"]; exit}' results/aggregate/results_master.tsv
grep -ri "validation" scripts/aggregate/concordance_m01_m02.R || echo "no 'validation' wording remains"
```
Expected: `TRUE`; and no "validation" framing remains in the concordance script.

- [ ] **Step 4: Commit.**
```bash
repo-commit "feat(aggregate): dose-confound flag on MBRT_vs_SBRT + reframe M01 cross-cohort as disagreement" \
  scripts/aggregate/assemble_results.R scripts/aggregate/concordance_m01_m02.R
```

---

## Task 9: Split UCell compute from plotting (cache against plot failures)

**Files:**
- Create: `scripts/aggregate/pathway_scores.R` (UCell + AddModuleScore + the limma test + cache)
- Create: `scripts/aggregate/pathway_plots.R` (read cache → the 3 plots)
- Modify: `workflows/aggregate_differential.smk` (split `pathway_summary` into `pathway_scores` + `pathway_plots`)

**Interfaces:**
- `pathway_scores.R` produces the existing tables (`pathway_scores_summary.tsv`, `pathway_test_m02day2.tsv`, `pathway_ucell_ams_concordance.tsv`) — the ~5h compute. `pathway_plots.R` consumes those tsvs and produces the three PNGs. A plot-layer failure can no longer roll back the compute.

- [ ] **Step 1: Carve `pathway_summary.R` at the plotting boundary.** Move the compute + the three `fwrite` table outputs into `pathway_scores.R` (args: `<rds> <labels> <sets> <summary> <test> <conc>`). Move the three `ggsave` blocks (incl. the already-fixed `tier=="primary"` heatmap, line 190) into `pathway_plots.R` (args: `<summary> <test> <conc> <heatmap> <timecourse> <scatter>`), reading the tsvs back with `fread`.

- [ ] **Step 2: Replace the rule with two rules.**
```python
rule pathway_scores:
    input:
        script = "scripts/aggregate/pathway_scores.R",
        rds    = f"{AGG}/merged.rds",
        labels = "results/aggregate/full_labels.parquet",
        sets   = "results/data_model/pathway_sets.tsv",
    output:
        summary = "results/aggregate/pathway_scores_summary.tsv",
        test    = "results/aggregate/pathway_test_m02day2.tsv",
        conc    = "results/aggregate/pathway_ucell_ams_concordance.tsv",
    threads: 4
    log: "logs/aggregate/pathway_scores.log",
    shell:
        "{RSCRIPT} {input.script} {input.rds} {input.labels} {input.sets} "
        "{output.summary} {output.test} {output.conc} > {log} 2>&1"

rule pathway_plots:
    input:
        script  = "scripts/aggregate/pathway_plots.R",
        summary = "results/aggregate/pathway_scores_summary.tsv",
        test    = "results/aggregate/pathway_test_m02day2.tsv",
        conc    = "results/aggregate/pathway_ucell_ams_concordance.tsv",
    output:
        heatmap    = "results/aggregate/plots/pathway_heatmap_m02day2.png",
        timecourse = "results/aggregate/plots/pathway_timecourse_m01.png",
        scatter    = "results/aggregate/plots/pathway_ucell_vs_ams_scatter.png",
    threads: 1
    log: "logs/aggregate/pathway_plots.log",
    shell:
        "{RSCRIPT} {input.script} {input.summary} {input.test} {input.conc} "
        "{output.heatmap} {output.timecourse} {output.scatter} > {log} 2>&1"
```

- [ ] **Step 3: Verify the cache survives a plot failure.** Build `pathway_scores` outputs, then temporarily break `pathway_plots.R` (e.g. reference a missing column), run `pathway_plots`, confirm it fails but the three score tsvs are untouched (mtimes unchanged), then revert the break.
```bash
TMPDIR=/mnt/data/projects/spatial-rads/tmp conda run -n basecamp snakemake -s workflows/aggregate_differential.smk \
  results/aggregate/pathway_scores_summary.tsv --cores 4
ls -l --time-style=+%s results/aggregate/pathway_scores_summary.tsv   # note mtime; must not change on plot failure
```
Expected: scores tsvs exist; after a deliberate plot break the compute tsvs keep their mtime (no recompute).

- [ ] **Step 4: Commit.**
```bash
repo-commit "refactor(aggregate): split pathway UCell compute from plotting so a plot bug can't burn the 5h recompute" \
  scripts/aggregate/pathway_scores.R scripts/aggregate/pathway_plots.R workflows/aggregate_differential.smk
```

---

## Task 10: Optional — Milo cluster-free DA on a 200k subsample

**Files:**
- Create: `scripts/aggregate/milo_da.py` (env `spatial-rads-scvi`, `milopy`)
- Modify: `workflows/aggregate_differential.smk` (rule **not** in `rule all`; run on demand)

**Interfaces:**
- Consumes: the scVI-latent AnnData (`{FULL}/cluster_checkpoint.h5ad` or the latent in `obs.parquet` + embedding) — verify the latent's location first.
- Produces: `results/aggregate/milo_da_m02day2.tsv` + a beeswarm PNG; design `~ condition` (NOT `~ slide_id + condition`).

- [ ] **Step 1: Confirm `milopy` is installed**, else add to the env:
```bash
conda run -n spatial-rads-scvi python -c "import milopy; print(milopy.__version__)" || \
  conda run -n spatial-rads-scvi pip install milopy
```
- [ ] **Step 2: Write `milo_da.py`** — load the latent AnnData, subset to Mutter_02 day-2, **stratified-subsample ~200k cells proportional by `cell_subtype`**, `milopy.core.make_nhoods` on the scVI latent, `count_nhoods`, `DA_nhoods(design="~condition")`; write the neighbourhood DA table + report the SpatialFDR distribution as the calibration check.
- [ ] **Step 3: Run + verify** the SpatialFDR is finite and the table is nonempty; this is a feasibility/calibration gate, not a `rule all` dependency.
- [ ] **Step 4: Commit** as optional tooling.
```bash
repo-commit "feat(aggregate): optional Milo cluster-free DA (200k subsample, ~condition) as label-free cross-check" \
  scripts/aggregate/milo_da.py workflows/aggregate_differential.smk
```

> Skip this task entirely if Task 6's propeller adjudication is unambiguous. It is the deferred cross-check, not on the critical path.

---

## Task 11: Consolidated rerun + acceptance validation

**Files:** none (orchestration only)

- [ ] **Step 1: Dry-run the full differential DAG.**
```bash
cd /home/jeszyman/repos/spatial-rads
TMPDIR=/mnt/data/projects/spatial-rads/tmp conda run -n basecamp snakemake -s workflows/aggregate_differential.smk \
  results/aggregate/results_master.tsv --dry-run
```
Expected: the DAG includes `pseudobulk_build → {deg_pseudobulk, detection_test}`, `composition` (+ sensitivity + substate), `pathway_scores → pathway_plots`, `geneset_overlap`, and `assemble_results`; no orphaned/typo paths.

- [ ] **Step 2: Execute the rerun.** (~5h pathway_scores dominates; run with adequate cores, BLAS pinned.)
```bash
TMPDIR=/mnt/data/projects/spatial-rads/tmp conda run -n basecamp snakemake -s workflows/aggregate_differential.smk \
  results/aggregate/results_master.tsv --cores 8 2>&1 | tee logs/aggregate/rerun_$(date +%Y%m%d).log
```
Run in the background and wait on completion (no polling); read the log on the completion notification.

- [ ] **Step 3: Acceptance checks (all must pass).**
```bash
M=results/aggregate/results_master.tsv
# (a) new DE columns present
head -1 $M | tr '\t' '\n' | grep -E "call_class|detection_padj|level_padj|interpretable|dose_confounded|independent"
# (b) p21 SBRT stromal rows classify fraction_shift
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i} $h["feature"]=="Cdkn1a"&&$h["contrast"]=="SBRT_vs_Ctrl"{print $h["unit"],$h["call_class"]}' $M
# (c) unassigned DE rows non-interpretable
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i} $h["unit"]=="unassigned"&&$h["readout_class"]=="DE"{print $h["interpretable"]}' $M | sort -u
# (d) sidecar tables exist
ls -la results/aggregate/{detection_test_m02day2,composition_unassigned_sensitivity,composition_substate_test_m02day2,geneset_overlap,substate_gate_report}.tsv
```
Expected: (a) all six columns listed; (b) stromal `Cdkn1a` rows = `fraction_shift`; (c) only `FALSE`; (d) all five sidecars present.

- [ ] **Step 4: Final commit** (results are git-ignored data; commit only any residual code/doc + mark the spec done).
```bash
repo-commit "chore(aggregate): consolidated differential rerun — robustness fixes live; results_master refreshed" \
  plan-differential-robustness.md plan-differential-robustness-impl.md
```

---

## Self-review (spec coverage)

- Fix 1 (unassigned) → Task 4. Fix 2 (detection) → Tasks 1–3. Fix 3 (state mixing) → Tasks 5, 6, optional 10. Fix 4 (claim scoping) → Task 8. Fix 5 (STING) → Task 7. Plumbing (pathway fix already done; UCell cache; one rerun) → Tasks 9, 11. All six fixes + plumbing covered.
- No-beta-binomial / no-naive-MAST constraint honored (Task 2 uses limma-eBayes). Pre-registration honored (Task 5 Step 1 commits markers before any sub-state DE). Specificity-anchor + detectability gates both present (Task 5). Milo `~condition` only, subsampled (Task 10).
- Column-name consistency: `call_class`, `detection_padj`, `level_padj`, `mean_among_expr_max`, `interpretable`, `dose_confounded`, `independent` used identically across Tasks 2/3/4/7/8/11.
