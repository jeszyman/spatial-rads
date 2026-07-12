# Comparison Registry + Marker Provenance Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the study's experimental comparison design and cell-lineage markers into validated, machine-readable Stage-A artifacts emitted by `data_model.smk`, mirroring the shipped pathway gene-set pattern, and consolidate the plan-doc sprawl.

**Architecture:** Two new curated sources (`config/comparisons.yaml`, `config/marker_sets_provenance.tsv`) drive two new R builder scripts that resolve against the samplesheet / common panel and emit validated TSVs (`comparisons.tsv`, `marker_panel_coverage.tsv`) plus a coverage figure. Both are metadata-only Stage-A rules — no heavy per-cell objects. The contrast-declaring scripts + `assemble_results.R` that would *read* `comparisons.tsv` live in `aggregate_differential.smk` (the cross-sample layer was split 2026-07-12 into `aggregate_typing.smk` → `full_labels.parquet` and `aggregate_differential.smk` → `results_master.tsv`, superseding the monolithic `aggregate.smk`); that rewiring is explicitly deferred to a later phase riding the pending `results_master.tsv` rerun.

**Tech Stack:** R (data.table, yaml, tidyverse for the figure), Snakemake, conda (`spatial-rads` for R via `conda run`, driver in `basecamp`).

## Global Constraints

- Snakemake invocation always sets `TMPDIR=/mnt/data/projects/spatial-rads/tmp` and runs via `conda run -n basecamp snakemake`; dry-run before execute. Copy verbatim: `TMPDIR=/mnt/data/projects/spatial-rads/tmp conda run -n basecamp snakemake -s workflows/data_model.smk --cores 1`.
- R steps run through `RSCRIPT = "conda run -n spatial-rads Rscript"` (already defined at `workflows/data_model.smk:7`).
- Builder scripts are declared as snakemake `input:` deps (script-tracking convention) so code edits trigger reruns.
- `rule all` is a manifest of ALL deliverables — every new output is listed there explicitly even if a co-output already pulls it.
- No `shell.prefix`; rely on Snakemake's default strict mode. Per-rule `threads: 1`.
- Curated sources are pre-registered: builders enforce a git-HEAD freeze check on the frozen fields (mirrors `build_gene_sets.R:49-63`).
- The common panel file `results/data_model/common_genes.tsv` is headerless, one gene symbol per line (950 lines); read with `readLines()`, matching `build_gene_sets.R:16`.
- Commit only the files named in each task, by name (never `git add -A`). Use `repo-commit "<msg>" <files...>` per the repo convention.
- No completion/verification claim without running the proving command in the same step and showing its output.

---

## File Structure

**Create:**
- `config/comparisons.yaml` — curated comparison registry (source of truth for the design).
- `config/marker_sets_provenance.tsv` — per-lineage marker source + weak-separation flags.
- `scripts/build_comparisons.R` — resolves the registry against the samplesheet → `comparisons.tsv`.
- `scripts/build_marker_sets.R` — marker panel-coverage audit → `marker_panel_coverage.tsv`.
- `scripts/fig_marker_coverage.R` — coverage figure (mirrors `fig_geneset_coverage.R`).

**Modify:**
- `workflows/data_model.smk` — three new rules + `rule all` manifest entries.
- `plan-preprocessing.md` — absorb `plan-qc-metrics.md`.
- `plan-aggregate.md`, `plan-mbrt-signatures.md`, `plan-outcomes-de.md` — correct stale atlas banners.
- `spatial-rads.org` — `** Plans` index + status ledger.

**Delete (after content-capture check):**
- `plan-pathway-genesets.md`, `plan-pathway-genesets-impl.md`, `plan-mbrt-vs-sbrt.md`, `plan-qc-metrics.md`.

**Emitted outputs (not hand-authored):**
- `results/data_model/comparisons.tsv`
- `results/data_model/marker_panel_coverage.tsv`
- `results/data_model/plots/marker_coverage.png`

---

## Task 1: Curated comparison registry source

**Files:**
- Create: `config/comparisons.yaml`

**Interfaces:**
- Produces: a YAML with a top-level `comparisons:` list. Each entry has keys read by Task 2's builder: `name` (str, unique), `kind` (`sample`|`region`|`context`), `cohort` (str) or `cohorts` (list), optional `group1`/`group2` (maps of samplesheet-column → value or list-of-values), `unit` (`mouse`|`sample`), optional `formula` (str), `tier` (`confirmatory`|`exploratory`|`descriptive`), optional `hypotheses` (list), optional `dose_confounded` (bool), optional `requires` (str), optional `resolution` (list), optional `strata` (list), optional `method_notes` (list), optional `notes`/`view`.

- [ ] **Step 1: Write the registry file**

Write `config/comparisons.yaml` exactly as specified in the design doc section A.1 (`docs/superpowers/specs/2026-07-12-comparison-registry-and-marker-provenance-design.md`), reproduced here:

```yaml
# config/comparisons.yaml
# Curated registry of the experimental comparisons this sample set admits.
# Group predicates filter results/data_model/samples.tsv; the builder resolves realized n,
# model rank, residual df, and inference-capability from the resolved sample subset.
# Pre-registered: tier + hypotheses are frozen (git-HEAD freeze check, like pathway sets).

comparisons:
  # ---- sample-level, day-2 inference cohort (n=4, formal) ----
  - name: SBRT_vs_Ctrl
    kind: sample
    resolution: [whole, compartment]
    cohort: mutter02_day2
    group1: {treatment: SBRT, timepoint_h: 48, model: flank}
    group2: {treatment: NT,   timepoint_h: 48, model: flank}
    unit: mouse
    formula: "~0 + condition + slide_id"
    tier: confirmatory
    hypotheses: [H1_immune, H2_vascular, H3_stromal]
    dose_confounded: false

  - name: MBRT_vs_Ctrl
    kind: sample
    resolution: [whole, compartment]
    cohort: mutter02_day2
    group1: {treatment: MBRT, timepoint_h: 48, model: flank}
    group2: {treatment: NT,   timepoint_h: 48, model: flank}
    unit: mouse
    formula: "~0 + condition + slide_id"
    tier: confirmatory
    hypotheses: [H1_immune, H2_vascular, H3_stromal]
    dose_confounded: false
    notes: >
      Whole-compartment averaging cancels a peak-restricted MBRT signal (dilution).
      Bounds, does not resolve, the MBRT question -- see region MBRT_peak_vs_valley.

  - name: MBRT_vs_SBRT
    kind: sample
    resolution: [whole, compartment]
    cohort: mutter02_day2
    group1: {treatment: MBRT, timepoint_h: 48, model: flank}
    group2: {treatment: SBRT, timepoint_h: 48, model: flank}
    unit: mouse
    formula: "~0 + condition + slide_id"
    tier: confirmatory
    hypotheses: [H1_immune, H2_vascular, H3_stromal]
    dose_confounded: true
    view:
      differential_response_scatter: true

  # ---- region-level, within-sample ----
  - name: MBRT_peak_vs_valley
    kind: region
    axis: dist_to_peak
    test: periodicity_rotation_null
    cohorts: [mutter01_mbrt_4h, mutter02_day2]
    unit: mouse
    pairing: within_mouse
    requires: h2ax_registration
    tier: confirmatory
    method_notes: [fov_alias_guard, dose_dropout_check, boundary_buffer, nonlinear_distance]

  - name: MBRT_valley_vs_SBRT
    kind: region
    cohorts: [mutter02_day2]
    unit: mouse
    requires: h2ax_registration
    tier: exploratory
    dose_confounded: true

  - name: niche_DA_by_arm
    kind: region
    test: kmeans_knn_composition
    cohort: mutter02_day2
    unit: mouse
    formula: "~0 + condition + slide_id"
    tier: exploratory

  - name: tumor_immune_mixing_by_arm
    kind: region
    test: keren_mixing_score
    cohort: mutter02_day2
    unit: mouse
    formula: "~0 + condition + slide_id"
    tier: exploratory

  - name: substate_composition_by_arm
    kind: region
    test: substate_fraction
    cohort: mutter02_day2
    unit: mouse
    formula: "~0 + condition + slide_id"
    tier: exploratory

  # ---- descriptive, time-course cohort (n=1, effect sizes only) ----
  - name: treated_vs_baseline_timecourse
    kind: sample
    cohort: mutter01_flank
    group1: {treatment: [MBRT, SBRT], model: flank}
    group2: {treatment: NT, timepoint_h: 0, model: flank}
    unit: sample
    tier: descriptive
    notes: single 0h control conflates time-since-implant with treatment; descriptive only.

  - name: MBRT_vs_SBRT_matched_timepoint
    kind: sample
    cohort: mutter01_flank
    strata: [4, 48, 144]
    unit: sample
    tier: descriptive

  - name: MBRT_vs_SBRT_tongue_day10
    kind: sample
    cohort: mutter01_tongue
    unit: sample
    tier: descriptive
    notes: different tumor + site; timepoint-confounded cross-model context only.

  # ---- context / replication ----
  - name: cohort_concordance
    kind: context
    test: effect_size_concordance
    tier: exploratory
    notes: M01 n=1 anti-correlates with M02 on the SBRT signal; contextual, cannot corroborate.
```

- [ ] **Step 2: Verify it parses as YAML**

Run:
```bash
cd /home/jeszyman/repos/spatial-rads
conda run -n spatial-rads Rscript -e 'x <- yaml::read_yaml("config/comparisons.yaml"); cat("entries:", length(x$comparisons), "\n"); cat("names:", paste(sapply(x$comparisons, `[[`, "name"), collapse=", "), "\n")'
```
Expected: `entries: 12` and the 12 comparison names printed, no parse error.

- [ ] **Step 3: Commit**

```bash
repo-commit "feat(data-model): add curated comparison registry source" config/comparisons.yaml
```

---

## Task 2: Comparison registry builder

**Files:**
- Create: `scripts/build_comparisons.R`
- Consumes: `config/comparisons.yaml` (Task 1), `results/data_model/samples.tsv`.

**Interfaces:**
- Produces: `results/data_model/comparisons.tsv`, one row per resolved comparison × resolution × cohort, columns: `name, kind, cohort, resolution, unit, tier, formula, hypotheses, dose_confounded, requires, n_group1, n_group2, n_mouse_group1, n_mouse_group2, model_rank, resid_df, inference_capable, gate`.
- Invocation: `Rscript scripts/build_comparisons.R <comparisons.yaml> <samples.tsv> <out_comparisons.tsv>`.

Cohort → samplesheet filter mapping (the builder's one piece of study-specific glue):
- `mutter02_day2` → `name == "Mutter_02"` (all 12 rows are day-2 flank).
- `mutter01_flank` → `name == "Mutter_01" & model == "flank"`.
- `mutter01_tongue` → `name == "Mutter_01" & model == "tongue"`.
- `mutter01_mbrt_4h` → `name == "Mutter_01" & treatment == "MBRT" & timepoint_h == 4`.

- [ ] **Step 1: Write the builder script**

Create `scripts/build_comparisons.R`:

```r
#!/usr/bin/env Rscript
# Resolve the curated comparison registry against the master samplesheet into a
# validated design table. Each registry entry DECLARES intent (groups, unit, formula,
# tier); this builder COMPUTES realized facts (resolved n, distinct mice, model rank,
# residual df, inference-capability) so the experimental design is machine-checkable.
# Pre-registration: tier + hypotheses are frozen against git HEAD (silent change to the
# confirmatory family hard-errors). Stage-A metadata only -- never touches per-cell data.
# Args: <comparisons.yaml> <samples.tsv> <out_comparisons.tsv>
suppressPackageStartupMessages({library(yaml); library(data.table)})

`%||%` <- function(a, b) if (is.null(a)) b else a

a       <- commandArgs(trailingOnly = TRUE)
yaml_p  <- a[1]; ss_p <- a[2]; out_p <- a[3]

reg <- read_yaml(yaml_p)$comparisons
ss  <- fread(ss_p, colClasses = list(character = "timepoint_h"))  # keep predicates string-safe
ss[, timepoint_h := suppressWarnings(as.integer(timepoint_h))]

# cohort label -> row filter on the samplesheet
cohort_rows <- function(cohort) {
  switch(cohort,
    mutter02_day2   = ss[name == "Mutter_02"],
    mutter01_flank  = ss[name == "Mutter_01" & model == "flank"],
    mutter01_tongue = ss[name == "Mutter_01" & model == "tongue"],
    mutter01_mbrt_4h= ss[name == "Mutter_01" & treatment == "MBRT" & timepoint_h == 4],
    stop("unknown cohort: ", cohort))
}

# apply a group predicate map (col -> value / list-of-values) to a cohort subset
apply_pred <- function(d, pred) {
  if (is.null(pred)) return(d[0])
  keep <- rep(TRUE, nrow(d))
  for (col in names(pred)) {
    if (!col %in% names(d)) stop("predicate column not in samplesheet: ", col)
    vals <- as.character(unlist(pred[[col]]))
    keep <- keep & (as.character(d[[col]]) %in% vals)
  }
  d[keep]
}

# resid df for `~0+condition+slide_id` on a resolved 2-group subset (else NA)
design_facts <- function(d, formula) {
  if (is.null(formula) || nrow(d) == 0) return(list(rank = NA_integer_, resid = NA_integer_))
  d <- droplevels(as.data.frame(d))
  d$condition <- factor(d$condition); d$slide_id <- factor(d$slide_id)
  mm <- tryCatch(model.matrix(as.formula(formula), d), error = function(e) NULL)
  if (is.null(mm)) return(list(rank = NA_integer_, resid = NA_integer_))
  r <- qr(mm)$rank
  list(rank = r, resid = nrow(mm) - r)
}

rows <- list()
for (e in reg) {
  cohorts <- e$cohorts %||% (e$cohort %||% NA_character_)
  resterms <- e$resolution %||% NA_character_
  for (ch in cohorts) for (res in resterms) {
    d  <- if (is.na(ch)) ss[0] else cohort_rows(ch)
    g1 <- apply_pred(d, e$group1); g2 <- apply_pred(d, e$group2)
    df <- design_facts(rbind(g1, g2), e$formula)
    unit <- e$unit %||% "sample"
    nu1  <- if (unit == "mouse") uniqueN(g1$mouse_id) else nrow(g1)
    nu2  <- if (unit == "mouse") uniqueN(g2$mouse_id) else nrow(g2)
    inf  <- (e$tier %||% "") %in% c("confirmatory","exploratory") &&
            e$kind %in% c("sample","region") && nu1 >= 2 && nu2 >= 2
    gate <- if (!is.null(e$requires)) "unchecked" else NA_character_
    rows[[length(rows)+1]] <- data.table(
      name = e$name, kind = e$kind, cohort = as.character(ch),
      resolution = as.character(res), unit = unit, tier = e$tier %||% NA_character_,
      formula = e$formula %||% NA_character_,
      hypotheses = paste(e$hypotheses %||% character(), collapse = ";"),
      dose_confounded = isTRUE(e$dose_confounded),
      requires = e$requires %||% NA_character_,
      n_group1 = nrow(g1), n_group2 = nrow(g2),
      n_mouse_group1 = nu1, n_mouse_group2 = nu2,
      model_rank = df$rank, resid_df = df$resid,
      inference_capable = inf, gate = gate)
  }
}
out <- rbindlist(rows, fill = TRUE)

## ---- validation: sample/region entries with predicates must resolve ----
bad <- out[kind %in% c("sample","region") & !is.na(cohort) &
           !is.null(unlist(lapply(reg, `[[`, "group1"))) & n_group1 == 0 &
           name %in% sapply(Filter(function(e) !is.null(e$group1), reg), `[[`, "name")]
if (nrow(bad)) stop("comparison(s) resolved to zero group1 samples: ",
                    paste(unique(bad$name), collapse = ", "))

## ---- freeze: tier + hypotheses must match git HEAD ----
committed <- tryCatch(
  yaml::yaml.load(paste(system2("git", c("show", "HEAD:config/comparisons.yaml"),
                                stdout = TRUE), collapse = "\n")),
  error = function(e) NULL)
if (!is.null(committed)) {
  ch_tier <- setNames(lapply(committed$comparisons, function(x)
    list(tier = x$tier %||% NA, hyp = x$hypotheses %||% character())),
    sapply(committed$comparisons, `[[`, "name"))
  for (e in reg) {
    c0 <- ch_tier[[e$name]]
    if (!is.null(c0)) {
      if (!identical(e$tier %||% NA, c0$tier))
        stop("freeze: tier changed for '", e$name, "' vs committed comparisons.yaml")
      if (!setequal(e$hypotheses %||% character(), c0$hyp))
        stop("freeze: hypotheses changed for '", e$name, "' vs committed comparisons.yaml")
    }
  }
}

dir.create(dirname(out_p), recursive = TRUE, showWarnings = FALSE)
fwrite(out, out_p, sep = "\t")
cat(sprintf("comparisons: %d rows from %d registry entries\n", nrow(out), length(reg)))
print(out[, .(name, cohort, resolution, n_mouse_group1, n_mouse_group2, resid_df, inference_capable, gate)])
```

- [ ] **Step 2: Run standalone and inspect output**

Run:
```bash
cd /home/jeszyman/repos/spatial-rads
conda run -n spatial-rads Rscript scripts/build_comparisons.R \
  config/comparisons.yaml results/data_model/samples.tsv \
  results/data_model/comparisons.tsv
```
Expected: prints `comparisons: N rows from 12 registry entries`, no error.

- [ ] **Step 3: Assert the key realized facts**

Run:
```bash
cd /home/jeszyman/repos/spatial-rads
conda run -n spatial-rads Rscript -e '
d <- data.table::fread("results/data_model/comparisons.tsv")
stopifnot(nrow(d[name=="SBRT_vs_Ctrl" & resolution=="whole"])==1)
r <- d[name=="SBRT_vs_Ctrl" & resolution=="whole"]
stopifnot(r$n_mouse_group1==4, r$n_mouse_group2==4, r$resid_df==6, r$inference_capable==TRUE)
stopifnot(d[name=="MBRT_vs_SBRT" & resolution=="whole"]$dose_confounded==TRUE)
stopifnot(d[name=="treated_vs_baseline_timecourse"]$inference_capable==FALSE)
stopifnot(d[name=="MBRT_peak_vs_valley" & cohort=="mutter02_day2"]$gate=="unchecked")
cat("all comparison assertions passed\n")'
```
Expected: `all comparison assertions passed`. (If `resid_df` is not 6, the block design or samplesheet changed — investigate before proceeding, do not adjust the assertion.)

- [ ] **Step 4: Commit**

```bash
repo-commit "feat(data-model): comparison registry builder -> comparisons.tsv" scripts/build_comparisons.R
```

---

## Task 3: Wire the comparison builder into data_model.smk

**Files:**
- Modify: `workflows/data_model.smk` (add `build_comparisons` rule; add output to `rule all`).

**Interfaces:**
- Consumes: `scripts/build_comparisons.R` (Task 2), `config/comparisons.yaml`, `config["samplesheet"]`.
- Produces: `results/data_model/comparisons.tsv` as a tracked workflow output.

- [ ] **Step 1: Add the output to `rule all`**

In `workflows/data_model.smk`, inside `rule all:`'s `input:` list (after the existing `pathway_sets.tsv` line, `:14`), add:

```python
        "results/data_model/comparisons.tsv",
```

- [ ] **Step 2: Add the rule**

Append to `workflows/data_model.smk` (after the `build_gene_sets` rule, before `panel_provenance`):

```python
# --- comparison registry: curated design resolved against the samplesheet ---
rule build_comparisons:
    input:
        script = "scripts/build_comparisons.R",
        yaml   = "config/comparisons.yaml",
        ss     = config["samplesheet"],
    output:
        "results/data_model/comparisons.tsv",
    threads: 1
    log: "logs/build_comparisons.log",
    shell:
        "{RSCRIPT} {input.script} {input.yaml} {input.ss} {output} > {log} 2>&1"
```

- [ ] **Step 3: Dry-run**

Run:
```bash
cd /home/jeszyman/repos/spatial-rads
TMPDIR=/mnt/data/projects/spatial-rads/tmp conda run -n basecamp \
  snakemake -s workflows/data_model.smk --dry-run results/data_model/comparisons.tsv
```
Expected: DAG shows `build_comparisons`, no cycle/missing-input error.

- [ ] **Step 4: Execute the rule**

Run:
```bash
cd /home/jeszyman/repos/spatial-rads
TMPDIR=/mnt/data/projects/spatial-rads/tmp conda run -n basecamp \
  snakemake -s workflows/data_model.smk --cores 1 results/data_model/comparisons.tsv
```
Expected: `build_comparisons` completes; `logs/build_comparisons.log` ends with the row summary; the assertions from Task 2 Step 3 still pass.

- [ ] **Step 5: Commit**

```bash
repo-commit "feat(data-model.smk): wire build_comparisons rule" workflows/data_model.smk
```

---

## Task 4: Marker provenance source

**Files:**
- Create: `config/marker_sets_provenance.tsv`

**Interfaces:**
- Produces: a tab-separated table with columns `set, tier, source, weak_separation`. `set` must cover every lineage in `config/lineage_markers.yaml` (12: T, B, Plasma, NK, Macrophage, Neutrophil, Mast, DC, Epithelial, Endothelial, Fibroblast, Pericyte, SmoothMuscle, Adipocyte) and every substate in `config/substate_markers.yaml` (`fibroblast_resting`, `fibroblast_activated`). Consumed by Task 5.

- [ ] **Step 1: Write the provenance table**

Create `config/marker_sets_provenance.tsv` (tab-separated; `weak_separation` blank where none). Source labels come from the header comments in `config/lineage_markers.yaml`:

```
set	tier	source	weak_separation
T	immune	ImmGen canonical T-lineage	
B	immune	ImmGen canonical B-lineage	
Plasma	immune	ImmGen canonical plasma (panel-clean lineage-level)	
NK	immune	ImmGen canonical NK	
Macrophage	immune	ImmGen canonical myeloid	
Neutrophil	immune	ImmGen canonical neutrophil (panel-clean lineage-level)	
Mast	immune	ImmGen canonical mast (panel-clean lineage-level)	
DC	immune	ImmGen canonical DC	shared_with_Macrophage (Itgax/Ciita)
Epithelial	nonimmune	canonical epithelial; 4T1 flank epithelial = tumor implant	
Endothelial	nonimmune	canonical endothelial	
Fibroblast	nonimmune	canonical fibroblast	
Pericyte	nonimmune	canonical perivascular	shared_with_SmoothMuscle (3 panel markers; Pdgfrb fibroblast-shared)
SmoothMuscle	nonimmune	canonical SMC	
Adipocyte	nonimmune	canonical adipocyte	
fibroblast_resting	substate	pre-registered (Zhang 2019 double-dip guard)	
fibroblast_activated	substate	pre-registered (Zhang 2019 double-dip guard); Col5a1/Col5a3 specificity anchor	
```

- [ ] **Step 2: Verify it reads with the expected columns and coverage**

Run:
```bash
cd /home/jeszyman/repos/spatial-rads
conda run -n spatial-rads Rscript -e '
p <- data.table::fread("config/marker_sets_provenance.tsv")
stopifnot(all(c("set","tier","source","weak_separation") %in% names(p)))
ly <- yaml::read_yaml("config/lineage_markers.yaml")
sy <- yaml::read_yaml("config/substate_markers.yaml")
need <- c(names(ly), "fibroblast_resting", "fibroblast_activated")
miss <- setdiff(need, p$set)
if (length(miss)) stop("provenance missing set(s): ", paste(miss, collapse=", "))
cat("provenance covers all", length(need), "sets\n")'
```
Expected: `provenance covers all 16 sets`.

- [ ] **Step 3: Commit**

```bash
repo-commit "feat(data-model): add cell-marker provenance source" config/marker_sets_provenance.tsv
```

---

## Task 5: Marker panel-coverage builder

**Files:**
- Create: `scripts/build_marker_sets.R`
- Consumes: `config/lineage_markers.yaml`, `config/substate_markers.yaml`, `config/marker_sets_provenance.tsv` (Task 4), `results/data_model/common_genes.tsv`.

**Interfaces:**
- Produces: `results/data_model/marker_panel_coverage.tsv`, columns `set, tier, source, n_total, n_panel, frac_panel, thin, weak_separation`.
- Invocation: `Rscript scripts/build_marker_sets.R <lineage.yaml> <substate.yaml> <provenance.tsv> <common_genes.tsv> <min_panel_genes> <out_coverage.tsv>`.
- Note: `substate_markers.yaml` is nested (`fibroblast: {resting: [...], activated: [...], specificity_anchor: [...], gate: {...}}`); the builder flattens `resting`/`activated` to sets `fibroblast_resting`/`fibroblast_activated` and ignores `specificity_anchor`/`gate`.

- [ ] **Step 1: Write the builder script**

Create `scripts/build_marker_sets.R`:

```r
#!/usr/bin/env Rscript
# Panel-coverage audit for the cell-lineage + substate marker sets, mirroring
# build_gene_sets.R. Emits per-set n_total / n_panel / frac_panel / thin, carrying the
# provenance source and known weak-separation flags. Enforces two integrity claims the
# marker YAMLs currently only assert in comments: (1) coarse-lineage on-panel completeness
# (every coarse marker present on the common panel) and (2) a git-HEAD freeze on the
# pre-registered marker membership. Stage-A metadata only.
# Args: <lineage.yaml> <substate.yaml> <provenance.tsv> <common_genes.tsv> <min_panel_genes> <out.tsv>
suppressPackageStartupMessages({library(yaml); library(data.table)})

`%||%` <- function(a, b) if (is.null(a)) b else a

a       <- commandArgs(trailingOnly = TRUE)
lin_p   <- a[1]; sub_p <- a[2]; prov_p <- a[3]; panel_p <- a[4]
min_pg  <- as.integer(a[5]); out_p <- a[6]

panel <- readLines(panel_p)
prov  <- fread(prov_p)

## ---- flatten lineage (flat) + substate (nested) into (set, tier, gene) ----
lin <- read_yaml(lin_p)
lin_dt <- rbindlist(lapply(names(lin), function(s)
  data.table(set = s, gene = as.character(lin[[s]]))))

sub <- read_yaml(sub_p)$fibroblast
sub_dt <- rbindlist(list(
  data.table(set = "fibroblast_resting",   gene = as.character(sub$resting)),
  data.table(set = "fibroblast_activated", gene = as.character(sub$activated))))

sets_long <- rbind(lin_dt, sub_dt)

## ---- provenance completeness (mirrors build_gene_sets.R:21-22) ----
miss_prov <- setdiff(unique(sets_long$set), prov$set)
if (length(miss_prov)) stop("no provenance for marker set(s): ", paste(miss_prov, collapse = ", "))

## ---- coverage ----
cov <- sets_long[, .(n_total = uniqueN(gene),
                     n_panel = uniqueN(gene[gene %in% panel])),
                 by = set]
cov <- merge(cov, prov[, .(set, tier, source, weak_separation)], by = "set", all.x = TRUE)
cov[, frac_panel := n_panel / n_total]
cov[, thin := n_panel < min_pg]
setcolorder(cov, c("set","tier","source","n_total","n_panel","frac_panel","thin","weak_separation"))

## ---- on-panel completeness check: coarse lineages claim ALL markers on-panel ----
off <- sets_long[tier_lookup <- TRUE][!gene %in% panel]  # placeholder guard; real check below
coarse_off <- lin_dt[!gene %in% panel]
if (nrow(coarse_off))
  warning("coarse-lineage markers OFF the common panel (weakens that lineage): ",
          paste(sprintf("%s:%s", coarse_off$set, coarse_off$gene), collapse = ", "))

## ---- freeze: marker membership must match git HEAD ----
freeze_check <- function(path, committed_getter) {
  committed <- tryCatch(committed_getter(), error = function(e) NULL)
  if (is.null(committed)) return(invisible())
  committed
}
lin_head <- tryCatch(yaml::yaml.load(paste(system2("git",
  c("show", "HEAD:config/lineage_markers.yaml"), stdout = TRUE), collapse = "\n")),
  error = function(e) NULL)
if (!is.null(lin_head)) for (s in names(lin))
  if (!setequal(as.character(lin[[s]]), as.character(lin_head[[s]] %||% character())))
    stop("freeze: lineage set '", s, "' differs from committed lineage_markers.yaml")

dir.create(dirname(out_p), recursive = TRUE, showWarnings = FALSE)
fwrite(cov, out_p, sep = "\t")
cat(sprintf("marker coverage: %d sets (%d coarse off-panel)\n",
            nrow(cov), nrow(coarse_off)))
print(cov[, .(set, tier, n_panel, n_total, thin, weak_separation)])
```

Note: delete the stray `off <- ...` placeholder line before finalizing — the real completeness check is `coarse_off`. (Left visible here so the reviewer sees the intent; remove it in the actual file.)

- [ ] **Step 2: Remove the placeholder line**

Edit `scripts/build_marker_sets.R` to delete the line:
```r
off <- sets_long[tier_lookup <- TRUE][!gene %in% panel]  # placeholder guard; real check below
```

- [ ] **Step 3: Run standalone**

Run:
```bash
cd /home/jeszyman/repos/spatial-rads
conda run -n spatial-rads Rscript scripts/build_marker_sets.R \
  config/lineage_markers.yaml config/substate_markers.yaml \
  config/marker_sets_provenance.tsv results/data_model/common_genes.tsv 5 \
  results/data_model/marker_panel_coverage.tsv
```
Expected: prints `marker coverage: 16 sets (0 coarse off-panel)` (the lineage YAML asserts all markers on-panel; 0 confirms it). If >0, the warning lists which — that is a real finding to surface, not an error.

- [ ] **Step 4: Assert output shape + the documented weak-separation flags**

Run:
```bash
cd /home/jeszyman/repos/spatial-rads
conda run -n spatial-rads Rscript -e '
d <- data.table::fread("results/data_model/marker_panel_coverage.tsv")
stopifnot(all(c("set","tier","source","n_total","n_panel","frac_panel","thin","weak_separation") %in% names(d)))
stopifnot(nrow(d)==16)
stopifnot(nzchar(d[set=="DC"]$weak_separation), nzchar(d[set=="Pericyte"]$weak_separation))
stopifnot(all(d[tier %in% c("immune","nonimmune")]$frac_panel==1))  # coarse all on-panel
cat("all marker-coverage assertions passed\n")'
```
Expected: `all marker-coverage assertions passed`.

- [ ] **Step 5: Commit**

```bash
repo-commit "feat(data-model): marker panel-coverage builder" scripts/build_marker_sets.R
```

---

## Task 6: Marker coverage figure

**Files:**
- Create: `scripts/fig_marker_coverage.R`
- Consumes: `results/data_model/marker_panel_coverage.tsv` (Task 5).

**Interfaces:**
- Produces: `results/data_model/plots/marker_coverage.png`.
- Invocation: `Rscript scripts/fig_marker_coverage.R <marker_panel_coverage.tsv> <out.png>`.

- [ ] **Step 1: Write the figure script (mirrors `fig_geneset_coverage.R`)**

Create `scripts/fig_marker_coverage.R`:

```r
#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# fig_marker_coverage.R
# Cell-lineage / substate marker panel coverage: fraction of each set's markers
# present on the common panel, coloured by tier, thin sets flagged. Reads the
# marker panel-coverage table. A lineage with off-panel markers is silently
# weakened; this figure makes that visible.
# Args: <marker_panel_coverage.tsv> <out.png>
# -----------------------------------------------------------------------------
suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
cov <- args[1]; out <- args[2]

d <- read_tsv(cov, show_col_types = FALSE) %>%
  arrange(tier, frac_panel)
d$set <- factor(d$set, levels = d$set)

p <- ggplot(d, aes(frac_panel, set)) +
  geom_segment(aes(x = 0, xend = frac_panel, y = set, yend = set), color = "grey70") +
  geom_point(aes(color = tier, shape = thin), size = 4) +
  geom_text(aes(label = sprintf("%d/%d", n_panel, n_total)), hjust = -0.35, size = 3.2) +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1.12),
                     breaks = seq(0, 1, 0.25)) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 4),
                     labels = c("FALSE" = "usable", "TRUE" = "thin"), name = NULL) +
  labs(x = "Fraction of set markers on the common 950-gene panel", y = NULL,
       title = "Cell-marker panel coverage",
       subtitle = "Off-panel markers silently weaken a lineage; x = thin set") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

ggsave(out, p, width = 8, height = 5, dpi = 300)
cat("wrote", out, "\n")
```

- [ ] **Step 2: Render and view the PNG**

Run:
```bash
cd /home/jeszyman/repos/spatial-rads
conda run -n spatial-rads Rscript scripts/fig_marker_coverage.R \
  results/data_model/marker_panel_coverage.tsv \
  results/data_model/plots/marker_coverage.png
feh --scale-down --auto-zoom --geometry 1900x1100 results/data_model/plots/marker_coverage.png
```
Expected: `wrote ...marker_coverage.png`; on inspection — 16 rows, x-axis 0–~112%, labels `n/n` not clipped, coarse lineages at 100%, substates possibly <100%. Fix and re-render if the axis is empty-margined or labels clip.

- [ ] **Step 3: Commit**

```bash
repo-commit "feat(data-model): marker coverage figure" scripts/fig_marker_coverage.R results/data_model/plots/marker_coverage.png
```

---

## Task 7: Wire marker builder + figure into data_model.smk

**Files:**
- Modify: `workflows/data_model.smk` (add `build_marker_sets` + `fig_marker_coverage` rules; add outputs to `rule all`).

**Interfaces:**
- Consumes: Task 5 + Task 6 scripts and their inputs; `config["pathway"]["min_panel_genes"]` (reuse the existing `5`).
- Produces: `results/data_model/marker_panel_coverage.tsv`, `results/data_model/plots/marker_coverage.png` as tracked outputs.

- [ ] **Step 1: Add outputs to `rule all`**

In `workflows/data_model.smk` `rule all:`'s `input:` list, add:

```python
        "results/data_model/marker_panel_coverage.tsv",
        "results/data_model/plots/marker_coverage.png",
```

- [ ] **Step 2: Add the two rules**

Append to `workflows/data_model.smk` (after `build_comparisons`):

```python
# --- cell-marker provenance + panel coverage (mirrors build_gene_sets) ---
rule build_marker_sets:
    input:
        script   = "scripts/build_marker_sets.R",
        lineage  = "config/lineage_markers.yaml",
        substate = "config/substate_markers.yaml",
        prov     = "config/marker_sets_provenance.tsv",
        panel    = "results/data_model/common_genes.tsv",
    output:
        coverage = "results/data_model/marker_panel_coverage.tsv",
    params:
        minpg = config["pathway"]["min_panel_genes"],
    threads: 1
    log: "logs/build_marker_sets.log",
    shell:
        "{RSCRIPT} {input.script} {input.lineage} {input.substate} {input.prov} "
        "{input.panel} {params.minpg} {output.coverage} > {log} 2>&1"
rule fig_marker_coverage:
    input:
        script = "scripts/fig_marker_coverage.R",
        cov    = "results/data_model/marker_panel_coverage.tsv",
    output:
        "results/data_model/plots/marker_coverage.png",
    threads: 1
    log: "logs/fig_marker_coverage.log",
    shell:
        "{RSCRIPT} {input.script} {input.cov} {output} > {log} 2>&1"
```

- [ ] **Step 3: Full dry-run of the workflow**

Run:
```bash
cd /home/jeszyman/repos/spatial-rads
TMPDIR=/mnt/data/projects/spatial-rads/tmp conda run -n basecamp \
  snakemake -s workflows/data_model.smk --dry-run
```
Expected: DAG includes `build_comparisons`, `build_marker_sets`, `fig_marker_coverage`; no error; `rule all` lists all new outputs.

- [ ] **Step 4: Execute the new rules**

Run:
```bash
cd /home/jeszyman/repos/spatial-rads
TMPDIR=/mnt/data/projects/spatial-rads/tmp conda run -n basecamp \
  snakemake -s workflows/data_model.smk --cores 1 \
  results/data_model/marker_panel_coverage.tsv \
  results/data_model/plots/marker_coverage.png
```
Expected: both rules complete green; assertions from Task 5 Step 4 still pass.

- [ ] **Step 5: Commit**

```bash
repo-commit "feat(data-model.smk): wire marker coverage rules" workflows/data_model.smk
```

---

## Task 8: Register this plan in the org index (additive only)

**Scope note (2026-07-12 — supersedes the original delete/fold/consolidate tasks):** Between
this plan being drafted and executed, the entire cross-sample plan-doc cluster was
reorganized by a concurrent session to reflect the `aggregate_typing.smk` /
`aggregate_differential.smk` split. Every deletion candidate from the original survey
(`plan-pathway-genesets{,-impl}.md`, `plan-mbrt-vs-sbrt.md`, `plan-qc-metrics.md`) was **kept and
actively re-edited** by that reorg — so they are owned, not abandoned, and deleting them here
would race a live owner. This plan therefore performs **only the additive, unambiguous** index
update and **defers all deletion, folding, and banner correction to the user** (see the Deferred
block below). The data-model artifacts (Tasks 1–7) are wholly independent of this and stand.

**Files:**
- Modify: `spatial-rads.org` (`** Plans` bullet index + `*** Plan status ledger` table) — ADD the
  new registry plan only; do not remove or rewrite any existing row.

- [ ] **Step 1: Locate the Plans section**

Run:
```bash
cd /home/jeszyman/repos/spatial-rads
grep -nE "^\*+ Plans|Plan status ledger|plan-comparison-registry" spatial-rads.org
```

- [ ] **Step 2: Add the registry-plan bullet + ledger row**

In `spatial-rads.org`, add one bullet to the `** Plans` list and one row to the
`*** Plan status ledger` table for `plan-comparison-registry.md`, describing the
comparison-registry + marker-provenance Stage-A artifacts (status: ACTIVE / this build). Match
the existing bullet + row shape exactly (flat list items per the org conventions; use the
org-edit skill). Do **not** touch any other plan's bullet or row — the concurrent reorg owns those.

- [ ] **Step 3: Verify the addition is well-formed and nothing else moved**

Run:
```bash
cd /home/jeszyman/repos/spatial-rads
grep -nE "plan-comparison-registry" spatial-rads.org
git diff --stat spatial-rads.org
```
Expected: the new bullet/row present; `git diff` shows only additions around the Plans section.

- [ ] **Step 4: Commit**

```bash
repo-commit "docs(org): index the comparison-registry plan" spatial-rads.org
```

---

## Deferred to the user (NOT executed by this plan)

The concurrent reorg left these open; each mutates a doc the reorg actively owns, so they are the
user's call, not this plan's inference:

1. **Stale atlas banners persist.** `plan-aggregate.md` (throughout) and `plan-outcomes-de.md`
   (lines 3, 18, 23–24) still assert "v2.0 external-atlas / NanoString CellProfileLibrary mouse
   mammary atlas" re-typing, which **contradicts** the actual shipped method in `CLAUDE.md` (scVI
   integration → Leiden cluster-then-annotate + tier-2 SingleR/UCell + marker rescue). These want
   correction, but `plan-aggregate.md` was just re-edited by the reorg — confirm intent before
   touching.
2. **Deletion candidates from the original survey are all still live** (`plan-pathway-genesets`
   pair, `plan-mbrt-vs-sbrt.md`, `plan-qc-metrics.md`). If any are genuinely superseded, the
   reorg owner should delete them; content-capture check + `CLAUDE.md` citation update required
   first.

---

## Task 11: Final end-to-end verification

**Files:** none (verification only).

- [ ] **Step 1: Clean full dry-run**

Run:
```bash
cd /home/jeszyman/repos/spatial-rads
TMPDIR=/mnt/data/projects/spatial-rads/tmp conda run -n basecamp \
  snakemake -s workflows/data_model.smk --dry-run
```
Expected: no error; the three new rules present; `rule all` lists `comparisons.tsv`, `marker_panel_coverage.tsv`, `marker_coverage.png`.

- [ ] **Step 2: Full workflow run**

Run:
```bash
cd /home/jeszyman/repos/spatial-rads
TMPDIR=/mnt/data/projects/spatial-rads/tmp conda run -n basecamp \
  snakemake -s workflows/data_model.smk --cores 1
```
Expected: all rules green (existing + 3 new). Confirms the registry/marker artifacts coexist with the shipped pathway artifacts and nothing regressed.

- [ ] **Step 3: Confirm the shipped pathway artifact is unchanged**

Run:
```bash
cd /home/jeszyman/repos/spatial-rads
git status --short results/data_model/pathway_sets.tsv results/data_model/gene_set_panel_coverage.tsv
```
Expected: no modification to the pathway outputs (the new work is additive; pillar 3 is confirm-only).

- [ ] **Step 4: Final assertions battery**

Run:
```bash
cd /home/jeszyman/repos/spatial-rads
conda run -n spatial-rads Rscript -e '
cmp <- data.table::fread("results/data_model/comparisons.tsv")
mk  <- data.table::fread("results/data_model/marker_panel_coverage.tsv")
stopifnot(cmp[name=="SBRT_vs_Ctrl" & resolution=="whole"]$resid_df==6)
stopifnot(cmp[tier=="descriptive", all(inference_capable==FALSE)])
stopifnot(nrow(mk)==16, all(mk[tier!="substate"]$frac_panel==1))
cat("end-to-end verification passed\n")'
```
Expected: `end-to-end verification passed`.

---

## Self-Review notes (addressed inline)

- **Spec coverage:** Artifact A → Tasks 1-3; Artifact B → Tasks 4-7; pathway confirm-only → Task 11 Step 3; plan consolidation → Tasks 8-10; method notes → carried in `comparisons.yaml` `method_notes` (Task 1) and the spec (not code). Deferred `aggregate_differential.smk` rewiring is explicitly out of scope (spec + this plan).
- **Split-awareness (2026-07-12):** the cross-sample layer is now two workflows (`aggregate_typing.smk`, `aggregate_differential.smk`), not the monolithic `aggregate.smk`. The `comparisons.tsv` consumers are all in the differential half. Plan-doc consolidation (Tasks 8-10) is narrowed to the docs NOT entangled with the in-flight aggregate plan-doc reorganization or with live `CLAUDE.md` citations — see the Task 8/9 scope notes.
- **Type consistency:** `comparisons.tsv` columns are defined once (Task 2 interface) and asserted against in Tasks 2/3/11 with the same names; `marker_panel_coverage.tsv` columns defined in Task 5 and asserted in Tasks 5/11. Builder invocation signatures match their smk `shell:` lines (Tasks 3, 7).
- **Placeholder scan:** the one intentional placeholder line in Task 5's script is called out and removed in Task 5 Step 2. No TBDs elsewhere.
- **Known risk:** if the samplesheet's `timepoint_h` for M02 rows is not all `48`, the `mutter02_day2` cohort filter still resolves (it keys on dataset name, not timepoint), but the `treated_vs_baseline_timecourse` group2 predicate needs a real `timepoint_h==0` control row — verified present (Mutter_01 Control, sam0001, timepoint_h=0). If Task 2 Step 3's `resid_df==6` assertion fails, the block structure changed; investigate rather than weaken the assertion.
```