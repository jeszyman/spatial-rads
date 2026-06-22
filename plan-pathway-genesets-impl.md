# Pathway Gene-Sets + Panel Generation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Generate the common gene panel and a tier-structured pathway gene-set artifact once, reproducibly, in `data_model.smk`; rewire every consumer to read them; collapse the redundant double pathway scoring.

**Architecture:** `data_model.smk` (Stage A) gains two rules emitting `results/data_model/common_genes.tsv` and `results/data_model/pathway_sets.tsv` (+ coverage). All gene-set/panel consumers read those artifacts instead of re-deriving (no more live `msigdbr` calls, no duplicated set definitions). Per-cell scoring collapses to the single merged-scale path. Delivered in two phases: Phase 1 touches nothing the locked aggregate must re-run; Phase 2 (the scored→typed rename) is sequenced with the pending aggregate rebuild.

**Tech Stack:** Snakemake (driver in `basecamp` env), R via `conda run -n spatial-rads Rscript` (data.table, yaml, arrow, Seurat, UCell, msigdbr, fgsea). Spec: `plan-pathway-genesets.md`.

> **STATUS (2026-06-22) — PHASE 1 COMPLETE, committed + pushed.** Tasks 1–5 done across commits `0ee5af6` (T1), `0c85dd0` (T2), `3b3f567` (T3), `0f6cbee` (T4), `12bc1d3` (T5), plus `7d7b201` (these plan docs). Verified: artifacts `common_genes.tsv` (950, invariance held) / `pathway_sets.tsv` (54 usable sets) / `gene_set_panel_coverage.tsv` (msigdbr 26.1.0, 58 logged) built + git-tracked; tier-1 freeze and panel-invariance guards fire; gene-set equivalence old-vs-new passes on all 54 shared sets (only the 4 gate-excluded thin Hallmark sets differ, by design); `msigdbr` confined to `scripts/build_gene_sets.R`; `data_model`/`processing` dry-runs clean, `aggregate` resolves with the audit `build_gene_sets` rule retired.
>
> **PHASE 2 (Task 6) NOT STARTED** — gated on the aggregate rebuild (`plan-aggregate-refactor.md`): the `scored.rds`→`typed.rds` merge-input change must run *with* that rebuild, not piecemeal against the locked aggregate. Code edits are stageable anytime; the run waits so the merge re-runs once. The plan's live numeric reproduction (execute `gsea.R`/`pathway_summary.R`, diff result tables) also rides that rebuild — Phase 1 validated only static gene-membership equivalence.

## Global Constraints

- Snakemake invocation: `TMPDIR=/mnt/data/projects/spatial-rads/tmp conda run -n basecamp snakemake -s workflows/<name>.smk --cores N`; dry-run before execute.
- R rules: `threads: 1` + the `shell.prefix` BLAS pinning already in the workflows; keep them.
- Tier-1 (curated) gene-set membership is **frozen** — pre-registered. The build hard-errors if generated tier-1 ≠ committed source.
- The 950-gene common panel content must not change — invariance is identity (`setequal`), not count.
- Nothing silently dropped: every excluded set is logged in the coverage table.
- Commit only the files each task changed, by name (never `git add -A`); use `repo-commit`. Do not push unless asked.
- Tier labels are `{primary, exploratory}` (matches existing scripts; `primary` = curated confirmatory).

---

## PHASE 1 — artifacts + rewire (no aggregate re-run)

### Task 1: Gene-set generator + provenance + config

**Files:**
- Create: `scripts/build_gene_sets.R` (promoted from `scripts/aggregate/build_gene_sets.R`, extended audit→generator)
- Create: `config/pathway_sets_provenance.tsv`
- Modify: `config/config.yaml` (add `pathway.min_panel_genes`)
- Modify: `config/spatial-rads-conda-env.yaml` (pin `r-msigdbr`)

**Interfaces:**
- Produces: `results/data_model/pathway_sets.tsv` (columns `set, tier, source, gene`), `results/data_model/gene_set_panel_coverage.tsv` (columns `set, tier, source, n_total, n_panel, usable, thin`), and a header `# msigdbr_version: X.Y.Z` line in the coverage file.
- Consumes: `config/pathway_gene_lists.yaml` (flat `name -> [genes]`, unchanged shape), `config/pathway_sets_provenance.tsv`, `results/data_model/common_genes.tsv` (from Task 2).

- [x] **Step 1: Create the provenance table** (`config/pathway_sets_provenance.tsv`)

```tsv
set	source
TypeI_interferon	NanoString Mouse UCC panel module annotation (Annotations sheet)
TypeII_interferon	NanoString Mouse UCC panel module annotation (Annotations sheet)
DNA_Damage_Repair	NanoString Mouse UCC panel module annotation (Annotations sheet)
STING	hand-curated cGAS-STING (panel-available genes)
Angiogenesis	hand-curated (Hallmark/canonical)
Hypoxia	hand-curated (Hallmark/canonical)
Fibrosis_remodeling	hand-curated (Hallmark/canonical)
Stromal_stress_senescence	hand-curated (Hallmark/canonical)
```

- [x] **Step 2: Add config knob** — append under the existing `normalize:` block in `config/config.yaml`:

```yaml
pathway:
  min_panel_genes: 5        # tier-2 (exploratory) set kept only if >= this many genes on panel; tier-1 always kept (flagged thin)
```

- [x] **Step 3: Pin msigdbr** — in `config/spatial-rads-conda-env.yaml`, change the `- r-msigdbr` line to pin the currently-installed version. First read it: `conda run -n spatial-rads Rscript -e 'cat(as.character(packageVersion("msigdbr")))'`, then set e.g. `- r-msigdbr=10.0.1` (use the printed value).

- [x] **Step 4: Write the generator** (`scripts/build_gene_sets.R`)

```r
#!/usr/bin/env Rscript
# Generate the canonical tier-structured pathway gene-set artifact + panel-coverage
# table. Tier-1 (primary) = curated sets from pathway_gene_lists.yaml (frozen, provenance
# from pathway_sets_provenance.tsv); Tier-2 (exploratory) = MSigDB mouse Hallmark via
# msigdbr (version recorded). Single source of truth for every downstream gene-set reader.
# Args: <pathway_gene_lists.yaml> <provenance.tsv> <common_genes.tsv> <min_panel_genes>
#       <out_pathway_sets.tsv> <out_coverage.tsv>
suppressPackageStartupMessages({library(yaml); library(data.table); library(msigdbr)})

a        <- commandArgs(trailingOnly = TRUE)
yaml_p   <- a[1]; prov_p <- a[2]; panel_p <- a[3]
min_pg   <- as.integer(a[4]); out_sets <- a[5]; out_cov <- a[6]

panel <- readLines(panel_p)
prov  <- fread(prov_p)                                  # set, source

## ---- tier-1: curated (frozen) ----
prim   <- lapply(read_yaml(yaml_p), as.character)
miss_prov <- setdiff(names(prim), prov$set)
if (length(miss_prov)) stop("no provenance for primary set(s): ", paste(miss_prov, collapse=", "))
prim_dt <- rbindlist(lapply(names(prim), function(s) data.table(
  set = s, tier = "primary",
  source = prov[set == s, source], gene = prim[[s]])))

## ---- tier-2: MSigDB mouse Hallmark (pinned, namespaced) ----
mver <- as.character(packageVersion("msigdbr"))
hm   <- as.data.table(msigdbr(species = "Mus musculus", collection = "H"))
hm_dt <- hm[, .(set = gs_name, tier = "exploratory",
                source = sprintf("MSigDB_Hallmark (msigdbr %s)", mver),
                gene = gene_symbol)]
hm_dt <- unique(hm_dt)
stopifnot(all(grepl("^HALLMARK_", hm_dt$set)))         # namespaced -> no collision with primary

sets_long <- rbind(prim_dt, hm_dt)

## ---- coverage + usability gate ----
cov <- sets_long[, .(n_total = uniqueN(gene),
                     n_panel = uniqueN(gene[gene %in% panel])),
                 by = .(set, tier, source)]
cov[, usable := tier == "primary" | n_panel >= min_pg] # tier-1 always usable
cov[, thin   := n_panel < min_pg]

## ---- completeness guard: every primary set present ----
miss_prim <- setdiff(names(prim), unique(sets_long[tier == "primary", set]))
if (length(miss_prim)) stop("primary set(s) missing from artifact: ", paste(miss_prim, collapse=", "))

## ---- write usable sets only; coverage logs ALL (incl. excluded) ----
keep <- cov[usable == TRUE, set]
out  <- sets_long[set %in% keep]
dir.create(dirname(out_sets), recursive = TRUE, showWarnings = FALSE)
fwrite(out, out_sets, sep = "\t")
writeLines(sprintf("# msigdbr_version: %s", mver), out_cov)
fwrite(cov, out_cov, sep = "\t", append = TRUE, col.names = TRUE)
cat(sprintf("pathway_sets: %d primary + %d hallmark; %d usable, %d excluded (logged)\n",
            uniqueN(prim_dt$set), uniqueN(hm_dt$set), length(keep), cov[usable==FALSE, .N]))
```

- [x] **Step 5: Smoke-test the generator standalone** (panel not yet at the new path — point at the current one for now):

Run:
```bash
conda run -n spatial-rads Rscript scripts/build_gene_sets.R \
  config/pathway_gene_lists.yaml config/pathway_sets_provenance.tsv \
  results/processing/common_genes.tsv 5 \
  /tmp/ps.tsv /tmp/cov.tsv && head -3 /tmp/cov.tsv && cut -f2 /tmp/ps.tsv | sort -u
```
Expected: prints `# msigdbr_version:` header, tier values `{primary, exploratory}`, ~58 sets total, no error.

- [x] **Step 6: Verify the tier-1 freeze guard fires** — temporarily add a junk gene to a scratch yaml copy and confirm... actually the freeze guard lives in Task 2's panel rule? No — re-read: tier-1 freeze means generated primary == committed yaml. Since the generator READS the yaml, primary always matches its own input. The freeze guard belongs at the **commit boundary**: the yaml is git-tracked, so any membership change shows in `git diff`. Add an explicit check instead: assert generated primary genes `setequal` the committed `config/pathway_gene_lists.yaml` via git:

Append to `scripts/build_gene_sets.R` after the completeness guard:
```r
## ---- tier-1 freeze: generated primary must match the git-committed yaml ----
committed <- tryCatch(yaml::yaml.load(system2("git", c("show", "HEAD:config/pathway_gene_lists.yaml"),
                                              stdout = TRUE)), error = function(e) NULL)
if (!is.null(committed)) {
  cg <- lapply(committed, as.character)
  for (s in names(prim))
    if (!setequal(prim[[s]], cg[[s]] %||% character()))
      stop("tier-1 freeze: primary set '", s, "' differs from committed config/pathway_gene_lists.yaml")
}
`%||%` <- function(a,b) if (is.null(a)) b else a
```
(Move the `%||%` definition to the top of the script.) Re-run Step 5; expect no error.

- [x] **Step 7: Commit**

```bash
repo-commit "feat(genesets): tier-structured gene-set generator + provenance + msigdbr pin" \
  scripts/build_gene_sets.R config/pathway_sets_provenance.tsv config/config.yaml config/spatial-rads-conda-env.yaml
```

---

### Task 2: Fix + relocate the common-panel generator

**Files:**
- Modify: `scripts/common_gene_panel.R` (read both datasets from RDS; add invariance guard)

**Interfaces:**
- Consumes: `results/data_model/samples.tsv` (cols `name, raw_input_path`).
- Produces: `results/data_model/common_genes.tsv` (950 lines).

- [x] **Step 1: Rewrite `scripts/common_gene_panel.R`**

```r
#!/usr/bin/env Rscript
# Common gene panel = Mutter_01 panel INT Mutter_02 panel, both read from the RNA-assay
# rownames of one per-slide raw RDS (both datasets are format=rds). Invariance guard:
# the freshly-built list must setequal the currently-committed common_genes.tsv before
# overwrite (protects the locked aggregate typed on this panel).
# Args: <samplesheet.tsv> <out.tsv>
suppressMessages({library(readr); library(dplyr); library(Seurat)})
args <- commandArgs(trailingOnly = TRUE); SS <- args[1]; OUT <- args[2]
ss <- read_tsv(SS, show_col_types = FALSE)

rds_genes <- function(name) {
  p <- ss %>% filter(name == !!name) %>% pull(raw_input_path) %>% unique()
  stopifnot(length(p) >= 1, file.exists(p[1]))
  o <- readRDS(p[1]); g <- rownames(o); rm(o); gc(); g
}
m01 <- rds_genes("Mutter_01"); m02 <- rds_genes("Mutter_02")
common <- sort(intersect(m01, m02))
cat(sprintf("Mutter_01: %d | Mutter_02: %d | common: %d\n", length(m01), length(m02), length(common)))
stopifnot(length(common) > 500)

## ---- invariance guard vs committed snapshot (identity, not count) ----
if (file.exists(OUT)) {
  old <- readLines(OUT)
  if (!setequal(common, old))
    stop(sprintf("panel invariance FAILED: %d added, %d removed vs committed %s",
                 length(setdiff(common, old)), length(setdiff(old, common)), OUT))
}
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
writeLines(common, OUT)
```

- [x] **Step 2: Move the committed snapshot to the new path (preserve git history)**

```bash
git mv results/processing/common_genes.tsv results/data_model/common_genes.tsv
```

- [x] **Step 3: Test the rewritten script reproduces the panel**

Run:
```bash
TMPDIR=/mnt/data/projects/spatial-rads/tmp conda run -n spatial-rads Rscript \
  scripts/common_gene_panel.R results/data_model/samples.tsv results/data_model/common_genes.tsv
```
Expected: `common: 950`, no invariance error, `git diff --stat results/data_model/common_genes.tsv` shows no content change.

- [x] **Step 4: Commit**

```bash
repo-commit "fix(panel): read M01/M02 genes from RDS + invariance guard; relocate to results/data_model" \
  scripts/common_gene_panel.R results/data_model/common_genes.tsv
```

---

### Task 3: Wire both rules into `data_model.smk`; drop the processing panel rule

**Files:**
- Modify: `workflows/data_model.smk`
- Modify: `workflows/processing.smk` (remove `common_gene_panel` rule)

**Interfaces:**
- Produces (DAG): `results/data_model/{common_genes.tsv, pathway_sets.tsv, gene_set_panel_coverage.tsv}`.

- [x] **Step 1: Add rules + targets to `workflows/data_model.smk`** (after `make_data_model`):

```python
rule all:
    input:
        config["samplesheet"],
        "results/data_model/common_genes.tsv",
        "results/data_model/pathway_sets.tsv",
        "results/data_model/gene_set_panel_coverage.tsv",

rule common_gene_panel:
    input:
        script = "scripts/common_gene_panel.R",
        ss     = config["samplesheet"],
    output:
        "results/data_model/common_genes.tsv",
    log: "logs/common_gene_panel.log",
    shell:
        "{RSCRIPT} {input.script} {input.ss} {output} > {log} 2>&1"

rule build_gene_sets:
    input:
        script = "scripts/build_gene_sets.R",
        yaml   = config["pathway"]["gene_lists"],
        prov   = "config/pathway_sets_provenance.tsv",
        panel  = "results/data_model/common_genes.tsv",
    output:
        sets     = "results/data_model/pathway_sets.tsv",
        coverage = "results/data_model/gene_set_panel_coverage.tsv",
    params:
        minpg = config["pathway"]["min_panel_genes"],
    log: "logs/build_gene_sets.log",
    shell:
        "{RSCRIPT} {input.script} {input.yaml} {input.prov} {input.panel} {params.minpg} "
        "{output.sets} {output.coverage} > {log} 2>&1"
```
Also add to `config/config.yaml` under `pathway:` a `gene_lists: config/pathway_gene_lists.yaml` key (so `config["pathway"]["gene_lists"]` resolves). The `rule all` `input:` replaces the existing single-line `all`.

- [x] **Step 2: Remove the panel rule from `processing.smk`** — delete the `rule common_gene_panel:` block (the `{RSCRIPT} scripts/common_gene_panel.R {input.ss} {output}` rule). Leave `GENES` for Task 4.

- [x] **Step 3: Dry-run data_model**

Run:
```bash
TMPDIR=/mnt/data/projects/spatial-rads/tmp conda run -n basecamp snakemake \
  -s workflows/data_model.smk --dry-run
```
Expected: lists `common_gene_panel`, `build_gene_sets`; no cycle/`MissingInputException`.

- [x] **Step 4: Execute data_model (produces the artifacts)**

Run the same command without `--dry-run --cores 1`. Expected: `common_genes.tsv` unchanged (invariance passes), `pathway_sets.tsv` + coverage written.

- [x] **Step 5: Commit**

```bash
repo-commit "feat(data_model): generate common panel + tiered gene sets as Stage-A artifacts" \
  workflows/data_model.smk workflows/processing.smk config/config.yaml \
  results/data_model/pathway_sets.tsv results/data_model/gene_set_panel_coverage.tsv
```

---

### Task 4: Repoint panel-path consumers

**Files (modify the panel path `results/processing/common_genes.tsv` → `results/data_model/common_genes.tsv`):**
- `workflows/processing.smk` (`GENES` var, L19)
- `workflows/aggregate.smk` (`PANEL` var, L21)
- `scripts/adapt_mutter01.R`, `scripts/adapt_mutter02.R`, `scripts/probe_qc.R`, `scripts/build_celltype_reference.R`, `scripts/aggregate/prepare_reference.R` (header arg comments only — the path is passed in; update any hardcoded default)
- `config/lineage_markers.yaml` (the verification comment, L2)

- [x] **Step 1: Find every reference**

Run: `grep -rn "results/processing/common_genes.tsv\|processing/common_genes" workflows/ scripts/ config/`

- [x] **Step 2: Replace each occurrence** with `results/data_model/common_genes.tsv` (the workflows pass `GENES`/`PANEL` as args, so this is the var definition + any comment/default).

- [x] **Step 3: Verify no stale refs**

Run: `grep -rn "processing/common_genes" workflows/ scripts/ config/` → expect no output.

- [x] **Step 4: Dry-run processing + aggregate to confirm the DAG still resolves**

Run:
```bash
TMPDIR=/mnt/data/projects/spatial-rads/tmp conda run -n basecamp snakemake -s workflows/processing.smk --dry-run
TMPDIR=/mnt/data/projects/spatial-rads/tmp conda run -n basecamp snakemake -s workflows/aggregate.smk --dry-run
```
Expected: no `MissingInputException`. (If aggregate dry-run flags unrelated pre-existing gaps per `plan-aggregate-refactor.md`, note them but they are out of scope.)

- [x] **Step 5: Commit**

```bash
repo-commit "refactor: point panel consumers at results/data_model/common_genes.tsv" \
  workflows/processing.smk workflows/aggregate.smk scripts/adapt_mutter01.R scripts/adapt_mutter02.R \
  scripts/probe_qc.R scripts/build_celltype_reference.R scripts/aggregate/prepare_reference.R config/lineage_markers.yaml
```

---

### Task 5: Rewire gene-set readers to the artifact

**Files:**
- Modify: `scripts/aggregate/gsea.R`, `scripts/aggregate/pathway_summary.R`, `scripts/aggregate/assemble_results.R`, `scripts/aggregate/panel_coverage.R`
- Modify: `workflows/aggregate.smk` (swap `yaml = "config/pathway_gene_lists.yaml"` → `sets = "results/data_model/pathway_sets.tsv"` in the `gsea`, `pathway_summary`, `panel_coverage`, `assemble_results` rules; retire the audit-only `build_gene_sets` rule + its `gene_set_panel_coverage.tsv` target now produced by data_model)
- Check (repoint only if they read pathway sets, not lineage markers): `scripts/contamination_qc.R`, `scripts/qc_plots.R`, `scripts/aggregate/tier2_stroma_ucell.R`

**Interfaces:**
- Consumes: `results/data_model/pathway_sets.tsv` (`set, tier, source, gene`).
- A shared reconstruction replaces the `read_yaml + msigdbr` block: build `all_sets` (named list) and `set_meta` (`pathway/pathway_name, pathway_source, tier`) from the TSV.

- [x] **Step 1: Define the replacement block** (used in gsea.R and pathway_summary.R). For `gsea.R`, replace lines 28-37 (the `prim_lists`…`set_size_full` construction and the `library(msigdbr)`) with:

```r
# --- gene sets: tier-structured artifact (single source of truth) ---
gs_long <- fread(sets_path)                      # sets_path = args[2]
all_sets <- split(gs_long$gene, gs_long$set)
all_sets <- lapply(all_sets, unique)
set_meta <- unique(gs_long[, .(pathway_name = set, pathway_source = source, tier)])
set_size_full <- vapply(all_sets, length, integer(1))
```
Remove `library(msigdbr)` from the `suppressPackageStartupMessages` block. Rename `yaml_path <- args[2]` to `sets_path <- args[2]`.

- [x] **Step 2: Apply the analogous edit to `pathway_summary.R`** — replace lines 51-62 (the `prim_lists`…`set_meta[, n_set_genes := ...]` block) with:

```r
gs_long  <- fread(sets_path)                      # sets_path = args[3]
all_sets <- lapply(split(gs_long$gene, gs_long$set), unique)
set_meta <- unique(gs_long[, .(pathway = set, pathway_source = source, tier)])
set_meta[, n_set_genes := vapply(all_sets[pathway], length, integer(1))]
```
Remove `library(msigdbr)`; rename `yaml_path <- args[3]` → `sets_path <- args[3]`.

- [x] **Step 3: Edit `assemble_results.R`** — it reads `yaml_p` for H1/H2/H3 gene membership (line ~38 `gs <- lapply(read_yaml(yaml_p), as.character)`). Replace with a read of the artifact filtered to primary:

```r
gs_dt <- fread(yaml_p)                             # now pathway_sets.tsv (renamed arg)
gs    <- split(gs_dt[tier == "primary", gene], gs_dt[tier == "primary", set])
```
Rename the CLI arg/comment from `<pathway_yaml>` to `<pathway_sets>`; update the `aggregate.smk` `assemble_results` rule input accordingly.

- [x] **Step 4: Edit `panel_coverage.R`** the same way (read `pathway_sets.tsv`, `split(gene, set)`); remove any `msigdbr`/`read_yaml`.

- [x] **Step 5: Update `workflows/aggregate.smk`** — in rules `gsea`, `pathway_summary`, `panel_coverage`, `assemble_results`: change the `yaml = "config/pathway_gene_lists.yaml"` input to `sets = "results/data_model/pathway_sets.tsv"` and the positional arg in the `shell:` line. Delete the `rule build_gene_sets:` block (now produced by data_model); update any rule that listed `results/aggregate/gene_set_panel_coverage.tsv` to consume `results/data_model/gene_set_panel_coverage.tsv` instead.

- [x] **Step 6: Check the three yaml-readers** — `grep -n "pathway_gene_lists\|read_yaml\|msigdbr" scripts/contamination_qc.R scripts/qc_plots.R scripts/aggregate/tier2_stroma_ucell.R`. If a script reads `pathway_gene_lists.yaml` for **pathway** sets, repoint to `pathway_sets.tsv` + `split`. If it reads `config/lineage_markers.yaml` (typing markers), leave it — out of scope.

- [x] **Step 7: Equivalence test (the real "did it change behavior" check)** — reconstruct the set list both ways and assert identical:

```bash
conda run -n spatial-rads Rscript -e '
suppressMessages({library(data.table); library(yaml); library(msigdbr)})
old_p <- lapply(read_yaml("config/pathway_gene_lists.yaml"), as.character)
hm <- as.data.table(msigdbr(species="Mus musculus", collection="H"))
old <- c(old_p, lapply(split(hm$gene_symbol, hm$gs_name), unique))
new_dt <- fread("results/data_model/pathway_sets.tsv")
new <- lapply(split(new_dt$gene, new_dt$set), unique)
# compare on the usable sets the artifact kept:
common <- intersect(names(old), names(new))
bad <- common[!vapply(common, function(s) setequal(old[[s]], new[[s]]), logical(1))]
cat(if (length(bad)==0) "EQUIVALENT\n" else paste("DIFFER:", paste(bad, collapse=", "), "\n"))'
```
Expected: `EQUIVALENT` (the only legitimate differences are tier-2 sets excluded by the `min_panel_genes` gate, which are absent from `new` by design — they won't be in `common`).

- [x] **Step 8: Dry-run aggregate**

Run: `TMPDIR=/mnt/data/projects/spatial-rads/tmp conda run -n basecamp snakemake -s workflows/aggregate.smk --dry-run`. Expected: no `msigdbr`/yaml-related MissingInput; `build_gene_sets` rule gone.

- [x] **Step 9: Confirm msigdbr fully removed from consumers**

Run: `grep -rn "msigdbr" scripts/` → only `scripts/build_gene_sets.R`.

- [x] **Step 10: Commit**

```bash
repo-commit "refactor(pathway): all gene-set readers consume results/data_model/pathway_sets.tsv (no live msigdbr)" \
  scripts/aggregate/gsea.R scripts/aggregate/pathway_summary.R scripts/aggregate/assemble_results.R \
  scripts/aggregate/panel_coverage.R workflows/aggregate.smk
```

---

## PHASE 2 — collapse double scoring (sequence with the aggregate rebuild)

> Do NOT run Phase 2 against the locked aggregate piecemeal. It changes the merge input (`scored.rds` → `typed.rds`), so it must run as part of the intentional aggregate rebuild in `plan-aggregate-refactor.md` (where `merged.rds` is regenerated anyway). Land the code edits, then run the rebuild once.

### Task 6: Remove per-sample scoring; merge from typed.rds

**Files:**
- Modify: `workflows/processing.smk` (delete `pathway_score` rule; `rule all` → `typed.rds`)
- Delete: `scripts/pathway_score.R`
- Modify: `workflows/aggregate.smk` (`SCORED` → typed dir; all `*.scored.rds` expands → `*.typed.rds`)
- Modify: `scripts/aggregate/merge.R`, `merge_pilot.R`, `coords_necrosis.R`, `recover_negprobes.R`, `cell_assignment_map.R` (the `sub("\\.scored\\.rds$", ...)` prefix-strip → `\\.typed\\.rds$`)

**Interfaces:**
- The per-sample terminal output becomes `processing/typed/{s}.typed.rds`; aggregate merges those. `pathway_summary.R` (merged scale, sourcing `pathway_sets.tsv`) is the single per-cell scoring path.

- [ ] **Step 1: Delete the per-sample scoring rule** — remove `rule pathway_score:` from `processing.smk`; change `rule all` `input:` from `…/scored/{s}.scored.rds` to `…/typed/{s}.typed.rds`. `git rm scripts/pathway_score.R`.

- [ ] **Step 2: Switch aggregate to typed inputs** — in `workflows/aggregate.smk` set `SCORED = f"{DATADIR}/processing/typed"` (rename var to `TYPED` for clarity) and replace every `{{s}}.scored.rds` with `{{s}}.typed.rds`.

- [ ] **Step 3: Update the prefix-strip in the 5 scripts** — replace `sub("\\.scored\\.rds$", "", ...)` with `sub("\\.typed\\.rds$", "", ...)` in `merge.R`, `merge_pilot.R`, `coords_necrosis.R`, `recover_negprobes.R`, `cell_assignment_map.R`.

- [ ] **Step 4: Grep clean**

Run: `grep -rn "scored\.rds\|/scored/\|pathway_score" workflows/ scripts/` → expect no output.

- [ ] **Step 5: Dry-run both workflows**

Run the processing + aggregate `--dry-run`. Expected: processing terminates at `typed.rds`; aggregate merges `typed.rds`; no `MissingInputException`.

- [ ] **Step 6: Commit (do not run the rebuild here — hand to plan-aggregate-refactor.md)**

```bash
repo-commit "refactor(scoring): collapse to single merged-scale path; per-sample ends at typed.rds" \
  workflows/processing.smk workflows/aggregate.smk scripts/aggregate/merge.R scripts/aggregate/merge_pilot.R \
  scripts/aggregate/coords_necrosis.R scripts/aggregate/recover_negprobes.R scripts/aggregate/cell_assignment_map.R
```

---

## Self-Review notes

- **Spec coverage:** artifacts (T1-T3), guards (panel invariance T2-S1, tier-1 freeze T1-S6, completeness T1-S4), consumer rewire (T4-T5), msigdbr pin (T1-S3), scoring collapse (T6), sequencing (Phase 2 banner) — all mapped.
- **Phase-1 reproducibility** is enforced by the T5 equivalence test (set lists identical) + the panel invariance guard, so no downstream numeric change before the intentional rebuild.
- **Out of scope:** `lineage_markers.yaml` typing markers; the aggregate typing method; adding new pathways.
