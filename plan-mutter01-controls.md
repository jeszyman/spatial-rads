# Plan: Mutter_01 raw-RDS adoption — staged control sidecar (falsecode + negprobe)

**Status:** control sidecar implemented 2026-06-15. **Adapter re-base + data-model `format=rds` flip DONE 2026-06-22** (see "Update 2026-06-22" below) — `adapt_mutter01.R` now ingests the 4 per-slide RDS, counts verified **byte-identical** to the parquet across all 11 M01 samples. Control-retirement (`recover_negprobes.R`) and the aggregate re-run remain deferred.
**Supersedes nothing.** Additive to the locked `aggregate.smk` typing + MBRT-vs-SBRT layer (2026-06-02/03) and the completed `processing.smk` per-sample pipeline (2026-05-30).

## Update 2026-06-22 — adapter re-base executed early, decoupled from the aggregate rebuild

Driven by adding CosMx-good-practice **contamination QC** to `processing.smk` (SpatialQM MECR, report-only, sample + per-FOV; new `contamination_qc` rule + `scripts/contamination_qc.R`), the raw-RDS adapter re-base was pulled forward and **decoupled** from the control-retirement work it was bundled with here:

- **Done:** `adapt_mutter01.R` rewritten onto the 4 per-slide RDS (RNA counts from the RDS; `Condition` + Yi labels + shared morphology joined from the metadata parquet by reconstructed `cell_id`; `setequal()` cohort guard). `format` flipped to `rds` + per-slide `raw_input_path` in `data/metadata.xlsx` + `scripts/build_metadata_xlsx.R`. **Invariance proven:** a scratch run reproduced all 11 M01 `raw.rds` byte-identical (cells, genes, counts).
- **Deliberately NOT done (kept lean / invariance-preserving):** negprobes+falsecode are **still dropped** from per-sample objects (controls stay in the `control_qc` sidecar / `cell_controls.parquet`, avoiding the ~194-feature memory cost); `recover_negprobes.R` + the aggregate's `cell_neg.tsv` **untouched**; `scored.rds` **not** regenerated; falsecode QC **still flag-only**.
- **Approach divergence:** no separate committed `config/mutter01_fov_sample_map.tsv` crosswalk was built — the metadata parquet is **retained as the M01 label/condition side-input** and joined directly (also keeps `probe_qc.R`'s parquet dependency satisfied). 100% `cell_id` reconstruction asserted, no NA-fill.
- **Armed state:** regenerating the samplesheet bumped its mtime, so `processing.smk` now *wants* to reprocess the M01 cascade (invariant outputs) and `aggregate.smk` its rebuild — both intentionally left for the deliberate frozen-seed rebuild below.

## Motivation

Mutter_01 has now been delivered as **raw per-slide Seurat RDS** (`/mnt/data/gdrive/seuratObject_0[1-4]_Mutter_01_CosMmR.RDS`, 327–584 MB), the same artifact type as Mutter_02. Each carries three assays the historical M01 parquet never held as raw vectors: `RNA` (1000), `negprobes` (10 individual `Negative1–10`), and `falsecode` (184 system-control barcodes), plus 64 metadata columns (morphology IF, `x/y_slide_mm`, `fov`, per-assay counts).

**Two platform negative-control types, distinct error sources:**
- **negprobes** (10) — real probes targeting nothing; estimate *hybridization* background. Drive the per-cell `neg` floor and `propNegative`.
- **falsecodes** (184) — codebook barcodes assigned to no gene; a count is an *optical/decoding* mis-call. A different failure mode, **never used anywhere in this project** (only `adapt_mutter02.R:24` drops the assay; `build_metadata_xlsx.R:25` records the count 184 as panel metadata).

The negprobe-floor good practice is already implemented (`scripts/probe_qc.R` flags genes at/below negprobe background in both datasets, signal-to-background ≤ 1; `scripts/aggregate/tier2_detectability.py` gates marker recall at ≥ 2× negprobe). **The falsecode side is the genuinely new capability** these RDS unlock.

## Key empirical finding — the aggregate is invariant to this data

The new RDS cover the **same cells** as the existing M01 parquet. Verified on slide 01 (all 324,129 cells joined): `nCount_NegativeProbes` (parquet) == `nCount_negprobes` (RDS) with **0 disagreements**, `nCount_RNA` bit-identical, per-slide cell counts match exactly across all 4 slides. The aggregate only ever consumes the *mean* negprobe per cell (`neg = nCount_negprobes/10`) as a scVI covariate + QC background — a scalar the parquet summary already provides exactly. **A full raw-RDS re-base would therefore reproduce the locked typing + differential results; there is no scientific reason to re-run the aggregate for this data.**

## Decision — staged sidecar, full re-base deferred

A blind devils-advocate panel (4 lenses, unanimous solid-with-caveats; logged) found that a full adapter re-base now would: (a) be result-invariant yet invasive (rewrite both adapters, flip the data model, re-run 23 objects); (b) entangle "retire `recover_negprobes.R`" with "don't re-run the locked aggregate" (the aggregate still consumes `recover_negprobes.R`'s `cell_neg.tsv`); (c) risk a memory-margin regression from retaining ~194 control features × 23 objects × 5 stages (CCA already OOM'd `sam0023`); (d) de-calibrate the tuned `max_prop_negative` gate if `propNegative` is redefined.

**Chosen path:** deliver the falsecode data + symmetric control handling **now, additively**, via a sidecar that touches neither adapter nor the locked aggregate inputs. **Defer the full raw-RDS re-base to the next aggregate rebuild**, when `scored.rds` regenerate anyway and `recover_negprobes.R` can genuinely be retired — getting the architectural win when it is free.

## Scope — NOW

1. **Relocate inputs.** `cp` (not symlink) the 4 raw RDS from `/mnt/data/gdrive/` → `/mnt/data/projects/spatial-rads/inputs/mutter01/`, mirroring `inputs/mutter02/`. Keep the parquets in place. The characterization/sidecar scripts take the M01 input dir directly (glob the 4 RDS); optionally record the per-slide RDS paths in the data-model `slides` sheet as informational provenance. **M01's `format` is NOT flipped** — it stays parquet-sourced for processing, so the adapter and `probe_qc.R` are untouched. gdrive is a slow FUSE mount — local copy first (standing rule).

2. **Control characterization (calibration before any threshold).** Standalone script: per-cell and per-FOV negprobe and falsecode distributions for **both** datasets (M01 from the new RDS assays, M02 from its raw assays). Output a report TSV + a cohort-vs-cohort summary. This is the reference distribution the falsecode flag is calibrated against — no QC cutoff is set before this exists.

3. **Per-cell control sidecar (additive, separate artifact).** Extend the `recover_negprobes.R` pattern into a standalone script emitting a per-cell control table for both datasets: `cell_id, dataset, slide, fov, neg, falsecode, falsecode_frac`, built from the raw RDS. **It is a NEW file (e.g. `results/processing/cell_controls.tsv`), NOT the aggregate's `cell_neg.tsv`** — `recover_negprobes.R` and `cell_neg.tsv` are left byte-for-byte untouched so the locked aggregate cannot re-fire. The sidecar's per-cell join asserts **100% coverage and errors on any unmatched cell — never silent NA-fill**.

4. **Flag-only falsecode QC metric.** A per-FOV (and per-cell) falsecode-fraction metric reported alongside existing QC, **report-only — excludes nothing**. `propNegative` stays **negprobe-only** (preserves the calibrated `config/config.yaml max_prop_negative: 0.5` gate and M02's exact behavior). Because no cell or FOV is removed, cell counts are unchanged and aggregate invariance holds by construction.

## Explicitly deferred — next aggregate rebuild

- Rewrite `adapt_mutter01.R` onto the raw RDS (mirroring `adapt_mutter02.R`) via a committed slide×FOV→sample crosswalk.
- Retain `negprobes` + `falsecode` assays on per-sample objects in both adapters.
- Retire `recover_negprobes.R`; collapse the parquet to a single (label-only) role.
- Flip M01 `format` to `rds` in the data model.
- Promote falsecode QC from flag-only to FOV exclusion (a deliberate decision that, by changing cell counts, **explicitly triggers** an aggregate re-run with frozen seeds).

A committed checklist for this deferred work lives in the "Deferred re-base" section below.

## Design detail

### Invariance + join guards (from panel triage)
- Extend the parquet↔RDS identity check from slide 01 to **all 4 slides + the full 987-FOV crosswalk**, and confirm the M01 RDS negprobe assay/column names match M02's (`negprobes` assay, `nCount_negprobes`). Run as a validation gate in the characterization step.
- The sidecar's barcode reconstruction (RDS local `c_1_<fov>_<cell>` → parquet `<run>_S<n>_<fov>_<cell>`) must assert exact, 100% join per slide and **fail loudly** on any miss.

### Do-not-touch list (protects the locked aggregate)
- `scripts/aggregate/recover_negprobes.R` and its `cell_neg.tsv` output — unchanged.
- `workflows/aggregate.smk` and all `scripts/aggregate/*` — not invoked.
- `data/.../scored/*.rds` — not regenerated (would re-touch timestamps and re-fire the aggregate by default).
- `propNegative` definition and `max_prop_negative` — unchanged.

### probe_qc.R (second parquet dependency, panel finding 6)
`scripts/probe_qc.R:21-27` reads M01 gene means + negprobe background directly from the parquet (`counts_path`, `metadata_path`), bypassing `raw.rds`. It is report-only and stays parquet-sourced; the data-model change must **keep `counts_path`/`metadata_path` populated** for M01 so it is unaffected. No edit required.

### Artifact homes
- Characterization + sidecar = standalone scripts first (fast iteration), then a committed `processing.smk` diagnostic rule adjacent to `probe_qc` (report-only, not on the aggregate critical path). Outputs are small TSVs under `results/processing/`.

## Risks / mitigations
- **Aggregate silently re-fires.** Mitigated by the do-not-touch list: no `scored.rds` regenerated, no aggregate rule invoked, `cell_neg.tsv` untouched.
- **Falsecode threshold premature.** Mitigated: characterization precedes any threshold; the metric ships flag-only.
- **M01 RDS column mismatch → all-NA.** Mitigated by the explicit name-match check + 100%-join assertion (error, not NA-fill).
- **Memory.** Sidecar reads one slide RDS at a time and writes a TSV; no object retains the 194 control features — the OOM margin is untouched.

## Documentation surfaces to sync at execution
(per the pilot→full doc-sync rule) — `CLAUDE.md` datasets/inputs + Analysis Architecture note; `spatial-rads.org` notebook entry; this plan's status; `plan-aggregate.md` "deferred negprobe/falsecode work" note. The org `* Plans` heading indexes this file.

## Deferred re-base — checklist (updated 2026-06-22)
- [x] ~~Build + commit `config/mutter01_fov_sample_map.tsv` crosswalk~~ — **superseded by approach:** the metadata parquet is retained as the M01 label/condition side-input and joined directly (1 condition/FOV verified: 987 FOVs, 0 multi); no separate crosswalk file needed.
- [x] Rewrite `adapt_mutter01.R` onto raw RDS — **done 2026-06-22**, Yi labels joined on reconstructed global `cell_id`, invariance proven byte-identical.
- [ ] Retain `negprobes`+`falsecode` on per-sample objects — **declined for now:** controls stay in the `control_qc` sidecar (memory + invariance). Revisit only if a per-object control need arises.
- [x] Flip M01 `format=rds` + per-slide `raw_input_path` — **done 2026-06-22**.
- [ ] Retire `recover_negprobes.R`; source `neg` from a retained assay — **still deferred** (aggregate still consumes `cell_neg.tsv`).
- [ ] Decide falsecode FOV exclusion; if adopted, re-run the aggregate with frozen scVI/Leiden seeds — **still deferred** (falsecode QC stays flag-only; aggregate not re-run).
