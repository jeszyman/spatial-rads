# Gene-set provenance: computational derivation + cited custom sets

**Date:** 2026-07-13
**Status:** design approved; steps 1-3 (build) to be reviewed before the rescore

## Problem

The eight "primary" pathway sets feeding pathway scoring and the coverage figure
(`fig_geneset_coverage`) are hand-typed in `config/pathway_gene_lists.yaml` and
labeled in `config/pathway_sets_provenance.tsv` as "NanoString module annotation"
or "hand-curated (Hallmark/canonical)". Both labels are false: the IFN-I/II/DDR
sets do NOT match the vendor Annotations-sheet modules (IFN-I 7/19, IFN-II 1/18,
DDR 7/19 overlap), and the stromal/vascular sets match neither the vendor modules
nor MSigDB Hallmark cleanly. No set carries a real citation, and the git-HEAD
"freeze" guards a mislabeled artifact.

## Goal

Every scored pathway set traces to one of three honest provenances:
- **nanostring**: derived computationally from the vendor Annotations sheet.
- **hallmark**: derived computationally from MSigDB (`msigdbr`, pinned version).
- **custom**: hand-curated for a concept no computational source covers, carrying
  a required citation.

The coverage figure color-codes rows by provenance and shows a concept from BOTH
computational sources where both exist (cross-source coverage audit).

## Scope revision (2026-07-13, post-design)

The design expanded during review from a relabel/rederive into a full pathway-tier
rebuild plus a pre-registration revision. Final decided state:

- **Confirmatory family: 7 concepts** (was 8). Per concept the highest-COVERAGE
  (not raw-count) source is confirmatory; the panel-designed NanoString module wins
  all 7 dual-source concepts (83-97% coverage). **STING is retired from the
  confirmatory family** (H1): its 3-gene set has only 2 on-panel genes, both shared
  with the IFN sets (already flagged weak in the 2026-06-23 robustness pass). The
  richer cGAS-STING/ICD biology moves to an exploratory custom set. `assemble_results.R`
  H1 definition updated to TypeI-IFN + TypeII-IFN only for the STING-gene DE/program rows.
- **msigdbr mouse-native**: `db_species="MM"`, collection `"MH"` (50 sets), not the
  human-ortholog default.
- **Exploratory tier**: all mouse Hallmark sets clearing the >=5 on-panel-gene gate
  (47 of 50) scored exploratory, restoring the full hypothesis-generating sweep the
  first registry cut had narrowed to 7 twins. Plus three custom-cited RT signatures
  (below). Tiered FDR keeps these isolated: exploratory rows carry their own
  within-analysis BH (`padj_exploratory`); the confirmatory family is pooled-BH over
  the 7 pre-registered concepts only (`padj_confirmatory`), so widening the
  exploratory scan cannot dilute H1/H2/H3.

### Custom-cited RT signatures (exploratory)

Three radiotherapy-response signatures MSigDB/NanoString do not cover, entering as
`source=custom` exploratory sets. **Gene memberships + DOIs must be verified against
the primary papers before build (no fabricated lists or DOIs).**

- **RSI** (radiosensitivity index): Eschrich et al., Int J Radiat Oncol Biol Phys 2009.
- **SASP** (senescence-associated secretory phenotype): Coppe 2010 / Basisty 2020.
- **ICD / cGAS-STING** (immunogenic cell death, cytosolic-DNA sensing): Galluzzi /
  Vanpouille-Box; folds in the retired STING genes.

## Architecture

### Concept registry (new single provenance surface)

`config/pathway_concept_registry.tsv`, one row per (concept, source):

| column | meaning |
|---|---|
| `concept` | study concept (`TypeI_interferon`, `DNA_Damage_Repair`, `Angiogenesis`, ...) |
| `source` | `nanostring` \| `hallmark` \| `custom` |
| `source_id` | `;`-separated vendor module name(s) / `HALLMARK_*` name(s); `NA` for custom |
| `citation` | DOI/PMID, required for `custom`, `NA` for derived |

`source_id` may name multiple modules/sets; the builder takes their union. This
keeps the union choice transparent and reviewable (recorded in the registry)
rather than buried in hand-typed genes.

### Members live where they are authoritative (never duplicated)

- `nanostring`: genes tagged `+` in the named module column(s) of the vendor
  Annotations sheet (`data/sources/2026-06-22-bruker-mouse-ucc-gene-list.xlsx`).
- `hallmark`: `msigdbr(species="Mus musculus", collection="H")` filtered to the
  named `HALLMARK_*` set(s).
- `custom`: `config/pathway_gene_lists.yaml`, which shrinks to only genuinely
  bespoke sets (e.g. STING), each now paired with a citation via the registry.

### Builder refactor (`scripts/build_gene_sets.R`)

Reads registry, derives nanostring/hallmark members by `source_id`, reads custom
members from YAML, emits one long `pathway_sets.tsv` tagged
`concept, source, source_id, provenance, gene`, plus the coverage table
(`gene_set_panel_coverage.tsv`) with per-(concept,source) `n_total`, `n_panel`,
`coverage`.

**Confirmatory selection:** per concept, the (concept, source) with the highest
panel overlap (`n_panel`) is flagged `confirmatory = TRUE`; ties go to nanostring.
All other sets are scored but tagged `exploratory`. Custom sets are confirmatory
for gap concepts. Result: exactly one confirmatory set per concept, giving a
stable, pre-registerable BH family.

**Freeze re-scope:** the git-HEAD freeze guards the two hand-authored surfaces,
`pathway_concept_registry.tsv` and the custom `pathway_gene_lists.yaml`, not the
machine-derived members (which legitimately change with the vendor file or msigdbr
version). Retire `config/pathway_sets_provenance.tsv`.

### Figure (`scripts/fig_geneset_coverage.R`)

Rows are (concept, source) pairs, grouped by concept; a dual-source concept renders
two rows. Fill encodes provenance (nanostring / hallmark / custom, categorical
colorblind-safe, no red-green). Confirmatory row marked (bold/asterisk or shape).
Custom rows carry a citation marker. `theme_scifig()`, caption baked into the PNG
per the scifig spec.

## Concept to source_id map (initial; the review surface)

| concept | nanostring source_id | hallmark source_id | custom |
|---|---|---|---|
| TypeI_interferon | Type I Interferon Signaling | HALLMARK_INTERFERON_ALPHA_RESPONSE | no |
| TypeII_interferon | Type II Interferon Signaling | HALLMARK_INTERFERON_GAMMA_RESPONSE | no |
| DNA_Damage_Repair | DNA Damage Repair | HALLMARK_DNA_REPAIR | no |
| Angiogenesis | Angiogenesis;VEGF Signaling | HALLMARK_ANGIOGENESIS | no |
| Hypoxia | HIF1 Signaling | HALLMARK_HYPOXIA | no |
| Fibrosis_remodeling | Collagen;TGF-beta Signaling;EMT | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION;HALLMARK_TGF_BETA_SIGNALING | no |
| Stromal_stress_senescence | Senescence;Cellular Stress;p53 Signaling | HALLMARK_P53_PATHWAY | no |
| STING | none | none | cite required |

## Sequencing

1. Build registry + refactor `build_gene_sets.R`, regenerate `pathway_sets.tsv` +
   coverage table (local, cheap).
2. Diff new confirmatory sets vs current frozen lists, report membership changes.
3. Rebuild figure; **user reviews** (checkpoint, stop here).
4. Rescore: `pathway_scores.R` (~5h UCell/AMS over `merged.rds`, confirmed present
   on beast), then GSEA, then `assemble_results.R` recomputes `padj_*`.
5. Re-validate (`validate_master.R`: SBRT fibrosis / MBRT null / p21-ambiguous
   hold) + report which confirmatory pathway rows moved.
6. Re-freeze: commit registry + custom YAML as new pre-registered state; update org
   notebook + CLAUDE.md.

Gene-level DE, composition, niches, mixing, myeloid, substate do NOT re-run
(pathway-set-independent).

## Files touched

- new `config/pathway_concept_registry.tsv`
- `config/pathway_gene_lists.yaml` (shrinks to custom + citations)
- `scripts/build_gene_sets.R`
- `scripts/fig_geneset_coverage.R`
- retire `config/pathway_sets_provenance.tsv`
- org notebook (Gene sets rule + figure block) + rule wiring; CLAUDE.md

## Open items

- STING citation: no fabricated DOI; placeholder flagged for the user to supply.
- Module-union choices in `source_id` are scientific judgments, recorded in the
  registry for review; the figure + registry are the review artifact at step 3.
