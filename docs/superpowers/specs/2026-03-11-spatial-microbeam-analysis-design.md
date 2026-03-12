# Spatial Microbeam Analysis Design

## Overview

Exploratory analysis of CosMx spatial transcriptomics data comparing microbeam radiotherapy (MBRT) vs single-beam radiotherapy (SBRT) vs control in murine xenograft models. Organized as four progressive layers, each building on the last. All analysis code lives in org-babel blocks in `spatial-rads.org`.

### Novelty Position

No published spatial transcriptomics studies of spatially fractionated radiation exist. This dataset enables first-in-field characterization of the peak/valley transcriptomic landscape and spatiotemporal kinetics of immune activation in MBRT vs SBRT.

### Key Literature Context

- MRT/MBRT produces stronger pro-inflammatory cytokine/IFN-gamma profiles vs broadbeam (Potez et al. 2019-2020, Bazyar et al. 2022)
- SFRT induces immunogenic cell death markers more potently (Fernandez-Palomo et al. 2021)
- SFRT shows enhanced abscopal responses with STING/IFN-I upregulation (Yang et al. 2021, Kanagavelu et al. 2014)
- Valley biology involves bystander signaling — valleys may be where immune priming occurs (Smyth et al. 2018)
- Spatial transcriptomics in radiation biology is nascent — no MRT/MBRT studies exist

## Experimental Design

### Samples (n=11)

| Condition | Model | Treatment | Timepoint |
|-----------|-------|-----------|-----------|
| Control | flank (TBD — confirm with lab) | NT | 0h |
| MBRT_1h | flank | MBRT | 1h |
| MBRT_4h | flank | MBRT | 4h |
| MBRT_day2 | flank | MBRT | 48h |
| MBRT_day6 | flank | MBRT | 144h |
| SBRT_4h | flank | SBRT | 4h |
| SBRT_day2 | flank | SBRT | 48h |
| SBRT_day6 | flank | SBRT | 144h |
| Tongue_MBRT_day10 | tongue | MBRT | 240h |
| Tongue_MBRT_day8 | tongue | MBRT | 192h |
| Tongue_SBRT_day10 | tongue | SBRT | 240h |

### Key Constraints

- **n=1 per condition** (more samples coming later). Analysis is descriptive/exploratory — characterize patterns and effect sizes, not formal significance testing.
- **Independent animals** at each timepoint. No spatial coordinate registration across timepoints.
- **Beam geometry is consistent** across animals (fixed collimator), but peak/valley spatial demarcation fades after 4h.
- **Peak/valley ground truth** is available at 4h via co-registered gamma-H2AX immunofluorescence (IF) from the CosMx run — not from transcript data (H2AX/H2afx is not in the CosMx gene panel).
- **Tongue samples** are a separate model; analyzed secondarily, not pooled with flank.
- **Platform:** CosMx spatial transcriptomics (~1000 targeted gene panel, not whole-transcriptome). Panel was likely designed for immune biology, so radiation-specific genes may be underrepresented.
- **No 1h SBRT sample** — earliest SBRT timepoint is 4h. The 1h timepoint is MBRT-only, limiting early kinetic comparison between modalities.

### Data Sources

- Expression counts: `projects_Mutter_01_CosMmR_app_data_raw_tx_counts_matrix.parquet`
- Cell metadata: `Analysis_Mutter_01_CosMmR_Mutter_updated_metadata.parquet`
- Co-registered IF images: gamma-H2AX channel (to be retrieved from lab)
- Stored at: `/mnt/gcs/jeszyman/projects/spatial-rads/inputs/`

## Analysis Layers

### Layer 1 — QC & Landscape

**Goal:** Process raw data into a clean, annotated Seurat object. Establish the cellular landscape across all samples.

**Processing:**
- Create Seurat object from CosMx count matrix + cell metadata (parquet files)
- QC filtering: `qcFlagsCell == "Pass"`, nCount_RNA > 20, nFeature_RNA > 10, propNegative < 0.5
- Normalization: LogNormalize, scale factor 1e4
- Variable features: 2000 genes, VST method
- Dimensionality reduction: PCA (20 PCs, validate with elbow plot), UMAP (30 neighbors, 0.3 min distance)
- Clustering: FindNeighbors + FindClusters, resolution 0.4
- Parallel processing: 48 workers on jeff-beast

**Cell type annotation — two-phase approach:**
- **Phase 1 (now):** Use Yi's pre-computed insitutype labels (`ImmuneAtlas_ImmGen_Main_cell_Types`) with validation:
  - Tabulate all cell type labels and frequencies
  - Spot-check top 5-6 populations against canonical markers in the CosMx panel (e.g., Cd3e/Cd4/Cd8a for T cells, Cd68/Adgre1 for macrophages, Epcam/Krt8 for tumor/epithelial)
  - Resolve the "a" cell type label (~32% of cells, likely tumor/epithelial) — classify using Epcam/Krt expression
  - If markers and labels broadly agree, proceed; if not, flag and consider re-annotation for that population
- **Phase 2 (later):** Re-annotate from raw counts using our own reference mapping once analyses stabilize and more samples arrive

**Batch assessment:**
- UMAP colored by Slide (4 slides: S1-S4) to check for slide-driven clustering
- If batch effects are present, note whether integration (e.g., Harmony) is needed

**Outputs:**
- Processed Seurat object saved as RDS
- QC summary table: cells per sample before/after filtering, gene detection rates
- UMAP plots: by cluster, by sample/condition, by cell type, by slide
- Spatial plots per sample showing cell type distributions
- Cell count summary table: sample x cell type
- Cell type validation summary (marker concordance with Yi's labels)

**Status:** Largely implemented in existing exploratory code. Formalize, add validation steps, and save intermediate objects.

### Layer 2 — MBRT vs SBRT vs Control Comparative Kinetics

**Goal:** Characterize how the two radiation modalities differ in transcriptomic response over time, across cell types.

#### 2a. Differential Expression by Treatment x Timepoint

- DEG analysis comparing:
  - MBRT vs Control at each timepoint
  - SBRT vs Control at each timepoint
  - MBRT vs SBRT at matched timepoints (4h, day 2, day 6)
- Run per cell type for major populations (macrophages, T cells, tumor/epithelial cells at minimum)
- **Pseudo-replication caveat:** With n=1 per condition, cell-level Wilcoxon p-values reflect within-animal variability, not treatment effects. Report effect sizes (log2FC) and percentage-expressed differences as primary metrics. Consider pseudobulk aggregation by FOV (~30-50 FOVs per sample) as a more conservative approach. Do not report p-values as evidence of differential expression until biological replicates arrive.
- Output: DEG tables (ranked by effect size, not p-value), volcano plots, summary heatmap of top DEGs across conditions

#### 2b. Pathway Score Kinetics

- Pathway scores: Type I IFN, Type II IFN, STING, DNA Damage Repair
- **Use Yi's pre-computed pathway score columns** from the metadata for initial exploration. Document their provenance (method, gene lists) — request from Yi/Mutter lab if not available. Recompute from defined gene sets in Phase 2 if needed.
- Summarize per sample x cell type x timepoint
- Kinetic line plots: pathway score over time, MBRT vs SBRT, faceted by cell type
- Descriptive comparisons (n=1); formal testing when replicates arrive

#### 2c. Immune Composition Shifts

- Cell type proportions by sample/condition/timepoint
- Stacked bar charts and composition heatmaps
- Focus: myeloid subsets (macrophage polarization if resolvable), T cell subsets, NK cells
- Composition change trajectories: MBRT vs SBRT over time

#### 2d. Spatial Neighborhood Analysis

- Cell-cell proximity using CosMx spatial coordinates
- Nearest-neighbor cell type counts or Ripley's K (per-FOV to manage compute on 1.4M cells)
- Question: are immune cells spatially clustering differently in MBRT vs SBRT?
- Kept simple — not full niche mapping at this stage

### Layer 3 — Peak/Valley at 4h

**Goal:** Define the transcriptional landscape of microbeam peaks vs valleys at the acute timepoint where spatial patterning is visible. Characterize both peak AND valley biology — valleys (bystander zone) may be where immune priming occurs.

#### 3a. Spatial Peak/Valley Classification

- **Primary method:** Use co-registered gamma-H2AX immunofluorescence intensity from the CosMx IF channel to classify cells as peak vs valley. This is protein-level ground truth, not transcript-based.
  - IF images to be retrieved from lab
  - Extract per-cell H2AX IF intensity, spatial smooth via KDE
  - Classify cells into: peak, valley, boundary/transition zones
  - Visual validation: do classified regions show stripe geometry consistent with microbeam spacing?
- **Fallback (if IF data unavailable per cell):** Composite DDR transcript score from available panel genes (Cdkn1a, Gadd45b, Ddit3, Chek1/2, Parp1, Atm). This is less reliable — these genes have different kinetics and are induced by multiple stresses. Would require careful validation that spatial pattern matches expected beam geometry.

#### 3b. Differential Expression: Peak vs Valley

- DEG analysis between peak and valley cells within the 4h MBRT sample
- Run globally and per cell type
- **Characterize both sides:**
  - Peak biology: DDR, apoptosis, direct damage signatures
  - Valley biology: bystander signaling, immune priming, stress response — this may be where the MBRT immunological advantage originates
- Explicitly exclude DDR classification genes from the signature to avoid circularity

#### 3c. Derive Peak and Valley Signatures

- From DEGs, define two compact gene signatures:
  - **Peak signature:** genes upregulated in peaks (excluding DDR markers used for classification)
  - **Valley signature:** genes upregulated in valleys (bystander/immune priming program)
- Validation: do the signatures spatially reconstruct the stripe pattern when scored across cells?
- **Critical caveat:** Both signatures are calibrated on n=1 animal. They are hypothesis-generating instruments, not validated biomarkers. Downstream projections inherit this uncertainty.

#### 3d. Compare to SBRT at 4h

- Score both peak and valley signatures in the 4h SBRT sample
- Key questions:
  - Does SBRT resemble MBRT peaks (high-dose biology), MBRT valleys (bystander biology), or neither?
  - Is the SBRT transcriptional program a blend of peak+valley, or something distinct?
- **Interpretation note:** If SBRT resembles peaks, this may reflect shared high-dose biology rather than peak-specific biology. Design comparisons to distinguish dose-level effects from spatial-pattern effects where possible (e.g., compare valley-specific genes that should NOT be dose-dependent).

### Layer 4 — Signature Projection Across Timepoints

**Goal:** Track the peak and valley transcriptional signatures beyond 4h to understand how the initial spatial heterogeneity evolves.

**Important framing:** At later timepoints, irradiated cells may have died and been replaced. "Signature persistence" could reflect either (a) cell-intrinsic transcriptional memory of prior irradiation, or (b) new cells recruited into or responding to the damage field. These are distinct mechanisms — analyses should be interpreted with this ambiguity in mind.

#### 4a. Score All Samples

- Apply both peak and valley signatures to every cell across all timepoints and treatments
- Per sample: compute signature score distributions, visualize spatially and on UMAP

#### 4b. Temporal Trajectory

- MBRT samples over time: do peak and valley signature scores decay? At different rates?
- Is decay uniform across cell types, or do some populations retain signatures longer?
- Kinetic plots: mean signature scores over timepoints, split by cell type
- **Valley signature trajectory is particularly interesting** — does early bystander/immune priming in valleys predict later immune infiltration patterns seen in Layer 2c?

#### 4c. Signatures in SBRT

- Does SBRT produce any component of the peak or valley signatures?
- Compare SBRT signature score distributions to MBRT peak cells vs MBRT valley cells
- Track SBRT signature scores over time — does the trajectory differ from MBRT?

#### 4d. Biological Interpretation

- Cross-reference with Layer 2 pathway kinetics: does signature persistence correlate with IFN/STING activation timing?
- If signatures resolve, what replaces them — homogenization or new programs?
- Do valleys show early immune priming signatures that presage later immune composition shifts?

**This is the most exploratory layer.** Findings depend on Layers 2-3. Some analyses may not be informative.

## Technical Stack

- **Literate programming:** org-babel blocks in `spatial-rads.org`
- **Language:** R (Seurat, tidyverse, ggplot2), Python (pandas/parquet for data loading)
- **Compute:** jeff-beast, 48 parallel workers
- **Conda environment:** spatial-rads
- **Data storage:** `/mnt/gcs/jeszyman/projects/spatial-rads/`
- **Version control:** git, code in org-babel blocks tangled as needed

## Design Decisions

1. **No Snakemake pipeline.** Analyses are exploratory and live in org-babel blocks. Pipeline formalization deferred until analyses stabilize and replicates arrive.
2. **Descriptive statistics only** for now (n=1 per condition). Effect sizes over p-values. Framework designed to scale to formal testing when more samples arrive.
3. **Flank is primary, tongue is secondary.** Different models not pooled.
4. **Peak/valley classification uses IF ground truth at 4h.** Transcript-based DDR score is a fallback only. At later timepoints, the derived signatures are projected as continuous scores, not used for binary classification.
5. **Progressive layers.** Each layer stands alone. Stop when data support is exhausted.
6. **Two-phase cell type annotation.** Phase 1: Yi's pre-computed labels with marker validation. Phase 2: re-annotate from raw counts once analyses stabilize.
7. **Valley biology gets equal emphasis.** The bystander zone may be where the MBRT immunological advantage originates — do not treat valleys as merely "not peaks."
8. **Signatures are hypothesis-generating.** Calibrated on n=1 animal. All downstream projections inherit that uncertainty. Validation requires biological replicates.
9. **1000-gene panel creates blind spots.** The most discriminating peak/valley genes may not be in the panel. Interpret negative results cautiously — absence of signal may reflect panel limitations, not biology.

## Open Items

- [ ] Confirm Control sample anatomical model (flank?) with lab
- [ ] Retrieve co-registered gamma-H2AX IF data from lab
- [ ] Request pathway score gene lists and insitutype parameters from Yi
- [ ] Verify microbeam spacing parameters (peak width, center-to-center distance) for spatial validation
