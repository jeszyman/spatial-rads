# M01 4-hour expression tables: summary

Descriptive differential expression from the Mutter_01 CosMx dataset at
4 hours post-irradiation, for R01 preliminary data. All analyses use
current merged-scale atlas labels (3,277,090-cell unified typing,
2026-06-02 lock).

## Samples

| Sample   | Condition   | Total cells | Tumor  | Stroma | Immune |
|----------|-------------|-------------|--------|--------|--------|
| sam0001  | Control 0h  | 127,399     | 56,946 | 57,408 | 13,045 |
| sam0003  | MBRT 4h     | 168,086     | 57,039 | 94,282 | 16,765 |
| sam0006  | SBRT 4h     | 166,163     | 61,788 | 89,799 | 14,575 |

Mutter_01 is n=1 per arm. No formal inference; all results are
descriptive effect sizes.

## Table index

### Sample-level descriptive DE (12 tables)

Pairwise pseudobulk log2FC via edgeR `exactTest` with fixed BCV=0.2
(the standard value for genetically identical inbred mice; edgeR User
Guide Section 2.9). TMM-normalized. Gene universe: 950-gene common panel.

| File | Contrast | Compartment |
|------|----------|-------------|
| `sample_level_MBRT_vs_Ctrl_all.tsv` | MBRT 4h vs Control | All cells |
| `sample_level_MBRT_vs_Ctrl_tumor.tsv` | MBRT 4h vs Control | Tumor |
| `sample_level_MBRT_vs_Ctrl_stroma.tsv` | MBRT 4h vs Control | Stroma |
| `sample_level_MBRT_vs_Ctrl_immune.tsv` | MBRT 4h vs Control | Immune |
| `sample_level_SBRT_vs_Ctrl_all.tsv` | SBRT 4h vs Control | All cells |
| `sample_level_SBRT_vs_Ctrl_tumor.tsv` | SBRT 4h vs Control | Tumor |
| `sample_level_SBRT_vs_Ctrl_stroma.tsv` | SBRT 4h vs Control | Stroma |
| `sample_level_SBRT_vs_Ctrl_immune.tsv` | SBRT 4h vs Control | Immune |
| `sample_level_MBRT_vs_SBRT_all.tsv` | MBRT 4h vs SBRT 4h | All cells |
| `sample_level_MBRT_vs_SBRT_tumor.tsv` | MBRT 4h vs SBRT 4h | Tumor |
| `sample_level_MBRT_vs_SBRT_stroma.tsv` | MBRT 4h vs SBRT 4h | Stroma |
| `sample_level_MBRT_vs_SBRT_immune.tsv` | MBRT 4h vs SBRT 4h | Immune |

Columns: `gene, log2FC, logCPM, pct_expressing_group1, pct_expressing_group2, hypothesis`

Positive log2FC = higher in the first-named arm. Genes belonging to
frozen hypothesis sets (immune_activation, vascular_hypoxia,
stromal_fibrosis, myofibroblast_expansion, normal_vessel_sparing,
myeloid_M2) are flagged in the `hypothesis` column.

**BCV assumption.** The fixed BCV of 0.2 assumes biological variability
between samples of the same condition is comparable to that of
genetically identical inbred mice (Balb/c host + 4T1 syngeneic tumor).
The p-values from `exactTest` are conditional on this assumed dispersion
and carry no inferential weight at n=1.

### Peak-vs-valley FOV-pseudobulk paired DE (8 tables)

Within the MBRT 4h sample (sam0003), cells assigned to peak or valley
zones by the fitted stripe model (15-deg tilt, 1.02 mm spacing, 4 beam
centers). Core bands: distance-to-nearest-peak < 0.10 mm = peak,
> 0.40 mm = valley.

| Zone | Tumor | Stroma | Immune | Total |
|------|-------|--------|--------|-------|
| Peak | 6,947 | 10,385 | 2,315 | 19,647 |
| Transition | 19,315 | 29,414 | 6,290 | 55,019 |
| Valley | 30,777 | 54,483 | 8,160 | 93,420 |

Model: `~ FOV + zone` (FOV as a blocking factor). Pseudobulk counts
summed per FOV x zone. DESeq2 is primary; limma-voom is the concordance
check. Gene universe: 940 genes (950-gene panel minus 10 on-panel DDR
genes excluded to avoid circularity with the stripe model fit).

| File | Method | Compartment | Paired FOVs |
|------|--------|-------------|-------------|
| `pv_deseq2_all.tsv` | DESeq2 | All cells | 26 |
| `pv_deseq2_tumor.tsv` | DESeq2 | Tumor | 25 |
| `pv_deseq2_stroma.tsv` | DESeq2 | Stroma | 26 |
| `pv_deseq2_immune.tsv` | DESeq2 | Immune | 24 |
| `pv_voom_all.tsv` | limma-voom | All cells | 26 |
| `pv_voom_tumor.tsv` | limma-voom | Tumor | 25 |
| `pv_voom_stroma.tsv` | limma-voom | Stroma | 26 |
| `pv_voom_immune.tsv` | limma-voom | Immune | 24 |

Columns: `gene, log2FC, SE, p, padj_BH, qvalue_storey, compartment`

**Multiple testing.** BH-adjusted p-values are primary. Storey q-values
computed as a concordance check. With pi0 = 1.0 across all compartments,
the two are identical: the estimated proportion of true nulls is 100%.

**KS uniformity test.** The p-value distribution was tested against the
uniform with a Kolmogorov-Smirnov test. KS p < 0.001 in all compartments,
but this reflects the heavy-null anti-conservative shape of the
distribution (slight excess of large p-values), not signal. Zero genes
reach padj < 0.10 in any compartment by either method.

### DDR genes excluded from peak-vs-valley

Cdkn1a, Gadd45a, Gadd45b, Mdm2, Bax, Pmaip1, Bbc3, Ddit3, Atm, Atr,
Chek1, Chek2, Rad51, Brca1, Brca2, Xrcc5, Xrcc6, Parp1 (18 genes;
10 present on the 950-gene panel). These genes were used to fit the
stripe model's beam registration; including them in the peak-vs-valley
DE test would be circular.

### Continuous distance-to-peak regression (4 tables)

Per-gene linear model of raw expression on continuous distance-to-peak,
with FOV as a fixed-effect blocking factor. All cells in the sample
(not restricted to core bands).

| File | Compartment |
|------|-------------|
| `pv_distance_all.tsv` | All cells |
| `pv_distance_tumor.tsv` | Tumor |
| `pv_distance_stroma.tsv` | Stroma |
| `pv_distance_immune.tsv` | Immune |

Columns: `gene, coefficient, SE, p, compartment`

### Validation tables

| File | Description |
|------|-------------|
| `validation_offset_scan.tsv` | Beam-center offset sweep (minus 0.50 to +0.50 mm, 0.05 mm steps) with per-module peak-vs-valley mean differences |
| `validation_pure_fov.tsv` | Whole-FOV classification by centroid distance, swept over offsets and thresholds |
| `validation_moran_tumor.tsv` | Moran's I spatial autocorrelation on distance-regression residuals (tumor, top 5 genes) |

## Key findings

### Sample-level: shared acute radiation response at 4h

Both MBRT and SBRT show the same acute radiation response relative to
unirradiated control: DNA damage response genes induced (Cdkn1a/p21
+0.86 log2FC in MBRT, +0.53 in SBRT), heat-shock proteins activated
(Hspa1a, Hspa1b), and early inflammatory cytokines elevated (Il1a,
Il6, Il11, Tnf). Structural and extracellular matrix genes (collagens,
Vim, Fn1, Lgals1) are lower in the irradiated samples, consistent with
acute tissue damage rather than fibrosis activation. SBRT shows larger
magnitude shifts in the structural component (36 genes with
|log2FC| > 1 in SBRT-vs-Ctrl tumor, vs 7 in MBRT-vs-Ctrl), but the
programs are qualitatively the same.

The direct MBRT-vs-SBRT contrast is nearly flat (1 gene with
|log2FC| > 1 in tumor). At 4h, the two modalities produce similar
acute responses; the divergence into distinct SBRT-fibrosis vs
MBRT-null programs appears later (M02 day-2 data).

### Peak-vs-valley: confirmed null on current labels

Zero genes reach BH-adjusted p < 0.10 in any compartment by either
DESeq2 or limma-voom. Storey pi0 = 1.0 across all compartments,
confirming the estimated proportion of true nulls is 100%. This result
holds on the current merged-scale atlas labels, consistent with the
prior audit on older per-sample labels.

### Registration does not rescue the null

The offset scan swept assumed beam-center positions across +/- 0.50 mm.
The p21 module shows a smooth periodic response peaking at offset
+0.20 mm (mean difference = 0.13 normalized-expression units), with
DDR and IFN modules peaking near offset -0.20 mm (differences
0.07-0.09). At +0.20 mm offset, 41% of cells change zone assignment,
confirming the offset is biologically meaningful. However, the signal
magnitude remains small (0.05-0.13 units) across the entire sweep:
no registration choice produces gene-level significance. The
peak/valley classification method is not the reason for the null result.

### Pure-FOV design is weaker, not stronger

The pure-FOV design (classifying whole FOVs by centroid position rather
than individual cells) produces similar or smaller module-mean
differences and t-statistics compared to the within-FOV paired design.
The paired `~ FOV + zone` approach is the more powerful test; the null
is not an artifact of within-FOV mixing.

### Distance regression: weak spatial gradient in structural genes

The continuous distance-to-peak regression identifies Lgals1, S100a6,
Lmna, S100a4, and Vim as the strongest spatial-gradient genes in tumor
(all with negative coefficients: expression decreases with distance
from peak). Moran's I on the FOV-blocked residuals is statistically
significant (p = 0.001 by permutation) but small in magnitude
(I = 0.03-0.08), indicating weak residual spatial autocorrelation
beyond the FOV blocking structure.

## Methods summary

**Cell labeling.** All analyses use unified per-cell labels from the
merged-scale atlas typing (scVI integration, Leiden clustering,
cluster-level annotation, locked 2026-06-02; 3,277,090 cells across
M01 + M02). Compartment assignments: tumor, stroma, immune.

**Sample-level DE.** Per-compartment pseudobulk (raw count sums across
all cells of one sample in one compartment). edgeR `exactTest` with
manually set BCV = 0.2 (genetically identical Balb/c + 4T1 syngeneic
system). TMM normalization. 950-gene common panel, all genes included.

**Peak-vs-valley DE.** Stripe model: 15-degree tilt, 1.02 mm beam
spacing, 4 beam centers at 1.14/2.16/3.18/4.20 mm (perpendicular
projection coordinates). Zone assignment by distance-to-nearest-peak:
< 0.10 mm = peak, > 0.40 mm = valley. FOV-pseudobulk paired design:
raw counts summed per FOV x zone, model `~ FOV + zone`. DESeq2
(primary) and limma-voom (concordance). 940 genes (10 on-panel DDR
genes excluded). BH correction primary; Storey q-values as concordance.

**Distance regression.** Per-gene `lm(expression ~ FOV + dist_to_peak)`
on raw counts, all cells (not restricted to core bands). Moran's I
computed on FOV-blocked residuals using k=15 nearest-neighbor spatial
weights (RANN), significance by 999 permutations.

**Software.** R 4.4, Seurat v5, edgeR, DESeq2, limma, qvalue, RANN.
Conda environment `spatial-rads` on jeff-beast.
