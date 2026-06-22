# spatial-rads -- cross-sample aggregation workflow (flank cohort, 20 samples).
# Consumes the 20 flank per-sample scored.rds from processing.smk and builds a
# merged object -> integrated embeddings -> three analysis tracks (composition,
# cell-type-resolved DE/pathway, spatial structure). See plan-aggregate.md.
# Run: TMPDIR=/mnt/data/projects/spatial-rads/tmp \
#        conda run -n basecamp snakemake -s workflows/aggregate.smk --cores N
# R steps run in the spatial-rads env via `conda run` (driver lives in basecamp).
import pandas as pd
configfile: "config/config.yaml"

# Pin BLAS to 1 thread per R process (see feedback_smk_thread_hygiene): snakemake
# schedules job-level parallelism, not R BLAS multi-threading. Every R rule is threads: 1.
shell.prefix("export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1; ")

RSCRIPT = "conda run -n spatial-rads Rscript"
DATADIR = config["datadir"]
MASTER  = config["samplesheet"]                  # results/data_model/samples.tsv
AGG     = f"{DATADIR}/aggregate"                  # heavy intermediates
SCORED  = f"{DATADIR}/processing/scored"          # per-sample inputs
CPL     = f"{DATADIR}/ref/CellProfileLibrary/Mouse/Adult"    # external atlas profiles
PANEL   = "results/processing/common_genes.tsv"   # 950 shared-panel gene list
# Raw inputs for negprobe recovery (dropped during 950-gene harmonization):
M01_META = f"{DATADIR}/inputs/mutter01/Analysis_Mutter_01_CosMmR_Mutter_updated_metadata.parquet"
M02_RAW  = f"{DATADIR}/inputs/mutter02"                       # 4 raw slide RDS

# --- cohort: flank only (tongue is separate biology, n=1; out of scope) ---
_s    = pd.read_csv(MASTER, sep="\t")
FLANK = _s.loc[_s["model"] == "flank", "sample_id"].tolist()    # 20: M01 x8, M02 x12

# Total post-QC flank cells, for the merge pilot's cell-based memory projection.
_qc         = pd.read_csv("results/processing/qc_summary.tsv", sep="\t")
FLANK_CELLS = int(_qc.loc[_qc["sample_id"].isin(FLANK), "post_filter"].sum())

# Memory pilot subset: 5 smallest M01 + 5 smallest M02 flank samples (~1.0M cells,
# ~30% of cohort) -- representative but safe to load before the full 20-sample merge.
PILOT = ["sam0002", "sam0007", "sam0005", "sam0008", "sam0001",
         "sam0018", "sam0020", "sam0021", "sam0017", "sam0019"]

rule all:
    input:
        "results/aggregate/merge_pilot_memory.tsv",
        f"{AGG}/merged.rds",
        "results/aggregate/merge_summary.tsv",
        "results/aggregate/composition_test_m02day2.tsv",
        "results/aggregate/plots/composition_m02day2_forest.png",
        "results/aggregate/pseudobulk_qc.tsv",
        "results/aggregate/degs_pseudobulk_m02day2.tsv",
        "results/aggregate/gsea_pseudobulk_m02day2.tsv",
        "results/aggregate/pathway_scores_summary.tsv",
        "results/aggregate/pathway_test_m02day2.tsv",
        "results/aggregate/plots/pathway_heatmap_m02day2.png",
        "results/aggregate/celltype_embed_qc.tsv",
        "results/aggregate/plots/celltype_umap.png",
        "results/aggregate/reference_coverage.tsv",
        f"{AGG}/merged_typed.rds",
        "results/aggregate/celltype_atlas_summary.tsv",
        "results/aggregate/celltype_atlas_validation.tsv",
        # --- MBRT-vs-SBRT downstream differential layer (plan-mbrt-vs-sbrt-impl.md) ---
        "results/aggregate/gene_set_panel_coverage.tsv",
        "results/aggregate/readout_detection_m02.tsv",
        "results/aggregate/celltype_qc_markers.tsv",
        "results/aggregate/concordance_m01_m02.tsv",
        "results/aggregate/results_master.tsv",   # terminal: tier-tagged master table

# --- Stage 0a: memory pilot (characterize peak RSS, choose merge strategy) ---
rule merge_pilot:
    input:
        script = "scripts/aggregate/merge_pilot.R",
        rds    = expand(f"{SCORED}/{{s}}.scored.rds", s=PILOT),
    output:
        "results/aggregate/merge_pilot_memory.tsv",
    params:
        flank_cells = FLANK_CELLS,
    threads: 1
    log:
        "logs/aggregate/merge_pilot.log",
    shell:
        "{RSCRIPT} {input.script} {output} {params.flank_cells} {input.rds} > {log} 2>&1"

# --- Stage 0b: single-pass merge of all 20 flank scored.rds (sparse, no densify) ---
rule merge:
    input:
        script = "scripts/aggregate/merge.R",
        ss     = MASTER,
        rds    = expand(f"{SCORED}/{{s}}.scored.rds", s=FLANK),
    output:
        rds     = f"{AGG}/merged.rds",
        summary = "results/aggregate/merge_summary.tsv",
        meta    = f"{AGG}/cell_metadata.tsv",
    params:
        scale_factor = config["normalize"]["scale_factor"],
    threads: 1
    log:
        "logs/aggregate/merge.log",
    shell:
        "{RSCRIPT} {input.script} {input.ss} {params.scale_factor} {output.rds} {output.summary} {output.meta} {input.rds} > {log} 2>&1"

# --- Stage 1a: cell-type integration embedding (Harmony) + QC clusters ---
# Heavy: ScaleData densifies 950 x 3.27M (~25 GB) + UMAP/SNN on all cells. Run
# after pathway_summary frees memory; not concurrent with other merged.rds jobs.
# lisi_min is the post-Harmony LISI(dataset) gate (rule fails below it). Set to
# 1.5: the in-script LISI is unweighted inverse-Simpson, so the imbalanced flank
# cohort (M01 29% / M02 71%) caps LISI(dataset) at 1/(.29^2+.71^2) ~= 1.70 even
# with perfect mixing -- the plan's 1.8 was unachievable. 1.5 ~= 88% of ceiling.
rule embed_celltype:
    input:
        script = "scripts/aggregate/embed_celltype.R",
        rds    = f"{AGG}/merged.rds",
    output:
        rds     = f"{AGG}/merged_celltype.rds",
        qc      = "results/aggregate/celltype_embed_qc.tsv",
        metrics = "results/aggregate/celltype_embed_metrics.tsv",
        umap    = "results/aggregate/plots/celltype_umap.png",
    params:
        lisi_min = 1.5,
    threads: 1
    log:
        "logs/aggregate/embed_celltype.log",
    shell:
        "{RSCRIPT} {input.script} {input.rds} {output.rds} {output.qc} "
        "{output.metrics} {output.umap} {params.lisi_min} > {log} 2>&1"

# --- Stage 1b inputs: external-atlas reference profile matrix (de-anchored) ---
# Mammary-centered per the agreed plan: NanoString CellProfileLibrary MammaryGland_Virgin
# (tissue-matched, 4T1 = mammary carcinoma) + ImmGen (immune depth). Lung + Muscle dropped
# 2026-06-01 (off-plan, diluted the match); SmoothMuscle left unanchored. Coarse-relabels to
# lineages, subsets to the 950 panel. No M01 data: identical knowledge to M01/M02. See prepare_reference.R.
rule prepare_reference:
    input:
        script = "scripts/aggregate/prepare_reference.R",
        panel  = PANEL,
    output:
        ref      = f"{AGG}/ref_profiles.rds",
        coverage = "results/aggregate/reference_coverage.tsv",
    params:
        cpl = CPL,
    threads: 1
    log:
        "logs/aggregate/prepare_reference.log",
    shell:
        "{RSCRIPT} {input.script} {params.cpl} {input.panel} {output.ref} "
        "{output.coverage} > {log} 2>&1"

# --- Stage 1b inputs: per-cell negative-probe background for InSituType ---
# Negprobes (the platform's background model) were dropped during 950-gene
# harmonization; recover from raw inputs. M01: metadata parquet nCount_NegativeProbes;
# M02: raw slide RDS negprobes assay. neg = nCount_neg / 10 (mouse 1k panel has 10
# negative probes -- one definition spans both datasets). Output barcodes match
# merge.R's {sample_id}_{raw_barcode} key.
rule recover_negprobes:
    input:
        script = "scripts/aggregate/recover_negprobes.R",
        m01    = M01_META,
        m02    = expand(f"{M02_RAW}/seuratObject_0{{i}}_Mutter_02_CosMmR.RDS", i=[1,2,3,4]),
        rds    = expand(f"{SCORED}/{{s}}.scored.rds", s=FLANK),
    output:
        neg = f"{AGG}/cell_neg.tsv",
    params:
        m02_dir = M02_RAW,
    threads: 1
    log:
        "logs/aggregate/recover_negprobes.log",
    shell:
        "{RSCRIPT} {input.script} {output.neg} {input.m01} {params.m02_dir} "
        "{input.rds} > {log} 2>&1"

# --- Stage 1b: field-standard CosMx cell typing via InSituType (Danaher 2022) ---
# Platform-native semi-supervised classifier: models the CosMx noise structure
# (per-cell negprobe background + Poisson counts) that flat correlation ignores, and
# types from RAW COUNTS (merged.rds) -- decoupled from the heavy embed. Runs ONCE on
# merged M01+M02 against the 12-lineage external atlas, so both cohorts get one
# identical rule -> structural comparability. De novo clusters absorb the 4T1 tumor
# (Epcam/Krt8 overlay -> tumor_epithelial); the rest map to nearest lineage.
rule typing_insitutype:
    input:
        script = "scripts/aggregate/typing_insitutype.R",
        rds    = f"{AGG}/merged.rds",
        ref    = f"{AGG}/ref_profiles.rds",
        neg    = f"{AGG}/cell_neg.tsv",
    output:
        typed      = f"{AGG}/merged_typed.rds",
        summary    = "results/aggregate/celltype_atlas_summary.tsv",
        validation = "results/aggregate/celltype_atlas_validation.tsv",
        labels     = f"{AGG}/cell_atlas_labels.tsv",
    threads: 1
    log:
        "logs/aggregate/typing_insitutype.log",
    shell:
        "{RSCRIPT} {input.script} {input.rds} {input.ref} {input.neg} "
        "{output.typed} {output.summary} {output.validation} {output.labels} "
        "> {log} 2>&1"

# --- Track 1: cell-type composition (M02 day2 propeller test + M01 descriptive) ---
rule composition:
    input:
        script = "scripts/aggregate/composition.R",
        obs    = f"{AGG}/full/obs.parquet",
        labels = "results/aggregate/full_labels.parquet",
    output:
        by_sample  = "results/aggregate/composition_by_sample.tsv",
        test       = "results/aggregate/composition_test_m02day2.tsv",
        dropped    = "results/aggregate/composition_dropped_celltypes.tsv",
        bars       = "results/aggregate/plots/composition_m02day2_bars.png",
        forest     = "results/aggregate/plots/composition_m02day2_forest.png",
        timecourse = "results/aggregate/plots/composition_m01_timecourse.png",
    threads: 1
    log:
        "logs/aggregate/composition.log",
    shell:
        "{RSCRIPT} {input.script} {input.obs} {input.labels} {output.by_sample} "
        "{output.test} {output.dropped} {output.bars} {output.forest} "
        "{output.timecourse} > {log} 2>&1"

# --- Track 2: pseudobulk construction (M02 day2, sample x cell_type count sums) ---
rule pseudobulk_build:
    input:
        script = "scripts/aggregate/pseudobulk_build.R",
        rds    = f"{AGG}/merged_typed.rds",
        labels = "results/aggregate/full_labels.parquet",
    output:
        se = f"{AGG}/pseudobulk_se.rds",
        qc = "results/aggregate/pseudobulk_qc.tsv",
    threads: 1
    log:
        "logs/aggregate/pseudobulk_build.log",
    shell:
        "{RSCRIPT} {input.script} {input.rds} {input.labels} {output.se} {output.qc} > {log} 2>&1"

# --- Track 2 inference: pseudobulk DESeq2 DE (M02 day2, abundance-floored) ---
rule deg_pseudobulk:
    input:
        script = "scripts/aggregate/deg_pseudobulk.R",
        se     = f"{AGG}/pseudobulk_se.rds",
    output:
        degs    = "results/aggregate/degs_pseudobulk_m02day2.tsv",
        summary = "results/aggregate/deg_summary_m02day2.tsv",
        skipped = "results/aggregate/degs_pseudobulk_skipped.tsv",
    threads: 1
    log:
        "logs/aggregate/deg_pseudobulk.log",
    shell:
        "{RSCRIPT} {input.script} {input.se} {output.degs} {output.summary} "
        "{output.skipped} > {log} 2>&1"

# --- Track 2 pathway: GSEA on pseudobulk stat-ranked genes (primary + Hallmark) ---
rule gsea:
    input:
        script = "scripts/aggregate/gsea.R",
        degs   = "results/aggregate/degs_pseudobulk_m02day2.tsv",
        yaml   = "config/pathway_gene_lists.yaml",
    output:
        gsea = "results/aggregate/gsea_pseudobulk_m02day2.tsv",
    threads: 1
    log:
        "logs/aggregate/gsea.log",
    shell:
        "{RSCRIPT} {input.script} {input.degs} {input.yaml} {output.gsea} > {log} 2>&1"

# --- Track 2 pathway: per-cell UCell + AddModuleScore scoring + M02 limma test ---
# threads: 4 -- UCell runs ncores=4 fork (BiocParallel) parallelism, not BLAS, so the
# threads:1 BLAS-hygiene convention is unaffected (BLAS stays pinned to 1 thread).
rule pathway_summary:
    input:
        script = "scripts/aggregate/pathway_summary.R",
        rds    = f"{AGG}/merged_typed.rds",
        labels = "results/aggregate/full_labels.parquet",
        yaml   = "config/pathway_gene_lists.yaml",
    output:
        summary    = "results/aggregate/pathway_scores_summary.tsv",
        test       = "results/aggregate/pathway_test_m02day2.tsv",
        conc       = "results/aggregate/pathway_ucell_ams_concordance.tsv",
        heatmap    = "results/aggregate/plots/pathway_heatmap_m02day2.png",
        timecourse = "results/aggregate/plots/pathway_timecourse_m01.png",
        scatter    = "results/aggregate/plots/pathway_ucell_vs_ams_scatter.png",
    threads: 4
    log:
        "logs/aggregate/pathway_summary.log",
    shell:
        "{RSCRIPT} {input.script} {input.rds} {input.labels} {input.yaml} {output.summary} "
        "{output.test} {output.conc} {output.heatmap} {output.timecourse} "
        "{output.scatter} > {log} 2>&1"

# ============================================================================
# MBRT-vs-SBRT downstream differential layer (plan-mbrt-vs-sbrt-impl.md T4-T13).
# All consume the unified per-cell labels (results/aggregate/full_labels.parquet)
# and the day-2 readout tables, terminating in the tier-tagged results_master.tsv.
# ============================================================================

# --- T4: panel coverage of the curated gene sets (which programs are scorable) ---
rule build_gene_sets:
    input:
        script = "scripts/aggregate/build_gene_sets.R",
        panel  = PANEL,
        yaml   = "config/pathway_gene_lists.yaml",
    output:
        coverage = "results/aggregate/gene_set_panel_coverage.tsv",
    threads: 1
    log:
        "logs/aggregate/build_gene_sets.log",
    shell:
        "{RSCRIPT} {input.script} {input.panel} {input.yaml} {output.coverage} > {log} 2>&1"

# --- T5: per-cell-type readout detection report (panel detectability in M02) ---
rule panel_coverage:
    input:
        script = "scripts/aggregate/panel_coverage.R",
        rds    = f"{AGG}/merged_typed.rds",
        labels = "results/aggregate/full_labels.parquet",
        yaml   = "config/pathway_gene_lists.yaml",
    output:
        detection = "results/aggregate/readout_detection_m02.tsv",
    threads: 1
    log:
        "logs/aggregate/panel_coverage.log",
    shell:
        "{RSCRIPT} {input.script} {input.rds} {input.labels} {input.yaml} "
        "{output.detection} > {log} 2>&1"

# --- T6: power / minimum-detectable-effect table at n=4 ---
rule power_mde:
    input:
        script = "scripts/aggregate/power_mde.R",
        comp   = "results/aggregate/composition_by_sample.tsv",
        se     = f"{AGG}/pseudobulk_se.rds",
        path   = "results/aggregate/pathway_scores_summary.tsv",
    output:
        mde = "results/aggregate/power_mde.tsv",
    threads: 1
    log:
        "logs/aggregate/power_mde.log",
    shell:
        "{RSCRIPT} {input.script} {input.comp} {input.se} {input.path} {output.mde} > {log} 2>&1"

# --- T7: per-cell slide coords + per-sample necrosis flag (spatial-track input) ---
rule coords_necrosis:
    input:
        script = "scripts/aggregate/coords_necrosis.R",
        ss     = MASTER,
        rds    = expand(f"{SCORED}/{{s}}.scored.rds", s=FLANK),
    output:
        coords = f"{AGG}/coords_necrosis.parquet",
    threads: 1
    log:
        "logs/aggregate/coords_necrosis.log",
    shell:
        "{RSCRIPT} {input.script} {input.ss} {output.coords} {input.rds} > {log} 2>&1"

# --- T8: cell-type-label QC marker dotplot (does each label express its markers) ---
rule celltype_qc:
    input:
        script = "scripts/aggregate/celltype_qc.R",
        rds    = f"{AGG}/merged.rds",
        labels = "results/aggregate/full_labels.parquet",
    output:
        markers = "results/aggregate/celltype_qc_markers.tsv",
        dotplot = "results/aggregate/plots/celltype_qc_dotplot.png",
    threads: 1
    log:
        "logs/aggregate/celltype_qc.log",
    shell:
        "{RSCRIPT} {input.script} {input.rds} {input.labels} {output.markers} "
        "{output.dotplot} > {log} 2>&1"

# --- T9: data-driven spatial niches (k=20 NN composition -> K=6 k-means) ---
rule niches:
    input:
        script = "scripts/aggregate/niches.R",
        labels = "results/aggregate/full_labels.parquet",
        coords = f"{AGG}/coords_necrosis.parquet",
        obs    = f"{AGG}/full/obs.parquet",
    output:
        per_cell  = f"{AGG}/niche_per_cell.parquet",
        centroids = "results/aggregate/niche_centroids.tsv",
        freq      = "results/aggregate/niche_frequency.tsv",
        test      = "results/aggregate/niche_test_m02day2.tsv",
        heatmap   = "results/aggregate/plots/niche_centroids_heatmap.png",
        freq_plot = "results/aggregate/plots/niche_frequency_m02.png",
    threads: 1
    log:
        "logs/aggregate/niches.log",
    shell:
        "{RSCRIPT} {input.script} {input.labels} {input.coords} {input.obs} "
        "{output.per_cell} {output.centroids} {output.freq} {output.test} "
        "{output.heatmap} {output.freq_plot} > {log} 2>&1"

# --- T10: tumor-immune spatial mixing (immune-neighbour fraction + Keren score) ---
rule spatial_mixing:
    input:
        script = "scripts/aggregate/spatial_mixing.R",
        labels = "results/aggregate/full_labels.parquet",
        coords = f"{AGG}/coords_necrosis.parquet",
        obs    = f"{AGG}/full/obs.parquet",
    output:
        per_sample = "results/aggregate/spatial_mixing_per_sample.tsv",
        test       = "results/aggregate/spatial_mixing_test_m02day2.tsv",
        per_cell   = f"{AGG}/spatial_mixing_per_cell.parquet",
        plot       = "results/aggregate/plots/mixing_m02.png",
    threads: 1
    log:
        "logs/aggregate/spatial_mixing.log",
    shell:
        "{RSCRIPT} {input.script} {input.labels} {input.coords} {input.obs} "
        "{output.per_sample} {output.test} {output.per_cell} {output.plot} > {log} 2>&1"

# --- T11: myeloid M1/M2 polarization (UCell M1/M2 panels -> per-sample ratio) ---
# threads: 8 -- AddModuleScore_UCell forks UCELL_CORES=8 workers (BiocParallel, not
# BLAS), so reserve 8 to avoid oversubscription; BLAS stays pinned to 1 via shell.prefix.
rule myeloid_polarization:
    input:
        script = "scripts/aggregate/myeloid_polarization.R",
        rds    = f"{AGG}/merged.rds",
        labels = "results/aggregate/full_labels.parquet",
    output:
        scores = "results/aggregate/myeloid_m1m2_scores.tsv",
        test   = "results/aggregate/myeloid_m1m2_test_m02day2.tsv",
        plot   = "results/aggregate/plots/myeloid_m1m2_m02.png",
    threads: 8
    log:
        "logs/aggregate/myeloid_polarization.log",
    shell:
        "{RSCRIPT} {input.script} {input.rds} {input.labels} {output.scores} "
        "{output.test} {output.plot} > {log} 2>&1"

# --- T12: M01<->M02 day-2 effect-size concordance on the unified labels ---
rule concordance_m01_m02:
    input:
        script = "scripts/aggregate/concordance_m01_m02.R",
        rds    = f"{AGG}/merged.rds",
        labels = "results/aggregate/full_labels.parquet",
        degs   = "results/aggregate/degs_pseudobulk_m02day2.tsv",
    output:
        tsv     = "results/aggregate/concordance_m01_m02.tsv",
        scatter = "results/aggregate/plots/concordance_scatter.png",
    threads: 1
    log:
        "logs/aggregate/concordance_m01_m02.log",
    shell:
        "{RSCRIPT} {input.script} {input.rds} {input.labels} {input.degs} "
        "{output.tsv} {output.scatter} > {log} 2>&1"

# --- T13: tier-tagged master results table with confirmatory-family FDR ---
rule assemble_results:
    input:
        script  = "scripts/aggregate/assemble_results.R",
        comp    = "results/aggregate/composition_test_m02day2.tsv",
        degs    = "results/aggregate/degs_pseudobulk_m02day2.tsv",
        gsea    = "results/aggregate/gsea_pseudobulk_m02day2.tsv",
        pathway = "results/aggregate/pathway_test_m02day2.tsv",
        niche   = "results/aggregate/niche_test_m02day2.tsv",
        mixing  = "results/aggregate/spatial_mixing_test_m02day2.tsv",
        myeloid = "results/aggregate/myeloid_m1m2_test_m02day2.tsv",
        mde     = "results/aggregate/power_mde.tsv",
        yaml    = "config/pathway_gene_lists.yaml",
    output:
        master = "results/aggregate/results_master.tsv",
    threads: 1
    log:
        "logs/aggregate/assemble_results.log",
    shell:
        "{RSCRIPT} {input.script} {input.comp} {input.degs} {input.gsea} "
        "{input.pathway} {input.niche} {input.mixing} {input.myeloid} "
        "{input.mde} {input.yaml} {output.master} > {log} 2>&1"
