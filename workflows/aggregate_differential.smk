# spatial-rads -- cross-sample differential workflow (MBRT vs SBRT vs Control, day 2).
# Consumes the typing workflow's unified per-cell labels + merged counts + obs table
# and runs the differential layer: composition (propeller), pseudobulk DE + GSEA, UCell
# pathway scores, spatial niches, tumor-immune mixing, myeloid polarization, M01<->M02
# concordance, and QC confound checks. Terminus = results_master.tsv (tier-tagged master).
# Requires aggregate_typing.smk to have produced its three handoff outputs first.
# Run: TMPDIR=/mnt/data/projects/spatial-rads/tmp \
#        conda run -n basecamp snakemake -s workflows/aggregate_differential.smk --cores N
# R steps run in the spatial-rads env via conda run (driver lives in basecamp).
import pandas as pd
configfile: "config/config.yaml"

# --- interpreter (driver runs in basecamp; R steps in spatial-rads) ---
RSCRIPT = "conda run -n spatial-rads Rscript"
# muscat is not conda/BiocManager-installable on this host; the differential-detection
# step runs in the prebuilt biocontainer via apptainer (bind the data disk).
MUSCAT_SIF = f"{config['datadir']}/containers/muscat.sif"

# --- paths: D_ data dirs, R_ repo locations ---
D_DATA    = config["datadir"]
D_AGG     = f"{D_DATA}/aggregate"                  # heavy intermediates
D_FULL    = f"{D_AGG}/full"                        # typing outputs consumed here
D_NORM    = f"{D_DATA}/processing/norm"            # per-sample inputs (coords/necrosis)
D_RES     = "results/aggregate"                    # small TSV/PNG reports
D_LOGS    = "logs/aggregate"
R_SCRIPTS = "scripts/aggregate"

MASTER = config["samplesheet"]                     # results/data_model/samples.tsv

# --- typing handoff: leaf inputs from aggregate_typing.smk ---
LABELS = "results/aggregate/full_labels.parquet"
MERGED = f"{D_AGG}/merged.rds"
OBS    = f"{D_FULL}/obs.parquet"
LATENT = f"{D_FULL}/scvi_latent.parquet"           # tier-1 scVI latent, reused for neighborhood-purity QC

# --- GPU-env Python interpreter (neighborhood-purity QC recomputes kNN on the latent) ---
PYSCVI = "conda run -n spatial-rads-scvi python"

# --- cohort: flank only ---
_s    = pd.read_csv(MASTER, sep="\t")
FLANK = _s.loc[_s["model"] == "flank", "sample_id"].tolist()
DE_COHORTS = ["mutter02_day2", "combined_4h_treated", "combined_4h", "mutter02_day2_pooledctrl", "combined_4h_pooledctrl"]
LM_COHORTS = ["mutter02_day2", "combined_4h_treated", "combined_4h", "mutter02_day2_pooledctrl", "combined_4h_pooledctrl"]
# smiDE model modes: "screen" is the nebula NB GLMM discovery pass, "spatial" is the
# per-spatial-unit GP_Matern fit plus inverse-variance meta-analysis.
SMIDE_MODES = ["screen", "spatial"]
rule all:
    input:
        "results/aggregate/composition_by_sample.tsv",
        "results/aggregate/composition_test_m02day2.tsv",
        "results/aggregate/composition_unassigned_sensitivity.tsv",
        expand("results/aggregate/engine/composition_engine_{coh}.tsv", coh=LM_COHORTS),
        expand("results/aggregate/engine/de_engine_{coh}.tsv", coh=DE_COHORTS),
        expand("results/aggregate/engine/niche_engine_{coh}.tsv", coh=LM_COHORTS),
        expand("results/aggregate/engine/mixing_engine_{coh}.tsv", coh=LM_COHORTS),
        expand("results/aggregate/engine/myeloid_engine_{coh}.tsv", coh=LM_COHORTS),
        expand("results/aggregate/engine/substate_engine_{coh}.tsv", coh=LM_COHORTS),
        "results/aggregate/plots/composition_m02day2_bars.png",
        "results/aggregate/plots/composition_m02day2_forest.png",
        "results/aggregate/plots/composition_m01_timecourse.png",
        "results/aggregate/pseudobulk_qc.tsv",
        "results/aggregate/differential_detection.tsv",
        expand("results/aggregate/gsea_pseudobulk_{coh}.tsv", coh=DE_COHORTS),
        "results/aggregate/pathway_scores_summary.tsv",
        expand("results/aggregate/pathway_test_{coh}.tsv", coh=LM_COHORTS),
        "results/aggregate/pathway_ucell_ams_concordance.tsv",
        "results/aggregate/plots/pathway_heatmap_m02day2.png",
        "results/aggregate/plots/pathway_timecourse_m01.png",
        "results/aggregate/plots/pathway_ucell_vs_ams_scatter.png",
        "results/aggregate/readout_detection_m02.tsv",
        "results/aggregate/power_mde.tsv",
        "results/aggregate/celltype_qc_markers.tsv",
        "results/aggregate/plots/celltype_qc_dotplot.png",
        "results/aggregate/niche_frequency.tsv",
        "results/aggregate/plots/niche_frequency_m02.png",
        "results/aggregate/plots/niche_centroids_heatmap.png",
        "results/aggregate/spatial_mixing_per_sample.tsv",
        "results/aggregate/plots/mixing_m02.png",
        "results/aggregate/myeloid_m1m2_scores.tsv",
        "results/aggregate/plots/myeloid_m1m2_m02.png",
        "results/aggregate/concordance_m01_m02.tsv",
        "results/aggregate/plots/concordance_scatter.png",
        "results/aggregate/geneset_overlap.tsv",
        "results/aggregate/fibroblast_substate.parquet",
        "results/aggregate/substate_gate_report.tsv",
        "results/aggregate/detectability_summary.tsv",
        "results/aggregate/results_master.tsv",
        "results/aggregate/overlap_ratio_qc.tsv",
        "results/aggregate/qc_arm_balance.tsv",
        "results/aggregate/qc_arm_balance_samples.tsv",
        "results/aggregate/qc_reproducibility.tsv",
        "results/aggregate/qc_panel_sparsity.tsv",
        "results/aggregate/qc_fov_signal.tsv",
        "results/aggregate/celltype_neighborhood_purity.tsv",
        expand("results/aggregate/engine/smide_de_{coh}_{mode}.tsv", coh=LM_COHORTS, mode=SMIDE_MODES),
        expand("results/aggregate/smide_concordance_{coh}.tsv", coh=LM_COHORTS),
        expand("results/aggregate/plots/smide_concordance_{coh}.png", coh=LM_COHORTS),
        "results/aggregate/insitucor_modules.tsv",
        "results/aggregate/insitucor_module_summary.tsv",
# --- Track 1: cell-type composition (M02 day2 propeller test + M01 descriptive) ---
rule composition:
    message: "composition: cell-type composition M02 day2 propeller test + M01 descriptive"
    input:
        script = f"{R_SCRIPTS}/composition.R",
        obs    = OBS,
        labels = LABELS,
    output:
        by_sample  = "results/aggregate/composition_by_sample.tsv",
        test       = "results/aggregate/composition_test_m02day2.tsv",
        dropped    = "results/aggregate/composition_dropped_celltypes.tsv",
        bars       = "results/aggregate/plots/composition_m02day2_bars.png",
        forest     = "results/aggregate/plots/composition_m02day2_forest.png",
        timecourse = "results/aggregate/plots/composition_m01_timecourse.png",
        sensitivity = "results/aggregate/composition_unassigned_sensitivity.tsv",
        lm_input   = "results/aggregate/engine_inputs/composition_cells.tsv",
    threads: 1
    log:
        f"{D_LOGS}/composition.log",
    shell:
        "{RSCRIPT} {input.script} {input.obs} {input.labels} {output.by_sample} "
        "{output.test} {output.dropped} {output.bars} {output.forest} "
        "{output.timecourse} {output.sensitivity} {output.lm_input} > {log} 2>&1"
# --- engine: composition arm test (propeller/logit path) on the per-cell labels ---
rule composition_engine:
    message: "composition_engine: lm_engine proportion test ({wildcards.cohort})"
    input:
        script = "scripts/engines/lm_engine.R",
        cells  = "results/aggregate/engine_inputs/composition_cells.tsv",
        comp   = "results/data_model/comparisons.tsv",
        params = "config/engine_params.yaml",
    output:
        stats = "results/aggregate/engine/composition_engine_{cohort}.tsv",
    threads: 1
    log:
        f"{D_LOGS}/composition_engine_{{cohort}}.log",
    shell:
        "{RSCRIPT} {input.script} {input.cells} proportion {wildcards.cohort} composition "
        "{input.comp} {input.params} {output.stats} > {log} 2>&1"
# --- Track 2: pseudobulk construction (all flank, sample x cell_type count sums) ---
rule pseudobulk_build:
    message: "pseudobulk_build: sample x cell_type pseudobulk count sums (all flank)"
    input:
        script = f"{R_SCRIPTS}/pseudobulk_build.R",
        rds    = MERGED,
        labels = LABELS,
        ss     = MASTER,
    output:
        se = f"{D_AGG}/pseudobulk_se.rds",
        qc = "results/aggregate/pseudobulk_qc.tsv",
    threads: 1
    log:
        f"{D_LOGS}/pseudobulk_build.log",
    shell:
        "{RSCRIPT} {input.script} {input.rds} {input.labels} {input.ss} {output.se} {output.qc} > {log} 2>&1"
# --- Track 2 inference: pseudobulk NB DE via the count engine (M02 day2, abundance-floored,
# apeglm two-fit). Sole pseudobulk-DE source: feeds results_master, gsea, and concordance.
# (Superseded the standalone deg_pseudobulk.R -- retired 2026-07-13, output bit-identical.) ---
rule de_engine:
    message: "de_engine: count_engine pseudobulk NB DE ({wildcards.cohort}) -> sufficient stats"
    input:
        script = "scripts/engines/count_engine.R",
        se     = f"{D_AGG}/pseudobulk_se.rds",
        comp   = "results/data_model/comparisons.tsv",
        params = "config/engine_params.yaml",
    output:
        stats   = "results/aggregate/engine/de_engine_{cohort}.tsv",
        skipped = "results/aggregate/engine/de_engine_{cohort}_skipped.tsv",
    threads: 1
    log:
        f"{D_LOGS}/de_engine_{{cohort}}.log",
    shell:
        "{RSCRIPT} {input.script} {input.se} {wildcards.cohort} {input.comp} {input.params} "
        "{output.stats} {output.skipped} > {log} 2>&1"
# --- Track 3 differential detection: fraction-of-cells-expressing change (muscat) ---
# dd_prep exports an SCE from merged.rds (main env has Seurat/arrow); dd_muscat runs
# muscat's num.detected aggregation + edgeR pbDS INSIDE the biocontainer (aptainer).
rule dd_prep:
    message: "dd_prep: merged.rds -> SCE for muscat differential detection (M02 day2)"
    input:
        script  = f"{R_SCRIPTS}/dd_prep.R",
        merged  = f"{D_FULL}/merged.rds",
        samples = MASTER,
    output:
        sce = f"{D_FULL}/dd_sce.rds",
    threads: 1
    log:
        f"{D_LOGS}/dd_prep.log",
    shell:
        "{RSCRIPT} {input.script} {input.merged} {input.samples} {output.sce} > {log} 2>&1"

rule dd_muscat:
    message: "dd_muscat: muscat differential detection (num.detected + edgeR pbDS) in container"
    input:
        script = f"{R_SCRIPTS}/dd_muscat.R",
        sce    = f"{D_FULL}/dd_sce.rds",
    output:
        tsv = "results/aggregate/differential_detection.tsv",
    threads: 1
    log:
        f"{D_LOGS}/dd_muscat.log",
    shell:
        "apptainer exec --bind {D_DATA} {MUSCAT_SIF} "
        "Rscript {input.script} {input.sce} {output.tsv} > {log} 2>&1"
# --- Track 2 pathway: GSEA on pseudobulk stat-ranked genes (primary + Hallmark) ---
rule gsea:
    message: "gsea: GSEA on pseudobulk stat-ranked genes ({wildcards.cohort})"
    input:
        script = f"{R_SCRIPTS}/gsea.R",
        degs   = "results/aggregate/engine/de_engine_{cohort}.tsv",
        sets   = "results/data_model/pathway_sets.tsv",
    output:
        gsea = "results/aggregate/gsea_pseudobulk_{cohort}.tsv",
    threads: 1
    log:
        f"{D_LOGS}/gsea_{{cohort}}.log",
    shell:
        "{RSCRIPT} {input.script} {input.degs} {input.sets} {output.gsea} > {log} 2>&1"
# --- Track 2 pathway: per-cell UCell + AddModuleScore scoring (cohort-agnostic compute) ---
# threads: 4 -- UCell runs ncores=4 fork (BiocParallel) parallelism, not BLAS.
rule pathway_scores:
    message: "pathway_scores: per-cell UCell + AddModuleScore scoring"
    input:
        script = f"{R_SCRIPTS}/pathway_scores.R",
        rds    = MERGED,   # count-only consumer; labels from full_labels.parquet
        labels = LABELS,
        sets   = "results/data_model/pathway_sets.tsv",
    output:
        summary = "results/aggregate/pathway_scores_summary.tsv",
        conc    = "results/aggregate/pathway_ucell_ams_concordance.tsv",
    threads: 4
    log:
        f"{D_LOGS}/pathway_scores.log",
    shell:
        "{RSCRIPT} {input.script} {input.rds} {input.labels} {input.sets} "
        "{output.summary} {output.conc} > {log} 2>&1"
# Per-cohort arm test: reads the cached summary TSV, cheap, parameterized by the comparison
# registry (no hardcoded day-2 literals) so it runs for every LM_COHORTS entry.
rule pathway_arm_test:
    message: "pathway_arm_test: UCell pathway arm test ({wildcards.cohort})"
    input:
        script  = f"{R_SCRIPTS}/pathway_arm_test.R",
        summary = "results/aggregate/pathway_scores_summary.tsv",
        comp    = "results/data_model/comparisons.tsv",
    output:
        test = "results/aggregate/pathway_test_{cohort}.tsv",
    threads: 1
    log:
        f"{D_LOGS}/pathway_arm_test_{{cohort}}.log",
    shell:
        "{RSCRIPT} {input.script} {input.summary} {input.comp} {wildcards.cohort} "
        "{output.test} > {log} 2>&1"
# Plot half: reads the cached score/test TSVs so a plot bug can't roll back the ~5h compute.
rule pathway_plots:
    message: "pathway_plots: pathway heatmap/timecourse/concordance plots from cached scores"
    input:
        script  = f"{R_SCRIPTS}/pathway_plots.R",
        summary = "results/aggregate/pathway_scores_summary.tsv",
        test    = "results/aggregate/pathway_test_mutter02_day2.tsv",
    output:
        heatmap    = "results/aggregate/plots/pathway_heatmap_m02day2.png",
        timecourse = "results/aggregate/plots/pathway_timecourse_m01.png",
        scatter    = "results/aggregate/plots/pathway_ucell_vs_ams_scatter.png",
    threads: 1
    log:
        f"{D_LOGS}/pathway_plots.log",
    shell:
        "{RSCRIPT} {input.script} {input.summary} {input.test} "
        "{output.heatmap} {output.timecourse} {output.scatter} > {log} 2>&1"
# ============================================================================
# MBRT-vs-SBRT downstream differential layer (plan-mbrt-vs-sbrt-impl.md T4-T13).
# All consume the unified per-cell labels (results/aggregate/full_labels.parquet)
# and the day-2 readout tables, terminating in the tier-tagged results_master.tsv.
# ============================================================================
# --- T5: per-cell-type readout detection report (panel detectability in M02) ---
rule panel_coverage:
    message: "panel_coverage: per-cell-type readout detection report (panel detectability in M02)"
    input:
        script = f"{R_SCRIPTS}/panel_coverage.R",
        rds    = MERGED,   # re-pointed off the deleted merged_typed.rds (count-only consumer; labels from full_labels.parquet)
        labels = LABELS,
        sets   = "results/data_model/pathway_sets.tsv",
    output:
        detection = "results/aggregate/readout_detection_m02.tsv",
    threads: 1
    log:
        f"{D_LOGS}/panel_coverage.log",
    shell:
        "{RSCRIPT} {input.script} {input.rds} {input.labels} {input.sets} "
        "{output.detection} > {log} 2>&1"
# --- T6: power / minimum-detectable-effect table at n=4 ---
rule power_mde:
    message: "power_mde: power / minimum-detectable-effect table at n=4"
    input:
        script = f"{R_SCRIPTS}/power_mde.R",
        comp   = "results/aggregate/composition_by_sample.tsv",
        se     = f"{D_AGG}/pseudobulk_se.rds",
        path   = "results/aggregate/pathway_scores_summary.tsv",
    output:
        mde = "results/aggregate/power_mde.tsv",
    threads: 1
    log:
        f"{D_LOGS}/power_mde.log",
    shell:
        "{RSCRIPT} {input.script} {input.comp} {input.se} {input.path} {output.mde} > {log} 2>&1"
# --- T7: per-cell slide coords + per-sample necrosis flag (spatial-track input) ---
rule coords_necrosis:
    message: "coords_necrosis: per-cell slide coords + per-sample necrosis flag"
    input:
        script = f"{R_SCRIPTS}/coords_necrosis.R",
        ss     = MASTER,
        rds    = expand(f"{D_NORM}/{{s}}.norm.rds", s=FLANK),
    output:
        coords = f"{D_AGG}/coords_necrosis.parquet",
    threads: 1
    log:
        f"{D_LOGS}/coords_necrosis.log",
    shell:
        "{RSCRIPT} {input.script} {input.ss} {output.coords} {input.rds} > {log} 2>&1"
# --- T8: cell-type-label QC marker dotplot (does each label express its markers) ---
rule celltype_qc:
    message: "celltype_qc: cell-type-label QC marker dotplot"
    input:
        script = f"{R_SCRIPTS}/celltype_qc.R",
        rds    = MERGED,
        labels = LABELS,
    output:
        markers = "results/aggregate/celltype_qc_markers.tsv",
        dotplot = "results/aggregate/plots/celltype_qc_dotplot.png",
    threads: 1
    log:
        f"{D_LOGS}/celltype_qc.log",
    shell:
        "{RSCRIPT} {input.script} {input.rds} {input.labels} {output.markers} "
        "{output.dotplot} > {log} 2>&1"
# --- T9: data-driven spatial niches (k=50 NN composition, SPIN 50%, Mclust BIC K selection) ---
rule niches:
    message: "niches: data-driven spatial niches (k=50 NN composition, SPIN 50%, Mclust BIC K selection)"
    input:
        script = f"{R_SCRIPTS}/niches.R",
        labels = LABELS,
        coords = f"{D_AGG}/coords_necrosis.parquet",
        obs    = OBS,
    output:
        per_cell  = f"{D_AGG}/niche_per_cell.parquet",
        centroids = "results/aggregate/niche_centroids.tsv",
        freq      = "results/aggregate/niche_frequency.tsv",
        lm_input  = "results/aggregate/engine_inputs/niche_cells.tsv",
        heatmap   = "results/aggregate/plots/niche_centroids_heatmap.png",
        freq_plot = "results/aggregate/plots/niche_frequency_m02.png",
    threads: 1
    log:
        f"{D_LOGS}/niches.log",
    shell:
        "{RSCRIPT} {input.script} {input.labels} {input.coords} {input.obs} "
        "{output.per_cell} {output.centroids} {output.freq} {output.lm_input} "
        "{output.heatmap} {output.freq_plot} > {log} 2>&1"
# --- engine: niche frequency arm test (propeller/logit path) ---
rule niche_engine:
    message: "niche_engine: lm_engine proportion test on niche frequency ({wildcards.cohort})"
    input:
        script = "scripts/engines/lm_engine.R",
        cells  = "results/aggregate/engine_inputs/niche_cells.tsv",
        comp   = "results/data_model/comparisons.tsv",
        params = "config/engine_params.yaml",
    output:
        stats = "results/aggregate/engine/niche_engine_{cohort}.tsv",
    threads: 1
    log:
        f"{D_LOGS}/niche_engine_{{cohort}}.log",
    shell:
        "{RSCRIPT} {input.script} {input.cells} proportion {wildcards.cohort} niche "
        "{input.comp} {input.params} {output.stats} > {log} 2>&1"
# --- T10: tumor-immune spatial mixing (k=50 NN, immune-neighbour fraction + Keren score) ---
rule spatial_mixing:
    message: "spatial_mixing: tumor-immune spatial mixing (k=50 NN)"
    input:
        script = f"{R_SCRIPTS}/spatial_mixing.R",
        labels = LABELS,
        coords = f"{D_AGG}/coords_necrosis.parquet",
        obs    = OBS,
    output:
        per_sample = "results/aggregate/spatial_mixing_per_sample.tsv",
        lm_input   = "results/aggregate/engine_inputs/mixing_metrics.tsv",
        per_cell   = f"{D_AGG}/spatial_mixing_per_cell.parquet",
        plot       = "results/aggregate/plots/mixing_m02.png",
    threads: 1
    log:
        f"{D_LOGS}/spatial_mixing.log",
    shell:
        "{RSCRIPT} {input.script} {input.labels} {input.coords} {input.obs} "
        "{output.per_sample} {output.lm_input} {output.per_cell} {output.plot} > {log} 2>&1"
# --- engine: tumor-immune mixing arm test (matrix/identity path, non-robust) ---
rule mixing_engine:
    message: "mixing_engine: lm_engine identity test on spatial mixing metrics ({wildcards.cohort})"
    input:
        script = "scripts/engines/lm_engine.R",
        matrix = "results/aggregate/engine_inputs/mixing_metrics.tsv",
        comp   = "results/data_model/comparisons.tsv",
        params = "config/engine_params.yaml",
    output:
        stats = "results/aggregate/engine/mixing_engine_{cohort}.tsv",
    threads: 1
    log:
        f"{D_LOGS}/mixing_engine_{{cohort}}.log",
    shell:
        "{RSCRIPT} {input.script} {input.matrix} matrix {wildcards.cohort} mixing "
        "{input.comp} {input.params} {output.stats} > {log} 2>&1"
# --- T11: myeloid M1/M2 polarization (UCell M1/M2 panels -> per-sample ratio) ---
# threads: 8 -- AddModuleScore_UCell forks UCELL_CORES=8 workers (BiocParallel, not
# BLAS), so reserve 8 to avoid oversubscription.
rule myeloid_polarization:
    message: "myeloid_polarization: myeloid M1/M2 polarization (UCell M1/M2 panels -> per-sample ratio)"
    input:
        script  = f"{R_SCRIPTS}/myeloid_polarization.R",
        rds     = MERGED,
        labels  = LABELS,
        samples = MASTER,
    output:
        scores = "results/aggregate/myeloid_m1m2_scores.tsv",
        lm_input = "results/aggregate/engine_inputs/myeloid_metrics.tsv",
        plot   = "results/aggregate/plots/myeloid_m1m2_m02.png",
    threads: 8
    log:
        f"{D_LOGS}/myeloid_polarization.log",
    shell:
        "{RSCRIPT} {input.script} {input.rds} {input.labels} {input.samples} "
        "{output.scores} {output.lm_input} {output.plot} > {log} 2>&1"
# --- engine: myeloid M1/M2 arm test (matrix/identity path, non-robust) ---
rule myeloid_engine:
    message: "myeloid_engine: lm_engine identity test on macrophage M1/M2 metrics ({wildcards.cohort})"
    input:
        script = "scripts/engines/lm_engine.R",
        matrix = "results/aggregate/engine_inputs/myeloid_metrics.tsv",
        comp   = "results/data_model/comparisons.tsv",
        params = "config/engine_params.yaml",
    output:
        stats = "results/aggregate/engine/myeloid_engine_{cohort}.tsv",
    threads: 1
    log:
        f"{D_LOGS}/myeloid_engine_{{cohort}}.log",
    shell:
        "{RSCRIPT} {input.script} {input.matrix} matrix {wildcards.cohort} myeloid "
        "{input.comp} {input.params} {output.stats} > {log} 2>&1"
# --- T12: M01<->M02 day-2 effect-size concordance on the unified labels ---
rule concordance_m01_m02:
    message: "concordance_m01_m02: M01<->M02 day-2 effect-size concordance on the unified labels"
    input:
        script = f"{R_SCRIPTS}/concordance_m01_m02.R",
        rds    = MERGED,
        labels = LABELS,
        degs   = "results/aggregate/engine/de_engine_mutter02_day2.tsv",
    output:
        tsv     = "results/aggregate/concordance_m01_m02.tsv",
        scatter = "results/aggregate/plots/concordance_scatter.png",
    threads: 1
    log:
        f"{D_LOGS}/concordance_m01_m02.log",
    shell:
        "{RSCRIPT} {input.script} {input.rds} {input.labels} {input.degs} "
        "{output.tsv} {output.scatter} > {log} 2>&1"
# --- T13: tier-tagged master results table with confirmatory-family FDR ---
rule substate_split:
    message: "substate_split: fibroblast resting/activated substate split with marker gate"
    input:
        script  = f"{R_SCRIPTS}/substate_split.R",
        rds     = MERGED,
        labels  = LABELS,
        markers = "config/substate_markers.yaml",
    output:
        parquet = "results/aggregate/fibroblast_substate.parquet",
        gate    = "results/aggregate/substate_gate_report.tsv",
    threads: 4
    log:
        f"{D_LOGS}/substate_split.log",
    shell:
        "{RSCRIPT} {input.script} {input.rds} {input.labels} {input.markers} "
        "{output.parquet} {output.gate} > {log} 2>&1"
rule substate_composition:
    message: "substate_composition: build fibroblast-split per-cell engine input (M02 day2)"
    input:
        script   = f"{R_SCRIPTS}/substate_composition.R",
        obs      = OBS,
        labels   = LABELS,
        substate = "results/aggregate/fibroblast_substate.parquet",
    output:
        lm_input = "results/aggregate/engine_inputs/substate_cells.tsv",
    threads: 1
    log:
        f"{D_LOGS}/substate_composition.log",
    shell:
        "{RSCRIPT} {input.script} {input.obs} {input.labels} {input.substate} {output.lm_input} > {log} 2>&1"
# --- engine: fibroblast sub-state composition arm test (propeller/logit path) ---
rule substate_engine:
    message: "substate_engine: lm_engine proportion test on fibroblast sub-state composition ({wildcards.cohort})"
    input:
        script = "scripts/engines/lm_engine.R",
        cells  = "results/aggregate/engine_inputs/substate_cells.tsv",
        comp   = "results/data_model/comparisons.tsv",
        params = "config/engine_params.yaml",
    output:
        stats = "results/aggregate/engine/substate_engine_{cohort}.tsv",
    threads: 1
    log:
        f"{D_LOGS}/substate_engine_{{cohort}}.log",
    shell:
        "{RSCRIPT} {input.script} {input.cells} proportion {wildcards.cohort} substate "
        "{input.comp} {input.params} {output.stats} > {log} 2>&1"
rule geneset_overlap:
    message: "geneset_overlap: pairwise gene-set overlap for the curated pathway sets"
    input:
        script = f"{R_SCRIPTS}/geneset_overlap.R",
        sets   = "results/data_model/pathway_sets.tsv",
    output:
        overlap = "results/aggregate/geneset_overlap.tsv",
    threads: 1
    log:
        f"{D_LOGS}/geneset_overlap.log",
    shell:
        "{RSCRIPT} {input.script} {input.sets} {output.overlap} > {log} 2>&1"
rule assemble_results:
    message: "assemble_results: tier-tagged master results table with confirmatory-family FDR"
    input:
        script  = f"{R_SCRIPTS}/assemble_results.R",
        # Declared expand()-ed lists give Snakemake DAG tracking over the full cohort
        # roster; assemble_results.R itself globs the per-cohort filename pattern out of
        # results/aggregate (and results/aggregate/engine) rather than taking one path
        # per cohort, so this list only needs to be complete, not individually threaded
        # through to the shell command.
        comp     = expand("results/aggregate/engine/composition_engine_{coh}.tsv", coh=LM_COHORTS),
        degs     = expand("results/aggregate/engine/de_engine_{coh}.tsv", coh=DE_COHORTS),
        gsea     = expand("results/aggregate/gsea_pseudobulk_{coh}.tsv", coh=DE_COHORTS),
        pathway  = expand("results/aggregate/pathway_test_{coh}.tsv", coh=LM_COHORTS),
        niche    = expand("results/aggregate/engine/niche_engine_{coh}.tsv", coh=LM_COHORTS),
        mixing   = expand("results/aggregate/engine/mixing_engine_{coh}.tsv", coh=LM_COHORTS),
        myeloid  = expand("results/aggregate/engine/myeloid_engine_{coh}.tsv", coh=LM_COHORTS),
        substate = expand("results/aggregate/engine/substate_engine_{coh}.tsv", coh=LM_COHORTS),
        smide    = expand("results/aggregate/engine/smide_de_{coh}_{mode}.tsv", coh=LM_COHORTS, mode=SMIDE_MODES),
        reg      = "results/data_model/comparisons.tsv",
        mde      = "results/aggregate/power_mde.tsv",
        sets     = "results/data_model/pathway_sets.tsv",
        cov      = "results/data_model/gene_set_panel_coverage.tsv",
        dd       = "results/aggregate/differential_detection.tsv",
        overlap  = "results/aggregate/overlap_ratio_qc.tsv",
    output:
        master = "results/aggregate/results_master.tsv",
        detect = "results/aggregate/detectability_summary.tsv",
    threads: 1
    log:
        f"{D_LOGS}/assemble_results.log",
    shell:
        "{RSCRIPT} {input.script} results/aggregate {input.reg} {input.mde} {input.sets} "
        "{input.cov} {input.dd} {input.overlap} {output.master} > {log} 2>&1"
# ============================================================================
# Spatial discovery: de novo spatial gene-gene co-expression modules, independent
# of the curated pathway sets. Discovery input for the planned peak/valley signature
# work. Runs on 4h timepoint cells only (all datasets, all arms).
# ============================================================================
rule insitucor_discovery:
    message: "insitucor_discovery: spatial gene-correlation modules on 4h cells (InSituCor k=100)"
    input:
        script = f"{R_SCRIPTS}/insitucor_discovery.R",
        rds    = MERGED,
        labels = LABELS,
        coords = f"{D_AGG}/coords_necrosis.parquet",
        obs    = OBS,
    output:
        modules = "results/aggregate/insitucor_modules.tsv",
        summary = "results/aggregate/insitucor_module_summary.tsv",
    threads: 1
    log:
        f"{D_LOGS}/insitucor_discovery.log",
    shell:
        "{RSCRIPT} {input.script} {input.rds} {input.labels} {input.coords} {input.obs} "
        "{output.modules} {output.summary} > {log} 2>&1"
# --- QC: cross-arm balance -- confound check on the day-2 composition result. Joins the per-sample
# technical metrics + MECR (both from preprocessing.smk) to the arm design; balanced arms => the
# fraction shift is not a sensitivity/contamination artifact. ---
rule agg_qc_arm_balance:
    message: "agg_qc_arm_balance: cross-arm technical/contamination balance confound check"
    input:
        script  = f"{R_SCRIPTS}/qc_arm_balance.R",
        tech    = "results/processing/sample_tech_metrics.tsv",
        contam  = "results/processing/contamination_qc.tsv",
        samples = MASTER,
    output:
        tsv     = "results/aggregate/qc_arm_balance.tsv",
        samples = "results/aggregate/qc_arm_balance_samples.tsv",
    threads: 1
    log:
        f"{D_LOGS}/qc_arm_balance.log",
    shell:
        "{RSCRIPT} {input.script} {input.tech} {input.contam} {input.samples} "
        "{output.tsv} {output.samples} > {log} 2>&1"
# --- QC: smiDE overlap ratio -- per-gene, per-cell_subtype segmentation-contamination
# metric (avg self expression vs avg neighbor-othertype expression within radius). Joined
# onto results_master as contamination_ratio so a DE hit's segmentation risk is legible
# alongside its effect size. ---
rule overlap_ratio_qc:
    message: "overlap_ratio_qc: smiDE per-gene per-cell_subtype segmentation contamination metric"
    input:
        script = f"{R_SCRIPTS}/overlap_ratio_qc.R",
        rds    = MERGED,
        labels = LABELS,
        coords = f"{D_AGG}/coords_necrosis.parquet",
    output:
        tsv = "results/aggregate/overlap_ratio_qc.tsv",
    threads: 1
    log:
        f"{D_LOGS}/overlap_ratio_qc.log",
    shell:
        "{RSCRIPT} {input.script} {input.rds} {input.labels} {input.coords} "
        "{output.tsv} > {log} 2>&1"
# --- genome-wide smiDE per-cell DE (cohort- and mode-parameterized) ---
# Follows the smiDE authors' published protocol: overlap-ratio gene prefilter per
# cell type, adjacency-only pre_de with on-the-fly neighbor expression, inverse-
# distance edge weighting, totalcount-normalized neighbor covariate. Two modes ---
# "screen" (nebula NB GLMM with a sample random intercept; discovery only) and
# "spatial" (per-spatial-unit GP_Matern fits combined by inverse-variance meta-
# analysis; the authors' fix for inflated type I error). Depends on the overlap
# ratio QC table for the prefilter. Runs BEFORE assemble_results, which maps
# screen -> readout_class "smiDE_screen" and spatial meta -> "smiDE".
rule smide_de:
    message: "smide_de: genome-wide per-cell smiDE DE ({wildcards.cohort}, {wildcards.mode})"
    wildcard_constraints:
        mode = "screen|spatial",
    input:
        script  = f"{R_SCRIPTS}/smide_de.R",
        rds     = MERGED,
        labels  = LABELS,
        coords  = f"{D_AGG}/coords_necrosis.parquet",
        obs     = OBS,
        comp    = "results/data_model/comparisons.tsv",
        samples = MASTER,
        overlap = "results/aggregate/overlap_ratio_qc.tsv",
    output:
        tsv = "results/aggregate/engine/smide_de_{cohort}_{mode}.tsv",
    threads: 16
    log:
        f"{D_LOGS}/smide_de_{{cohort}}_{{mode}}.log",
    shell:
        "{RSCRIPT} {input.script} {input.rds} {input.labels} {input.coords} {input.obs} "
        "{input.comp} {input.samples} {input.overlap} {wildcards.cohort} {wildcards.mode} "
        "{output.tsv} > {log} 2>&1"
# --- QC: pseudobulk vs smiDE concordance -- effect-size correlation, hit overlap,
# contamination-ratio enrichment in discordant hits, confirmatory hit survival. Runs
# after both assemble_results (results_master.tsv) and overlap_ratio_qc, scoped to
# one cohort at a time (results_master.tsv rows are filtered by `comparison`). ---
rule smide_concordance:
    message: "smide_concordance: pseudobulk vs smiDE effect-size concordance ({wildcards.cohort})"
    input:
        script  = f"{R_SCRIPTS}/smide_concordance.R",
        master  = "results/aggregate/results_master.tsv",
        overlap = "results/aggregate/overlap_ratio_qc.tsv",
    output:
        tsv  = "results/aggregate/smide_concordance_{cohort}.tsv",
        plot = "results/aggregate/plots/smide_concordance_{cohort}.png",
    threads: 1
    log:
        f"{D_LOGS}/smide_concordance_{{cohort}}.log",
    shell:
        "{RSCRIPT} {input.script} {input.master} {input.overlap} {wildcards.cohort} "
        "{output.tsv} {output.plot} > {log} 2>&1"
# --- QC: replicate reproducibility -- per-arm pseudobulk concordance (SpatialQM getCorrelation) +
# technical-metric PCA over the n=4/arm M02 day-2 cohort; flags an outlier slide driving an arm. ---
rule agg_qc_reproducibility:
    message: "agg_qc_reproducibility: per-arm pseudobulk concordance + technical-metric PCA reproducibility check"
    input:
        script = f"{R_SCRIPTS}/qc_reproducibility.R",
        se     = f"{D_AGG}/pseudobulk_se.rds",
        tech   = "results/processing/sample_tech_metrics.tsv",
    output:
        tsv = "results/aggregate/qc_reproducibility.tsv",
    threads: 1
    log:
        f"{D_LOGS}/qc_reproducibility.log",
    shell:
        "{RSCRIPT} {input.script} {input.se} {input.tech} {output.tsv} > {log} 2>&1"
# --- QC: per-cell panel detection -- the platform-sparsity fact motivating cluster-then-annotate
# typing (per-cell InSituType unreachable at ~6% panel detection). Joins per-cell nFeature (merged.rds)
# to the unified compartment labels + arm/dataset design -> box quantiles by compartment. ---
rule agg_qc_panel_sparsity:
    message: "agg_qc_panel_sparsity: per-cell panel detection by compartment/dataset/arm"
    input:
        script = f"{R_SCRIPTS}/qc_panel_sparsity.R",
        merged = MERGED,
        labels = LABELS,
        samples = MASTER,
    output:
        tsv = "results/aggregate/qc_panel_sparsity.tsv",
    threads: 1
    log:
        f"{D_LOGS}/qc_panel_sparsity.log",
    shell:
        "{RSCRIPT} {input.script} {input.merged} {input.labels} {input.samples} "
        "{output.tsv} > {log} 2>&1"

# --- QC: per-FOV signal-loss screening -- flags FOVs whose mean RNA count falls below
# 40% of the slide's median per-FOV mean (>60% signal loss), a CosMx best-practice check
# for degraded-signal regions (e.g. focal plane / staining dropout). Report-only: does not
# exclude cells. Verifies flagged FOVs do not concentrate in one treatment arm. ---
rule qc_fov_signal:
    message: "qc_fov_signal: per-FOV signal-loss screening (report-only)"
    input:
        script = f"{R_SCRIPTS}/qc_fov_signal.R",
        rds    = MERGED,
        labels = LABELS,
        obs    = OBS,
    output:
        tsv = "results/aggregate/qc_fov_signal.tsv",
    threads: 1
    log:
        f"{D_LOGS}/qc_fov_signal.log",
    shell:
        "{RSCRIPT} {input.script} {input.rds} {input.labels} {input.obs} "
        "{output.tsv} > {log} 2>&1"

# --- QC: transcriptional neighborhood purity (Plummer et al. Nat Biotechnol 2025, Fig 4d) --
# report-only coherence check on the locked typing. Recomputes kNN on the reused tier-1 scVI
# latent (no retyping, no GPU training) and scores per-cell purity by final cell_subtype. Runs
# in the scVI env for scanpy; plotting is the standalone fig_qc_neighborhood_purity.R block.
rule agg_celltype_purity:
    message: "agg_celltype_purity: transcriptional neighborhood purity by cell subtype (typing coherence QC)"
    input:
        script = f"{R_SCRIPTS}/celltype_neighborhood_purity.py",
        latent = LATENT,
        labels = LABELS,
    output:
        tsv = "results/aggregate/celltype_neighborhood_purity.tsv",
    threads: 1
    log:
        f"{D_LOGS}/celltype_neighborhood_purity.log",
    shell:
        "{PYSCVI} {input.script} {input.latent} {input.labels} {output.tsv} > {log} 2>&1"
