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

# --- paths: D_ data dirs, R_ repo locations ---
D_DATA    = config["datadir"]
D_AGG     = f"{D_DATA}/aggregate"                  # heavy intermediates
D_FULL    = f"{D_AGG}/full"                        # typing outputs consumed here
D_NORM    = f"{D_DATA}/processing/norm"            # per-sample inputs (coords/necrosis)
D_RES     = "results/aggregate"                    # small TSV/PNG reports
D_LOGS    = "logs/aggregate"
R_SCRIPTS = "scripts/aggregate"

MASTER = config["samplesheet"]                     # results/data_model/samples.tsv

# --- typing handoff: three pre-existing leaf inputs from aggregate_typing.smk ---
LABELS = "results/aggregate/full_labels.parquet"
MERGED = f"{D_AGG}/merged.rds"
OBS    = f"{D_FULL}/obs.parquet"

# --- cohort: flank only ---
_s    = pd.read_csv(MASTER, sep="\t")
FLANK = _s.loc[_s["model"] == "flank", "sample_id"].tolist()

rule all:
    input:
        "results/aggregate/composition_by_sample.tsv",
        "results/aggregate/composition_test_m02day2.tsv",
        "results/aggregate/composition_unassigned_sensitivity.tsv",
        "results/aggregate/engine/composition_engine.tsv",
        "results/aggregate/engine/de_engine.tsv",
        "results/aggregate/engine/niche_engine.tsv",
        "results/aggregate/engine/mixing_engine.tsv",
        "results/aggregate/engine/myeloid_engine.tsv",
        "results/aggregate/engine/substate_engine.tsv",
        "results/aggregate/plots/composition_m02day2_bars.png",
        "results/aggregate/plots/composition_m02day2_forest.png",
        "results/aggregate/plots/composition_m01_timecourse.png",
        "results/aggregate/pseudobulk_qc.tsv",
        "results/aggregate/engine/de_engine.tsv",
        "results/aggregate/detection_test_m02day2.tsv",
        "results/aggregate/gsea_pseudobulk_m02day2.tsv",
        "results/aggregate/pathway_scores_summary.tsv",
        "results/aggregate/pathway_test_m02day2.tsv",
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
        "results/aggregate/qc_arm_balance.tsv",
        "results/aggregate/plots/qc_arm_balance.png",
        "results/aggregate/qc_reproducibility.tsv",
        "results/aggregate/plots/qc_reproducibility.png",
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
    message: "composition_engine: lm_engine proportion test on cell-type composition (M02 day2)"
    input:
        script = "scripts/engines/lm_engine.R",
        cells  = "results/aggregate/engine_inputs/composition_cells.tsv",
        comp   = "results/data_model/comparisons.tsv",
        params = "config/engine_params.yaml",
    output:
        stats = "results/aggregate/engine/composition_engine.tsv",
    threads: 1
    log:
        f"{D_LOGS}/composition_engine.log",
    shell:
        "{RSCRIPT} {input.script} {input.cells} proportion mutter02_day2 composition "
        "{input.comp} {input.params} {output.stats} > {log} 2>&1"
# --- Track 2: pseudobulk construction (M02 day2, sample x cell_type count sums) ---
rule pseudobulk_build:
    message: "pseudobulk_build: sample x cell_type pseudobulk count sums (M02 day2)"
    input:
        script = f"{R_SCRIPTS}/pseudobulk_build.R",
        rds    = MERGED,   # re-pointed off the deleted merged_typed.rds (count-only consumer; labels from full_labels.parquet)
        labels = LABELS,
    output:
        se = f"{D_AGG}/pseudobulk_se.rds",
        qc = "results/aggregate/pseudobulk_qc.tsv",
    threads: 1
    log:
        f"{D_LOGS}/pseudobulk_build.log",
    shell:
        "{RSCRIPT} {input.script} {input.rds} {input.labels} {output.se} {output.qc} > {log} 2>&1"
# --- Track 2 inference: pseudobulk NB DE via the count engine (M02 day2, abundance-floored,
# apeglm two-fit). Sole pseudobulk-DE source: feeds results_master, gsea, and concordance.
# (Superseded the standalone deg_pseudobulk.R -- retired 2026-07-13, output bit-identical.) ---
rule de_engine:
    message: "de_engine: count_engine pseudobulk NB DE (apeglm two-fit) -> sufficient stats"
    input:
        script = "scripts/engines/count_engine.R",
        se     = f"{D_AGG}/pseudobulk_se.rds",
        comp   = "results/data_model/comparisons.tsv",
        params = "config/engine_params.yaml",
    output:
        stats   = "results/aggregate/engine/de_engine.tsv",
        skipped = "results/aggregate/engine/de_engine_skipped.tsv",
    threads: 1
    log:
        f"{D_LOGS}/de_engine.log",
    shell:
        "{RSCRIPT} {input.script} {input.se} mutter02_day2 {input.comp} {input.params} "
        "{output.stats} {output.skipped} > {log} 2>&1"
# --- Track 2: detection-vs-level decomposition + sparsity-aware call_class (Fix 2) ---
rule detection_test:
    message: "detection_test: detection-vs-level decomposition + sparsity-aware call_class"
    input:
        script = f"{R_SCRIPTS}/detection_test.R",
        se     = f"{D_AGG}/pseudobulk_se.rds",
    output:
        det = "results/aggregate/detection_test_m02day2.tsv",
    threads: 1
    log:
        f"{D_LOGS}/detection_test.log",
    shell:
        "{RSCRIPT} {input.script} {input.se} {output.det} > {log} 2>&1"
# --- Track 2 pathway: GSEA on pseudobulk stat-ranked genes (primary + Hallmark) ---
rule gsea:
    message: "gsea: GSEA on pseudobulk stat-ranked genes (primary + Hallmark)"
    input:
        script = f"{R_SCRIPTS}/gsea.R",
        degs   = "results/aggregate/engine/de_engine.tsv",
        sets   = "results/data_model/pathway_sets.tsv",
    output:
        gsea = "results/aggregate/gsea_pseudobulk_m02day2.tsv",
    threads: 1
    log:
        f"{D_LOGS}/gsea.log",
    shell:
        "{RSCRIPT} {input.script} {input.degs} {input.sets} {output.gsea} > {log} 2>&1"
# --- Track 2 pathway: per-cell UCell + AddModuleScore scoring + M02 limma test ---
# threads: 4 -- UCell runs ncores=4 fork (BiocParallel) parallelism, not BLAS.
rule pathway_scores:
    message: "pathway_scores: per-cell UCell + AddModuleScore scoring + M02 limma test"
    input:
        script = f"{R_SCRIPTS}/pathway_scores.R",
        rds    = MERGED,   # count-only consumer; labels from full_labels.parquet
        labels = LABELS,
        sets   = "results/data_model/pathway_sets.tsv",
    output:
        summary = "results/aggregate/pathway_scores_summary.tsv",
        test    = "results/aggregate/pathway_test_m02day2.tsv",
        conc    = "results/aggregate/pathway_ucell_ams_concordance.tsv",
    threads: 4
    log:
        f"{D_LOGS}/pathway_scores.log",
    shell:
        "{RSCRIPT} {input.script} {input.rds} {input.labels} {input.sets} "
        "{output.summary} {output.test} {output.conc} > {log} 2>&1"
# Plot half: reads the cached score/test TSVs so a plot bug can't roll back the ~5h compute.
rule pathway_plots:
    message: "pathway_plots: pathway heatmap/timecourse/concordance plots from cached scores"
    input:
        script  = f"{R_SCRIPTS}/pathway_plots.R",
        summary = "results/aggregate/pathway_scores_summary.tsv",
        test    = "results/aggregate/pathway_test_m02day2.tsv",
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
# --- T9: data-driven spatial niches (k=20 NN composition -> K=6 k-means) ---
rule niches:
    message: "niches: data-driven spatial niches (k=20 NN composition -> K=6 k-means)"
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
    message: "niche_engine: lm_engine proportion test on niche frequency (M02 day2)"
    input:
        script = "scripts/engines/lm_engine.R",
        cells  = "results/aggregate/engine_inputs/niche_cells.tsv",
        comp   = "results/data_model/comparisons.tsv",
        params = "config/engine_params.yaml",
    output:
        stats = "results/aggregate/engine/niche_engine.tsv",
    threads: 1
    log:
        f"{D_LOGS}/niche_engine.log",
    shell:
        "{RSCRIPT} {input.script} {input.cells} proportion mutter02_day2 niche "
        "{input.comp} {input.params} {output.stats} > {log} 2>&1"
# --- T10: tumor-immune spatial mixing (immune-neighbour fraction + Keren score) ---
rule spatial_mixing:
    message: "spatial_mixing: tumor-immune spatial mixing (immune-neighbour fraction + Keren score)"
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
    message: "mixing_engine: lm_engine identity test on spatial mixing metrics (M02 day2)"
    input:
        script = "scripts/engines/lm_engine.R",
        matrix = "results/aggregate/engine_inputs/mixing_metrics.tsv",
        comp   = "results/data_model/comparisons.tsv",
        params = "config/engine_params.yaml",
    output:
        stats = "results/aggregate/engine/mixing_engine.tsv",
    threads: 1
    log:
        f"{D_LOGS}/mixing_engine.log",
    shell:
        "{RSCRIPT} {input.script} {input.matrix} matrix mutter02_day2 mixing "
        "{input.comp} {input.params} {output.stats} > {log} 2>&1"
# --- T11: myeloid M1/M2 polarization (UCell M1/M2 panels -> per-sample ratio) ---
# threads: 8 -- AddModuleScore_UCell forks UCELL_CORES=8 workers (BiocParallel, not
# BLAS), so reserve 8 to avoid oversubscription.
rule myeloid_polarization:
    message: "myeloid_polarization: myeloid M1/M2 polarization (UCell M1/M2 panels -> per-sample ratio)"
    input:
        script = f"{R_SCRIPTS}/myeloid_polarization.R",
        rds    = MERGED,
        labels = LABELS,
    output:
        scores = "results/aggregate/myeloid_m1m2_scores.tsv",
        lm_input = "results/aggregate/engine_inputs/myeloid_metrics.tsv",
        plot   = "results/aggregate/plots/myeloid_m1m2_m02.png",
    threads: 8
    log:
        f"{D_LOGS}/myeloid_polarization.log",
    shell:
        "{RSCRIPT} {input.script} {input.rds} {input.labels} {output.scores} "
        "{output.lm_input} {output.plot} > {log} 2>&1"
# --- engine: myeloid M1/M2 arm test (matrix/identity path, non-robust) ---
rule myeloid_engine:
    message: "myeloid_engine: lm_engine identity test on macrophage M1/M2 metrics (M02 day2)"
    input:
        script = "scripts/engines/lm_engine.R",
        matrix = "results/aggregate/engine_inputs/myeloid_metrics.tsv",
        comp   = "results/data_model/comparisons.tsv",
        params = "config/engine_params.yaml",
    output:
        stats = "results/aggregate/engine/myeloid_engine.tsv",
    threads: 1
    log:
        f"{D_LOGS}/myeloid_engine.log",
    shell:
        "{RSCRIPT} {input.script} {input.matrix} matrix mutter02_day2 myeloid "
        "{input.comp} {input.params} {output.stats} > {log} 2>&1"
# --- T12: M01<->M02 day-2 effect-size concordance on the unified labels ---
rule concordance_m01_m02:
    message: "concordance_m01_m02: M01<->M02 day-2 effect-size concordance on the unified labels"
    input:
        script = f"{R_SCRIPTS}/concordance_m01_m02.R",
        rds    = MERGED,
        labels = LABELS,
        degs   = "results/aggregate/engine/de_engine.tsv",
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
    message: "substate_engine: lm_engine proportion test on fibroblast sub-state composition (M02 day2)"
    input:
        script = "scripts/engines/lm_engine.R",
        cells  = "results/aggregate/engine_inputs/substate_cells.tsv",
        comp   = "results/data_model/comparisons.tsv",
        params = "config/engine_params.yaml",
    output:
        stats = "results/aggregate/engine/substate_engine.tsv",
    threads: 1
    log:
        f"{D_LOGS}/substate_engine.log",
    shell:
        "{RSCRIPT} {input.script} {input.cells} proportion mutter02_day2 substate "
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
        comp    = "results/aggregate/engine/composition_engine.tsv",
        degs    = "results/aggregate/engine/de_engine.tsv",
        gsea    = "results/aggregate/gsea_pseudobulk_m02day2.tsv",
        pathway = "results/aggregate/pathway_test_m02day2.tsv",
        niche   = "results/aggregate/engine/niche_engine.tsv",
        mixing  = "results/aggregate/engine/mixing_engine.tsv",
        myeloid = "results/aggregate/engine/myeloid_engine.tsv",
        mde     = "results/aggregate/power_mde.tsv",
        sets    = "results/data_model/pathway_sets.tsv",
        cov     = "results/data_model/gene_set_panel_coverage.tsv",
        det     = "results/aggregate/detection_test_m02day2.tsv",
    output:
        master = "results/aggregate/results_master.tsv",
        detect = "results/aggregate/detectability_summary.tsv",
    threads: 1
    log:
        f"{D_LOGS}/assemble_results.log",
    shell:
        "{RSCRIPT} {input.script} {input.comp} {input.degs} {input.gsea} "
        "{input.pathway} {input.niche} {input.mixing} {input.myeloid} "
        "{input.mde} {input.sets} {input.cov} {input.det} {output.master} > {log} 2>&1"
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
        tsv = "results/aggregate/qc_arm_balance.tsv",
        png = "results/aggregate/plots/qc_arm_balance.png",
    threads: 1
    log:
        f"{D_LOGS}/qc_arm_balance.log",
    shell:
        "{RSCRIPT} {input.script} {input.tech} {input.contam} {input.samples} "
        "{output.tsv} {output.png} > {log} 2>&1"
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
        png = "results/aggregate/plots/qc_reproducibility.png",
    threads: 1
    log:
        f"{D_LOGS}/qc_reproducibility.log",
    shell:
        "{RSCRIPT} {input.script} {input.se} {input.tech} {output.tsv} {output.png} > {log} 2>&1"
