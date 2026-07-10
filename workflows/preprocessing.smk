# spatial-rads -- invariant per-sample preprocessing workflow (pre-aggregate).
# Deterministic function of (raw input, pinned config): dataset adapter (raw per-slide
# RDS -> common 950-gene Seurat) -> 4-criteria QC filter -> LogNormalize, plus report-only
# QC (probe QC, control characterization, SpatialQM contamination MECR, QC summaries/plots).
# Per-sample terminus = norm.rds; canonical cell typing happens at merged scale in the
# aggregate workflow. Run: conda run -n basecamp snakemake -s workflows/preprocessing.smk --cores 8
# R steps run in the spatial-rads env via `conda run` (driver lives in basecamp).

import os
import pandas as pd

configfile: "config/config.yaml"

# Strict mode (Snakemake default, preserved) + BLAS pinned to 1 thread per R process:
# snakemake schedules job-level parallelism, not R BLAS multi-threading. Without the pin,
# --cores N x ~48 BLAS threads each swamps the box.
shell.prefix("set -euo pipefail; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1; ")

# --- interpreter (driver runs in basecamp; R steps in spatial-rads) ---
RSCRIPT = "conda run -n spatial-rads Rscript"

# --- paths: D_ data dirs, R_ repo locations ---
D_DATA    = config["datadir"]
D_PROC    = os.path.join(D_DATA, "processing")    # heavy per-sample RDS
D_RES     = "results/processing"                  # small TSV/PNG reports
D_LOGS    = "logs"
R_SCRIPTS = "scripts"

MASTER  = config["samplesheet"]                   # master sheet from the data_model workflow
SCOPED  = f"{D_RES}/samplesheet.tsv"
GENES   = "results/data_model/common_genes.tsv"   # common panel (M01 n M02; built in data_model.smk)
YIREF   = f"{D_RES}/yi_reference.tsv"             # M01 vendor-label side-output (retired with per-sample typing)
QC      = config["qc"]
NRM     = config["normalize"]
M01_RDS_DIR = f"{D_DATA}/inputs/mutter01"
M02_RDS_DIR = f"{D_DATA}/inputs/mutter02"
M01_META_PQ = f"{D_DATA}/inputs/mutter01/Analysis_Mutter_01_CosMmR_Mutter_updated_metadata.parquet"

# Parse-time sample lists from the master sheet.
_s = pd.read_csv(MASTER, sep="\t")
M01 = _s.loc[_s["name"] == "Mutter_01", "sample_id"].tolist()
M02 = _s.loc[_s["name"] == "Mutter_02", "sample_id"].tolist()
ALL = _s["sample_id"].tolist()
rule all:
    input:
        expand(f"{D_PROC}/norm/{{s}}.norm.rds", s=ALL),
        f"{D_RES}/qc_summary.tsv",
        f"{D_RES}/probe_qc_report.tsv",
        f"{D_RES}/m01_rds_validation.tsv",
        f"{D_RES}/control_characterization.tsv",
        f"{D_RES}/fov_falsecode_qc.tsv",
        f"{D_RES}/contamination_qc.tsv",
        f"{D_RES}/contamination_fov_qc.tsv",
        f"{D_RES}/plots/qc_removal_attribution.png",
        f"{D_RES}/plots/qc_metric_distributions.png",
# Workflow-linked sample sheet (scoped view of the data model).
rule preproc_samplesheet:
    message: "[samplesheet] scoped view of the data model"
    input:
        script = f"{R_SCRIPTS}/processing_samplesheet.R",
        rda    = "data/data_model.rda",
    log:
        f"{D_LOGS}/preproc_samplesheet.log",
    threads: 1
    output:
        ss = SCOPED,
    shell:
        "{RSCRIPT} {input.script} {input.rda} {output.ss} > {log} 2>&1"
# Probe-vs-negative-control QC DIAGNOSTIC (report-only; does NOT drop genes).
rule preproc_probe_qc:
    message: "[probe_qc] background-level probe flags"
    input:
        script = f"{R_SCRIPTS}/probe_qc.R",
        ss     = SCOPED,
        panel  = GENES,
    log:
        f"{D_LOGS}/preproc_probe_qc.log",
    threads: 1
    output:
        report = f"{D_RES}/probe_qc_report.tsv",
    shell:
        "{RSCRIPT} {input.script} {input.ss} {input.panel} {output.report} > {log} 2>&1"
# Control QC (report-only): negprobe + falsecode characterization & per-cell sidecar.
# Additive; depends on raw RDS + M01 parquet, never scored objects.
rule preproc_control_qc:
    message: "[control_qc] negprobe/falsecode characterization + sidecar"
    input:
        validate     = f"{R_SCRIPTS}/validate_m01_rds.R",
        characterize = f"{R_SCRIPTS}/characterize_controls.R",
        sidecar      = f"{R_SCRIPTS}/control_sidecar.R",
        m01_meta     = M01_META_PQ,
    log:
        f"{D_LOGS}/preproc_control_qc.log",
    threads: 1
    output:
        valid   = f"{D_RES}/m01_rds_validation.tsv",
        charac  = f"{D_RES}/control_characterization.tsv",
        density = f"{D_RES}/plots/fov_falsecode_density.png",
        cells   = f"{D_PROC}/cell_controls.parquet",
        fovqc   = f"{D_RES}/fov_falsecode_qc.tsv",
    shell:
        "{RSCRIPT} {input.validate} " + M01_RDS_DIR + " {input.m01_meta} {output.valid} > {log} 2>&1 && "
        "{RSCRIPT} {input.characterize} " + M01_RDS_DIR + " " + M02_RDS_DIR + " {output.charac} {output.density} >> {log} 2>&1 && "
        "{RSCRIPT} {input.sidecar} " + M01_RDS_DIR + " {input.m01_meta} " + M02_RDS_DIR + " {output.cells} {output.fovqc} >> {log} 2>&1"
# Contamination QC (report-only): SpatialQM MECR marker-bleed, sample + FOV grain.
rule preproc_contamination_qc:
    message: "[contamination_qc] SpatialQM MECR marker-bleed"
    input:
        script  = f"{R_SCRIPTS}/contamination_qc.R",
        qc      = expand(f"{D_PROC}/qc/{{s}}.qc.rds", s=ALL),
        sqm     = f"{R_SCRIPTS}/aggregate/spatialqm_metrics.R",
        markers = "config/lineage_markers.yaml",
    log:
        f"{D_LOGS}/preproc_contamination_qc.log",
    threads: 1
    output:
        sample = f"{D_RES}/contamination_qc.tsv",
        fov    = f"{D_RES}/contamination_fov_qc.tsv",
    shell:
        "{RSCRIPT} {input.script} {D_PROC}/qc {input.sqm} "
        "{input.markers} {output.sample} {output.fov} > {log} 2>&1"
rule preproc_adapt_mutter01:
    message: "[adapt_mutter01] raw per-slide RDS -> common Seurat"
    input:
        script = f"{R_SCRIPTS}/adapt_mutter01.R",
        ss     = SCOPED,
        genes  = GENES,
    log:
        f"{D_LOGS}/preproc_adapt_mutter01.log",
    threads: 1
    output:
        rds = expand(f"{D_PROC}/raw/{{s}}.raw.rds", s=M01),
        yi  = YIREF,
    shell:
        "{RSCRIPT} {input.script} {input.ss} {input.genes} {D_DATA} {output.yi} > {log} 2>&1"
rule preproc_adapt_mutter02:
    message: "[adapt_mutter02] raw per-slide RDS -> common Seurat"
    input:
        script = f"{R_SCRIPTS}/adapt_mutter02.R",
        ss     = SCOPED,
        genes  = GENES,
    log:
        f"{D_LOGS}/preproc_adapt_mutter02.log",
    threads: 1
    output:
        rds = expand(f"{D_PROC}/raw/{{s}}.raw.rds", s=M02),
    shell:
        "{RSCRIPT} {input.script} {input.ss} {input.genes} {D_DATA} > {log} 2>&1"
rule preproc_qc_filter:
    message: "[qc_filter] {wildcards.s}"
    input:
        script = f"{R_SCRIPTS}/qc_filter.R",
        rds    = f"{D_PROC}/raw/{{s}}.raw.rds",
    log:
        f"{D_LOGS}/preproc_qc_filter/{{s}}.log",
    params:
        min_counts   = QC["min_counts"],
        min_features = QC["min_features"],
        max_propneg  = QC["max_prop_negative"],
        area_nmads   = QC["area_nmads"],
    threads: 1
    output:
        rds     = f"{D_PROC}/qc/{{s}}.qc.rds",
        summary = f"{D_RES}/qc/{{s}}.qcsummary.tsv",
    shell:
        "{RSCRIPT} {input.script} {input.rds} {output.rds} {output.summary} "
        "{params.min_counts} {params.min_features} {params.max_propneg} {params.area_nmads} > {log} 2>&1"
rule preproc_normalize:
    message: "[normalize] {wildcards.s}"
    input:
        script = f"{R_SCRIPTS}/normalize.R",
        rds    = f"{D_PROC}/qc/{{s}}.qc.rds",
    log:
        f"{D_LOGS}/preproc_normalize/{{s}}.log",
    params:
        scale_factor = NRM["scale_factor"],
        nfeatures    = NRM["n_variable_features"],
    threads: 1
    output:
        rds = f"{D_PROC}/norm/{{s}}.norm.rds",
    shell:
        "{RSCRIPT} {input.script} {input.rds} {output.rds} {params.scale_factor} {params.nfeatures} > {log} 2>&1"
rule preproc_qc_report:
    message: "[qc_report] pooled per-sample QC summaries"
    input:
        script = f"{R_SCRIPTS}/qc_report.R",
        sums   = expand(f"{D_RES}/qc/{{s}}.qcsummary.tsv", s=ALL),
    log:
        f"{D_LOGS}/preproc_qc_report.log",
    threads: 1
    output:
        report = f"{D_RES}/qc_summary.tsv",
    shell:
        "{RSCRIPT} {input.script} {D_RES}/qc {output.report} > {log} 2>&1"
rule preproc_qc_plots:
    message: "[qc_plots] QC removal attribution + metric distributions"
    input:
        script = f"{R_SCRIPTS}/qc_plots.R",
        raw    = expand(f"{D_PROC}/raw/{{s}}.raw.rds", s=ALL),
        cfg    = "config/config.yaml",
    log:
        f"{D_LOGS}/preproc_qc_plots.log",
    threads: 1
    output:
        attribution   = f"{D_RES}/plots/qc_removal_attribution.png",
        distributions = f"{D_RES}/plots/qc_metric_distributions.png",
        attr_tsv      = f"{D_RES}/plots/qc_removal_attribution.tsv",
    shell:
        "{RSCRIPT} {input.script} {input.cfg} {D_RES}/plots > {log} 2>&1"
