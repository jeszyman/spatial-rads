# spatial-rads -- initial processing & QC workflow (per-sample, pre-aggregate).
# adapter (raw -> common-format Seurat) -> QC filter -> LogNormalize -> cell typing
# (Yi labels for M01; Seurat anchor TransferData for M02 against pooled M01 cell reference)
# -> pathway scoring (UCell + AddModuleScore). probe-QC diagnostic flags background-level
# probes (report-only). Both datasets traverse identical Stage C onward; Yi's Mutter_01
# layers held as a check. Run: conda run -n basecamp snakemake -s workflows/processing.smk --cores 4
# R steps run in the spatial-rads env via `conda run` (driver lives in basecamp).
import pandas as pd
configfile: "config/config.yaml"

# Pin BLAS to 1 thread per R process: snakemake schedules job-level parallelism, not
# R BLAS-level multi-threading. Without this, --cores N x ~48 BLAS threads each swamps the box.
shell.prefix("export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1; ")

RSCRIPT = "conda run -n spatial-rads Rscript"
DATADIR = config["datadir"]
MASTER  = config["samplesheet"]            # master sheet from the data_model workflow
SCOPED  = "results/processing/samplesheet.tsv"
GENES   = "results/data_model/common_genes.tsv"  # common panel (Mutter_01 ∩ Mutter_02; built in data_model.smk); NOT gene-filtered
YIREF   = "results/processing/yi_reference.tsv"
QC      = config["qc"]
NRM     = config["normalize"]
M01_RDS_DIR = f"{DATADIR}/inputs/mutter01"
M02_RDS_DIR = f"{DATADIR}/inputs/mutter02"
M01_META_PQ = f"{DATADIR}/inputs/mutter01/Analysis_Mutter_01_CosMmR_Mutter_updated_metadata.parquet"

# Parse-time sample lists come from the existing master sheet.
_s = pd.read_csv(MASTER, sep="\t")
M01 = _s.loc[_s["name"] == "Mutter_01", "sample_id"].tolist()
M02 = _s.loc[_s["name"] == "Mutter_02", "sample_id"].tolist()
ALL = _s["sample_id"].tolist()

rule all:
    input:
        expand(f"{DATADIR}/processing/scored/{{s}}.scored.rds", s=ALL),
        "results/processing/qc_summary.tsv",
        "results/processing/yi_concordance.tsv",
        "results/processing/celltype_summary.tsv",
        "results/processing/probe_qc_report.tsv",
        "results/processing/plots/qc_removal_attribution.png",
        "results/processing/plots/qc_metric_distributions.png",
        "results/processing/m01_rds_validation.tsv",
        "results/processing/control_characterization.tsv",
        "results/processing/fov_falsecode_qc.tsv",
        "results/processing/contamination_qc.tsv",
        "results/processing/contamination_fov_qc.tsv",

# --- workflow-linked sample sheet (scoped view of the data model) ---
rule processing_samplesheet:
    input:
        rda = "data/data_model.rda",
    output:
        SCOPED,
    threads: 1
    log:
        "logs/processing_samplesheet.log",
    shell:
        "{RSCRIPT} scripts/processing_samplesheet.R {input.rda} {output} > {log} 2>&1"

# --- probe-vs-negative-control QC DIAGNOSTIC (report-only; does NOT drop genes -- see scripts/probe_qc.R) ---
rule probe_qc:
    input:
        ss    = SCOPED,
        panel = GENES,
    output:
        "results/processing/probe_qc_report.tsv",
    threads: 1
    log:
        "logs/probe_qc.log",
    shell:
        "{RSCRIPT} scripts/probe_qc.R {input.ss} {input.panel} {output} > {log} 2>&1"

# --- Control QC (report-only): negprobe + falsecode characterization & per-cell sidecar ---
# Additive diagnostic. Depends on raw RDS + M01 parquet, NOT scored.rds -> cannot re-fire the
# locked aggregate. propNegative/scored objects are untouched.
rule control_qc:
    input:
        validate     = "scripts/validate_m01_rds.R",
        characterize = "scripts/characterize_controls.R",
        sidecar      = "scripts/control_sidecar.R",
        m01_meta     = M01_META_PQ,
    output:
        valid   = "results/processing/m01_rds_validation.tsv",
        charac  = "results/processing/control_characterization.tsv",
        falsecode_density = "results/processing/plots/fov_falsecode_density.png",
        cells   = f"{DATADIR}/processing/cell_controls.parquet",
        fovqc   = "results/processing/fov_falsecode_qc.tsv",
    log:
        "logs/control_qc.log",
    threads: 1
    shell:
        "{RSCRIPT} {input.validate} " + M01_RDS_DIR + " {input.m01_meta} {output.valid} > {log} 2>&1 && "
        "{RSCRIPT} {input.characterize} " + M01_RDS_DIR + " " + M02_RDS_DIR + " {output.charac} {output.falsecode_density} >> {log} 2>&1 && "
        "{RSCRIPT} {input.sidecar} " + M01_RDS_DIR + " {input.m01_meta} " + M02_RDS_DIR + " {output.cells} {output.fovqc} >> {log} 2>&1"

# --- Contamination QC (report-only): SpatialQM MECR marker-bleed metric, sample + FOV grain ---
# Additive: reads qc.rds, writes NEW TSVs -> cannot re-fire the raw->qc->norm->typed->scored chain.
# Per-FOV table carries orthogonal triage cols (counts/propNeg/area) to separate segmentation
# contamination from necrosis/low-quality tissue vs genuine mixed biology. NEVER excludes (report-only).
rule contamination_qc:
    input:
        qc      = expand(f"{DATADIR}/processing/qc/{{s}}.qc.rds", s=ALL),
        sqm     = "scripts/aggregate/spatialqm_metrics.R",
        markers = "config/lineage_markers.yaml",
        script  = "scripts/contamination_qc.R",
    output:
        sample = "results/processing/contamination_qc.tsv",
        fov    = "results/processing/contamination_fov_qc.tsv",
    threads: 1
    log:
        "logs/contamination_qc.log",
    shell:
        "{RSCRIPT} scripts/contamination_qc.R {DATADIR}/processing/qc {input.sqm} "
        "{input.markers} {output.sample} {output.fov} > {log} 2>&1"

# --- Stage B adapters: raw -> per-sample common-format Seurat ---
rule adapt_mutter01:
    input:
        ss    = SCOPED,
        genes = GENES,
    output:
        rds = expand(f"{DATADIR}/processing/raw/{{s}}.raw.rds", s=M01),
        yi  = YIREF,
    threads: 1
    log:
        "logs/adapt_mutter01.log",
    shell:
        "{RSCRIPT} scripts/adapt_mutter01.R {input.ss} {input.genes} {DATADIR} {output.yi} > {log} 2>&1"

rule adapt_mutter02:
    input:
        ss    = SCOPED,
        genes = GENES,
    output:
        rds = expand(f"{DATADIR}/processing/raw/{{s}}.raw.rds", s=M02),
    threads: 1
    log:
        "logs/adapt_mutter02.log",
    shell:
        "{RSCRIPT} scripts/adapt_mutter02.R {input.ss} {input.genes} {DATADIR} > {log} 2>&1"

# --- Stage C (shared): QC filter -> normalize ---
rule qc_filter:
    input:
        f"{DATADIR}/processing/raw/{{s}}.raw.rds",
    output:
        rds     = f"{DATADIR}/processing/qc/{{s}}.qc.rds",
        summary = "results/processing/qc/{s}.qcsummary.tsv",
    params:
        min_counts   = QC["min_counts"],
        min_features = QC["min_features"],
        max_propneg  = QC["max_prop_negative"],
        area_nmads   = QC["area_nmads"],
    threads: 1
    log:
        "logs/qc_filter.{s}.log",
    shell:
        "{RSCRIPT} scripts/qc_filter.R {input} {output.rds} {output.summary} "
        "{params.min_counts} {params.min_features} {params.max_propneg} {params.area_nmads} > {log} 2>&1"

rule normalize:
    input:
        f"{DATADIR}/processing/qc/{{s}}.qc.rds",
    output:
        f"{DATADIR}/processing/norm/{{s}}.norm.rds",
    params:
        scale_factor = NRM["scale_factor"],
        nfeatures    = NRM["n_variable_features"],
    threads: 1
    log:
        "logs/normalize.{s}.log",
    shell:
        "{RSCRIPT} scripts/normalize.R {input} {output} {params.scale_factor} {params.nfeatures} > {log} 2>&1"

# --- cell typing (label transfer against M01-derived 42-type reference) ---
rule build_celltype_reference:
    input:
        yi    = YIREF,
        genes = GENES,
        m01   = expand(f"{DATADIR}/processing/norm/{{s}}.norm.rds", s=M01),
    output:
        f"{DATADIR}/processing/celltype_reference.rds",
    threads: 1
    log:
        "logs/build_celltype_reference.log",
    shell:
        "{RSCRIPT} scripts/build_celltype_reference.R {input.yi} {input.genes} {output} {input.m01} > {log} 2>&1"

rule celltype:
    input:
        norm = f"{DATADIR}/processing/norm/{{s}}.norm.rds",
        ref  = f"{DATADIR}/processing/celltype_reference.rds",
        yi   = YIREF,
    output:
        rds     = f"{DATADIR}/processing/typed/{{s}}.typed.rds",
        summary = "results/processing/celltype/{s}.celltype.tsv",
    threads: 1
    log:
        "logs/celltype.{s}.log",
    shell:
        "{RSCRIPT} scripts/celltype.R {input.norm} {input.ref} {input.yi} {output.rds} {output.summary} > {log} 2>&1"

rule celltype_report:
    input:
        expand("results/processing/celltype/{s}.celltype.tsv", s=ALL),
    output:
        "results/processing/celltype_summary.tsv",
    threads: 1
    log:
        "logs/celltype_report.log",
    shell:
        "{RSCRIPT} scripts/celltype_report.R results/processing/celltype {output} > {log} 2>&1"

# --- pathway scoring (UCell primary + AddModuleScore seed=42 secondary) ---
rule pathway_score:
    input:
        typed = f"{DATADIR}/processing/typed/{{s}}.typed.rds",
        genes = "config/pathway_gene_lists.yaml",
    output:
        f"{DATADIR}/processing/scored/{{s}}.scored.rds",
    threads: 1
    log:
        "logs/pathway_score.{s}.log",
    shell:
        "{RSCRIPT} scripts/pathway_score.R {input.typed} {input.genes} {output} > {log} 2>&1"

# --- reports / checks ---
rule qc_report:
    input:
        expand("results/processing/qc/{s}.qcsummary.tsv", s=ALL),
    output:
        "results/processing/qc_summary.tsv",
    threads: 1
    log:
        "logs/qc_report.log",
    shell:
        "{RSCRIPT} scripts/qc_report.R results/processing/qc {output} > {log} 2>&1"

rule yi_concordance:
    input:
        yi = YIREF,
        qc = expand(f"{DATADIR}/processing/qc/{{s}}.qc.rds", s=M01),
    output:
        "results/processing/yi_concordance.tsv",
    threads: 1
    log:
        "logs/yi_concordance.log",
    shell:
        "{RSCRIPT} scripts/yi_concordance.R {input.yi} {DATADIR} {output} > {log} 2>&1"

rule qc_plots:
    input:
        raw = expand(f"{DATADIR}/processing/raw/{{s}}.raw.rds", s=ALL),
        cfg = "config/config.yaml",
    output:
        "results/processing/plots/qc_removal_attribution.png",
        "results/processing/plots/qc_metric_distributions.png",
        "results/processing/plots/qc_removal_attribution.tsv",
    threads: 1
    log:
        "logs/qc_plots.log",
    shell:
        "{RSCRIPT} scripts/qc_plots.R {input.cfg} results/processing/plots > {log} 2>&1"
