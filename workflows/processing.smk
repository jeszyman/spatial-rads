# spatial-rads -- initial processing & QC workflow (Stage B + C)
# Adapter (raw -> common-format per-sample Seurat) -> QC filter -> LogNormalize.
# Both datasets traverse identical Stage C code; Yi's Mutter_01 layers held as a check.
# Run: conda run -n basecamp snakemake -s workflows/processing.smk --cores 8
# R steps run in the spatial-rads env via `conda run` (driver lives in basecamp).
import pandas as pd
configfile: "config/config.yaml"

RSCRIPT = "conda run -n spatial-rads Rscript"
DATADIR = config["datadir"]
MASTER  = config["samplesheet"]            # master sheet from the data_model workflow
SCOPED  = "results/processing/samplesheet.tsv"
GENES   = "results/processing/common_genes.tsv"  # common panel (Mutter_01 ∩ Mutter_02); NOT gene-filtered
YIREF   = "results/processing/yi_reference.tsv"
QC      = config["qc"]
NRM     = config["normalize"]

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

# --- workflow-linked sample sheet (scoped view of the data model) ---
rule processing_samplesheet:
    input:
        rda = "data/data_model.rda",
    output:
        SCOPED,
    log:
        "logs/processing_samplesheet.log",
    shell:
        "{RSCRIPT} scripts/processing_samplesheet.R {input.rda} {output} > {log} 2>&1"

# --- common gene panel (Mutter_01 panel intersect Mutter_02 panel) ---
rule common_gene_panel:
    input:
        ss = SCOPED,
    output:
        GENES,
    log:
        "logs/common_gene_panel.log",
    shell:
        "{RSCRIPT} scripts/common_gene_panel.R {input.ss} {output} > {log} 2>&1"

# --- probe-vs-negative-control QC DIAGNOSTIC (report-only; does NOT drop genes -- see scripts/probe_qc.R) ---
rule probe_qc:
    input:
        ss    = SCOPED,
        panel = GENES,
    output:
        "results/processing/probe_qc_report.tsv",
    log:
        "logs/probe_qc.log",
    shell:
        "{RSCRIPT} scripts/probe_qc.R {input.ss} {input.panel} {output} > {log} 2>&1"

# --- Stage B adapters: raw -> per-sample common-format Seurat ---
rule adapt_mutter01:
    input:
        ss    = SCOPED,
        genes = GENES,
    output:
        rds = expand(f"{DATADIR}/processing/raw/{{s}}.raw.rds", s=M01),
        yi  = YIREF,
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
    log:
        "logs/celltype.{s}.log",
    shell:
        "{RSCRIPT} scripts/celltype.R {input.norm} {input.ref} {input.yi} {output.rds} {output.summary} > {log} 2>&1"

rule celltype_report:
    input:
        expand("results/processing/celltype/{s}.celltype.tsv", s=ALL),
    output:
        "results/processing/celltype_summary.tsv",
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
    log:
        "logs/qc_plots.log",
    shell:
        "{RSCRIPT} scripts/qc_plots.R {input.cfg} results/processing/plots > {log} 2>&1"
