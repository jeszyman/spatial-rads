# spatial-rads -- data model workflow (Stage A)
# Builds the relational metadata.xlsx into data_model.rda + the master sample sheet.
# Run: conda run -n basecamp snakemake -s workflows/data_model.smk --cores 1
# R steps run in the spatial-rads env via `conda run` (snakemake driver lives in basecamp).
configfile: "config/config.yaml"

RSCRIPT = "conda run -n spatial-rads Rscript"
rule all:
    input:
        config["samplesheet"],
        "data/data_model.rda",
        "results/data_model/panel_provenance_summary.tsv",
        "results/data_model/common_genes.tsv",
        "results/data_model/pathway_sets.tsv",
        "results/data_model/gene_set_panel_coverage.tsv",
        "results/data_model/comparisons.tsv",
        "results/data_model/marker_panel_coverage.tsv",
        "results/data_model/panel_provenance.tsv",
        "results/data_model/plots/design_grid.png",
        "results/data_model/plots/panel_provenance.png",
        "results/data_model/plots/geneset_coverage.png",
        "results/data_model/plots/marker_coverage.png",
rule make_data_model:
    input:
        xlsx   = config["metadata"]["xlsx"],
        schema = config["metadata"]["schema"],
    output:
        rda         = "data/data_model.rda",
        samplesheet = config["samplesheet"],
    log:
        "logs/make_data_model.log",
    shell:
        "{RSCRIPT} scripts/make_data_model.R {input.xlsx} {input.schema} "
        "{output.rda} {output.samplesheet} > {log} 2>&1"
# --- common gene panel (Mutter_01 panel INT Mutter_02 panel; read from per-slide RDS) ---
rule common_gene_panel:
    input:
        script = "scripts/common_gene_panel.R",
        ss     = config["samplesheet"],
    output:
        "results/data_model/common_genes.tsv",
    threads: 1
    log: "logs/common_gene_panel.log",
    shell:
        "{RSCRIPT} {input.script} {input.ss} {output} > {log} 2>&1"
# --- tiered pathway gene-set artifact + panel coverage (single source of truth) ---
rule build_gene_sets:
    input:
        script = "scripts/build_gene_sets.R",
        yaml   = config["pathway"]["gene_lists"],
        prov   = "config/pathway_sets_provenance.tsv",
        panel  = "results/data_model/common_genes.tsv",
    output:
        sets     = "results/data_model/pathway_sets.tsv",
        coverage = "results/data_model/gene_set_panel_coverage.tsv",
    params:
        minpg = config["pathway"]["min_panel_genes"],
    threads: 1
    log: "logs/build_gene_sets.log",
    shell:
        "{RSCRIPT} {input.script} {input.yaml} {input.prov} {input.panel} {params.minpg} "
        "{output.sets} {output.coverage} > {log} 2>&1"
# --- comparison registry: curated design resolved against the samplesheet ---
rule build_comparisons:
    input:
        script = "scripts/build_comparisons.R",
        yaml   = "config/comparisons.yaml",
        ss     = config["samplesheet"],
    output:
        "results/data_model/comparisons.tsv",
    threads: 1
    log: "logs/build_comparisons.log",
    shell:
        "{RSCRIPT} {input.script} {input.yaml} {input.ss} {output} > {log} 2>&1"
# --- cell-marker provenance + panel coverage (mirrors build_gene_sets) ---
rule build_marker_sets:
    input:
        script   = "scripts/build_marker_sets.R",
        lineage  = "config/lineage_markers.yaml",
        substate = "config/substate_markers.yaml",
        prov     = "config/marker_sets_provenance.tsv",
        panel    = "results/data_model/common_genes.tsv",
    output:
        coverage = "results/data_model/marker_panel_coverage.tsv",
    params:
        minpg = config["pathway"]["min_panel_genes"],
    threads: 1
    log: "logs/build_marker_sets.log",
    shell:
        "{RSCRIPT} {input.script} {input.lineage} {input.substate} {input.prov} "
        "{input.panel} {params.minpg} {output.coverage} > {log} 2>&1"
rule fig_marker_coverage:
    input:
        script = "scripts/fig_marker_coverage.R",
        cov    = "results/data_model/marker_panel_coverage.tsv",
    output:
        "results/data_model/plots/marker_coverage.png",
    threads: 1
    log: "logs/fig_marker_coverage.log",
    shell:
        "{RSCRIPT} {input.script} {input.cov} {output} > {log} 2>&1"
# --- panel provenance: UCC-standard vs per-experiment custom (from the delivered objects) ---
rule panel_provenance:
    input:
        script = "scripts/panel_provenance.R",
        vendor = "data/sources/2026-06-22-bruker-mouse-ucc-gene-list.xlsx",
    output:
        membership = "results/data_model/panel_provenance.tsv",
        summary    = "results/data_model/panel_provenance_summary.tsv",
    params:
        m01 = config["datadir"] + "/inputs/mutter01",
        m02 = config["datadir"] + "/inputs/mutter02",
    threads: 1
    log: "logs/panel_provenance.log",
    shell:
        "{RSCRIPT} {input.script} {input.vendor} {params.m01} {params.m02} "
        "{output.membership} {output.summary} > {log} 2>&1"
rule fig_design_grid:
    input:
        script = "scripts/fig_design_grid.R",
        ss = config["samplesheet"],
    output:
        "results/data_model/plots/design_grid.png",
    threads: 1
    log: "logs/fig_design_grid.log",
    shell:
        "{RSCRIPT} {input.script} {input.ss} {output} > {log} 2>&1"
rule fig_panel_provenance:
    input:
        script = "scripts/fig_panel_provenance.R",
        pv = "results/data_model/panel_provenance.tsv",
    output:
        "results/data_model/plots/panel_provenance.png",
    threads: 1
    log: "logs/fig_panel_provenance.log",
    shell:
        "{RSCRIPT} {input.script} {input.pv} {output} > {log} 2>&1"
rule fig_geneset_coverage:
    input:
        script = "scripts/fig_geneset_coverage.R",
        cov = "results/data_model/gene_set_panel_coverage.tsv",
    output:
        "results/data_model/plots/geneset_coverage.png",
    threads: 1
    log: "logs/fig_geneset_coverage.log",
    shell:
        "{RSCRIPT} {input.script} {input.cov} {output} > {log} 2>&1"
