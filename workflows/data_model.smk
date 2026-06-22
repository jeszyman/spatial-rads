# spatial-rads -- data model workflow (Stage A)
# Builds the relational metadata.xlsx into data_model.rda + the master sample sheet.
# Run: conda run -n basecamp snakemake -s workflows/data_model.smk --cores 1
# R steps run in the spatial-rads env via `conda run` (snakemake driver lives in basecamp).
configfile: "config/config.yaml"

RSCRIPT = "conda run -n spatial-rads Rscript"

rule all:
    input:
        config["samplesheet"],
        "results/data_model/common_genes.tsv",
        "results/data_model/pathway_sets.tsv",
        "results/data_model/gene_set_panel_coverage.tsv",

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
