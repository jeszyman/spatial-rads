# spatial-rads -- data model workflow (Stage A)
# Builds the relational metadata.xlsx into data_model.rda + the master sample sheet.
# Run: conda run -n basecamp snakemake -s workflows/data_model.smk --cores 1
# R steps run in the spatial-rads env via `conda run` (snakemake driver lives in basecamp).
configfile: "config/config.yaml"

RSCRIPT = "conda run -n spatial-rads Rscript"

rule all:
    input:
        config["samplesheet"],

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
