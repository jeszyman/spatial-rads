# spatial-rads -- cross-sample cell-typing workflow (flank cohort, 20 samples).
# Merges the 20 flank per-sample norm.rds into one 3.27M-cell object, then assigns
# identity once at merged scale: scVI integration (GPU) -> Leiden cluster-then-annotate
# (tier-1 coarse) -> tier-2 immune (SingleR/ImmGen + marker rescue) + stroma (UCell +
# marker rescue) -> unified per-cell label table. Terminus = full_labels.parquet.
# The differential workflow consumes merged.rds + full/obs.parquet + full_labels.parquet.
# Run: TMPDIR=/mnt/data/projects/spatial-rads/tmp \
#        conda run -n basecamp snakemake -s workflows/aggregate_typing.smk --cores N --resources gpu=1
# R steps run in spatial-rads; GPU Python steps in spatial-rads-scvi (driver in basecamp).
import pandas as pd
configfile: "config/config.yaml"

# --- interpreters: R (spatial-rads), GPU Python (spatial-rads-scvi) ---
RSCRIPT = "conda run -n spatial-rads Rscript"
PYSCVI  = "conda run -n spatial-rads-scvi python"

# --- paths: D_ data dirs, R_ repo locations ---
D_DATA    = config["datadir"]
D_AGG     = f"{D_DATA}/aggregate"                  # heavy intermediates
D_FULL    = f"{D_AGG}/full"                        # tier-1/2 typing intermediate dir
D_NORM    = f"{D_DATA}/processing/norm"            # per-sample inputs (normalized)
D_RES     = "results/aggregate"                    # small TSV/PNG reports
D_LOGS    = "logs/aggregate"
R_SCRIPTS = "scripts/aggregate"

MASTER  = config["samplesheet"]                   # results/data_model/samples.tsv
PANEL   = "results/data_model/common_genes.tsv"   # 950 shared-panel gene list
MARKERS = "config/lineage_markers.yaml"           # coarse + tier-2 lineage marker sets
M01_META = f"{D_DATA}/inputs/mutter01/Analysis_Mutter_01_CosMmR_Mutter_updated_metadata.parquet"
M02_RAW  = f"{D_DATA}/inputs/mutter02"

# --- cohort: flank only (tongue is separate biology, n=1; out of scope) ---
_s    = pd.read_csv(MASTER, sep="\t")
FLANK = _s.loc[_s["model"] == "flank", "sample_id"].tolist()    # 20: M01 x8, M02 x12

# Total post-QC flank cells, for the merge pilot's cell-based memory projection.
_qc         = pd.read_csv("results/processing/qc_summary.tsv", sep="\t")
FLANK_CELLS = int(_qc.loc[_qc["sample_id"].isin(FLANK), "post_filter"].sum())

# Memory pilot subset: 5 smallest M01 + 5 smallest M02 flank samples.
PILOT = ["sam0002", "sam0007", "sam0005", "sam0008", "sam0001",
         "sam0018", "sam0020", "sam0021", "sam0017", "sam0019"]

rule all:
    input:
        "results/aggregate/merge_pilot_memory.tsv",
        f"{D_AGG}/merged.rds",
        "results/aggregate/merge_summary.tsv",
        "results/aggregate/full_coarse_labels.parquet",
        "results/aggregate/full_labels.parquet",
        f"{D_FULL}/immune_subtypes_rescued.parquet",
        f"{D_FULL}/stroma_subtypes_rescued.parquet",
        "results/aggregate/full_labels_summary.tsv",
        "results/aggregate/full_gates.json",
        "results/aggregate/full_coarse_summary.tsv",
        "results/aggregate/full_sweep_summary.tsv",
        "results/aggregate/stroma_sweep_summary.tsv",
# --- Stage 0a: memory pilot (characterize peak RSS, choose merge strategy) ---
rule merge_pilot:
    message: "merge_pilot: characterize peak RSS on a 10-sample subset, choose merge strategy"
    input:
        script = f"{R_SCRIPTS}/merge_pilot.R",
        rds    = expand(f"{D_NORM}/{{s}}.norm.rds", s=PILOT),
    output:
        "results/aggregate/merge_pilot_memory.tsv",
    params:
        flank_cells = FLANK_CELLS,
    threads: 1
    log:
        f"{D_LOGS}/merge_pilot.log",
    shell:
        "{RSCRIPT} {input.script} {output} {params.flank_cells} {input.rds} > {log} 2>&1"
# --- Stage 0b: single-pass merge of all 20 flank scored.rds (sparse, no densify) ---
rule merge:
    message: "merge: single-pass merge of all 20 flank norm.rds (sparse, no densify)"
    input:
        script = f"{R_SCRIPTS}/merge.R",
        ss     = MASTER,
        rds    = expand(f"{D_NORM}/{{s}}.norm.rds", s=FLANK),
    output:
        rds     = f"{D_AGG}/merged.rds",
        summary = "results/aggregate/merge_summary.tsv",
        meta    = f"{D_AGG}/cell_metadata.tsv",
    params:
        scale_factor = config["normalize"]["scale_factor"],
    threads: 1
    log:
        f"{D_LOGS}/merge.log",
    shell:
        "{RSCRIPT} {input.script} {input.ss} {params.scale_factor} {output.rds} {output.summary} {output.meta} {input.rds} > {log} 2>&1"
# --- Tier-1/2 input: per-cell negative-probe background (feeds the scVI export) ---
# Negprobes (the platform background model) were dropped during 950-gene harmonization;
# recover from raw inputs. M01: metadata parquet nCount_NegativeProbes -- the M01 raw
# per-slide RDS DO carry negprobes, but their object-local barcodes (c_1_*) do NOT match
# the adapter's scored barcodes (0% overlap), so the raw-RDS retire stays DEFERRED; the
# parquet neg is verified byte-identical to the RDS, so this is invariant (see
# typing_provenance.md S2 / plan-mutter01-controls.md). M02: raw slide RDS negprobes assay.
# neg = nCount_neg / 10. Output barcodes match merge.R's {sample_id}_{raw_barcode} key;
# consumed by full_export.R as the scVI continuous covariate.
rule recover_negprobes:
    message: "recover_negprobes: recover per-cell negative-probe background from raw inputs"
    input:
        script = f"{R_SCRIPTS}/recover_negprobes.R",
        m01    = M01_META,
        m02    = expand(f"{M02_RAW}/seuratObject_0{{i}}_Mutter_02_CosMmR.RDS", i=[1,2,3,4]),
        rds    = expand(f"{D_NORM}/{{s}}.norm.rds", s=FLANK),
    output:
        neg = f"{D_AGG}/cell_neg.tsv",
    params:
        m02_dir = M02_RAW,
    threads: 1
    log:
        f"{D_LOGS}/recover_negprobes.log",
    shell:
        "{RSCRIPT} {input.script} {output.neg} {input.m01} {params.m02_dir} "
        "{input.rds} > {log} 2>&1"
# --- T1.0: export merged counts (MTX) + obs (+ neg covariate) for scVI ---
rule full_export:
    message: "full_export: export merged counts (MTX) + obs (+ neg covariate) for scVI"
    input:
        script = f"{R_SCRIPTS}/full_export.R",
        rds    = f"{D_AGG}/merged.rds",
        neg    = f"{D_AGG}/cell_neg.tsv",
    output:
        mtx   = f"{D_FULL}/mtx/counts.mtx",
        feats = f"{D_FULL}/mtx/features.tsv",
        bars  = f"{D_FULL}/mtx/barcodes.tsv",
        obs   = f"{D_FULL}/obs.parquet",
    params:
        outdir = D_FULL,
    threads: 1
    log:
        f"{D_LOGS}/full_export.log",
    shell:
        "{RSCRIPT} {input.script} {input.rds} {input.neg} {params.outdir} > {log} 2>&1"
# --- T1.1: scVI integration (GPU-direct on the RTX A4000) -> 30-dim latent ---
# Params pinned to the locked scvi_model/model.pt (typing_provenance.md S3). GPU resource
# serializes this rule; run the workflow with --resources gpu=1.
rule full_scvi:
    message: "full_scvi: scVI integration (GPU-direct) -> 30-dim latent"
    input:
        script = f"{R_SCRIPTS}/full_scvi.py",
        mtx    = f"{D_FULL}/mtx/counts.mtx",
        obs    = f"{D_FULL}/obs.parquet",
    output:
        latent = f"{D_FULL}/scvi_latent.parquet",
        model  = directory(f"{D_FULL}/scvi_model"),
    params:
        full = D_FULL,
    threads: 4
    resources:
        gpu = 1,
    log:
        f"{D_LOGS}/full_scvi.log",
    shell:
        "{PYSCVI} {input.script} {params.full} > {log} 2>&1"
# --- T1.2: Leiden sweep on the scVI latent + cluster-level coarse annotation + Q4 gates ---
# out prefix "results/aggregate/full" -> per-res full_r{tag}_*; finalize consumes res=3.0.
rule full_cluster:
    message: "full_cluster: Leiden resolution sweep + cluster-level coarse annotation + Q4 gates"
    input:
        script  = f"{R_SCRIPTS}/full_cluster.py",
        latent  = f"{D_FULL}/scvi_latent.parquet",
        mtx     = f"{D_FULL}/mtx/counts.mtx",
        obs     = f"{D_FULL}/obs.parquet",
        markers = MARKERS,
    output:
        ckpt   = f"{D_FULL}/cluster_checkpoint.h5ad",
        labels = "results/aggregate/full_r3_coarse_labels.parquet",
        qc     = "results/aggregate/full_r3_cluster_qc.tsv",
        gates  = "results/aggregate/full_r3_gates.json",
        sweep  = "results/aggregate/full_sweep_summary.tsv",
    params:
        full   = D_FULL,
        prefix = "results/aggregate/full",
        res    = "1.0,2.0,3.0",
    threads: 4
    log:
        f"{D_LOGS}/full_cluster.log",
    shell:
        "{PYSCVI} {input.script} {params.full} {input.markers} {params.prefix} "
        "{params.res} > {log} 2>&1"
# --- T1.3: finalize tier-1 at the adopted res=3.0 -> canonical coarse labels ---
rule finalize_tier1:
    message: "finalize_tier1: finalize tier-1 at adopted res=3.0 -> canonical coarse labels"
    input:
        script = f"{R_SCRIPTS}/finalize_tier1.py",
        labels = "results/aggregate/full_r3_coarse_labels.parquet",
        qc     = "results/aggregate/full_r3_cluster_qc.tsv",
        gates  = "results/aggregate/full_r3_gates.json",
    output:
        coarse  = "results/aggregate/full_coarse_labels.parquet",
        gates   = "results/aggregate/full_gates.json",
        summary = "results/aggregate/full_coarse_summary.tsv",
    params:
        dir = "results/aggregate",
    threads: 1
    log:
        f"{D_LOGS}/finalize_tier1.log",
    shell:
        "{PYSCVI} {input.script} {params.dir} > {log} 2>&1"
# --- T2.immune.a: subcluster the immune compartment (res=3.0 substrate for SingleR) ---
rule tier2_immune_subcluster:
    message: "tier2_immune_subcluster: subcluster immune compartment (res=3.0 substrate for SingleR)"
    input:
        script = f"{R_SCRIPTS}/tier2_immune_subcluster.py",
        ckpt   = f"{D_FULL}/cluster_checkpoint.h5ad",
        labels = "results/aggregate/full_coarse_labels.parquet",
    output:
        sub   = f"{D_FULL}/immune_subclusters.parquet",
        mtx   = f"{D_FULL}/immune_mtx/counts.mtx",
        feats = f"{D_FULL}/immune_mtx/features.tsv",
        bars  = f"{D_FULL}/immune_mtx/barcodes.tsv",
        ckpt  = f"{D_FULL}/immune_checkpoint.h5ad",
    params:
        outdir = D_FULL,
        res    = "3.0",
    threads: 4
    log:
        f"{D_LOGS}/tier2_immune_subcluster.log",
    shell:
        "{PYSCVI} {input.script} {input.ckpt} {input.labels} {params.outdir} "
        "{params.res} > {log} 2>&1"
# --- T2.immune.b: finer res=6.0 subclusters (rare-lineage rescue substrate) ---
rule tier2_immune_subcluster_r6:
    message: "tier2_immune_subcluster_r6: finer res=6.0 immune subclusters (rare-lineage rescue substrate)"
    input:
        script = f"{R_SCRIPTS}/tier2_immune_subcluster.py",
        ckpt   = f"{D_FULL}/cluster_checkpoint.h5ad",
        labels = "results/aggregate/full_coarse_labels.parquet",
    output:
        sub = f"{D_FULL}/immune_r6/immune_subclusters.parquet",
    params:
        outdir = f"{D_FULL}/immune_r6",
        res    = "6.0",
    threads: 4
    log:
        f"{D_LOGS}/tier2_immune_subcluster_r6.log",
    shell:
        "{PYSCVI} {input.script} {input.ckpt} {input.labels} {params.outdir} "
        "{params.res} > {log} 2>&1"
# --- T2.immune.c: cluster-level SingleR vs celldex/ImmGen on the res=3.0 subclusters ---
rule tier2_singler:
    message: "tier2_singler: cluster-level SingleR vs celldex/ImmGen on res=3.0 subclusters"
    input:
        script = f"{R_SCRIPTS}/tier2_singler.R",
        mtx    = f"{D_FULL}/immune_mtx/counts.mtx",
        sub    = f"{D_FULL}/immune_subclusters.parquet",
    output:
        singler  = f"{D_FULL}/immune_subtype_singler.tsv",
        subtypes = f"{D_FULL}/immune_subtypes.parquet",
    params:
        mtxdir = f"{D_FULL}/immune_mtx",
        outdir = D_FULL,
    threads: 1
    log:
        f"{D_LOGS}/tier2_singler.log",
    shell:
        "{RSCRIPT} {input.script} {params.mtxdir} {input.sub} {params.outdir} > {log} 2>&1"
# --- T2.immune.d: marker-rescue (B/Plasma/Neutrophil) -- res=6.0 identity-marker override ---
rule tier2_immune_rescue:
    message: "tier2_immune_rescue: identity-marker rescue (B/Plasma/Neutrophil) via res=6.0 override"
    input:
        script  = f"{R_SCRIPTS}/tier2_immune_rescue.py",
        mtx     = f"{D_FULL}/immune_mtx/counts.mtx",
        obs     = f"{D_FULL}/obs.parquet",
        sub6    = f"{D_FULL}/immune_r6/immune_subclusters.parquet",
        singler = f"{D_FULL}/immune_subtypes.parquet",
    output:
        rescued = f"{D_FULL}/immune_subtypes_rescued.parquet",
        audit   = f"{D_FULL}/immune_rescue_audit.tsv",
    params:
        mtxdir = f"{D_FULL}/immune_mtx",
    threads: 1
    log:
        f"{D_LOGS}/tier2_immune_rescue.log",
    shell:
        "{PYSCVI} {input.script} {params.mtxdir} {input.obs} {input.sub6} "
        "{input.singler} {output.rescued} {output.audit} > {log} 2>&1"
# --- T2.stroma.a: subcluster stroma; sweep res (r2 base for UCell, r6 for rescue) ---
rule tier2_stroma_subcluster:
    message: "tier2_stroma_subcluster: subcluster stroma; sweep res (r2 base for UCell, r6 for rescue)"
    input:
        script = f"{R_SCRIPTS}/tier2_stroma_subcluster.py",
        ckpt   = f"{D_FULL}/cluster_checkpoint.h5ad",
        labels = "results/aggregate/full_coarse_labels.parquet",
    output:
        mtx   = f"{D_FULL}/stroma_mtx/counts.mtx",
        r2    = f"{D_FULL}/stroma_r2_subclusters.parquet",
        r6    = f"{D_FULL}/stroma_r6_subclusters.parquet",
        ckpt  = f"{D_FULL}/stroma_checkpoint.h5ad",
        sweep = "results/aggregate/stroma_sweep_summary.tsv",
    params:
        ddir   = D_FULL,
        prefix = "results/aggregate/stroma",
        res    = "0.5,1,2,3,4,6",
    threads: 4
    log:
        f"{D_LOGS}/tier2_stroma_subcluster.log",
    shell:
        "{PYSCVI} {input.script} {input.ckpt} {input.labels} {params.ddir} "
        "{params.prefix} {params.res} > {log} 2>&1"
# --- T2.stroma.b: UCell cluster-level argmax+margin on the r2 base subclusters ---
rule tier2_stroma_ucell:
    message: "tier2_stroma_ucell: UCell cluster-level argmax+margin on r2 base subclusters"
    input:
        script  = f"{R_SCRIPTS}/tier2_stroma_ucell.R",
        sub     = f"{D_FULL}/stroma_r2_subclusters.parquet",
        mtx     = f"{D_FULL}/stroma_mtx/counts.mtx",
        markers = MARKERS,
    output:
        subtypes = f"{D_FULL}/stroma_subtypes.parquet",
        summary  = f"{D_FULL}/stroma_subtype_summary.tsv",
        ucell    = f"{D_FULL}/stroma_subcluster_ucell.tsv",
    params:
        mtxdir = f"{D_FULL}/stroma_mtx",
        outdir = D_FULL,
    threads: 4
    log:
        f"{D_LOGS}/tier2_stroma_ucell.log",
    shell:
        "{RSCRIPT} {input.script} {params.mtxdir} {input.sub} {input.markers} "
        "{params.outdir} > {log} 2>&1"
# --- T2.stroma.c: identity-marker rescue (Endothelial/Adipocyte) -- r6 subclusters ---
rule tier2_stroma_rescue:
    message: "tier2_stroma_rescue: identity-marker rescue (Endothelial/Adipocyte) via r6 subclusters"
    input:
        script = f"{R_SCRIPTS}/tier2_stroma_rescue.py",
        mtx    = f"{D_FULL}/stroma_mtx/counts.mtx",
        obs    = f"{D_FULL}/obs.parquet",
        sub6   = f"{D_FULL}/stroma_r6_subclusters.parquet",
        base   = f"{D_FULL}/stroma_subtypes.parquet",
    output:
        rescued = f"{D_FULL}/stroma_subtypes_rescued.parquet",
        audit   = f"{D_FULL}/stroma_rescue_audit.tsv",
    params:
        mtxdir = f"{D_FULL}/stroma_mtx",
    threads: 1
    log:
        f"{D_LOGS}/tier2_stroma_rescue.log",
    shell:
        "{PYSCVI} {input.script} {params.mtxdir} {input.obs} {input.sub6} "
        "{input.base} {output.rescued} {output.audit} > {log} 2>&1"
# --- T2.unify: join tier-1 coarse + tier-2 immune + tier-2 stroma -> canonical labels ---
rule unify_labels:
    message: "unify_labels: join tier-1 coarse + tier-2 immune + tier-2 stroma -> canonical labels"
    input:
        script = f"{R_SCRIPTS}/unify_labels.py",
        coarse = "results/aggregate/full_coarse_labels.parquet",
        immune = f"{D_FULL}/immune_subtypes_rescued.parquet",
        stroma = f"{D_FULL}/stroma_subtypes_rescued.parquet",
    output:
        labels  = "results/aggregate/full_labels.parquet",
        summary = "results/aggregate/full_labels_summary.tsv",
    threads: 1
    log:
        f"{D_LOGS}/unify_labels.log",
    shell:
        "{PYSCVI} {input.script} {input.coarse} {input.immune} {input.stroma} "
        "{output.labels} {output.summary} > {log} 2>&1"
