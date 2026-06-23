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
PYSCVI  = "conda run -n spatial-rads-scvi python"   # tier-1/2 Python (scVI/scanpy, GPU)
DATADIR = config["datadir"]
MASTER  = config["samplesheet"]                  # results/data_model/samples.tsv
AGG     = f"{DATADIR}/aggregate"                  # heavy intermediates
FULL    = f"{AGG}/full"                           # tier-1/2 typing intermediate dir
SCORED  = f"{DATADIR}/processing/scored"          # per-sample inputs
PANEL   = "results/data_model/common_genes.tsv"   # 950 shared-panel gene list (built in data_model.smk)
MARKERS = "config/lineage_markers.yaml"           # coarse + tier-2 lineage marker sets
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
        "results/aggregate/full_coarse_labels.parquet",   # tier-1 coarse (finalized)
        "results/aggregate/full_labels.parquet",          # canonical unified per-cell labels
        # --- MBRT-vs-SBRT downstream differential layer (plan-mbrt-vs-sbrt-impl.md) ---
        "results/data_model/gene_set_panel_coverage.tsv",   # built in data_model.smk
        "results/aggregate/readout_detection_m02.tsv",
        "results/aggregate/celltype_qc_markers.tsv",
        "results/aggregate/concordance_m01_m02.tsv",
        "results/aggregate/geneset_overlap.tsv",            # Fix 5: set-overlap diagnostic
        "results/aggregate/fibroblast_substate.parquet",    # Fix 3: fibroblast resting/activated split
        "results/aggregate/composition_substate_test_m02day2.tsv",  # Fix 3: sub-state propeller
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

# ============================================================================
# TIER-1 + TIER-2 cell typing: the real scVI cluster-then-annotate chain
# (Step 3 of plan-aggregate-refactor.md -- replaces the dead per-cell InSituType
# chain that validated at 85% mistyped tumor). Multi-env: R (spatial-rads) ->
# GPU Python (spatial-rads-scvi) -> R -> Python. Method LOCKED to the current
# stack by the Workstream B review (plan-aggregate.md v2.7: keep scVI / swept
# res=3.0 / per-sample MECR). results/aggregate/full_labels.parquet is the
# canonical per-cell output; acceptance = re-validation of the 4 Q4 gates, not
# byte-identity (scVI-GPU + Leiden-igraph are not byte-reproducible).
# ============================================================================

# --- T1.0: export merged counts (MTX) + obs (+ neg covariate) for scVI ---
rule full_export:
    input:
        script = "scripts/aggregate/full_export.R",
        rds    = f"{AGG}/merged.rds",
        neg    = f"{AGG}/cell_neg.tsv",
    output:
        mtx   = f"{FULL}/mtx/counts.mtx",
        feats = f"{FULL}/mtx/features.tsv",
        bars  = f"{FULL}/mtx/barcodes.tsv",
        obs   = f"{FULL}/obs.parquet",
    params:
        outdir = FULL,
    threads: 1
    log:
        "logs/aggregate/full_export.log",
    shell:
        "{RSCRIPT} {input.script} {input.rds} {input.neg} {params.outdir} > {log} 2>&1"

# --- T1.1: scVI integration (GPU-direct on the RTX A4000) -> 30-dim latent ---
# Params pinned to the locked scvi_model/model.pt (typing_provenance.md S3). GPU resource
# serializes this rule; run the workflow with --resources gpu=1.
rule full_scvi:
    input:
        script = "scripts/aggregate/full_scvi.py",
        mtx    = f"{FULL}/mtx/counts.mtx",
        obs    = f"{FULL}/obs.parquet",
    output:
        latent = f"{FULL}/scvi_latent.parquet",
        model  = directory(f"{FULL}/scvi_model"),
    params:
        full = FULL,
    threads: 4
    resources:
        gpu = 1,
    log:
        "logs/aggregate/full_scvi.log",
    shell:
        "{PYSCVI} {input.script} {params.full} > {log} 2>&1"

# --- T1.2: Leiden sweep on the scVI latent + cluster-level coarse annotation + Q4 gates ---
# out prefix "results/aggregate/full" -> per-res full_r{tag}_*; finalize consumes res=3.0.
rule full_cluster:
    input:
        script  = "scripts/aggregate/full_cluster.py",
        latent  = f"{FULL}/scvi_latent.parquet",
        mtx     = f"{FULL}/mtx/counts.mtx",
        obs     = f"{FULL}/obs.parquet",
        markers = MARKERS,
    output:
        ckpt   = f"{FULL}/cluster_checkpoint.h5ad",
        labels = "results/aggregate/full_r3_coarse_labels.parquet",
        qc     = "results/aggregate/full_r3_cluster_qc.tsv",
        gates  = "results/aggregate/full_r3_gates.json",
        sweep  = "results/aggregate/full_sweep_summary.tsv",
    params:
        full   = FULL,
        prefix = "results/aggregate/full",
        res    = "1.0,2.0,3.0",
    threads: 4
    log:
        "logs/aggregate/full_cluster.log",
    shell:
        "{PYSCVI} {input.script} {params.full} {input.markers} {params.prefix} "
        "{params.res} > {log} 2>&1"

# --- T1.3: finalize tier-1 at the adopted res=3.0 -> canonical coarse labels ---
rule finalize_tier1:
    input:
        script = "scripts/aggregate/finalize_tier1.py",
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
        "logs/aggregate/finalize_tier1.log",
    shell:
        "{PYSCVI} {input.script} {params.dir} > {log} 2>&1"

# --- T2.immune.a: subcluster the immune compartment (res=3.0 substrate for SingleR) ---
rule tier2_immune_subcluster:
    input:
        script = "scripts/aggregate/tier2_immune_subcluster.py",
        ckpt   = f"{FULL}/cluster_checkpoint.h5ad",
        labels = "results/aggregate/full_coarse_labels.parquet",
    output:
        sub   = f"{FULL}/immune_subclusters.parquet",
        mtx   = f"{FULL}/immune_mtx/counts.mtx",
        feats = f"{FULL}/immune_mtx/features.tsv",
        bars  = f"{FULL}/immune_mtx/barcodes.tsv",
        ckpt  = f"{FULL}/immune_checkpoint.h5ad",
    params:
        outdir = FULL,
        res    = "3.0",
    threads: 4
    log:
        "logs/aggregate/tier2_immune_subcluster.log",
    shell:
        "{PYSCVI} {input.script} {input.ckpt} {input.labels} {params.outdir} "
        "{params.res} > {log} 2>&1"

# --- T2.immune.b: finer res=6.0 subclusters (rare-lineage rescue substrate) ---
rule tier2_immune_subcluster_r6:
    input:
        script = "scripts/aggregate/tier2_immune_subcluster.py",
        ckpt   = f"{FULL}/cluster_checkpoint.h5ad",
        labels = "results/aggregate/full_coarse_labels.parquet",
    output:
        sub = f"{FULL}/immune_r6/immune_subclusters.parquet",
    params:
        outdir = f"{FULL}/immune_r6",
        res    = "6.0",
    threads: 4
    log:
        "logs/aggregate/tier2_immune_subcluster_r6.log",
    shell:
        "{PYSCVI} {input.script} {input.ckpt} {input.labels} {params.outdir} "
        "{params.res} > {log} 2>&1"

# --- T2.immune.c: cluster-level SingleR vs celldex/ImmGen on the res=3.0 subclusters ---
rule tier2_singler:
    input:
        script = "scripts/aggregate/tier2_singler.R",
        mtx    = f"{FULL}/immune_mtx/counts.mtx",
        sub    = f"{FULL}/immune_subclusters.parquet",
    output:
        singler  = f"{FULL}/immune_subtype_singler.tsv",
        subtypes = f"{FULL}/immune_subtypes.parquet",
    params:
        mtxdir = f"{FULL}/immune_mtx",
        outdir = FULL,
    threads: 1
    log:
        "logs/aggregate/tier2_singler.log",
    shell:
        "{RSCRIPT} {input.script} {params.mtxdir} {input.sub} {params.outdir} > {log} 2>&1"

# --- T2.immune.d: marker-rescue (B/Plasma/Neutrophil) -- res=6.0 identity-marker override ---
rule tier2_immune_rescue:
    input:
        script  = "scripts/aggregate/tier2_immune_rescue.py",
        mtx     = f"{FULL}/immune_mtx/counts.mtx",
        obs     = f"{FULL}/obs.parquet",
        sub6    = f"{FULL}/immune_r6/immune_subclusters.parquet",
        singler = f"{FULL}/immune_subtypes.parquet",
    output:
        rescued = f"{FULL}/immune_subtypes_rescued.parquet",
        audit   = f"{FULL}/immune_rescue_audit.tsv",
    params:
        mtxdir = f"{FULL}/immune_mtx",
    threads: 1
    log:
        "logs/aggregate/tier2_immune_rescue.log",
    shell:
        "{PYSCVI} {input.script} {params.mtxdir} {input.obs} {input.sub6} "
        "{input.singler} {output.rescued} {output.audit} > {log} 2>&1"

# --- T2.stroma.a: subcluster stroma; sweep res (r2 base for UCell, r6 for rescue) ---
rule tier2_stroma_subcluster:
    input:
        script = "scripts/aggregate/tier2_stroma_subcluster.py",
        ckpt   = f"{FULL}/cluster_checkpoint.h5ad",
        labels = "results/aggregate/full_coarse_labels.parquet",
    output:
        mtx   = f"{FULL}/stroma_mtx/counts.mtx",
        r2    = f"{FULL}/stroma_r2_subclusters.parquet",
        r6    = f"{FULL}/stroma_r6_subclusters.parquet",
        ckpt  = f"{FULL}/stroma_checkpoint.h5ad",
        sweep = "results/aggregate/stroma_sweep_summary.tsv",
    params:
        ddir   = FULL,
        prefix = "results/aggregate/stroma",
        res    = "0.5,1,2,3,4,6",
    threads: 4
    log:
        "logs/aggregate/tier2_stroma_subcluster.log",
    shell:
        "{PYSCVI} {input.script} {input.ckpt} {input.labels} {params.ddir} "
        "{params.prefix} {params.res} > {log} 2>&1"

# --- T2.stroma.b: UCell cluster-level argmax+margin on the r2 base subclusters ---
rule tier2_stroma_ucell:
    input:
        script  = "scripts/aggregate/tier2_stroma_ucell.R",
        sub     = f"{FULL}/stroma_r2_subclusters.parquet",
        mtx     = f"{FULL}/stroma_mtx/counts.mtx",
        markers = MARKERS,
    output:
        subtypes = f"{FULL}/stroma_subtypes.parquet",
        summary  = f"{FULL}/stroma_subtype_summary.tsv",
        ucell    = f"{FULL}/stroma_subcluster_ucell.tsv",
    params:
        mtxdir = f"{FULL}/stroma_mtx",
        outdir = FULL,
    threads: 4
    log:
        "logs/aggregate/tier2_stroma_ucell.log",
    shell:
        "{RSCRIPT} {input.script} {params.mtxdir} {input.sub} {input.markers} "
        "{params.outdir} > {log} 2>&1"

# --- T2.stroma.c: identity-marker rescue (Endothelial/Adipocyte) -- r6 subclusters ---
rule tier2_stroma_rescue:
    input:
        script = "scripts/aggregate/tier2_stroma_rescue.py",
        mtx    = f"{FULL}/stroma_mtx/counts.mtx",
        obs    = f"{FULL}/obs.parquet",
        sub6   = f"{FULL}/stroma_r6_subclusters.parquet",
        base   = f"{FULL}/stroma_subtypes.parquet",
    output:
        rescued = f"{FULL}/stroma_subtypes_rescued.parquet",
        audit   = f"{FULL}/stroma_rescue_audit.tsv",
    params:
        mtxdir = f"{FULL}/stroma_mtx",
    threads: 1
    log:
        "logs/aggregate/tier2_stroma_rescue.log",
    shell:
        "{PYSCVI} {input.script} {params.mtxdir} {input.obs} {input.sub6} "
        "{input.base} {output.rescued} {output.audit} > {log} 2>&1"

# --- T2.unify: join tier-1 coarse + tier-2 immune + tier-2 stroma -> canonical labels ---
rule unify_labels:
    input:
        script = "scripts/aggregate/unify_labels.py",
        coarse = "results/aggregate/full_coarse_labels.parquet",
        immune = f"{FULL}/immune_subtypes_rescued.parquet",
        stroma = f"{FULL}/stroma_subtypes_rescued.parquet",
    output:
        labels  = "results/aggregate/full_labels.parquet",
        summary = "results/aggregate/full_labels_summary.tsv",
    threads: 1
    log:
        "logs/aggregate/unify_labels.log",
    shell:
        "{PYSCVI} {input.script} {input.coarse} {input.immune} {input.stroma} "
        "{output.labels} {output.summary} > {log} 2>&1"

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
        sensitivity = "results/aggregate/composition_unassigned_sensitivity.tsv",
    threads: 1
    log:
        "logs/aggregate/composition.log",
    shell:
        "{RSCRIPT} {input.script} {input.obs} {input.labels} {output.by_sample} "
        "{output.test} {output.dropped} {output.bars} {output.forest} "
        "{output.timecourse} {output.sensitivity} > {log} 2>&1"

# --- Track 2: pseudobulk construction (M02 day2, sample x cell_type count sums) ---
rule pseudobulk_build:
    input:
        script = "scripts/aggregate/pseudobulk_build.R",
        rds    = f"{AGG}/merged.rds",   # re-pointed off the deleted merged_typed.rds (count-only consumer; labels from full_labels.parquet)
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

# --- Track 2: detection-vs-level decomposition + sparsity-aware call_class (Fix 2) ---
rule detection_test:
    input:
        script = "scripts/aggregate/detection_test.R",
        se     = f"{AGG}/pseudobulk_se.rds",
    output:
        det = "results/aggregate/detection_test_m02day2.tsv",
    threads: 1
    log:
        "logs/aggregate/detection_test.log",
    shell:
        "{RSCRIPT} {input.script} {input.se} {output.det} > {log} 2>&1"

# --- Track 2 pathway: GSEA on pseudobulk stat-ranked genes (primary + Hallmark) ---
rule gsea:
    input:
        script = "scripts/aggregate/gsea.R",
        degs   = "results/aggregate/degs_pseudobulk_m02day2.tsv",
        sets   = "results/data_model/pathway_sets.tsv",
    output:
        gsea = "results/aggregate/gsea_pseudobulk_m02day2.tsv",
    threads: 1
    log:
        "logs/aggregate/gsea.log",
    shell:
        "{RSCRIPT} {input.script} {input.degs} {input.sets} {output.gsea} > {log} 2>&1"

# --- Track 2 pathway: per-cell UCell + AddModuleScore scoring + M02 limma test ---
# threads: 4 -- UCell runs ncores=4 fork (BiocParallel) parallelism, not BLAS, so the
# threads:1 BLAS-hygiene convention is unaffected (BLAS stays pinned to 1 thread).
rule pathway_summary:
    input:
        script = "scripts/aggregate/pathway_summary.R",
        rds    = f"{AGG}/merged.rds",   # re-pointed off the deleted merged_typed.rds (count-only consumer; labels from full_labels.parquet)
        labels = "results/aggregate/full_labels.parquet",
        sets   = "results/data_model/pathway_sets.tsv",
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
        "{RSCRIPT} {input.script} {input.rds} {input.labels} {input.sets} {output.summary} "
        "{output.test} {output.conc} {output.heatmap} {output.timecourse} "
        "{output.scatter} > {log} 2>&1"

# ============================================================================
# MBRT-vs-SBRT downstream differential layer (plan-mbrt-vs-sbrt-impl.md T4-T13).
# All consume the unified per-cell labels (results/aggregate/full_labels.parquet)
# and the day-2 readout tables, terminating in the tier-tagged results_master.tsv.
# ============================================================================

# --- T5: per-cell-type readout detection report (panel detectability in M02) ---
rule panel_coverage:
    input:
        script = "scripts/aggregate/panel_coverage.R",
        rds    = f"{AGG}/merged.rds",   # re-pointed off the deleted merged_typed.rds (count-only consumer; labels from full_labels.parquet)
        labels = "results/aggregate/full_labels.parquet",
        sets   = "results/data_model/pathway_sets.tsv",
    output:
        detection = "results/aggregate/readout_detection_m02.tsv",
    threads: 1
    log:
        "logs/aggregate/panel_coverage.log",
    shell:
        "{RSCRIPT} {input.script} {input.rds} {input.labels} {input.sets} "
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
rule substate_split:
    input:
        script  = "scripts/aggregate/substate_split.R",
        rds     = f"{AGG}/merged.rds",
        labels  = "results/aggregate/full_labels.parquet",
        markers = "config/substate_markers.yaml",
    output:
        parquet = "results/aggregate/fibroblast_substate.parquet",
        gate    = "results/aggregate/substate_gate_report.tsv",
    threads: 4
    log:
        "logs/aggregate/substate_split.log",
    shell:
        "{RSCRIPT} {input.script} {input.rds} {input.labels} {input.markers} "
        "{output.parquet} {output.gate} > {log} 2>&1"

rule substate_composition:
    input:
        script   = "scripts/aggregate/substate_composition.R",
        obs      = f"{AGG}/full/obs.parquet",
        labels   = "results/aggregate/full_labels.parquet",
        substate = "results/aggregate/fibroblast_substate.parquet",
    output:
        test = "results/aggregate/composition_substate_test_m02day2.tsv",
    threads: 1
    log:
        "logs/aggregate/substate_composition.log",
    shell:
        "{RSCRIPT} {input.script} {input.obs} {input.labels} {input.substate} {output.test} > {log} 2>&1"

rule geneset_overlap:
    input:
        script = "scripts/aggregate/geneset_overlap.R",
        sets   = "results/data_model/pathway_sets.tsv",
    output:
        overlap = "results/aggregate/geneset_overlap.tsv",
    threads: 1
    log:
        "logs/aggregate/geneset_overlap.log",
    shell:
        "{RSCRIPT} {input.script} {input.sets} {output.overlap} > {log} 2>&1"

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
        sets    = "results/data_model/pathway_sets.tsv",
        cov     = "results/data_model/gene_set_panel_coverage.tsv",
        det     = "results/aggregate/detection_test_m02day2.tsv",
    output:
        master = "results/aggregate/results_master.tsv",
        detect = "results/aggregate/detectability_summary.tsv",
    threads: 1
    log:
        "logs/aggregate/assemble_results.log",
    shell:
        "{RSCRIPT} {input.script} {input.comp} {input.degs} {input.gsea} "
        "{input.pathway} {input.niche} {input.mixing} {input.myeloid} "
        "{input.mde} {input.sets} {input.cov} {input.det} {output.master} > {log} 2>&1"
