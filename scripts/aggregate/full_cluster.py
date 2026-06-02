#!/usr/bin/env python3
"""Full-cohort coarse typing on the scVI latent (cluster-then-annotate, v2.3/v2.4).
Restores (or builds + checkpoints) the 3.27M-cell joint neighbor graph once, then sweeps
Leiden resolutions, annotating each cluster (not each cell) by argmax over coarse lineage
marker means and running the four pre-registered Q4 gates (plan-aggregate.md) at every res:
  Q4.1 premise   -- >3 major clusters, Immune >5%, Stroma >10%, no label >85%
  Q4.2 batch-mix -- scaled iLISI(dataset) >=0.5 AND minority dataset >=5% in every cluster >1%
  Q4.3 recall    -- per-lineage gold-marker recall at cluster mean, negprobe-relative, split M01/M02
  Q4.4 negprobe  -- per-cluster negprobe fraction (standing QC)
The neighbor graph (the ~10-min step) is cached to cluster_checkpoint.h5ad so the sweep,
and any later re-tries, cost only the bounded Leiden + annotation per resolution. scaled
iLISI is embedding-level (res-independent) so it is computed once. Emits per-resolution
QC/recall/gates/labels plus a sweep summary. Cluster-level annotation only; tier-2 and
downstream tracks stay parked pending review of the sweep.
Usage: full_cluster.py <full_dir> <lineage_markers.yaml> <out_prefix> [res_csv=1.0,2.0,3.0]
"""
import sys, os, json
import numpy as np
import pandas as pd
import scipy.io
import anndata as ad
import scanpy as sc
import yaml

FDIR, MARKERS_YAML, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
RES_LIST = [float(x) for x in (sys.argv[4].split(",") if len(sys.argv) > 4
                               else ["1.0", "2.0", "3.0"])]
SEED  = 0
N_SUB = 50000                      # subsample for the iLISI batch-mixing metric
CKPT  = os.path.join(FDIR, "cluster_checkpoint.h5ad")
np.random.seed(SEED)

COMPARTMENT = {
    "Epithelial": "tumor",
    "T": "immune", "B": "immune", "Plasma": "immune", "NK": "immune",
    "Macrophage": "immune", "Neutrophil": "immune", "Mast": "immune", "DC": "immune",
    "Endothelial": "stroma", "Fibroblast": "stroma", "Pericyte": "stroma",
    "SmoothMuscle": "stroma", "Adipocyte": "stroma",
}

# ---- restore checkpoint, or build it (raw counts X + obs + scVI latent + neighbors) ----
if os.path.exists(CKPT):
    print(f"restoring checkpoint {CKPT} ...", flush=True)
    adata = ad.read_h5ad(CKPT)
else:
    feats = [l.strip() for l in open(os.path.join(FDIR, "mtx", "features.tsv"))]
    bcs   = [l.strip() for l in open(os.path.join(FDIR, "mtx", "barcodes.tsv"))]
    X = scipy.io.mmread(os.path.join(FDIR, "mtx", "counts.mtx")).T.tocsr()
    adata = ad.AnnData(X=X); adata.var_names = feats; adata.obs_names = bcs
    obs = pd.read_parquet(os.path.join(FDIR, "obs.parquet")).set_index("cell").loc[adata.obs_names]
    adata.obs["dataset"]   = obs["dataset"].astype(str).values
    adata.obs["slide_id"]  = obs["slide_id"].astype(str).values
    adata.obs["cell_type"] = obs["cell_type"].fillna("unassigned").astype(str).values
    adata.obs["neg"]       = obs["neg"].astype(float).values
    lat = pd.read_parquet(os.path.join(FDIR, "scvi_latent.parquet")).set_index("cell").loc[adata.obs_names]
    adata.obsm["X_scvi"] = lat.values.astype(np.float32)
    print(f"loaded {adata.n_obs} cells x {adata.n_vars} genes | scVI {lat.shape[1]}d | "
          f"datasets {dict(adata.obs['dataset'].value_counts())}", flush=True)
    sc.pp.neighbors(adata, use_rep="X_scvi", n_neighbors=15, random_state=SEED)
    print("neighbors done; writing checkpoint ...", flush=True)
    adata.write_h5ad(CKPT)
print(f"adata: {adata.n_obs} cells x {adata.n_vars} genes", flush=True)
tot = adata.n_obs

# ---- per-cell lineage mean raw count (resolution-independent) ----
markers = yaml.safe_load(open(MARKERS_YAML))
LM = {}
for lin, genes in markers.items():
    gi = [adata.var_names.get_loc(g) for g in genes if g in adata.var_names]
    LM[lin] = np.asarray(adata[:, gi].X.mean(axis=1)).ravel() if gi else np.zeros(adata.n_obs)
LM = pd.DataFrame(LM, index=adata.obs_names)
neg = adata.obs["neg"].values.astype(float)
ds  = adata.obs["dataset"].astype(str).values

# ---- scaled iLISI once (embedding-level, res-independent) ----
idx = (adata.obs.groupby("dataset", group_keys=False)
       .apply(lambda d: d.sample(min(len(d), int(round(N_SUB * len(d) / tot))),
                                 random_state=SEED)).index)
sub = adata[idx].copy()
p = sub.obs["dataset"].value_counts(normalize=True)
ceiling = 1.0 / float(np.sum(p.values ** 2))
from scib_metrics.benchmark import Benchmarker
bm = Benchmarker(sub, batch_key="dataset", label_key="cell_type",
                 embedding_obsm_keys=["X_scvi"], n_jobs=-1)
bm.benchmark()
scib = bm.get_results(min_max_scale=False)
def col(d, needle):
    for c in d.columns:
        if needle.lower() in str(c).lower():
            return c
    return None
raw_ilisi = float(scib.loc["X_scvi", col(scib, "ilisi")]) + 1.0   # scib iLISI=(raw-1)/(nb-1), nb=2
scaled_ilisi = float(np.clip((raw_ilisi - 1) / (ceiling - 1), 0, 1))
print(f"iLISI: subsample {sub.n_obs} frac {dict(p.round(3))} ceiling {ceiling:.3f} "
      f"raw {raw_ilisi:.3f} scaled {scaled_ilisi:.3f}", flush=True)


def annotate_and_gate(res):
    tag = f"{res:g}"
    key = f"leiden_r{tag}"
    # n_iterations=2: leidenalg default -1 (run-to-convergence) is intractable here (8h+);
    # 2 passes is the standard practical Leiden setting.
    sc.tl.leiden(adata, resolution=res, flavor="leidenalg", n_iterations=2,
                 random_state=SEED, key_added=key)
    cl = adata.obs[key]
    n_clusters = cl.nunique()
    print(f"\n[res {tag}] leiden: {n_clusters} clusters", flush=True)

    clu_lin = LM.groupby(cl.values).mean()
    top_lin = clu_lin.idxmax(axis=1)
    clu_comp = top_lin.map(lambda l: COMPARTMENT.get(l, "other"))
    coarse_label = cl.map(top_lin).astype(str)
    compartment  = cl.map(clu_comp).astype(str)

    rows = []
    for c in clu_lin.index:
        m = (cl == c).values
        n = int(m.sum()); frac = n / tot
        dsv = pd.Series(ds[m]).value_counts(normalize=True)
        m01 = float(dsv.get("Mutter_01", 0.0)); m02 = float(dsv.get("Mutter_02", 0.0))
        negmean = float(neg[m].mean()); tl = top_lin[c]
        s2b = float(clu_lin.loc[c, tl] / max(negmean, 1e-6))
        rows.append(dict(cluster=c, n_cells=n, frac=round(frac, 4),
                         top_lineage=tl, compartment=clu_comp[c],
                         frac_M01=round(m01, 4), frac_M02=round(m02, 4),
                         minority_frac=round(min(m01, m02), 4),
                         neg_mean=round(negmean, 4), signal_to_bg=round(s2b, 2)))
    qc = pd.DataFrame(rows).sort_values("n_cells", ascending=False)
    qc.to_csv(f"{OUT}_r{tag}_cluster_qc.tsv", sep="\t", index=False)

    rec_rows = []
    for lin in sorted(set(top_lin.values)):
        clusters = top_lin[top_lin == lin].index
        m = cl.isin(clusters).values
        for dsname in ["Mutter_01", "Mutter_02"]:
            md = m & (ds == dsname)
            if md.sum() == 0:
                rec_rows.append(dict(lineage=lin, dataset=dsname, n_cells=0,
                                     marker_mean=np.nan, neg_mean=np.nan, recall_s2b=np.nan))
                continue
            mm = float(LM.loc[md, lin].mean()); nn = float(neg[md].mean())
            rec_rows.append(dict(lineage=lin, dataset=dsname, n_cells=int(md.sum()),
                                 marker_mean=round(mm, 4), neg_mean=round(nn, 4),
                                 recall_s2b=round(mm / max(nn, 1e-6), 2)))
    recall = pd.DataFrame(rec_rows)
    recall.to_csv(f"{OUT}_r{tag}_marker_recall.tsv", sep="\t", index=False)
    rwide = recall.pivot_table(index="lineage", columns="dataset", values="recall_s2b")

    major = qc[qc["frac"] > 0.01]
    n_major = int(len(major))
    immune_frac = float(compartment.eq("immune").mean())
    stroma_frac = float(compartment.eq("stroma").mean())
    lab_share = coarse_label.value_counts(normalize=True)
    max_label = str(lab_share.index[0]); max_label_share = float(lab_share.iloc[0])
    comps_present = sorted(set(clu_comp.values) & {"tumor", "immune", "stroma"})

    q41 = bool(n_major > 3 and immune_frac > 0.05 and stroma_frac > 0.10 and max_label_share < 0.85)
    minority_ok = bool((major["minority_frac"] >= 0.05).all())
    q42 = bool(scaled_ilisi >= 0.5 and minority_ok)
    q43 = bool(((rwide.fillna(0) > 1.0).all(axis=1)).all())

    gates = {
        "resolution": res, "n_clusters": n_clusters, "compartments_present": comps_present,
        "Q4.1_premise": {"pass": q41, "n_major_clusters": n_major,
                         "immune_frac": round(immune_frac, 4), "stroma_frac": round(stroma_frac, 4),
                         "max_label": max_label, "max_label_share": round(max_label_share, 4)},
        "Q4.2_batch_mixing": {"pass": q42, "scaled_iLISI": round(scaled_ilisi, 4),
                              "iLISI_ceiling": round(ceiling, 4),
                              "minority>=5%_all_major_clusters": minority_ok},
        "Q4.3_marker_recall": {"pass": q43, "lineages_failing_negbg": [
            str(l) for l in rwide.index[~(rwide.fillna(0) > 1.0).all(axis=1)]]},
        "Q4.4_negprobe": {"max_cluster_neg_mean": round(float(qc["neg_mean"].max()), 4),
                          "min_signal_to_bg": round(float(qc["signal_to_bg"].min()), 2)},
        "ALL_PASS": bool(q41 and q42 and q43),
    }
    json.dump(gates, open(f"{OUT}_r{tag}_gates.json", "w"), indent=2)

    out = pd.DataFrame({
        "cell": adata.obs_names, "leiden": cl.values,
        "coarse_label": coarse_label.values, "compartment": compartment.values,
        "dataset": ds, "slide_id": adata.obs["slide_id"].astype(str).values})
    out.to_parquet(f"{OUT}_r{tag}_coarse_labels.parquet", index=False)

    print(qc.to_string(index=False), flush=True)
    print(f"[res {tag}] Q4.1={'PASS' if q41 else 'FAIL'} "
          f"(major={n_major} immune={immune_frac:.3f} stroma={stroma_frac:.3f} "
          f"maxlabel={max_label_share:.3f}) Q4.2={'PASS' if q42 else 'FAIL'} "
          f"Q4.3={'PASS' if q43 else 'FAIL'} ALL={'PASS' if gates['ALL_PASS'] else 'FAIL'}", flush=True)
    return dict(resolution=res, n_clusters=n_clusters, n_major=n_major,
                immune_frac=round(immune_frac, 4), stroma_frac=round(stroma_frac, 4),
                tumor_frac=round(float(compartment.eq("tumor").mean()), 4),
                max_label=max_label, max_label_share=round(max_label_share, 4),
                scaled_iLISI=round(scaled_ilisi, 4),
                Q41=q41, Q42=q42, Q43=q43, ALL_PASS=gates["ALL_PASS"])


summary = pd.DataFrame([annotate_and_gate(r) for r in RES_LIST])
summary.to_csv(f"{OUT}_sweep_summary.tsv", sep="\t", index=False)
print("\n=== SWEEP SUMMARY ===\n" + summary.to_string(index=False), flush=True)
print(f"\nwrote: {OUT}_sweep_summary.tsv + per-res _r*_{{cluster_qc,marker_recall,gates,coarse_labels}}", flush=True)
