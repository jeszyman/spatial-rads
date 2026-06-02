#!/usr/bin/env python3
"""Tier-2 stroma subtyping, step 1/2 (mirror of tier2_immune_subcluster.py).
Subsets the stroma compartment from the tier-1 res=3.0 coarse labels, rebuilds a neighbor
graph on the stroma cells' scVI latent only (the global graph is not a valid subgraph),
exports the raw-count matrix (genes x cells MTX) for the R/UCell annotation step, then sweeps
Leiden resolutions writing per-cell subcluster assignments + method-agnostic per-subcluster
QC (n_cells, dataset mix, negprobe). Lineage annotation is NOT done here: stromal markers
differ by orders of magnitude in raw counts (Fibroblast ECM genes Col1a1/Col1a2/Col3a1
dominate), so a raw-count argmax collapses every subcluster to Fibroblast. Lineage calls are
made downstream by tier2_stroma_ucell.R using rank-based UCell scoring, which is
magnitude-invariant. The five stromal lineages live in config/lineage_markers.yaml.
Usage: tier2_stroma_subcluster.py <full_checkpoint.h5ad> <coarse_labels.parquet> <ddir>
       <out_prefix> [res_csv=0.5,1.0,2.0]
"""
import sys, os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy.io import mmwrite

CKPT, LABELS_PQ, DDIR, OUT = sys.argv[1:5]
RES_LIST = [float(x) for x in (sys.argv[5].split(",") if len(sys.argv) > 5
                               else ["0.5", "1.0", "2.0"])]
SEED = 0
STR_CKPT = os.path.join(DDIR, "stroma_checkpoint.h5ad")
MTXDIR = os.path.join(DDIR, "stroma_mtx")
os.makedirs(MTXDIR, exist_ok=True)

# ---- restore stroma checkpoint, or build it (subset + neighbors on scVI latent) ----
if os.path.exists(STR_CKPT):
    print(f"restoring stroma checkpoint {STR_CKPT} ...", flush=True)
    st = ad.read_h5ad(STR_CKPT)
    print(f"stroma cells: {st.n_obs}", flush=True)
else:
    print("reading full checkpoint ...", flush=True)
    adata = ad.read_h5ad(CKPT)
    lab = pd.read_parquet(LABELS_PQ, columns=["cell", "compartment"]).set_index("cell")
    adata.obs["compartment"] = lab["compartment"].reindex(adata.obs_names).values
    st = adata[adata.obs["compartment"] == "stroma"].copy()
    del adata
    print(f"stroma cells: {st.n_obs}", flush=True)
    sc.pp.neighbors(st, use_rep="X_scvi", n_neighbors=15, random_state=SEED)
    st.write_h5ad(STR_CKPT)
    print(f"wrote {STR_CKPT}", flush=True)

tot = st.n_obs

# ---- export raw counts (genes x cells) for the R/UCell annotation step ----
mmwrite(os.path.join(MTXDIR, "counts.mtx"), st.X.T.tocsr())
pd.Series(st.var_names).to_csv(os.path.join(MTXDIR, "features.tsv"), index=False, header=False)
pd.Series(st.obs_names).to_csv(os.path.join(MTXDIR, "barcodes.tsv"), index=False, header=False)
print(f"wrote stroma_mtx/ ({st.n_vars} x {st.n_obs})", flush=True)

neg = st.obs["neg"].values.astype(float)
ds = st.obs["dataset"].astype(str).values


def subcluster(res):
    tag = f"{res:g}"
    key = f"stroma_leiden_r{tag}"
    # n_iterations=2: leidenalg run-to-convergence default is intractable at this scale.
    sc.tl.leiden(st, resolution=res, flavor="leidenalg", n_iterations=2,
                 random_state=SEED, key_added=key)
    cl = st.obs[key]
    n_sub = cl.nunique()
    print(f"\n[res {tag}] leiden: {n_sub} subclusters", flush=True)

    rows = []
    for c in sorted(cl.cat.categories, key=int):
        m = (cl == c).values
        n = int(m.sum())
        dsv = pd.Series(ds[m]).value_counts(normalize=True)
        m01 = float(dsv.get("Mutter_01", 0.0)); m02 = float(dsv.get("Mutter_02", 0.0))
        rows.append(dict(subcluster=c, n_cells=n, frac=round(n / tot, 4),
                         frac_M01=round(m01, 4), frac_M02=round(m02, 4),
                         minority_frac=round(min(m01, m02), 4),
                         neg_mean=round(float(neg[m].mean()), 4)))
    qc = pd.DataFrame(rows).sort_values("n_cells", ascending=False)
    qc.to_csv(f"{OUT}_r{tag}_subcluster_qc.tsv", sep="\t", index=False)

    out = pd.DataFrame({
        "cell": st.obs_names, "stroma_subcluster": cl.values,
        "dataset": ds, "slide_id": st.obs["slide_id"].astype(str).values})
    out.to_parquet(os.path.join(DDIR, f"stroma_r{tag}_subclusters.parquet"), index=False)
    print(qc.to_string(index=False), flush=True)
    return dict(resolution=res, n_subclusters=n_sub)


summary = pd.DataFrame([subcluster(r) for r in RES_LIST])
summary.to_csv(f"{OUT}_sweep_summary.tsv", sep="\t", index=False)
print("\n=== STROMA SWEEP SUMMARY ===\n" + summary.to_string(index=False), flush=True)
print(f"\nwrote {OUT}_sweep_summary.tsv + per-res _r*_subcluster_qc.tsv "
      f"+ {DDIR}/stroma_r*_subclusters.parquet + stroma_mtx/", flush=True)
