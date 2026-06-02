#!/usr/bin/env python
# Tier-2 immune subtyping, step 1/2. Subsets the immune compartment from the tier-1
# res=3.0 coarse labels, rebuilds a neighbor graph on the immune cells' scVI latent only
# (the global graph is not a valid subgraph), and Leiden-subclusters them at a conservative
# resolution. Writes the immune raw-counts subset (genes x cells MTX) + barcodes/features
# and a per-cell subcluster parquet for the cluster-level SingleR/ImmGen step.
# Args: <checkpoint.h5ad> <full_coarse_labels.parquet> <outdir> [resolution=1.0]
import sys, os
import numpy as np
import scanpy as sc
import anndata as ad
import pandas as pd
from scipy.io import mmwrite

ckpt, labels_pq, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
res = float(sys.argv[4]) if len(sys.argv) > 4 else 1.0
os.makedirs(os.path.join(outdir, "immune_mtx"), exist_ok=True)
IMM_CKPT = os.path.join(outdir, "immune_checkpoint.h5ad")

if os.path.exists(IMM_CKPT):
    print(f"restoring immune checkpoint {IMM_CKPT} ...", flush=True)
    imm = ad.read_h5ad(IMM_CKPT)
    print(f"immune cells: {imm.n_obs}", flush=True)
else:
    print("reading full checkpoint ...", flush=True)
    adata = ad.read_h5ad(ckpt)
    lab = pd.read_parquet(labels_pq, columns=["cell", "compartment"]).set_index("cell")
    adata.obs["compartment"] = lab["compartment"].reindex(adata.obs_names).values
    imm = adata[adata.obs["compartment"] == "immune"].copy()
    del adata
    print(f"immune cells: {imm.n_obs}", flush=True)
    sc.pp.neighbors(imm, use_rep="X_scvi", n_neighbors=15)
    imm.write_h5ad(IMM_CKPT)
    print(f"wrote {IMM_CKPT}", flush=True)
# n_iterations=2: leidenalg default -1 (run-to-convergence) is intractable at scale.
sc.tl.leiden(imm, resolution=res, flavor="leidenalg", n_iterations=2,
             random_state=0, key_added="immune_subcluster")
nsub = imm.obs["immune_subcluster"].nunique()
print(f"immune subclusters @ res={res}: {nsub}", flush=True)
print(imm.obs["immune_subcluster"].value_counts().sort_index().to_string(), flush=True)

# raw counts subset, genes x cells, for SingleR
mmwrite(os.path.join(outdir, "immune_mtx", "counts.mtx"), imm.X.T.tocsr())
pd.Series(imm.var_names).to_csv(os.path.join(outdir, "immune_mtx", "features.tsv"),
                               index=False, header=False)
pd.Series(imm.obs_names).to_csv(os.path.join(outdir, "immune_mtx", "barcodes.tsv"),
                               index=False, header=False)

out = pd.DataFrame({
    "cell": imm.obs_names,
    "immune_subcluster": imm.obs["immune_subcluster"].astype(str).values,
    "dataset": imm.obs["dataset"].values,
    "slide_id": imm.obs["slide_id"].values,
})
out.to_parquet(os.path.join(outdir, "immune_subclusters.parquet"), index=False)
print(f"wrote immune_mtx/ ({imm.n_vars} x {imm.n_obs}) + immune_subclusters.parquet", flush=True)
