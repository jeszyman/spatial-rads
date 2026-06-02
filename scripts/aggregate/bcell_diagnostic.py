#!/usr/bin/env python
# What happened to B cells? Look cohort-wide (not just the immune compartment) for any
# pocket of B-marker-high cells. Three explanations to distinguish: (a) genuinely rare in
# this 4T1 model, (b) B markers absent/dead on the panel, (c) present but mis-compartmented
# at tier-1. Reports per-panel B marker prevalence, cohort-wide expressing fraction, and the
# top coarse clusters by B-marker signal with their compartment.
# Args: <full_checkpoint.h5ad> <full_coarse_labels.parquet>
import sys
import numpy as np
import scanpy as sc
import anndata as ad
import pandas as pd

ckpt, labels_pq = sys.argv[1], sys.argv[2]
B = ["Cd19", "Ms4a1", "Cd79a", "Cd79b", "Cd22", "Ighm", "Iglc1"]

a = ad.read_h5ad(ckpt)
lab = pd.read_parquet(labels_pq, columns=["cell", "leiden", "coarse_label", "compartment"]).set_index("cell")
for c in ["leiden", "coarse_label", "compartment"]:
    a.obs[c] = lab[c].reindex(a.obs_names).values

present = [g for g in B if g in a.var_names]
absent  = [g for g in B if g not in a.var_names]
print(f"B markers ON panel:  {present}")
print(f"B markers OFF panel: {absent}")

# raw-count prevalence (fraction of all cells with >=1 count) before normalization
raw = a[:, present].X
expr_frac = np.asarray((raw > 0).mean(axis=0)).ravel()
print("\n--- cohort-wide fraction of cells expressing (raw count >=1) ---")
for g, f in zip(present, expr_frac):
    print(f"  {g:8s} {f*100:6.2f}%")

sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
gi = [a.var_names.get_loc(g) for g in present]
a.obs["B_score"] = np.asarray(a[:, gi].X.mean(axis=1)).ravel()

print("\n--- B score by compartment (mean lognorm of B markers) ---")
print(a.obs.groupby("compartment")["B_score"].agg(["mean", "size"]).to_string())

clu = a.obs.groupby("leiden").agg(
    n=("B_score", "size"), B_score=("B_score", "mean"),
    compartment=("compartment", lambda s: s.iloc[0]),
    coarse=("coarse_label", lambda s: s.iloc[0])).reset_index()
clu = clu.sort_values("B_score", ascending=False)
print("\n--- top 12 coarse clusters by B score (any compartment) ---")
print(clu.head(12).to_string(index=False))
print(f"\nmax cluster B score = {clu.B_score.max():.3f}; "
      f"for reference a real lineage cluster scores ~1-1.8 in its own markers.")
