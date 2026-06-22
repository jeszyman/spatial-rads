#!/usr/bin/env python3
"""Integration-QC + atlas visualization for the M01/M02 scVI integration (task #8,
plan-aggregate-refactor.md). Three figures from a dataset-stratified ~200k subsample of
the integrated cohort:

  1. integration_umap.png  -- paired UMAP: UNINTEGRATED (PCA on log-norm counts) vs
     scVI-INTEGRATED latent, each colored BY DATASET (M01/M02: split -> blended = batch
     removed) and BY COMPARTMENT (stays structured = biology preserved).
  2. integration_cluster_composition.png -- per-Leiden-cluster M01/M02 fraction
     (computed on ALL cells, not the subsample = the honest quantitative view); the
     Q4.2 minority-presence gate drawn, scaled iLISI annotated.
  3. umap_atlas.png -- the scVI-integrated UMAP colored by final cell_subtype (the atlas).

Reads cluster_checkpoint.h5ad (already holds raw counts in X + the scVI latent in
obsm['X_scvi']), so no MTX re-read. Labels join from full_labels.parquet.
Usage: full_umap.py <cluster_checkpoint.h5ad> <full_labels.parquet> <full_gates.json>
       <out_integration.png> <out_composition.png> <out_atlas.png> <out_coords.parquet>
"""
import sys, json
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

CKPT, LABELS, GATES, OUT_INT, OUT_COMP, OUT_ATLAS, OUT_COORDS = sys.argv[1:8]
SEED, N_SUB = 0, 200_000
np.random.seed(SEED)

# ---- load checkpoint (X=counts, obsm['X_scvi']=latent) + join final labels ----
print("reading checkpoint ...", flush=True)
adata = ad.read_h5ad(CKPT)
lab = pd.read_parquet(LABELS).set_index("cell")
for c in ["dataset", "leiden", "compartment", "cell_subtype"]:
    adata.obs[c] = lab[c].reindex(adata.obs_names).values
print(f"cohort {adata.n_obs} cells | datasets {dict(pd.Series(adata.obs['dataset']).value_counts())}", flush=True)

# ---- dataset-stratified subsample to ~200k (proportional = honest M01/M02 ratio) ----
obs = adata.obs.reset_index(names="cell")
frac = min(1.0, N_SUB / len(obs))
idx = (obs.groupby("dataset", group_keys=False)
          .apply(lambda g: g.sample(frac=frac, random_state=SEED)).index.values)
sub = adata[np.sort(idx)].copy()
print(f"subsample {sub.n_obs} | {dict(pd.Series(sub.obs['dataset']).value_counts())}", flush=True)

# ---- integrated UMAP from the scVI latent ----
print("integrated UMAP (scVI latent) ...", flush=True)
sc.pp.neighbors(sub, use_rep="X_scvi", n_neighbors=15, random_state=SEED)
sc.tl.umap(sub, random_state=SEED)
umap_int = sub.obsm["X_umap"].copy()

# ---- unintegrated UMAP from PCA on log-norm counts (the 'before' panel) ----
print("unintegrated UMAP (PCA on counts) ...", flush=True)
un = sub.copy()
sc.pp.normalize_total(un, target_sum=1e4); sc.pp.log1p(un); sc.pp.scale(un, max_value=10)
sc.tl.pca(un, n_comps=30, random_state=SEED)
sc.pp.neighbors(un, use_rep="X_pca", n_neighbors=15, random_state=SEED)
sc.tl.umap(un, random_state=SEED)
umap_un = un.obsm["X_umap"].copy()

ds   = sub.obs["dataset"].astype(str).values
comp = sub.obs["compartment"].astype(str).values
subt = sub.obs["cell_subtype"].astype(str).values

def scatter(ax, xy, key, palette, title, legend=False, s=2):
    # per-point color in a SHUFFLED z-order so no category systematically overplots
    # another (else the 71% M02 hides M01 and the integrated panel looks all-orange).
    cols  = np.array([palette.get(k, "#dddddd") for k in key])
    order = np.random.permutation(len(key))
    ax.scatter(xy[order, 0], xy[order, 1], s=s, c=cols[order], alpha=0.45,
               linewidths=0, rasterized=True)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title, fontsize=11)
    for sp in ax.spines.values(): sp.set_visible(False)
    if legend:
        handles = [Line2D([0], [0], marker="o", linestyle="", markersize=6,
                          markerfacecolor=c, markeredgewidth=0, label=k)
                   for k, c in palette.items()]
        ax.legend(handles=handles, fontsize=7, frameon=False,
                  loc="center left", bbox_to_anchor=(1.0, 0.5))

DS_PAL   = {"Mutter_01": "#1f77b4", "Mutter_02": "#ff7f0e"}
COMP_PAL = {"tumor": "#d62728", "stroma": "#2ca02c", "immune": "#1f77b4", "unassigned": "#cccccc"}

# ===== Figure 1: paired UMAP (unintegrated vs scVI) x (by dataset / by compartment) =====
fig, axes = plt.subplots(2, 2, figsize=(11, 10))
scatter(axes[0, 0], umap_un,  ds,   DS_PAL,   "Unintegrated (PCA)  -  by dataset")
scatter(axes[0, 1], umap_un,  comp, COMP_PAL, "Unintegrated (PCA)  -  by compartment")
scatter(axes[1, 0], umap_int, ds,   DS_PAL,   "scVI-integrated  -  by dataset", legend=True)
scatter(axes[1, 1], umap_int, comp, COMP_PAL, "scVI-integrated  -  by compartment", legend=True)
fig.suptitle("M01/M02 integration: datasets blend, compartments stay separated", fontsize=13)
fig.tight_layout(rect=[0, 0, 0.88, 0.97])
fig.savefig(OUT_INT, dpi=200); plt.close(fig)
print(f"wrote {OUT_INT}", flush=True)

# ===== Figure 2: per-cluster M01/M02 composition (ALL cells) + iLISI =====
ilisi = json.load(open(GATES)).get("Q4.2_batch_mixing", {}).get("scaled_iLISI", float("nan"))
comp_tab = (lab.reset_index().groupby(["leiden", "dataset"]).size().unstack(fill_value=0))
comp_tab = comp_tab.reindex(columns=["Mutter_01", "Mutter_02"], fill_value=0)
tot = comp_tab.sum(axis=1)
comp_tab = comp_tab.loc[tot.sort_values(ascending=False).index]      # largest clusters first
n_tot = tot.sum(); m01_overall = comp_tab["Mutter_01"].sum() / n_tot
fr = comp_tab.div(comp_tab.sum(axis=1), axis=0)
fig, ax = plt.subplots(figsize=(14, 4.5))
x = np.arange(len(fr))
ax.bar(x, fr["Mutter_01"], width=1.0, color=DS_PAL["Mutter_01"], label="Mutter_01")
ax.bar(x, fr["Mutter_02"], width=1.0, bottom=fr["Mutter_01"], color=DS_PAL["Mutter_02"], label="Mutter_02")
ax.axhline(m01_overall, color="k", lw=1, ls="--", label=f"cohort M01 frac = {m01_overall:.2f}")
ax.set_xlim(-0.5, len(fr) - 0.5); ax.set_ylim(0, 1)
ax.set_xlabel(f"Leiden cluster (res=3.0), {len(fr)} clusters, largest first")
ax.set_ylabel("dataset fraction")
ax.set_title(f"Per-cluster M01/M02 composition  (scaled iLISI = {ilisi:.3f}, gate >=0.5)")
ax.legend(fontsize=8, frameon=False, ncol=3, loc="upper right")
fig.tight_layout(); fig.savefig(OUT_COMP, dpi=200); plt.close(fig)
print(f"wrote {OUT_COMP}", flush=True)

# ===== Figure 3: integrated atlas UMAP by cell_subtype =====
order = pd.Series(subt).value_counts().index.tolist()
cmap = plt.get_cmap("tab20")
SUB_PAL = {k: cmap(i % 20) for i, k in enumerate(order)}
fig, ax = plt.subplots(figsize=(9, 7.5))
scatter(ax, umap_int, subt, SUB_PAL, "scVI-integrated atlas  -  by cell subtype", legend=True, s=2)
fig.tight_layout(rect=[0, 0, 0.78, 1]); fig.savefig(OUT_ATLAS, dpi=200); plt.close(fig)
print(f"wrote {OUT_ATLAS}", flush=True)

# ---- coords for reuse ----
pd.DataFrame({
    "cell": sub.obs_names, "dataset": ds, "compartment": comp, "cell_subtype": subt,
    "leiden": sub.obs["leiden"].astype(str).values,
    "umap_int_1": umap_int[:, 0], "umap_int_2": umap_int[:, 1],
    "umap_unint_1": umap_un[:, 0], "umap_unint_2": umap_un[:, 1],
}).to_parquet(OUT_COORDS, index=False)
print(f"wrote {OUT_COORDS} ({sub.n_obs} cells)", flush=True)
