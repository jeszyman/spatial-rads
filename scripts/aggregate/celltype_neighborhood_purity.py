#!/usr/bin/env python3
# Transcriptional neighborhood purity (Plummer et al., Nat Biotechnol 2025, Fig 4d): per cell,
# the fraction of its k nearest neighbors in the scVI integrated latent that share its final
# cell_subtype. A typing-QC metric complementing marker recall -- it measures whether the locked
# cross-dataset labels are transcriptionally coherent (do like-typed cells sit together in the
# integrated space). k=30 matches the paper's neighborhood size; a k=15 pass (the tier-1 Leiden
# graph's own neighborhood size) checks the per-subtype ordering is not a k artifact. Reads the
# existing latent + locked labels; recomputes neighbors on CPU (no scVI, no GPU). Emits per
# compartment x cell_subtype box quantiles for the figure. Plotting lives in
# scripts/fig_qc_neighborhood_purity.R.
# Args: <scvi_latent.parquet> <full_labels.parquet> <out.tsv>
import sys
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from scipy.stats import spearmanr

LATENT, LABELS, OUT_TSV = sys.argv[1], sys.argv[2], sys.argv[3]
COMPARTMENTS = ["tumor", "stroma", "immune"]
QS = [0.05, 0.25, 0.50, 0.75, 0.95]


def per_cell_purity(adata, labels, k):
    """Fraction of each cell's k-1 graph neighbors sharing its label."""
    sc.pp.neighbors(adata, n_neighbors=k, use_rep="X_scvi", random_state=0)
    g = adata.obsp["distances"]  # CSR, k-1 explicit neighbors per row
    indptr, indices = g.indptr, g.indices
    counts = np.diff(indptr)
    assert counts.min() == counts.max(), (
        f"kNN graph not regular ({counts.min()}..{counts.max()} nbrs/row)")
    codes = pd.Categorical(labels).codes
    self_code = codes
    nbr_codes = codes[indices]
    same = (nbr_codes == np.repeat(self_code, counts))
    return same.reshape(len(labels), counts[0]).mean(axis=1)


lat = pd.read_parquet(LATENT)
scvi_cols = sorted([c for c in lat.columns if c.startswith("scVI_")],
                   key=lambda c: int(c.split("_")[1]))
lab = pd.read_parquet(LABELS)[["cell", "compartment", "cell_subtype"]]

d = lat.merge(lab, on="cell", how="inner")
d = d[d["compartment"].isin(COMPARTMENTS)].reset_index(drop=True)
assert len(d) > 0

ad = AnnData(X=np.zeros((len(d), 1), dtype=np.float32))
ad.obsm["X_scvi"] = d[scvi_cols].to_numpy(dtype=np.float32)
subtype = d["cell_subtype"].to_numpy()

pur30 = per_cell_purity(ad, subtype, 30)
pur15 = per_cell_purity(ad, subtype, 15)

d = d.assign(purity=pur30)


def box(x):
    q = np.quantile(x["purity"], QS)
    return pd.Series({"n": len(x), "ymin": q[0], "lower": q[1], "middle": q[2],
                      "upper": q[3], "ymax": q[4], "mean": x["purity"].mean()})


summ = (d.groupby(["compartment", "cell_subtype"], observed=True)
          .apply(box, include_groups=False).reset_index())

# k-robustness: Spearman of per-subtype medians, k=15 vs k=30 (constant column for the legend).
med30 = pd.Series(pur30).groupby(subtype).median()
med15 = pd.Series(pur15).groupby(subtype).median()
rho = spearmanr(med30.loc[med15.index], med15).correlation

summ["k_main"] = 30
summ["k_robust_spearman"] = round(float(rho), 4)
summ = summ.sort_values(["compartment", "middle"], ascending=[True, False])
summ.to_csv(OUT_TSV, sep="\t", index=False)

print(f"cells scored: {len(d):,}; subtypes: {summ['cell_subtype'].nunique()}")
print(f"k=15 vs k=30 per-subtype median Spearman: {rho:.4f}")
print(summ[["compartment", "cell_subtype", "n", "middle", "mean"]]
      .to_string(index=False))
print(f"wrote {OUT_TSV}")
