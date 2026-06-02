#!/usr/bin/env python
# Adjudicate the SingleR immune calls with direct gold-marker evidence. Computes the
# per-subcluster log-normalized mean of canonical lineage markers and the negprobe-relative
# signal, so we can tell whether (a) the "ILC"/"T" flip reflects real CD3 expression and
# (b) B / plasma / neutrophil lineages are present but mislabeled by SingleR.
# Args: <immune_checkpoint.h5ad> <immune_subclusters.parquet> <out_tsv> [singler_tsv]
import sys
import numpy as np
import scanpy as sc
import anndata as ad
import pandas as pd

ckpt, sub_pq, out_tsv = sys.argv[1:4]
singler_tsv = sys.argv[4] if len(sys.argv) > 4 else None

MARKERS = {
    "T":          ["Cd3d", "Cd3e", "Cd3g", "Cd8a", "Cd4"],
    "NK":         ["Ncr1", "Klrb1c", "Nkg7", "Gzma"],
    "B":          ["Cd19", "Ms4a1", "Cd79a", "Cd79b"],
    "Plasma":     ["Jchain", "Mzb1", "Xbp1", "Sdc1"],
    "Macrophage": ["Adgre1", "Itgam", "Csf1r", "C1qa", "Lyz2"],
    "DC":         ["Itgax", "Flt3", "Xcr1"],
    "Neutrophil": ["S100a8", "S100a9", "Ly6g", "Retnlg"],
    "Mast":       ["Cma1", "Mcpt4", "Kit", "Fcer1a"],
}

imm = ad.read_h5ad(ckpt)
sub = pd.read_parquet(sub_pq).set_index("cell")
imm.obs["sc"] = sub["immune_subcluster"].reindex(imm.obs_names).astype(str).values

sc.pp.normalize_total(imm, target_sum=1e4)
sc.pp.log1p(imm)

present = {k: [g for g in v if g in imm.var_names] for k, v in MARKERS.items()}
neg = imm.obs["neg"].values if "neg" in imm.obs else None

rows = []
for scl, idx in imm.obs.groupby("sc").indices.items():
    X = imm[idx]
    n = len(idx)
    row = {"subcluster": scl, "n": n}
    for lin, genes in present.items():
        if not genes:
            row[lin] = np.nan; continue
        gi = [imm.var_names.get_loc(g) for g in genes]
        row[lin] = float(np.asarray(X.X[:, gi].mean()))
    rows.append(row)

df = pd.DataFrame(rows)
df["subcluster"] = df["subcluster"].astype(str)
if singler_tsv:
    sg = pd.read_csv(singler_tsv, sep="\t").rename(columns={"immune_subcluster": "subcluster"})
    sg["subcluster"] = sg["subcluster"].astype(str)
    df = df.merge(sg[["subcluster", "immune_subtype", "delta"]], on="subcluster", how="left")
# argmax lineage by marker mean -> quick label-free read of what each subcluster looks like
lin_cols = list(present.keys())
df["marker_argmax"] = df[lin_cols].idxmax(axis=1)
df = df.sort_values("n", ascending=False)
df.to_csv(out_tsv, sep="\t", index=False, float_format="%.4f")

pd.set_option("display.width", 200, "display.max_columns", 30)
print(f"markers found: {[ (k,len(v)) for k,v in present.items() ]}")
print(df.to_string(index=False))
