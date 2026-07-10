#!/usr/bin/env python3
"""Top differential genes x the three contrasts -- log2FC heatmap so you can read which
genes move, which way, in which comparison. Rows = the union of the top-N |log2FC| genes
per contrast (gene @ cell_type); columns = MBRT_vs_Ctrl / SBRT_vs_Ctrl / MBRT_vs_SBRT.
Args: <results_master.tsv> <out.png>"""
import sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

M, OUT = sys.argv[1], sys.argv[2]
N_TOP = 14
CONTRASTS = ["MBRT_vs_Ctrl", "SBRT_vs_Ctrl", "MBRT_vs_SBRT"]
d = pd.read_csv(M, sep="\t", low_memory=False)
de = d[d.readout_class == "DE"].dropna(subset=["effect"]).copy()
de["row"] = de.feature.astype(str) + " @ " + de.unit.astype(str)
# union of top-N rows per contrast by |effect|
keep = set()
for c in CONTRASTS:
    s = de[de.contrast == c]
    keep |= set(s.reindex(s.effect.abs().sort_values(ascending=False).index).head(N_TOP).row)
mat = (de[de.row.isin(keep)]
       .pivot_table(index="row", columns="contrast", values="effect", aggfunc="first")
       .reindex(columns=CONTRASTS))
mat = mat.reindex(mat["SBRT_vs_Ctrl"].fillna(0).sort_values().index)   # order by SBRT effect
A = mat.values.astype(float)
vmax = np.nanmax(np.abs(A))
fig, ax = plt.subplots(figsize=(6.2, max(6, 0.26 * len(mat))))
im = ax.imshow(A, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
ax.set_xticks(range(len(CONTRASTS))); ax.set_xticklabels(CONTRASTS, rotation=20, ha="right", fontsize=9)
ax.set_yticks(range(len(mat))); ax.set_yticklabels(mat.index, fontsize=7)
for i in range(A.shape[0]):
    for j in range(A.shape[1]):
        v = A[i, j]
        if np.isfinite(v):
            ax.text(j, i, f"{v:+.1f}", ha="center", va="center", fontsize=6,
                    color="white" if abs(v) > 0.6 * vmax else "black")
ax.set_title("Top differential genes (log2FC) across the three contrasts", fontsize=11)
cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02); cb.set_label("log2FC", fontsize=8)
fig.tight_layout(); fig.savefig(OUT, dpi=175); print("wrote", OUT, "|", len(mat), "genes")
