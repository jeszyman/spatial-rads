#!/usr/bin/env python3
"""Top pathways x the three contrasts -- one clean heatmap. GSEA NES averaged across cell
types (one value per pathway per contrast), top sets by |NES|. Args: <results_master.tsv> <out.png>"""
import sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

M, OUT = sys.argv[1], sys.argv[2]
CONTRASTS = ["MBRT_vs_Ctrl", "SBRT_vs_Ctrl", "MBRT_vs_SBRT"]
g = pd.read_csv(M, sep="\t", low_memory=False)
g = g[g.readout_class == "gsea"].dropna(subset=["effect"]).copy()
g["set"] = g.feature.str.replace("HALLMARK_", "", regex=False)
mat = (g.pivot_table(index="set", columns="contrast", values="effect", aggfunc="mean")
        .reindex(columns=CONTRASTS))
mat = mat.loc[mat.abs().max(axis=1).sort_values(ascending=False).head(15).index]
mat = mat.reindex(mat["SBRT_vs_Ctrl"].fillna(0).sort_values().index)
A = mat.values.astype(float); vmax = np.nanmax(np.abs(A))
fig, ax = plt.subplots(figsize=(6.4, 0.42 * len(mat) + 1.6))
im = ax.imshow(A, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
ax.set_xticks(range(3)); ax.set_xticklabels(CONTRASTS, rotation=20, ha="right", fontsize=9)
ax.set_yticks(range(len(mat)))
ax.set_yticklabels([s.replace("_", " ").title() for s in mat.index], fontsize=8.5)
for i in range(A.shape[0]):
    for j in range(3):
        v = A[i, j]
        if np.isfinite(v):
            ax.text(j, i, f"{v:+.1f}", ha="center", va="center", fontsize=7.5,
                    color="white" if abs(v) > 0.6 * vmax else "black")
ax.set_title("Top pathways across the three contrasts (mean GSEA NES)", fontsize=11)
fig.colorbar(im, ax=ax, fraction=0.045, pad=0.02).set_label("NES (red=up, blue=down)", fontsize=8)
fig.tight_layout(); fig.savefig(OUT, dpi=175); print("wrote", OUT, "|", len(mat), "pathways")
