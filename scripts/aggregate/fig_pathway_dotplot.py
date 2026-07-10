#!/usr/bin/env python3
"""Top pathways as a single-cell DOT PLOT, faceted by contrast (M02 day-2, n=4).
y = pathway (GSEA set), x = cell type. COLOR + SIZE = NES (direction + magnitude).
BLACK RING = FDR-significant (padj < 0.05); no ring = trend only.
Args: <results_master.tsv> <out.png>"""
import sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

M, OUT = sys.argv[1], sys.argv[2]
CONTRASTS = ["MBRT_vs_Ctrl", "SBRT_vs_Ctrl", "MBRT_vs_SBRT"]
CT = ["Tumor","Fibroblast","SmoothMuscle","Adipocyte","Endothelial",
      "Macrophages","DC","T cells","NK cells","ILC","Plasma cells","Mast cells","Neutrophils"]
N_TOP = 7
d = pd.read_csv(M, sep="\t", low_memory=False)
g = d[(d.readout_class == "gsea") & (d.unit.isin(CT))].dropna(subset=["effect"]).copy()
g["set"] = g.feature.str.replace("HALLMARK_", "", regex=False).str.replace("_", " ").str.title()
g["sig"] = (g.padj_confirmatory.fillna(1) < 0.05) | (g.padj_exploratory.fillna(1) < 0.05)
sets = []
for c in CONTRASTS:
    s = g[g.contrast == c]
    sets += list(s.reindex(s.effect.abs().sort_values(ascending=False).index).head(N_TOP).set)
sets = list(dict.fromkeys(sets))
g = g[g.set.isin(sets)]
ord_s = (g[g.contrast == "SBRT_vs_Ctrl"].groupby("set").effect
         .apply(lambda x: x.loc[x.abs().idxmax()]).reindex(sets).fillna(0).sort_values().index.tolist())
si = {s: i for i, s in enumerate(ord_s)}; ci = {c: i for i, c in enumerate(CT)}
vmax = float(np.nanpercentile(g.effect.abs(), 98))
smap = lambda e: 18 + 64 * np.minimum(np.abs(e), 2.5) / 2.5
fig, axes = plt.subplots(1, 3, figsize=(14.5, 0.34 * len(ord_s) + 2.4), sharey=True)
sc = None
for ax, c in zip(axes, CONTRASTS):
    s = g[g.contrast == c]
    ec = np.where(s.sig.values, "black", "#cfcfcf"); lw = np.where(s.sig.values, 1.5, 0.3)
    sc = ax.scatter(s.unit.map(ci), s.set.map(si), s=smap(s.effect.values),
                    c=s.effect.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                    edgecolors=ec, linewidths=lw)
    ax.set_xticks(range(len(CT))); ax.set_xticklabels(CT, rotation=90, fontsize=7)
    ax.set_title(c, fontsize=10.5); ax.set_xlim(-0.5, len(CT) - 0.5)
    ax.grid(alpha=0.15, linewidth=0.4)
axes[0].set_yticks(range(len(ord_s))); axes[0].set_yticklabels(ord_s, fontsize=8)
cb = fig.colorbar(sc, ax=axes, fraction=0.014, pad=0.01); cb.set_label("NES", fontsize=8)
sig_leg = [
    Line2D([0],[0], marker="o", color="w", markerfacecolor="#dddddd", markeredgecolor="black",
           markeredgewidth=1.6, markersize=12, label="FDR-significant  (padj < 0.05)"),
    Line2D([0],[0], marker="o", color="w", markerfacecolor="#dddddd", markeredgecolor="#cfcfcf",
           markeredgewidth=0.6, markersize=12, label="not significant  (trend only)"),
]
axes[1].legend(handles=sig_leg, loc="upper center", bbox_to_anchor=(0.5, -0.20),
               ncol=2, fontsize=9, frameon=False, handletextpad=0.4)
fig.suptitle("Top pathways (GSEA NES) by cell type, faceted by contrast  --  Mutter_02 day 2 (48h, n=4)\n"
             "color+size = NES; black ring = FDR-significant", fontsize=11.5)
fig.savefig(OUT, dpi=160, bbox_inches="tight"); print("wrote", OUT, "|", len(ord_s), "pathways")
