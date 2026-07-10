#!/usr/bin/env python3
"""Top differential genes as a single-cell DOT PLOT, faceted by contrast.
y = gene, x = cell type. dot COLOR = log2FC (direction+magnitude), SIZE = |log2FC|.
RING encodes the detection-vs-level call_class (the robustness-pass headline):
  black ring  = regulation     (per-expresser level moved; genuine transcriptional change)
  gold ring   = fraction_shift (more/fewer cells express; a cell-state composition effect)
  no ring     = ambiguous      (significant but per-expresser level too sparse to separate)
Args: <results_master.tsv> <out.png>"""
import sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

M, OUT = sys.argv[1], sys.argv[2]
CONTRASTS = ["MBRT_vs_Ctrl", "SBRT_vs_Ctrl", "MBRT_vs_SBRT"]
CT = ["Tumor", "Fibroblast", "SmoothMuscle", "Adipocyte", "Endothelial",
      "Macrophages", "DC", "T cells", "NK cells", "ILC", "Plasma cells", "Mast cells", "Neutrophils"]
N_TOP = 8
REG, FRAC = "#000000", "#e8a13a"          # ring colours: regulation / fraction_shift
d = pd.read_csv(M, sep="\t", low_memory=False)
de = d[(d.readout_class == "DE") & (d.unit.isin(CT))].dropna(subset=["effect"]).copy()
if "call_class" not in de.columns:
    de["call_class"] = np.nan
de["cc"] = de.call_class.fillna("ambiguous")
genes = []
for c in CONTRASTS:
    s = de[de.contrast == c]
    genes += list(s.reindex(s.effect.abs().sort_values(ascending=False).index).head(N_TOP).feature)
genes = list(dict.fromkeys(genes))
de = de[de.feature.isin(genes)]
ord_g = (de[de.contrast == "SBRT_vs_Ctrl"].groupby("feature").effect
         .apply(lambda x: x.loc[x.abs().idxmax()]).reindex(genes).fillna(0).sort_values().index.tolist())
gi = {g: i for i, g in enumerate(ord_g)}; ci = {c: i for i, c in enumerate(CT)}
vmax = float(np.nanpercentile(de.effect.abs(), 98))
smap = lambda e: 18 + 64 * np.minimum(np.abs(e), 1.5) / 1.5

def ring(cc):
    ec = np.select([cc == "regulation", cc == "fraction_shift"], [REG, FRAC], default="#cfcfcf")
    lw = np.select([cc == "regulation", cc == "fraction_shift"], [1.7, 1.7], default=0.3)
    return ec, lw

fig, axes = plt.subplots(1, 3, figsize=(14, 0.33 * len(ord_g) + 2.4), sharey=True)
sc = None
for ax, c in zip(axes, CONTRASTS):
    s = de[de.contrast == c]
    ec, lw = ring(s.cc.values)
    sc = ax.scatter(s.unit.map(ci), s.feature.map(gi), s=smap(s.effect.values),
                    c=s.effect.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                    edgecolors=ec, linewidths=lw)
    ax.set_xticks(range(len(CT))); ax.set_xticklabels(CT, rotation=90, fontsize=7)
    ax.set_title(c, fontsize=10.5); ax.set_xlim(-0.5, len(CT) - 0.5)
    ax.grid(alpha=0.15, linewidth=0.4)
axes[0].set_yticks(range(len(ord_g))); axes[0].set_yticklabels(ord_g, fontsize=8)
cb = fig.colorbar(sc, ax=axes, fraction=0.014, pad=0.01); cb.set_label("log2FC", fontsize=8)
leg = [
    Line2D([0], [0], marker="o", color="w", markerfacecolor="#dddddd", markeredgecolor=REG,
           markeredgewidth=1.8, markersize=12, label="regulation  (per-cell level moved)"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor="#dddddd", markeredgecolor=FRAC,
           markeredgewidth=1.8, markersize=12, label="fraction_shift  (more/fewer cells express)"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor="#dddddd", markeredgecolor="#cfcfcf",
           markeredgewidth=0.6, markersize=12, label="ambiguous  (too sparse to separate)"),
]
axes[1].legend(handles=leg, loc="upper center", bbox_to_anchor=(0.5, -0.20),
               ncol=3, fontsize=8.5, frameon=False, handletextpad=0.4)
fig.suptitle("Top differential genes by cell type, faceted by contrast  "
             "(color+size = log2FC; ring = detection-vs-level class)", fontsize=12)
fig.savefig(OUT, dpi=160, bbox_inches="tight"); print("wrote", OUT, "|", len(ord_g), "genes")
