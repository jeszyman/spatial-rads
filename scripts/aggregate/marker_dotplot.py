#!/usr/bin/env python3
"""Marker-specificity dot-plot: does each label's cells express that cell type's
canonical markers? Left = ORIGINAL transfer labels, right = scVI integration labels,
SAME marker columns. A clean diagonal (each label lights up only its own markers) =
the better labeling. No interpretive mapping -- just label vs marker biology.

Reads cluster_checkpoint.h5ad (raw counts in X) + obs.parquet (old cell_type) +
full_labels.parquet (new cell_subtype). CPM-log normalizes on the fly over marker cols.
"""
import numpy as np, pandas as pd, scipy.sparse as sp
import anndata as ad
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

plt.rcParams.update({"font.family": "Arial", "font.size": 10, "pdf.fonttype": 42,
                     "ps.fonttype": 42, "svg.fonttype": "none"})
AGG = "/mnt/data/projects/spatial-rads/aggregate"; OUT = "results/aggregate/plots"

# ---- marker panel: (gene, lineage block) in column order ----
MARKERS = [
    ("Epcam","Tumor"),("Krt8","Tumor"),
    ("Col1a1","Fibroblast"),("Pdgfra","Fibroblast"),
    ("Pecam1","Endothelial"),("Cdh5","Endothelial"),
    ("Acta2","SmoothMuscle"),("Rgs5","SmoothMuscle"),
    ("Cidea","Adipocyte"),("Adipoq","Adipocyte"),
    ("Ptprc","Immune (pan)"),
    ("Cd3e","T / NK"),("Nkg7","T / NK"),
    ("Cd79a","B / Plasma"),("Mzb1","B / Plasma"),
    ("C1qa","Myeloid"),("Cd68","Myeloid"),("Itgax","Myeloid"),
    ("S100a8","Neutrophil"),("Elane","Neutrophil"),
    ("Cpa3","Mast"),
]
GENES = [g for g, _ in MARKERS]

NEW_ORDER = ["Tumor","Epithelial cells","Fibroblast","Endothelial","SmoothMuscle",
             "Adipocyte","T cells","NK cells","ILC","Plasma cells","Macrophages","DC",
             "Neutrophils","Mast cells","unassigned"]
OLD_ORDER = ["tumor_epithelial","a","Stem.Prog","Fibroblastic.reticular","Pericyte",
             "Blood.endothelial","CD8.T.cell","gdT","NKT","Ly6Clo.blood.monocytes",
             "spleen.red.pulp.macs","Dendritic"]

# lineage colour + short code for the column group bar; row labels are tinted to match
# so a cell label (row) points at the marker block (columns) it should light up.
LIN = {"Tumor":("#8c564b","Tumor"), "Fibroblast":("#2ca02c","Fibro"),
       "Endothelial":("#17becf","Endo"), "SmoothMuscle":("#bcbd22","SMC"),
       "Adipocyte":("#e377c2","Adipo"), "Immune (pan)":("#7f7f7f","Imm"),
       "T / NK":("#1f77b4","T/NK"), "B / Plasma":("#9467bd","B/Pl"),
       "Myeloid":("#ff7f0e","Mye"), "Neutrophil":("#d62728","Neut"), "Mast":("#d6276a","Mast")}
NEW_LIN = {"Tumor":"Tumor","Epithelial cells":"Tumor","Fibroblast":"Fibroblast",
           "Endothelial":"Endothelial","SmoothMuscle":"SmoothMuscle","Adipocyte":"Adipocyte",
           "T cells":"T / NK","NK cells":"T / NK","ILC":"T / NK","Plasma cells":"B / Plasma",
           "Macrophages":"Myeloid","DC":"Myeloid","Neutrophils":"Neutrophil","Mast cells":"Mast",
           "unassigned":None}
OLD_LIN = {"tumor_epithelial":"Tumor","a":None,"Stem.Prog":None,
           "Fibroblastic.reticular":"Fibroblast","Pericyte":"SmoothMuscle",
           "Blood.endothelial":"Endothelial","CD8.T.cell":"T / NK","gdT":"T / NK","NKT":"T / NK",
           "Ly6Clo.blood.monocytes":"Myeloid","spleen.red.pulp.macs":"Myeloid","Dendritic":"Myeloid"}

# ---- load counts (subset to marker genes) + labels ----
print("reading checkpoint ...", flush=True)
adata = ad.read_h5ad(f"{AGG}/full/cluster_checkpoint.h5ad")
old = pd.read_parquet(f"{AGG}/full/obs.parquet", columns=["cell","cell_type"]).set_index("cell")
new = pd.read_parquet("results/aggregate/full_labels.parquet",
                      columns=["cell","cell_subtype"]).set_index("cell")
adata.obs["old"] = old["cell_type"].reindex(adata.obs_names).values
adata.obs["new"] = new["cell_subtype"].reindex(adata.obs_names).values

X = adata[:, GENES].X
X = sp.csr_matrix(X) if not sp.isspmatrix_csr(X) else X
lib = np.asarray(adata.X.sum(1)).ravel(); lib[lib == 0] = 1     # library size over ALL genes
norm = X.multiply(1e4 / lib[:, None]).tocsr(); norm.data = np.log1p(norm.data)  # CPM-log
expr = pd.DataFrame(norm.toarray(), columns=GENES, index=adata.obs_names)
expr["old"] = adata.obs["old"].values; expr["new"] = adata.obs["new"].values
print("expr ready", expr.shape, flush=True)

def summarise(group_col, order):
    g = expr.groupby(group_col)
    mean = g[GENES].mean().reindex(order)                       # mean log-expr per group
    pct = g[GENES].apply(lambda d: (d > 0).mean()).reindex(order) * 100   # % expressing
    return mean, pct

m_new, p_new = summarise("new", NEW_ORDER)
m_old, p_old = summarise("old", OLD_ORDER)

# per-gene 0-1 colour scale shared across BOTH panels (comparable)
gmin = pd.concat([m_new, m_old]).min(); gmax = pd.concat([m_new, m_old]).max()
def scale(m): return ((m - gmin) / (gmax - gmin).replace(0, 1)).clip(0, 1)
c_new, c_old = scale(m_new), scale(m_old)

# ---- draw paired dot-plot ----
blocks = [b for _, b in MARKERS]
# contiguous column spans per lineage block, in order
SPANS, _s = [], None
for j, b in enumerate(blocks):
    if b != _s: SPANS.append([b, j, j]); _s = b
    else: SPANS[-1][2] = j

def panel(ax, color, pct, rows, row_lin, title):
    nr, nc = len(rows), len(GENES)
    for i in range(nr):
        for j, g in enumerate(GENES):
            ax.scatter(j, nr - 1 - i, s=(pct.iloc[i, j] / 100) * 170 + 4,
                       c=[plt.cm.Reds(color.iloc[i, j])], edgecolors="0.4", linewidths=0.3)
    ax.set_xticks(range(nc)); ax.set_xticklabels(GENES, rotation=90, fontsize=8, style="italic")
    ax.set_yticks(range(nr)); ax.set_yticklabels(rows[::-1], fontsize=9)
    # tint each cell-label row by the lineage it CLAIMS (points to its marker block)
    for t in ax.get_yticklabels():
        lin = row_lin.get(t.get_text())
        if lin: t.set_color(LIN[lin][0]); t.set_fontweight("bold")
        else: t.set_color("0.55")
    ax.set_xlim(-0.6, nc - 0.4); ax.set_ylim(-0.6, nr - 0.4 + 1.0)
    ax.set_title(title, fontsize=12, pad=16)
    # full gridlines so every dot sits in its own cell (trace row label <-> marker)
    ax.vlines(np.arange(-0.5, nc), -0.5, nr - 0.5, color="0.9", lw=0.5, zorder=0)
    ax.hlines(np.arange(-0.5, nr), -0.5, nc - 0.5, color="0.9", lw=0.5, zorder=0)
    # darker separators between lineage blocks
    for x in [s[1] - 0.5 for s in SPANS[1:]]:
        ax.vlines(x, -0.5, nr - 0.5, color="0.7", lw=1.0, zorder=0)
    # lineage group bar across the top, coloured + short code, matching row tints
    ytop = nr - 0.4 + 0.18
    for lin, c0, c1 in SPANS:
        col, code = LIN[lin]
        r, g, b = matplotlib.colors.to_rgb(col)
        txtc = "black" if (0.299*r + 0.587*g + 0.114*b) > 0.6 else "white"
        ax.add_patch(plt.Rectangle((c0 - 0.45, ytop), (c1 - c0) + 0.9, 0.58,
                                   facecolor=col, edgecolor="none", clip_on=False))
        ax.text((c0 + c1) / 2, ytop + 0.29, code, ha="center", va="center",
                fontsize=7.5, color=txtc, fontweight="bold", clip_on=False)
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.tick_params(length=0)

fig, axes = plt.subplots(1, 2, figsize=(15.5, 6.4),
                         gridspec_kw={"height_ratios":[1], "width_ratios":[1,1]})
panel(axes[0], c_old, p_old, OLD_ORDER, OLD_LIN, "Original transfer labels")
panel(axes[1], c_new, p_new, NEW_ORDER, NEW_LIN, "scVI integration labels")

fig.subplots_adjust(left=0.07, right=0.99, top=0.84, bottom=0.20, wspace=0.28)

# legends in a clean reserved bottom band: size = % expressing, color = scaled mean expr
size_handles = [Line2D([0],[0],marker="o",linestyle="",markersize=np.sqrt((p/100)*170+4),
                       markerfacecolor="0.6",markeredgecolor="0.4",label=f"{p}%") for p in [10,25,50,100]]
leg1 = fig.legend(handles=size_handles, title="% cells expressing", loc="lower left",
                  bbox_to_anchor=(0.07,0.005), ncol=4, frameon=False, fontsize=9, title_fontsize=9,
                  handletextpad=0.1, columnspacing=1.2)
fig.add_artist(leg1)
sm = plt.cm.ScalarMappable(cmap="Reds", norm=plt.Normalize(0,1)); sm.set_array([])
cax = fig.add_axes([0.55, 0.045, 0.20, 0.02])
cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
cb.set_label("mean expression (scaled per gene)", fontsize=9); cb.set_ticks([0,1]); cb.set_ticklabels(["low","high"])

fig.suptitle("Do cells express the canonical markers for what they are called?\n"
             "scVI labels track their markers (clean diagonal); transfer labels do not — "
             "immune calls lack Cd3e/Cd79a and the 1.59M-cell ‘a’ bucket is unresolved",
             fontsize=12, y=0.985)
fig.savefig(f"{OUT}/marker_dotplot_compare.png", dpi=300); plt.close(fig)
print(f"wrote {OUT}/marker_dotplot_compare.png", flush=True)
