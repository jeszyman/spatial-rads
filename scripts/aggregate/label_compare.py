#!/usr/bin/env python3
"""Two views of 'original transfer labels vs full-integration labels'.

  1. label_compare_umap.png  -- paired UMAP on the SAME scVI-integrated embedding:
     left  = each cell colored by the compartment its ORIGINAL transfer label implies
             (ImmGen/TransferData name -> nominal tumor/stroma/immune),
     right = each cell colored by its FULL-INTEGRATION compartment.
     Same coords both panels, so the visual difference IS the reassignment.

  2. label_reassignment_heatmap.png -- rows = original transfer-label lineage families,
     cols = new integration compartment, color = row-normalized % of that family's cells.
     A clean diagonal would mean the labels agree; off-diagonal mass = scrambling.

Reuses aggregate/full/umap_coords.parquet (200k integrated subsample, carries new
compartment) + aggregate/full/obs.parquet (original `cell_type`) for the UMAP, and the
full 3.27M-cell join for the honest heatmap fractions.
"""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

plt.rcParams.update({
    "font.family": "Arial", "font.size": 11, "figure.dpi": 300,
    "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
})

AGG = "/mnt/data/projects/spatial-rads/aggregate"
OUT = "results/aggregate/plots"
SEED = 0
np.random.seed(SEED)

COMP_PAL = {"tumor": "#d62728", "stroma": "#2ca02c", "immune": "#1f77b4", "unassigned": "#cccccc"}
COMP_ORDER = ["tumor", "stroma", "immune"]

# ---- original transfer label -> lineage family -> nominal compartment ----
# tumor-anchor InSituType bucket + epithelial; stroma = endothelial/fibroblast/pericyte;
# everything else is an immune-named ImmGen / TransferData class.
FAMILY = {
    "Tumor bucket (a / epithelial)": ["a", "Stem.Prog", "tumor_epithelial"],
    "Endothelial":                   ["Blood.endothelial", "Lymphatic.endothelial"],
    "Fibroblast / Pericyte":         ["Fibroblastic.reticular", "Pericyte"],
    "T / NK / NKT":                  ["gdT", "CD8.T.cell", "NKT", "NK", "ILC",
                                      "Thymic.4.8.CD3lo.DP", "Thymic.preT.DN1", "Thymic.CD4SP",
                                      "Thymic.CD8SP", "DN2a", "DN2b", "DN3", "DN4", "ISP",
                                      "Colon.Treg.Nrplo", "Spleen.Naive.CD4", "Spleen.Naive.CD8",
                                      "Spleen.Treg", "Spleen.CD4Act.48hrs", "Spleen.LN.Naive.CD4"],
    "B / Plasma":                    ["b", "B.cell", "Memory.B", "Plasma", "Plasmablast",
                                      "GC_centrocyes", "GC_centroblasts", "Spleen.CD19"],
    "Myeloid (mono / mac / DC)":     ["Ly6Clo.blood.monocytes", "Ly6Chi.blood.monocytes",
                                      "spleen.red.pulp.macs", "PerC.macrophage", "Macrophage",
                                      "small.peritoneal.macs", "Microglia", "Dendritic"],
    "Neutrophil":                    ["BM.Neutrophil", "Thio.induced.Peritoneal.Neutrophil"],
}
FAM_COMPARTMENT = {  # nominal compartment for the UMAP recolor
    "Tumor bucket (a / epithelial)": "tumor",
    "Endothelial": "stroma", "Fibroblast / Pericyte": "stroma",
    "T / NK / NKT": "immune", "B / Plasma": "immune",
    "Myeloid (mono / mac / DC)": "immune", "Neutrophil": "immune",
}
LABEL2FAM = {lab: fam for fam, labs in FAMILY.items() for lab in labs}
FAM_ORDER = list(FAMILY.keys())  # tumor, stroma-ish, then immune-ish

# =====================================================================
# Figure 1: paired UMAP on the integrated embedding
# =====================================================================
uc = pd.read_parquet(f"{AGG}/full/umap_coords.parquet",
                     columns=["cell", "compartment", "umap_int_1", "umap_int_2"])
obs = pd.read_parquet(f"{AGG}/full/obs.parquet", columns=["cell", "cell_type"])
u = uc.merge(obs, on="cell", how="left")
unmapped = sorted(set(u["cell_type"].dropna()) - set(LABEL2FAM))
assert not unmapped, f"old labels with no family: {unmapped}"
u["old_nominal"] = u["cell_type"].map(LABEL2FAM).map(FAM_COMPARTMENT)
print(f"UMAP subsample {len(u)} cells | old labels mapped, {u.old_nominal.isna().sum()} NA", flush=True)

xy = u[["umap_int_1", "umap_int_2"]].values
order = np.random.permutation(len(u))  # shuffle z-order so no class systematically overplots

def umap_panel(ax, color_key, title):
    cols = np.array([COMP_PAL.get(k, "#dddddd") for k in color_key])
    ax.scatter(xy[order, 0], xy[order, 1], s=2, c=cols[order], alpha=0.45,
               linewidths=0, rasterized=True)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title, fontsize=12)
    for sp in ax.spines.values(): sp.set_visible(False)

fig, axes = plt.subplots(1, 2, figsize=(12, 6.2))
umap_panel(axes[0], u["old_nominal"].values, "Original transfer labels\n(ImmGen / TransferData → nominal compartment)")
umap_panel(axes[1], u["compartment"].values, "Full-integration labels\n(scVI cluster-then-annotate)")
handles = [Line2D([0], [0], marker="o", linestyle="", markersize=8, markerfacecolor=COMP_PAL[k],
                  markeredgewidth=0, label=k.capitalize()) for k in COMP_ORDER]
fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=11,
           bbox_to_anchor=(0.5, 0.02))
fig.suptitle("Same integrated embedding; difference between panels = cells re-typed by integration",
             fontsize=12.5, y=0.99)
fig.tight_layout(rect=[0, 0.06, 1, 0.97])
fig.savefig(f"{OUT}/label_compare_umap.png", dpi=300); plt.close(fig)
print(f"wrote {OUT}/label_compare_umap.png", flush=True)

# =====================================================================
# Figure 2: reassignment heatmap (full 3.27M cells, honest fractions)
# =====================================================================
fo = pd.read_parquet(f"{AGG}/full/obs.parquet", columns=["cell", "cell_type"]).set_index("cell")
fn = pd.read_parquet("results/aggregate/full_labels.parquet",
                     columns=["cell", "compartment"]).set_index("cell")
d = fo.join(fn, how="inner")
d["family"] = d["cell_type"].map(LABEL2FAM)
assert d["family"].isna().sum() == 0, "unmapped old labels in full join"

n_by_fam = d.groupby("family").size().reindex(FAM_ORDER)
ct = (pd.crosstab(d["family"], d["compartment"], normalize="index") * 100)
ct = ct.reindex(index=FAM_ORDER, columns=COMP_ORDER).fillna(0)

fig, ax = plt.subplots(figsize=(6.6, 5.2))
im = ax.imshow(ct.values, cmap="Blues", vmin=0, vmax=100, aspect="auto")
ax.set_xticks(range(len(COMP_ORDER)))
ax.set_xticklabels([c.capitalize() for c in COMP_ORDER], fontsize=11)
ax.set_yticks(range(len(FAM_ORDER)))
ax.set_yticklabels([f"{f}\n(n={n_by_fam[f]:,})" for f in FAM_ORDER], fontsize=9.5)
ax.set_xlabel("Full-integration compartment", fontsize=11, labelpad=8)
ax.set_title("Where each original transfer-label family's cells\nland after full integration (row = 100%)",
             fontsize=11.5, pad=10)
for i in range(len(FAM_ORDER)):
    for j in range(len(COMP_ORDER)):
        v = ct.values[i, j]
        ax.text(j, i, f"{v:.0f}%", ha="center", va="center", fontsize=10,
                color="white" if v >= 55 else "black")
for sp in ax.spines.values(): sp.set_visible(False)
ax.set_xticks(np.arange(-.5, len(COMP_ORDER), 1), minor=True)
ax.set_yticks(np.arange(-.5, len(FAM_ORDER), 1), minor=True)
ax.grid(which="minor", color="white", linewidth=1.5)
ax.tick_params(which="both", length=0)
cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cb.set_label("% of family's cells", fontsize=10)
fig.tight_layout()
fig.savefig(f"{OUT}/label_reassignment_heatmap.png", dpi=300); plt.close(fig)
print(f"wrote {OUT}/label_reassignment_heatmap.png", flush=True)

print("\nfamily x compartment (%):"); print(ct.round(1).to_string())

# =====================================================================
# Figure 3: concordance bars -- ONE number per original lineage:
# % of its cells that integration confirms in the EXPECTED compartment.
# =====================================================================
conf = pd.DataFrame({
    "family": FAM_ORDER,
    "expected": [FAM_COMPARTMENT[f] for f in FAM_ORDER],
    "pct": [ct.loc[f, FAM_COMPARTMENT[f]] for f in FAM_ORDER],
    "n": [int(n_by_fam[f]) for f in FAM_ORDER],
}).sort_values("pct", ascending=True).reset_index(drop=True)  # worst at bottom of barh = top visually

fig, ax = plt.subplots(figsize=(8.2, 4.6))
y = np.arange(len(conf))
ax.barh(y, conf["pct"], color=[COMP_PAL[e] for e in conf["expected"]],
        edgecolor="black", linewidth=0.6, height=0.72)
for i, r in conf.iterrows():
    ax.text(r["pct"] + 1.5, i, f"{r['pct']:.0f}%", va="center", ha="left", fontsize=10.5)
ax.axvline(50, color="black", ls="--", lw=1)
ax.text(50, len(conf) - 0.35, "majority confirmed →", fontsize=8.5, color="gray30" if False else "#555555",
        ha="left", va="center")
ax.set_yticks(y)
ax.set_yticklabels([f"{r['family']}\n→ expected {r['expected']}  (n={r['n']:,})"
                    for _, r in conf.iterrows()], fontsize=9.5)
ax.set_xlim(0, 100); ax.set_xticks(range(0, 101, 25))
ax.set_xlabel("% of cells confirmed in the expected compartment by full integration", fontsize=11, labelpad=8)
ax.set_title("Original transfer-label calls confirmed by integration:\nimmune calls were least reliable",
             fontsize=12, pad=10)
handles = [Line2D([0], [0], marker="s", linestyle="", markersize=9, markerfacecolor=COMP_PAL[k],
                  markeredgecolor="black", markeredgewidth=0.5, label=k.capitalize()) for k in COMP_ORDER]
ax.legend(handles=handles, title="Expected compartment", fontsize=9, title_fontsize=9.5,
          frameon=False, loc="lower right")
for sp in ["top", "right"]: ax.spines[sp].set_visible(False)
fig.tight_layout()
fig.savefig(f"{OUT}/label_concordance_bars.png", dpi=300); plt.close(fig)
print(f"wrote {OUT}/label_concordance_bars.png", flush=True)
print("\nconcordance (expected-compartment %):"); print(conf.to_string(index=False))
