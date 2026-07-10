#!/usr/bin/env python3
"""Concrete audit: how many cells did the original transfer labels get right?
Bars sized by CELL COUNT (so magnitude shows), split into confirmed vs overturned
by the full integration. The 'a' bucket (InSituType catch-all, no cell type) is its
own bar. Reference = integration compartment (cluster-level, more reliable at this
panel sparsity), so 'overturned' = the better method disagrees.
"""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
plt.rcParams.update({"font.family":"Arial","font.size":11,"pdf.fonttype":42,
                     "ps.fonttype":42,"svg.fonttype":"none",
                     "axes.spines.top":False,"axes.spines.right":False})
OUT="results/aggregate/plots"

STROMA={"Blood.endothelial","Lymphatic.endothelial","Fibroblastic.reticular","Pericyte"}
TUMOR_NAMED={"tumor_epithelial","Stem.Prog"}   # named tumor classes (NOT the 'a' catch-all)
def claimed(l):
    if l=="a": return "a"
    if l in STROMA: return "stroma"
    if l in TUMOR_NAMED: return "tumor"
    return "immune"                              # all remaining ImmGen/TransferData names

old=pd.read_parquet("/mnt/data/projects/spatial-rads/aggregate/full/obs.parquet",
                    columns=["cell","cell_type"]).set_index("cell")
new=pd.read_parquet("results/aggregate/full_labels.parquet",
                    columns=["cell","compartment"]).set_index("cell")
d=old.join(new,how="inner")
d["claim"]=d["cell_type"].map(claimed)
d["ok"]=np.where(d["claim"]=="a", d["compartment"].eq("tumor"), d["claim"]==d["compartment"])
N=len(d)

GROUPS=[("immune","Immune-named\n(T/NK, B, myeloid, ...)"),
        ("a",     "'a' catch-all\n(InSituType, no cell type)"),
        ("stroma","Stroma-named\n(endothelial, fibroblast, pericyte)"),
        ("tumor", "Tumor-named\n(epithelial, Stem.Prog)")]
rows=[]
for key,lab in GROUPS:
    sub=d[d["claim"]==key]; n=len(sub); conf=int(sub["ok"].sum()); over=n-conf
    rows.append(dict(key=key,lab=lab,n=n,conf=conf,over=over,pct_conf=conf/n*100))
df=pd.DataFrame(rows)

GREEN,RED,GRAY="#2ca02c","#d62728","#9e9e9e"
fig,ax=plt.subplots(figsize=(10,4.6))
y=np.arange(len(df))[::-1]                        # first group at top
for yi,r in zip(y,df.itertuples()):
    if r.key=="a":                                # catch-all: single gray bar
        ax.barh(yi,r.n,color=GRAY,edgecolor="white")
        ax.text(r.n+25000,yi,f"{r.n/1e6:.2f}M  (65% resolve to tumor)",va="center",fontsize=9.5,color="0.3")
    else:
        ax.barh(yi,r.conf,color=GREEN,edgecolor="white")
        ax.barh(yi,r.over,left=r.conf,color=RED,edgecolor="white")
        if r.conf>=80000: ax.text(r.conf/2,yi,f"{r.pct_conf:.0f}%",va="center",ha="center",color="white",fontsize=10,fontweight="bold")
        if r.over>=80000: ax.text(r.conf+r.over/2,yi,f"{100-r.pct_conf:.0f}%",va="center",ha="center",color="white",fontsize=10,fontweight="bold")
        ax.text(r.n+25000,yi,f"{r.n/1e6:.2f}M",va="center",fontsize=9.5,color="0.3")
ax.set_yticks(y); ax.set_yticklabels(df["lab"],fontsize=9.5)
ax.set_xlim(0,2.0e6); ax.set_xticks(np.arange(0,2.01e6,0.5e6))
ax.set_xticklabels(["0","0.5M","1.0M","1.5M","2.0M"])
ax.set_xlabel("Number of cells",fontsize=11,labelpad=8)
ax.set_title("Only 19% of cells carried a transfer label that integration confirmed\n"
             "Immune transfer labels were 75% overturned; half the cohort was an uninterpretable bucket",
             fontsize=11.5,pad=12)
handles=[Line2D([0],[0],marker="s",linestyle="",markersize=10,markerfacecolor=GREEN,markeredgecolor="none",label="Confirmed by integration"),
         Line2D([0],[0],marker="s",linestyle="",markersize=10,markerfacecolor=RED,markeredgecolor="none",label="Overturned by integration"),
         Line2D([0],[0],marker="s",linestyle="",markersize=10,markerfacecolor=GRAY,markeredgecolor="none",label="'a' catch-all (no cell type)")]
ax.legend(handles=handles,fontsize=9,frameon=False,loc="lower right",bbox_to_anchor=(1.0,-0.02))
fig.tight_layout()
fig.savefig(f"{OUT}/transfer_label_verdict.png",dpi=300); plt.close(fig)
print(f"wrote {OUT}/transfer_label_verdict.png")
print(df[["key","n","conf","over","pct_conf"]].to_string(index=False))
print(f"named-confirmed / all = {df[df.key!='a']['conf'].sum()/N*100:.1f}%")
