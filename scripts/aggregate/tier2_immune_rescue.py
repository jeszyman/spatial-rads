#!/usr/bin/env python3
# Tier-2 immune RESCUE -- marker evidence overrides reference (SingleR) confidence for the
# rare lineages SingleR/ImmGen systematically fails to assign at panel sparsity.
#
# Precedent (the "Mast rule"): in the BMC Bioinformatics 2025 Xenium cell-typing benchmark
# (Cheng et al.), reference-based methods missed rare populations (e.g. dendritic cells) that a
# marker-panel check recovered. Cluster-level SingleR (tier2_singler.R, res=3.0) labeled the
# 367,800 immune cells with NO B / Plasma / Neutrophil class at all; this step re-examines the
# finer res=6.0 subclusters for those three and overrides SingleR where the marker evidence is
# unambiguous.
#
# WHY a plain set-mean signal-to-background argmax is NOT used here (it is for the tier-1/tier-2
# detectability GATE, which only asks "does ANY lineage's panel clear background per subcluster"):
# for ASSIGNING identity it is corrupted by two panel artifacts --
#   1. AMBIENT secreted transcripts: plasma immunoglobulin (Jchain/Igkc/Igha/Ighm) bleeds into
#      neighboring cells; on a T-cell subcluster (Cd3e ~44x bg) the Ig genes alone push the full
#      Plasma panel mean above the T panel mean -> a false T->Plasma flip of 8k cells.
#   2. NON-SPECIFIC panel members: Cd37 (pan-leukocyte) inflates the B panel mean on T cells;
#      Xbp1 (UPR, any secretory cell) inflates the Plasma mean on epithelial-contamination cells.
# So identity here uses IDENTITY-marker discriminators (below): the lineage-defining, non-secreted,
# non-shared subset of each panel. A subcluster's call = argmax over identity-marker mean s2b, and a
# rescue additionally REQUIRES >=2 of the rescue lineage's identity markers to individually clear
# the gate (corroboration -- blocks single-gene argmax wins like Xbp1-only epithelial cells).
#
# Args: <immune_mtx_dir> <obs.parquet> <res6_subclusters.parquet> <singler_subtypes.parquet>
#       <out_parquet> <audit_tsv>
import sys
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

N_MIN = 100   # min cells to consider a subcluster (tier-1/tier-2 gate)
S_MIN = 2.0   # signal-to-background floor (tier-1/tier-2 gate)
MIN_CORROBORATE = 2   # >= this many identity markers must individually pass for a rescue

# IDENTITY (lineage-defining) marker subsets -- the non-secreted, non-shared core of each panel.
# Excluded vs config/lineage_markers.yaml and why: Plasma drops Ig (Jchain/Igkc/Igha/Ighm: ambient
# secreted) keeping the program genes Mzb1/Xbp1/Tnfrsf17; B drops Cd37 (pan-leukocyte) + the weak
# Blk/Fcrla/Tcl1/Vpreb3, keeping canonical Cd19/Cd79a/Ms4a1; Neutrophil drops the S100a8/S100a9
# alarmins (ambient) keeping the granule genes; DC drops Itgax/Ciita (shared with macrophage); NK
# drops Gzma (shared with cytotoxic T).
IDENTITY = {
    "T":          ["Cd3d", "Cd3e", "Cd3g"],
    "B":          ["Cd19", "Cd79a", "Ms4a1"],
    "Plasma":     ["Mzb1", "Xbp1", "Tnfrsf17"],
    "NK":         ["Ncr1", "Klrb1c", "Nkg7", "Klrk1"],
    "Macrophage": ["Cd68", "Csf1r", "C1qa", "C1qb", "C1qc", "Mrc1", "Cd163"],
    "Neutrophil": ["Mpo", "Elane", "Prtn3", "Ltf", "Camp", "Csf3r"],
    "Mast":       ["Cpa3", "Tpsab1", "Tpsb2", "Mcpt8"],
    "DC":         ["Clec10a", "Cd209a", "Cd209e", "Lamp3", "Ly75"],
}
RESCUE = {"B", "Plasma", "Neutrophil"}                 # lineages SingleR res=3.0 never assigns
RENAME = {"B": "B cells", "Plasma": "Plasma cells", "Neutrophil": "Neutrophils"}

mtxdir, obs_pq, sub_pq, singler_pq, out_pq, audit_tsv = sys.argv[1:7]

feats = [l.strip() for l in open(f"{mtxdir}/features.tsv")]
bcs   = [l.strip() for l in open(f"{mtxdir}/barcodes.tsv")]
cts   = scipy.io.mmread(f"{mtxdir}/counts.mtx").tocsr()      # genes x cells
assert cts.shape == (len(feats), len(bcs)), (cts.shape, len(feats), len(bcs))
gidx = {g: i for i, g in enumerate(feats)}
LIN = list(IDENTITY)
for L in LIN:
    IDENTITY[L] = [g for g in IDENTITY[L] if g in gidx]     # panel-filter

neg_map = pd.read_parquet(obs_pq, columns=["cell", "neg"]).set_index("cell")["neg"]
neg = neg_map.reindex(bcs).to_numpy(dtype=float)
assert not np.isnan(neg).any(), "missing neg for some immune cells"

sub = pd.read_parquet(sub_pq).set_index("cell")["immune_subcluster"].astype(str)
grp = sub.reindex(bcs)
assert not grp.isna().any(), "missing res6 subcluster for some immune cells"
levels = sorted(grp.unique(), key=lambda x: (len(x), x))
code = pd.Categorical(grp, categories=levels).codes
ind = sp.csr_matrix((np.ones(len(code)), (np.arange(len(code)), code)),
                    shape=(len(code), len(levels)))                # cells x subclusters
ncell = np.asarray(ind.sum(axis=0)).ravel().astype(int)
gene_mean = np.asarray((cts @ ind).todense()) / np.maximum(ncell, 1)   # genes x subclusters
negf = np.maximum(np.asarray(neg @ ind).ravel() / np.maximum(ncell, 1), 1e-6)

# identity-marker set-level s2b per lineage x subcluster, and per-marker pass counts
s2b = np.full((len(LIN), len(levels)), -np.inf)
npass = np.zeros((len(LIN), len(levels)), dtype=int)
for li, L in enumerate(LIN):
    gi = [gidx[g] for g in IDENTITY[L]]
    if gi:
        per_marker = gene_mean[gi, :] / negf                       # markers x subclusters
        s2b[li, :] = per_marker.mean(axis=0)
        npass[li, :] = (per_marker >= S_MIN).sum(axis=0)

best_li = np.argmax(s2b, axis=0)
best_s2b = s2b[best_li, np.arange(len(levels))]
marker_call = np.array(LIN)[best_li]
best_npass = npass[best_li, np.arange(len(levels))]
# rescue: eligible + dominant identity is a missed lineage + corroborated by >=2 identity markers
rescue_sub = ((ncell >= N_MIN) & (best_s2b >= S_MIN) &
              np.isin(marker_call, list(RESCUE)) & (best_npass >= MIN_CORROBORATE))

# BASE per-cell SingleR label (res=3.0)
base = pd.read_parquet(singler_pq)
lab = base.set_index("cell")["immune_subtype"].reindex(bcs)
assert not lab.isna().any(), "missing SingleR label for some immune cells"

cell_is_rescue = rescue_sub[code]
cell_call = marker_call[code]
rescued = lab.to_numpy().copy()
override_to = pd.Series(cell_call).map(RENAME).to_numpy()
rescued[cell_is_rescue] = override_to[cell_is_rescue]

out = base.set_index("cell").reindex(bcs).reset_index()
out["immune_subtype_singler"] = lab.to_numpy()
out["immune_subtype"] = rescued
out["res6_subcluster"] = grp.to_numpy()
out["rescued"] = cell_is_rescue
out[["cell", "immune_subcluster", "res6_subcluster",
     "immune_subtype_singler", "immune_subtype", "rescued",
     "dataset", "slide_id"]].to_parquet(out_pq, index=False)

# ---- audit: one row per rescued res6 subcluster ----
rows = []
for j in np.where(rescue_sub)[0]:
    L = marker_call[j]
    prior = lab.to_numpy()[code == j]
    prior_mix = ", ".join(f"{k}:{v}" for k, v in pd.Series(prior).value_counts().items())
    rows.append(dict(res6_subcluster=levels[j], marker_call=L, new_label=RENAME[L],
                     n_cells=int(ncell[j]), id_s2b=round(float(best_s2b[j]), 2),
                     id_markers_pass=int(best_npass[j]), prior_singler=prior_mix))
audit = pd.DataFrame(rows).sort_values(["marker_call", "id_s2b"], ascending=[True, False])
audit.to_csv(audit_tsv, sep="\t", index=False)

print("=== immune rescue (identity-marker override of SingleR for {}) ===".format(", ".join(sorted(RESCUE))))
print(f"res6 subclusters: {len(levels)} | rescued: {int(rescue_sub.sum())} "
      f"| cells overridden: {int(cell_is_rescue.sum())}")
nB = int((marker_call[(ncell >= N_MIN)] == "B").sum())
print(f"B-dominant eligible subclusters: {nB}  (0 => B near-absent, confirmed)")
print("\n--- rescued subclusters ---")
print(audit.to_string(index=False) if rows else "(none)")
print("\n--- composition BEFORE (SingleR) ---")
print(lab.value_counts().to_string())
print("\n--- composition AFTER (rescued) ---")
print(pd.Series(rescued).value_counts().to_string())
print(f"\nwrote {out_pq}\nwrote {audit_tsv}")
