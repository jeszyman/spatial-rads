#!/usr/bin/env python3
# Tier-2 stroma RESCUE -- marker evidence overrides the UCell rank-argmax for the minor stromal
# lineages (Endothelial, Adipocyte, Pericyte) that UCell systematically fails to assign because
# the fibroblast-dominated background washes out its within-cell rank score. Direct analog of
# tier2_immune_rescue.py (same Mast-rule precedent: Cheng et al. BMC Bioinformatics 2025); the
# detectability gate (tier2_detectability.py) already PROVED these three clear background in
# specific r6 subclusters (Endothelial s2b 27 / Adipocyte 56 / Pericyte 14), yet UCell assigned
# 0 of them at r2 -- present-but-unresolved. This step re-examines the finer r6 subclusters and
# overrides the UCell label where the marker evidence is unambiguous.
#
# WHY identity markers, not a plain set-mean s2b argmax (which stays the detectability GATE):
# for ASSIGNING identity the full panel is corrupted by two artifacts --
#   1. AMBIENT secreted matrix/adipokine transcripts: fibrillar collagens (Col1a1/Col1a2/Col3a1/
#      Col5a1, Fn1, Mfap5, Vcan) deposit as ECM and bleed across cells, inflating the Fibroblast
#      mean on EVERY subcluster -- the exact mechanism that lets fibroblast background mask a true
#      endothelial/adipocyte subcluster; adipokines Adipoq/Lep secrete and bleed likewise.
#   2. NON-SPECIFIC / SHARED panel members: Acta2 (pericyte + myofibroblast) inflates SmoothMuscle;
#      Pdgfrb is broadly mesenchymal (fibroblast-shared); Vwf is stored/secreted (Weibel-Palade).
# So identity here uses IDENTITY-marker discriminators (below): the lineage-defining, non-secreted,
# non-shared subset of each panel. A subcluster's call = argmax over identity-marker mean s2b, and a
# rescue additionally REQUIRES >=2 of the rescue lineage's identity markers to individually clear
# the gate (corroboration -- blocks single-gene argmax wins).
#
# NOTE (panel-weak Pericyte): Pericyte and Endothelial share their best r6 subcluster (8), and
# Endothelial s2b (27) > Pericyte (14) there, so subcluster 8 calls Endothelial; Pericyte (only
# Rgs5/Pdgfrb/Notch3 on panel, Pdgfrb fibroblast-shared) may win argmax in no subcluster. A zero
# Pericyte recovery is the honest outcome, not a bug -- reported in the audit.
#
# Args: <stroma_mtx_dir> <obs.parquet> <r6_subclusters.parquet> <stroma_subtypes.parquet>
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
# Excluded vs config/lineage_markers.yaml and why: Fibroblast drops fibrillar collagens +
# Fn1/Mfap5/Vcan (ambient ECM) keeping cell-associated Pdgfra/Dcn/Lum; SmoothMuscle drops shared
# Acta2 keeping the SMC-specific contractile genes; Endothelial drops stored/secreted Vwf + the
# mesenchyme-shared Eng keeping the surface RTK/junction core; Adipocyte drops secreted adipokines
# Adipoq/Lep keeping intrinsic Cidea/Pparg/Fabp4; Pericyte keeps all 3 (panel-floor) but Pdgfrb is
# the fibroblast-shared soft member, so the >=2 corroboration leans on Rgs5/Notch3.
IDENTITY = {
    "Fibroblast":   ["Pdgfra", "Dcn", "Lum"],
    "SmoothMuscle": ["Myh11", "Tagln", "Myl9", "Actg2", "Carmn", "Tpm2"],
    "Endothelial":  ["Pecam1", "Cdh5", "Kdr", "Flt1", "Tek", "Tie1", "Clec14a", "Esam"],
    "Adipocyte":    ["Cidea", "Pparg", "Fabp4"],
    "Pericyte":     ["Rgs5", "Pdgfrb", "Notch3"],
}
# SPECIFICITY ANCHOR -- for a rescue lineage, >=1 of these lineage-EXCLUSIVE markers must
# individually clear the gate. Blocks shared-gene corroboration: Fabp4 is adipocyte+endothelial+
# macrophage (s2b 254 on the adipocyte cluster but also 39 on the endothelial one), so a
# Fabp4-driven "adipocyte" call with no lipid-droplet Cidea is false. Cidea (lipid droplet) is the
# adipocyte anchor; Pecam1/Cdh5 (CD31/VE-cadherin) the endothelial; Rgs5 the pericyte. This is the
# stroma analog of the immune rescue's >=2-corroboration guard against ambient-Ig/Xbp1 false flips.
ANCHOR = {
    "Endothelial": ["Pecam1", "Cdh5"],
    "Adipocyte":   ["Cidea"],
    "Pericyte":    ["Rgs5"],
}
RESCUE = {"Endothelial", "Adipocyte", "Pericyte"}      # lineages UCell @ r2 never assigns
RENAME = {L: L for L in RESCUE}                         # stroma names already final vocabulary

mtxdir, obs_pq, sub_pq, base_pq, out_pq, audit_tsv = sys.argv[1:7]

feats = [l.strip() for l in open(f"{mtxdir}/features.tsv")]
bcs   = [l.strip() for l in open(f"{mtxdir}/barcodes.tsv")]
cts   = scipy.io.mmread(f"{mtxdir}/counts.mtx").tocsr()      # genes x cells
assert cts.shape == (len(feats), len(bcs)), (cts.shape, len(feats), len(bcs))
gidx = {g: i for i, g in enumerate(feats)}
LIN = list(IDENTITY)
for L in LIN:
    IDENTITY[L] = [g for g in IDENTITY[L] if g in gidx]     # panel-filter
for L in ANCHOR:
    ANCHOR[L] = [g for g in ANCHOR[L] if g in gidx]

neg_map = pd.read_parquet(obs_pq, columns=["cell", "neg"]).set_index("cell")["neg"]
neg = neg_map.reindex(bcs).to_numpy(dtype=float)
assert not np.isnan(neg).any(), "missing neg for some stroma cells"

sub = pd.read_parquet(sub_pq).set_index("cell")["stroma_subcluster"].astype(str)
grp = sub.reindex(bcs)
assert not grp.isna().any(), "missing r6 subcluster for some stroma cells"
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

# specificity-anchor pass per lineage x subcluster (>=1 exclusive anchor marker clears the gate)
anchor_pass = np.ones((len(LIN), len(levels)), dtype=bool)   # non-anchored lineages unconstrained
for li, L in enumerate(LIN):
    am = [gidx[g] for g in ANCHOR.get(L, [])]
    if am:
        anchor_pass[li, :] = (gene_mean[am, :] / negf >= S_MIN).any(axis=0)

best_li = np.argmax(s2b, axis=0)
best_s2b = s2b[best_li, np.arange(len(levels))]
marker_call = np.array(LIN)[best_li]
best_npass = npass[best_li, np.arange(len(levels))]
best_anchor = anchor_pass[best_li, np.arange(len(levels))]
# rescue: eligible + dominant identity is a missed lineage + >=2 identity markers corroborate
#         + the lineage-exclusive specificity anchor clears the gate
rescue_sub = ((ncell >= N_MIN) & (best_s2b >= S_MIN) &
              np.isin(marker_call, list(RESCUE)) & (best_npass >= MIN_CORROBORATE) & best_anchor)

# BASE per-cell UCell label (r2)
base = pd.read_parquet(base_pq)
lab = base.set_index("cell")["stroma_subtype"].reindex(bcs)
assert not lab.isna().any(), "missing UCell label for some stroma cells"

cell_is_rescue = rescue_sub[code]
cell_call = marker_call[code]
rescued = lab.to_numpy().copy()
override_to = pd.Series(cell_call).map(RENAME).to_numpy()
rescued[cell_is_rescue] = override_to[cell_is_rescue]

out = base.set_index("cell").reindex(bcs).reset_index()
out["stroma_subtype_ucell"] = lab.to_numpy()
out["stroma_subtype"] = rescued
out["r6_subcluster"] = grp.to_numpy()
out["rescued"] = cell_is_rescue
out[["cell", "stroma_subcluster", "r6_subcluster",
     "stroma_subtype_ucell", "stroma_subtype", "rescued",
     "dataset", "slide_id"]].to_parquet(out_pq, index=False)

# ---- audit: one row per rescued r6 subcluster ----
rows = []
for j in np.where(rescue_sub)[0]:
    L = marker_call[j]
    prior = lab.to_numpy()[code == j]
    prior_mix = ", ".join(f"{k}:{v}" for k, v in pd.Series(prior).value_counts().items())
    anchor_s2b = {g: round(float(gene_mean[gidx[g], j] / negf[j]), 2) for g in ANCHOR.get(L, [])}
    rows.append(dict(r6_subcluster=levels[j], marker_call=L, new_label=RENAME[L],
                     n_cells=int(ncell[j]), id_s2b=round(float(best_s2b[j]), 2),
                     id_markers_pass=int(best_npass[j]),
                     anchor=";".join(f"{g}={v}" for g, v in anchor_s2b.items()),
                     prior_ucell=prior_mix))
audit = pd.DataFrame(rows).sort_values(["marker_call", "id_s2b"], ascending=[True, False])
audit.to_csv(audit_tsv, sep="\t", index=False)

print("=== stroma rescue (identity-marker override of UCell for {}) ===".format(", ".join(sorted(RESCUE))))
print(f"r6 subclusters: {len(levels)} | rescued: {int(rescue_sub.sum())} "
      f"| cells overridden: {int(cell_is_rescue.sum())}")
for L in sorted(RESCUE):
    nL = int((marker_call[(ncell >= N_MIN)] == L).sum())
    print(f"  {L}-dominant eligible subclusters (argmax): {nL}")
print("\n--- rescued subclusters ---")
print(audit.to_string(index=False) if rows else "(none)")
print("\n--- composition BEFORE (UCell) ---")
print(lab.value_counts().to_string())
print("\n--- composition AFTER (rescued) ---")
print(pd.Series(rescued).value_counts().to_string())
print(f"\nwrote {out_pq}\nwrote {audit_tsv}")
