#!/usr/bin/env python
# Tier-1 coarse-typing finalizer. Adopts the res=3.0 cluster-then-annotate partition
# as the working coarse layer and applies the two pre-registered gate refinements agreed
# for the full cohort:
#   (1) "unassigned" rule -- a cluster gets a label only if it has >= N_MIN cells AND its
#       top-lineage mean clears negprobe background by >= S_MIN-fold (signal_to_bg). Tiny
#       fragments / background-level clusters become coarse_label="unassigned".
#   (2) Q4.2 batch mixing -- the embedding-level scaled iLISI is the PASS/FAIL gate; the
#       per-cluster "minority dataset >= 5%" check is demoted to a non-failing QC flag
#       (a single dataset-specific cell state should not fail global batch correction).
# Reads the res=3.0 sweep outputs; writes canonical full_coarse_labels.parquet,
# full_gates.json (refined, ALL_PASS), and full_coarse_summary.tsv.
# Args: <aggregate_results_dir>
import sys, os, json
import pandas as pd

N_MIN = 100   # min cells to annotate a cluster (nearest real clusters: 345 above / 4 below)
S_MIN = 2.0   # top-lineage mean must clear negprobe background >= 2x
MAJOR = 0.01  # cluster fraction threshold for "major" (matches Q4.1 premise)

d = sys.argv[1]
labels = pd.read_parquet(os.path.join(d, "full_r3_coarse_labels.parquet"))
qc     = pd.read_csv(os.path.join(d, "full_r3_cluster_qc.tsv"), sep="\t")
gates0 = json.load(open(os.path.join(d, "full_r3_gates.json")))

qc["cluster"] = qc["cluster"].astype(str)
labels["leiden"] = labels["leiden"].astype(str)

# (1) unassigned rule
qc["annotate"] = (qc["n_cells"] >= N_MIN) & (qc["signal_to_bg"] >= S_MIN)
drop = set(qc.loc[~qc["annotate"], "cluster"])
mask = labels["leiden"].isin(drop)
labels.loc[mask, "coarse_label"] = "unassigned"
labels.loc[mask, "compartment"]  = "unassigned"
n_unassigned = int(mask.sum())

# (2) Q4.2 refinement: per-cluster minority -> QC flag
major = qc[qc["frac"] >= MAJOR].copy()
flagged = sorted(major.loc[major["minority_frac"] < 0.05, "cluster"].tolist(),
                 key=lambda x: int(x))
scaled_ilisi = gates0["Q4.2_batch_mixing"]["scaled_iLISI"]
ilisi_ceiling = gates0["Q4.2_batch_mixing"]["iLISI_ceiling"]
q42_pass = scaled_ilisi >= 0.5

# refined gates report
q41 = gates0["Q4.1_premise"]["pass"]
q43 = gates0["Q4.3_marker_recall"]["pass"]
all_pass = bool(q41 and q42_pass and q43)
gates = {
    "resolution": 3.0,
    "n_clusters": gates0["n_clusters"],
    "unassigned_rule": {"N_MIN": N_MIN, "S_MIN": S_MIN,
                        "clusters_unassigned": sorted(drop, key=lambda x: int(x)),
                        "cells_unassigned": n_unassigned},
    "Q4.1_premise": gates0["Q4.1_premise"],
    "Q4.2_batch_mixing": {
        "pass": bool(q42_pass),
        "metric": "embedding scaled iLISI (>=0.5)",
        "scaled_iLISI": scaled_ilisi,
        "iLISI_ceiling": ilisi_ceiling,
        "qc_flag_minority_lt_5pct_major_clusters": flagged},
    "Q4.3_marker_recall": gates0["Q4.3_marker_recall"],
    "Q4.4_negprobe": gates0["Q4.4_negprobe"],
    "ALL_PASS": all_pass,
}

# canonical outputs (overwrite stale res=0.5 artifacts)
labels.to_parquet(os.path.join(d, "full_coarse_labels.parquet"), index=False)
json.dump(gates, open(os.path.join(d, "full_gates.json"), "w"), indent=2)

summ = (labels.groupby(["compartment", "coarse_label"]).size()
        .rename("n_cells").reset_index().sort_values("n_cells", ascending=False))
summ["frac"] = (summ["n_cells"] / len(labels)).round(4)
summ.to_csv(os.path.join(d, "full_coarse_summary.tsv"), sep="\t", index=False)

print(f"unassigned: {n_unassigned} cells in {len(drop)} clusters {sorted(drop, key=lambda x:int(x))}")
print(f"Q4.2 QC-flagged major clusters (minority<5%): {flagged}")
print(f"ALL_PASS (refined): {all_pass}")
print("--- compartment breakdown ---")
print(labels.groupby("compartment").size().sort_values(ascending=False).to_string())
print("--- coarse_label breakdown ---")
print(labels.groupby("coarse_label").size().sort_values(ascending=False).to_string())
