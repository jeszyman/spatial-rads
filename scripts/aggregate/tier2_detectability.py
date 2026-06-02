#!/usr/bin/env python3
# Tier-2 lineage detectability gate -- the SAME signal_to_bg criterion that gated tier-1
# coarse typing (full_cluster.py / finalize_tier1.py), generalized from the cluster's
# argmax lineage to EVERY candidate lineage, applied identically to immune and stroma.
#
# Published precedent: this marker-set-signal-over-negprobe-background ratio is SpatialQM's
# getMeanSignalRatio (Plummer/Segato Dezem et al., "Standardized metrics for assessment and
# reproducibility of imaging-based spatial transcriptomics datasets," Nature Biotechnology
# 2025, doi:10.1038/s41587-025-02811-9). We reuse the already-validated, pre-registered
# tier-1 gate rather than importing redundant code; SpatialQM is the field-standard citation
# for the metric, not a second implementation of it.
#
# Metric (identical to full_cluster.py:65-72,121-127):
#   per-cell lineage value  = mean over the lineage's panel markers of that cell's raw counts
#   per-subcluster s2b(L)   = mean_cells(value_L) / max(mean_cells(neg), 1e-6)
#   where neg = per-cell mean of the 10 panel negative probes (obs.parquet column `neg`).
# Because mean-over-genes-then-over-cells == mean-over-cells-then-over-genes, computing the
# per-gene per-subcluster mean and averaging over a lineage's markers reproduces the tier-1
# set-level s2b exactly, while also exposing the per-marker breakdown.
#
# Decision: a lineage is RETAINED if its best subcluster (>= N_MIN cells) has set-level
# s2b >= S_MIN (the tier-1 thresholds). n_markers_pass (individual markers clearing S_MIN)
# is reported as a single-gene-artifact flag, NOT a second gate -- the gate is getMeanSignalRatio.
# Adjudicates contested lineages (e.g. B cells) on the same citable basis as tier-1.
# Args: <tag> <mtx_dir> <subclusters.parquet> <subcl_col> <obs.parquet> <markers.yaml> <lineages_csv> <out_tsv>
import sys
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import yaml

N_MIN = 100   # min cells to consider a subcluster (tier-1 gate)
S_MIN = 2.0   # set-level signal-to-background floor (tier-1 gate)

tag, mtxdir, sub_pq, subcol, obs_pq, markers_yaml, lineages_csv, out_tsv = sys.argv[1:9]
lineages = lineages_csv.split(",")

feats = [l.strip() for l in open(f"{mtxdir}/features.tsv")]
bcs   = [l.strip() for l in open(f"{mtxdir}/barcodes.tsv")]
cts   = scipy.io.mmread(f"{mtxdir}/counts.mtx").tocsr()      # genes x cells
assert cts.shape == (len(feats), len(bcs)), (cts.shape, len(feats), len(bcs))
gidx = {g: i for i, g in enumerate(feats)}

neg_map = pd.read_parquet(obs_pq, columns=["cell", "neg"]).set_index("cell")["neg"]
neg = neg_map.reindex(bcs).to_numpy(dtype=float)
assert not np.isnan(neg).any(), "missing neg for some compartment cells"

sub = pd.read_parquet(sub_pq).set_index("cell")[subcol].astype(str)
grp = sub.reindex(bcs)
assert not grp.isna().any(), "missing subcluster for some compartment cells"
levels = sorted(grp.unique(), key=lambda x: (len(x), x))
code = pd.Categorical(grp, categories=levels).codes
ind = sp.csr_matrix((np.ones(len(code)), (np.arange(len(code)), code)),
                    shape=(len(code), len(levels)))                # cells x subclusters
ncell = np.asarray(ind.sum(axis=0)).ravel().astype(int)

gene_sum = np.asarray((cts @ ind).todense())                      # genes x subclusters (count sums)
gene_mean = gene_sum / np.maximum(ncell, 1)                       # genes x subclusters mean count
neg_mean = np.asarray(neg @ ind).ravel() / np.maximum(ncell, 1)  # per-subcluster background
negf = np.maximum(neg_mean, 1e-6)

markers = yaml.safe_load(open(markers_yaml))
print(f"\n================ {tag.upper()} : tier-2 lineage detectability "
      f"(N_MIN={N_MIN}, S_MIN={S_MIN}) ================")

rows = []
for L in lineages:
    genes = [g for g in markers.get(L, []) if g in gidx]
    if not genes:
        print(f"\n--- {L}: no panel markers ---")
        rows.append(dict(lineage=L, best_subcl="", n_cells=0, set_s2b=np.nan,
                         n_markers=0, n_markers_pass=0, retain=False, markers_pass=""))
        continue
    gi = [gidx[g] for g in genes]
    set_mean = gene_mean[gi, :].mean(axis=0)                      # per-subcluster set-level mean
    set_s2b = set_mean / negf
    elig = ncell >= N_MIN
    s2b_masked = np.where(elig, set_s2b, -np.inf)
    best = int(np.argmax(s2b_masked))
    per_marker_s2b = gene_mean[gi, best] / negf[best]
    passed = [g for g, v in zip(genes, per_marker_s2b) if v >= S_MIN]
    best_s2b = float(set_s2b[best])
    retain = bool(elig[best] and best_s2b >= S_MIN)
    rows.append(dict(lineage=L, best_subcl=levels[best], n_cells=int(ncell[best]),
                     set_s2b=round(best_s2b, 2), n_markers=len(genes),
                     n_markers_pass=len(passed), retain=retain,
                     markers_pass=",".join(passed)))
    print(f"\n--- {L}: best subcluster {levels[best]} (n={ncell[best]}) | "
          f"set s2b={best_s2b:.2f} | {len(passed)}/{len(genes)} markers >= {S_MIN}x | "
          f"{'RETAIN' if retain else 'DROP'} ---")
    det = pd.DataFrame({"marker": genes,
                        "mean_count": np.round(gene_mean[gi, best], 4),
                        "s2b": np.round(per_marker_s2b, 2)}).sort_values("s2b", ascending=False)
    print(det.to_string(index=False))

summ = pd.DataFrame(rows).sort_values(["retain", "set_s2b"], ascending=[False, False])
summ.to_csv(out_tsv, sep="\t", index=False)
print(f"\n=== SUMMARY ({tag}) | neg_mean range [{neg_mean.min():.4f}, {neg_mean.max():.4f}] ===")
print(summ.to_string(index=False))
print(f"\nwrote {out_tsv}")
