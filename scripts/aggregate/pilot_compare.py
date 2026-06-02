#!/usr/bin/env python3
"""Integration bake-off pilot, decision step. Loads the three latents (unintegrated
PCA baseline, Harmony, scVI), and on a common dataset-stratified subsample scores each
on: (1) standard scIB metrics via scib-metrics (iLISI, silhouette batch/label, ...);
(2) an imbalance-scaled iLISI that accounts for the 27/73 dataset split (raw iLISI
ceiling is 1/sum(p^2) ~= 1.66, NOT the n_batch value, so the naive scib normalization
under-credits mixing); (3) a label-independent premise gate -- leidenalg clusters,
per-cluster lineage marker recall over the per-cell negprobe floor, and coarse
compartment coverage (tumor/immune/stroma all recovered). Emits a comparison TSV and a
go/no-go decision following plan-aggregate.md v2.4.
Usage: pilot_compare.py <pilot_dir> <lineage_markers.yaml> <out_prefix>
"""
import sys, os, json
import numpy as np
import pandas as pd
import scipy.io, scipy.sparse
import anndata as ad
import scanpy as sc
import yaml

PDIR, MARKERS_YAML, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
N_SUB   = 50000
RES     = 0.5
SEED    = 0
np.random.seed(SEED)

COMPARTMENT = {
    "Epithelial": "tumor",
    "T": "immune", "B": "immune", "Plasma": "immune", "NK": "immune",
    "Macrophage": "immune", "Neutrophil": "immune", "Mast": "immune", "DC": "immune",
    "Endothelial": "stroma", "Fibroblast": "stroma", "Pericyte": "stroma",
    "SmoothMuscle": "stroma", "Adipocyte": "stroma",
}
EMB = {"X_pca": "pca_latent.parquet",
       "X_harmony": "harmony_latent.parquet",
       "X_scvi": "scvi_latent.parquet"}

# ---- assemble AnnData: raw counts X + obs + three latents ----
feats = [l.strip() for l in open(os.path.join(PDIR, "mtx", "features.tsv"))]
bcs   = [l.strip() for l in open(os.path.join(PDIR, "mtx", "barcodes.tsv"))]
X = scipy.io.mmread(os.path.join(PDIR, "mtx", "counts.mtx")).T.tocsr()
adata = ad.AnnData(X=X); adata.var_names = feats; adata.obs_names = bcs

obs = pd.read_parquet(os.path.join(PDIR, "obs.parquet")).set_index("cell").loc[adata.obs_names]
adata.obs["dataset"]   = obs["dataset"].astype(str).values
adata.obs["cell_type"] = obs["cell_type"].fillna("unassigned").astype(str).values
adata.obs["neg"]       = obs["neg"].astype(float).values

for key, fn in EMB.items():
    lat = pd.read_parquet(os.path.join(PDIR, fn)).set_index("cell").loc[adata.obs_names]
    adata.obsm[key] = lat.values.astype(np.float32)
print(f"loaded {adata.n_obs} cells | latents {list(EMB)}", flush=True)

# ---- common dataset-stratified subsample ----
idx = (adata.obs.groupby("dataset", group_keys=False)
       .apply(lambda d: d.sample(min(len(d), int(round(N_SUB * len(d) / adata.n_obs))),
                                 random_state=SEED)).index)
sub = adata[idx].copy()
p = sub.obs["dataset"].value_counts(normalize=True)
ceiling = 1.0 / np.sum(p.values ** 2)            # inverse-Simpson mixing ceiling
print(f"subsample {sub.n_obs} cells | dataset frac {dict(p.round(3))} | iLISI ceiling {ceiling:.3f}", flush=True)

# ---- scIB benchmark across the three embeddings ----
from scib_metrics.benchmark import Benchmarker
bm = Benchmarker(sub, batch_key="dataset", label_key="cell_type",
                 embedding_obsm_keys=list(EMB), n_jobs=-1)
bm.benchmark()
scib = bm.get_results(min_max_scale=False)
scib.to_csv(f"{OUT}_scib_full.tsv", sep="\t")

def col(df, needle):                              # fuzzy column pull
    for c in df.columns:
        if needle.lower() in str(c).lower():
            return c
    return None
c_ilisi = col(scib, "ilisi"); c_sl = col(scib, "silhouette label")
c_sb = col(scib, "silhouette batch"); c_bio = col(scib, "bio conservation")
c_bat = col(scib, "batch correction"); c_tot = col(scib, "total")

# ---- per-cell lineage mean raw count (for the negprobe-floor premise gate) ----
markers = yaml.safe_load(open(MARKERS_YAML))
lineage_mean = {}
for lin, genes in markers.items():
    gi = [sub.var_names.get_loc(g) for g in genes if g in sub.var_names]
    lineage_mean[lin] = np.asarray(sub[:, gi].X.mean(axis=1)).ravel()
LM = pd.DataFrame(lineage_mean, index=sub.obs_names)
neg = sub.obs["neg"].values

rows = []
for key in EMB:
    sc.pp.neighbors(sub, use_rep=key, n_neighbors=15, random_state=SEED)
    sc.tl.leiden(sub, resolution=RES, flavor="leidenalg", random_state=SEED,
                 key_added=f"cl_{key}")
    cl = sub.obs[f"cl_{key}"]
    n_clusters = cl.nunique()

    both = []                                     # per-cluster: both datasets present >=10%
    above = []; comps = set(); s2b_min = np.inf
    for c in cl.cat.categories:
        m = (cl == c).values
        ds = sub.obs.loc[m, "dataset"].value_counts(normalize=True)
        both.append(int((ds >= 0.10).sum() >= 2))
        lin_means = LM.loc[m].mean(axis=0)
        top_lin = lin_means.idxmax()
        s2b = lin_means[top_lin] / max(neg[m].mean(), 1e-6)
        above.append(int(s2b > 1.0)); s2b_min = min(s2b_min, s2b)
        comps.add(COMPARTMENT.get(top_lin, "other"))
    comps &= {"tumor", "immune", "stroma"}

    raw_ilisi = float(scib.loc[key, c_ilisi]) + 1.0          # scib iLISI is (raw-1)/(nb-1); nb=2
    scaled_ilisi = float(np.clip((raw_ilisi - 1) / (ceiling - 1), 0, 1))
    frac_both = float(np.mean(both)); frac_above = float(np.mean(above))
    sl = float(scib.loc[key, c_sl]) if c_sl else np.nan
    n_comp = len(comps)
    bio = np.nanmean([sl, n_comp / 3.0, frac_above])
    total = 0.4 * scaled_ilisi + 0.6 * bio
    premise = (n_comp == 3) and (frac_above >= 0.8) and (scaled_ilisi >= 0.5)

    rows.append(dict(
        embedding=key, n_clusters=n_clusters,
        scib_ilisi=round(float(scib.loc[key, c_ilisi]), 4),
        scaled_ilisi=round(scaled_ilisi, 4),
        pct_clusters_both_ds=round(frac_both, 3),
        silhouette_label=round(sl, 4) if c_sl else np.nan,
        silhouette_batch=round(float(scib.loc[key, c_sb]), 4) if c_sb else np.nan,
        compartments=n_comp, compartment_set="+".join(sorted(comps)) or "none",
        frac_clusters_above_negbg=round(frac_above, 3), min_s2b=round(float(s2b_min), 2),
        bio=round(float(bio), 4), pilot_total=round(float(total), 4),
        premise_pass=bool(premise)))
    print(f"{key}: clusters={n_clusters} scaled_iLISI={scaled_ilisi:.3f} "
          f"comps={n_comp}/3({'+'.join(sorted(comps))}) above_negbg={frac_above:.2f} "
          f"sil_label={sl:.3f} total={total:.3f} premise={'PASS' if premise else 'FAIL'}", flush=True)

res = pd.DataFrame(rows).set_index("embedding")
res.to_csv(f"{OUT}_compare.tsv", sep="\t")

# ---- decision (plan-aggregate.md v2.4) ----
h, s, pca = res.loc["X_harmony"], res.loc["X_scvi"], res.loc["X_pca"]
if not h.premise_pass and not s.premise_pass:
    if pca.compartments == 3:
        decision = "ESCALATE: both integrators fail premise but PCA recovers 3 compartments -> integration is over-correcting; revisit integration params/covariates before full run."
    else:
        decision = "ESCALATE: no embedding (incl. unintegrated PCA) recovers 3 compartments -> upstream sparsity/embedding problem; the cluster-then-annotate premise itself is unmet."
elif s.premise_pass and (not h.premise_pass or s.pilot_total > 1.05 * h.pilot_total):
    decision = f"USE scVI (total {s.pilot_total:.3f} vs Harmony {h.pilot_total:.3f}; >5% better or Harmony failed premise)."
else:
    decision = f"USE Harmony (total {h.pilot_total:.3f} vs scVI {s.pilot_total:.3f}; scVI not >5% better -> prefer the simpler R-native integrator)."

print("\n=== PILOT DECISION ===\n" + decision, flush=True)
json.dump({"decision": decision, "ceiling": ceiling,
           "results": res.reset_index().to_dict("records")},
          open(f"{OUT}_decision.json", "w"), indent=2)
print(res.to_string())
