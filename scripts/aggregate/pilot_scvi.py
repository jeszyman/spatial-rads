#!/usr/bin/env python3
"""Integration bake-off pilot, scVI arm. Reads the raw-count MTX + obs parquet from
pilot_export.R, trains scVI with slide_id as the batch key and the per-cell negprobe
background as a continuous covariate (the principled analog of InSituType's negprobe
model), and writes a 30-dim latent parquet. n_latent matches the Harmony arm's 30 PCs
so the only difference downstream is the integrator. GPU if available.
Usage: pilot_scvi.py <pilot_dir>
"""
import sys, os
import numpy as np
import pandas as pd
import scipy.io, scipy.sparse
import anndata as ad
import scvi
import torch

PDIR    = sys.argv[1]
N_LATENT = 30
MAX_EPOCHS = 50

scvi.settings.seed = 0
print(f"torch {torch.__version__} | cuda {torch.cuda.is_available()} "
      f"| device {'gpu' if torch.cuda.is_available() else 'cpu'}", flush=True)

# ---- load MTX (genes x cells from R writeMM) -> AnnData (cells x genes) ----
mtx   = os.path.join(PDIR, "mtx", "counts.mtx")
feats = [l.strip() for l in open(os.path.join(PDIR, "mtx", "features.tsv"))]
bcs   = [l.strip() for l in open(os.path.join(PDIR, "mtx", "barcodes.tsv"))]
X = scipy.io.mmread(mtx).T.tocsr()                  # -> cells x genes
adata = ad.AnnData(X=X)
adata.var_names = feats
adata.obs_names = bcs

obs = pd.read_parquet(os.path.join(PDIR, "obs.parquet")).set_index("cell")
obs = obs.loc[adata.obs_names]                       # align to matrix order
for c in ["slide_id", "dataset", "sample_id", "cell_type", "condition"]:
    adata.obs[c] = obs[c].astype("category").values
adata.obs["neg"] = obs["neg"].astype(float).values
print(f"adata: {adata.n_obs} cells x {adata.n_vars} genes | "
      f"datasets {dict(obs['dataset'].value_counts())}", flush=True)

# ---- scVI: slide_id batch (nests dataset), neg as continuous covariate ----
scvi.model.SCVI.setup_anndata(
    adata, batch_key="slide_id", continuous_covariate_keys=["neg"])
model = scvi.model.SCVI(adata, n_latent=N_LATENT, n_layers=2, gene_likelihood="nb")
model.train(max_epochs=MAX_EPOCHS, early_stopping=True,
            early_stopping_patience=10, batch_size=1024,
            accelerator="gpu" if torch.cuda.is_available() else "cpu")

n_ep = len(model.history["elbo_train"])
print(f"trained {n_ep} epochs", flush=True)

latent = model.get_latent_representation()
lat = pd.DataFrame(latent, index=adata.obs_names,
                   columns=[f"scVI_{i+1}" for i in range(latent.shape[1])])
lat.index.name = "cell"
lat.reset_index().to_parquet(os.path.join(PDIR, "scvi_latent.parquet"))
model.save(os.path.join(PDIR, "scvi_model"), overwrite=True)
print(f"wrote scvi_latent.parquet ({lat.shape}) + scvi_model/", flush=True)
