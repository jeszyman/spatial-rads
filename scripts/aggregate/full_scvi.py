#!/usr/bin/env python3
"""Full-cohort scVI integration (tier-1 integrator) -- the DAG-tracked successor to the
off-DAG run that produced the locked scvi_latent.parquet (Step 3 of plan-aggregate-refactor.md;
the `full_scvi.py` gap recorded in typing_provenance.md). Reads the raw-count MTX + obs parquet
emitted by full_export.R, trains scVI with slide_id as the batch key (nests dataset) and the
per-cell negprobe background as a continuous covariate -- the principled analog of InSituType's
negprobe model -- and writes the 30-dim latent parquet + saved model into the same full dir.

Hyperparameters are pinned to the locked run recovered from scvi_model/model.pt
(typing_provenance.md S3): n_hidden=128, n_latent=30, n_layers=2, dropout_rate=0.1,
dispersion="gene", gene_likelihood="nb", use_observed_lib_size=True, latent_distribution="normal";
setup_anndata batch_key="slide_id", continuous_covariate_keys=["neg"]; 50 epochs, seed 0.
Acceptance is re-validation of the Q4 gates, not byte-identity (scVI-GPU is not byte-reproducible).
Usage: full_scvi.py <full_dir>
"""
import sys, os
import numpy as np
import pandas as pd
import scipy.io
import anndata as ad
import scvi
import torch

FDIR = sys.argv[1]
N_LATENT, N_LAYERS, N_HIDDEN = 30, 2, 128
DROPOUT, MAX_EPOCHS, BATCH = 0.1, 50, 1024

scvi.settings.seed = 0
print(f"torch {torch.__version__} | cuda {torch.cuda.is_available()} "
      f"| device {'gpu' if torch.cuda.is_available() else 'cpu'}", flush=True)

# ---- load MTX (genes x cells from R writeMM) -> AnnData (cells x genes) ----
mtx   = os.path.join(FDIR, "mtx", "counts.mtx")
feats = [l.strip() for l in open(os.path.join(FDIR, "mtx", "features.tsv"))]
bcs   = [l.strip() for l in open(os.path.join(FDIR, "mtx", "barcodes.tsv"))]
X = scipy.io.mmread(mtx).T.tocsr()                       # -> cells x genes
adata = ad.AnnData(X=X)
adata.var_names = feats
adata.obs_names = bcs

obs = pd.read_parquet(os.path.join(FDIR, "obs.parquet")).set_index("cell")
obs = obs.loc[adata.obs_names]                            # align to matrix order
adata.obs["slide_id"] = obs["slide_id"].astype("category").values
adata.obs["dataset"]  = obs["dataset"].astype("category").values
adata.obs["neg"]      = obs["neg"].astype(float).values
print(f"adata: {adata.n_obs} cells x {adata.n_vars} genes | "
      f"datasets {dict(obs['dataset'].value_counts())}", flush=True)

# ---- scVI: slide_id batch (nests dataset), neg as continuous covariate ----
scvi.model.SCVI.setup_anndata(
    adata, batch_key="slide_id", continuous_covariate_keys=["neg"])
model = scvi.model.SCVI(adata, n_latent=N_LATENT, n_layers=N_LAYERS, n_hidden=N_HIDDEN,
                        dropout_rate=DROPOUT, dispersion="gene", gene_likelihood="nb")
model.train(max_epochs=MAX_EPOCHS, early_stopping=True, early_stopping_patience=10,
            batch_size=BATCH, accelerator="gpu" if torch.cuda.is_available() else "cpu")
print(f"trained {len(model.history['elbo_train'])} epochs", flush=True)

latent = model.get_latent_representation()
lat = pd.DataFrame(latent, index=adata.obs_names,
                   columns=[f"scVI_{i+1}" for i in range(latent.shape[1])])
lat.index.name = "cell"
lat.reset_index().to_parquet(os.path.join(FDIR, "scvi_latent.parquet"))
model.save(os.path.join(FDIR, "scvi_model"), overwrite=True)
print(f"wrote scvi_latent.parquet ({lat.shape}) + scvi_model/", flush=True)
