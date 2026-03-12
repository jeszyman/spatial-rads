library(Seurat)
library(future)

obj <- readRDS("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/objects/seurat_qc.rds")
cat(sprintf("Loaded: %d genes x %d cells\n", nrow(obj), ncol(obj)))

# --- Normalize ---
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 1e4)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
cat(sprintf("Variable features: %d (requested 2000, panel has %d genes)\n",
            length(VariableFeatures(obj)), nrow(obj)))

obj <- ScaleData(obj, features = VariableFeatures(obj))
obj <- RunPCA(obj, features = VariableFeatures(obj))

pdf("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/figures/elbow_plot.pdf", width = 6, height = 4)
ElbowPlot(obj, ndims = 30)
dev.off()
cat("Elbow plot saved.\n")

# --- Cluster (parallel) ---
options(future.globals.maxSize = 50 * 1024^3)
plan(multicore, workers = 48)
obj <- FindNeighbors(obj, dims = 1:20)
obj <- FindClusters(obj, resolution = 0.4)
plan(sequential)

# --- UMAP ---
obj <- RunUMAP(obj, dims = 1:20, n.neighbors = 30, min.dist = 0.3)
cat(sprintf("Clusters: %d | Cells: %d\n", length(levels(obj$seurat_clusters)), ncol(obj)))

# --- Save ---
saveRDS(obj, "/mnt/gcs/jeszyman/projects/spatial-rads/analysis/objects/seurat_clustered.rds")
cat("Clustered object saved.\n")
