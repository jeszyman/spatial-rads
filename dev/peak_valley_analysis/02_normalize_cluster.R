if (!exists("obj")) source("00_load_data.R")
library(future)

# Load QC'd object if not already in memory from 01
if (!"qcFlagsCell" %in% colnames(obj@meta.data) || any(obj$qcFlagsCell != "Pass", na.rm = TRUE)) {
  obj <- readRDS(file.path(ANALYSIS_DIR, "objects", "seurat_qc.rds"))
  cat(sprintf("Loaded QC'd object: %d cells\n", ncol(obj)))
}

# --- Normalize ---
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 1e4)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
cat(sprintf("Variable features: %d\n", length(VariableFeatures(obj))))
obj <- ScaleData(obj, features = VariableFeatures(obj))
obj <- RunPCA(obj, features = VariableFeatures(obj))

# --- Cluster ---
options(future.globals.maxSize = 50 * 1024^3)
plan(multicore, workers = 48)
obj <- FindNeighbors(obj, dims = 1:20)
obj <- FindClusters(obj, resolution = 0.4)
plan(sequential)

# --- UMAP ---
obj <- RunUMAP(obj, dims = 1:20, n.neighbors = 30, min.dist = 0.3)
cat(sprintf("Clusters: %d | Cells: %d\n", length(levels(obj$seurat_clusters)), ncol(obj)))

# --- Plots ---
library(patchwork)
p1 <- DimPlot(obj, group.by = "seurat_clusters", raster = TRUE, label = TRUE) + labs(title = "Clusters")
p2 <- DimPlot(obj, group.by = "Condition", raster = TRUE) + labs(title = "Condition")
p3 <- DimPlot(obj, group.by = "Slide", raster = TRUE) + labs(title = "Slide (batch)")
p4 <- DimPlot(obj, group.by = "treatment", raster = TRUE) + labs(title = "Treatment")
p_umap <- (p1 | p2) / (p3 | p4)
ggsave(file.path(PLOT_DIR, "umap_landscape.png"), plot = p_umap, width = 16, height = 12, dpi = 150)

# --- Save ---
saveRDS(obj, file.path(ANALYSIS_DIR, "objects", "seurat_clustered.rds"))
cat("Clustered object saved.\n")
