library(Seurat)
library(tidyverse)
library(patchwork)

obj <- readRDS("/mnt/data/projects/spatial-rads/analysis/objects/seurat_clustered.rds")

# --- Batch assessment ---
p_batch <- DimPlot(obj, group.by = "Slide", raster = TRUE) +
  labs(title = "UMAP by Slide — batch check")
ggsave("/mnt/data/projects/spatial-rads/analysis/figures/umap_by_slide.pdf",
       plot = p_batch, width = 8, height = 6)
cat("Batch assessment plot saved.\n")

# --- Landscape UMAPs ---
p1 <- DimPlot(obj, group.by = "seurat_clusters", raster = TRUE, label = TRUE) + labs(title = "Clusters")
p2 <- DimPlot(obj, group.by = "Condition", raster = TRUE) + labs(title = "Condition")
p3 <- DimPlot(obj, group.by = "ImmuneAtlas_ImmGen_Main_cell_Types",
              raster = TRUE, label = TRUE, repel = TRUE) + labs(title = "Cell type")
p4 <- DimPlot(obj, group.by = "treatment", raster = TRUE) + labs(title = "Treatment")
p_landscape <- (p1 | p2) / (p3 | p4)
ggsave("/mnt/data/projects/spatial-rads/analysis/figures/umap_landscape.pdf",
       plot = p_landscape, width = 16, height = 12)
cat("Landscape UMAP saved.\n")

# --- Spatial cell type maps (flank) ---
p_spatial <- obj@meta.data %>%
  as.data.frame() %>%
  filter(model == "flank") %>%
  ggplot(aes(x = x_slide_mm, y = y_slide_mm, color = ImmuneAtlas_ImmGen_Main_cell_Types)) +
  geom_point(size = 0.1, alpha = 0.3) +
  facet_wrap(~Condition, scales = "free") +
  labs(title = "Spatial cell type distribution (flank)", color = "Cell type") +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))
ggsave("/mnt/data/projects/spatial-rads/analysis/figures/spatial_celltype_flank.pdf",
       plot = p_spatial, width = 16, height = 12)
cat("Spatial cell type map saved.\n")

# --- Cell count summary ---
cell_counts <- obj@meta.data %>%
  as.data.frame() %>%
  count(Condition, ImmuneAtlas_ImmGen_Main_cell_Types) %>%
  pivot_wider(names_from = ImmuneAtlas_ImmGen_Main_cell_Types, values_from = n, values_fill = 0)
write.csv(cell_counts, "/mnt/data/projects/spatial-rads/analysis/tables/cell_counts_by_condition_celltype.csv", row.names = FALSE)
cat("Cell count table saved.\n")
print(cell_counts)
