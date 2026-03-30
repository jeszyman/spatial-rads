if (!exists("PATHWAY_COLS")) source("00_load_data.R")
obj <- readRDS(file.path(ANALYSIS_DIR, "objects", "seurat_clustered.rds"))
cat(sprintf("Loaded: %d genes x %d cells\n", nrow(obj), ncol(obj)))

# --- Tabulate Yi's insitutype labels ---
ct_table <- sort(table(obj$ImmuneAtlas_ImmGen_Main_cell_Types), decreasing = TRUE)
cat("Cell type frequencies:\n")
print(ct_table)

ct_family <- sort(table(obj$ImmuneAtlas_ImmGen_cellFamily_Main_cell_Types), decreasing = TRUE)
cat("\nCell family frequencies:\n")
print(ct_family)

# --- Canonical marker validation ---
marker_candidates <- c(
  "Cd3e", "Cd3d", "Cd4", "Cd8a", "Cd8b1",
  "Cd68", "Adgre1", "Csf1r", "Itgam",
  "Epcam", "Krt8", "Krt18", "Krt14",
  "Cd19", "Ms4a1",
  "Ncr1", "Klrb1c",
  "Pecam1", "Cdh5"
)
available_markers <- marker_candidates[marker_candidates %in% rownames(obj)]
cat(sprintf("\nAvailable markers: %d / %d\n", length(available_markers), length(marker_candidates)))

p_dot <- DotPlot(obj, features = available_markers,
                 group.by = "ImmuneAtlas_ImmGen_Main_cell_Types") +
  RotatedAxis() +
  labs(title = "Marker validation: Yi's cell type labels")
ggsave(file.path(PLOT_DIR, "celltype_marker_dotplot.png"), plot = p_dot, width = 14, height = 8, dpi = 150)

# --- Characterize "a" population as tumor_epithelial ---
epi_markers <- available_markers[available_markers %in% c("Epcam", "Krt8", "Krt18", "Krt14")]
norm_data <- GetAssayData(obj, layer = "data")
a_cells <- which(obj$ImmuneAtlas_ImmGen_Main_cell_Types == "a")
non_a <- which(obj$ImmuneAtlas_ImmGen_Main_cell_Types != "a")
cat("\nEpithelial marker enrichment in 'a' population:\n")
for (m in epi_markers) {
  a_mean <- mean(norm_data[m, a_cells])
  other_mean <- mean(norm_data[m, non_a])
  cat(sprintf("  %s: 'a' mean=%.3f, others mean=%.3f, ratio=%.1f\n",
              m, a_mean, other_mean, a_mean / max(other_mean, 0.001)))
}

# --- Create validated cell type label ---
obj$cell_type_validated <- obj$ImmuneAtlas_ImmGen_Main_cell_Types
obj$cell_type_validated[obj$cell_type_validated == "a"] <- "tumor_epithelial"

# --- UMAP by validated cell type ---
p_ct <- DimPlot(obj, group.by = "cell_type_validated", raster = TRUE, label = TRUE, repel = TRUE) +
  labs(title = "Cell types (validated)")
ggsave(file.path(PLOT_DIR, "umap_celltype_validated.png"), plot = p_ct, width = 10, height = 8, dpi = 150)

# --- Spatial cell type maps (flank) ---
p_spatial <- obj@meta.data %>%
  as.data.frame() %>%
  filter(model == "flank") %>%
  ggplot(aes(x = x_slide_mm, y = y_slide_mm, color = cell_type_validated)) +
  geom_point(size = 0.1, alpha = 0.3) +
  facet_wrap(~Condition, scales = "free") +
  labs(title = "Spatial cell type distribution (flank)", color = "Cell type") +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))
ggsave(file.path(PLOT_DIR, "spatial_celltype_flank.png"), plot = p_spatial, width = 16, height = 12, dpi = 150)

# --- Save updated object ---
saveRDS(obj, file.path(ANALYSIS_DIR, "objects", "seurat_clustered.rds"))
cat("Updated object saved with validated cell types.\n")
