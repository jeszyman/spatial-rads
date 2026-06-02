library(Seurat)
library(tidyverse)

obj <- readRDS("/mnt/data/projects/spatial-rads/analysis/objects/seurat_clustered.rds")

# --- Tabulate Yi's labels ---
ct_table <- sort(table(obj$ImmuneAtlas_ImmGen_Main_cell_Types), decreasing = TRUE)
cat("Cell type frequencies:\n")
print(ct_table)

ct_family <- sort(table(obj$ImmuneAtlas_ImmGen_cellFamily_Main_cell_Types), decreasing = TRUE)
cat("\nCell family frequencies:\n")
print(ct_family)

# --- Check available canonical markers ---
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
cat(paste(available_markers, collapse = ", "), "\n")

# --- Dot plot: markers vs cell types ---
p <- DotPlot(obj, features = available_markers,
             group.by = "ImmuneAtlas_ImmGen_Main_cell_Types") +
  RotatedAxis() +
  labs(title = "Marker validation: Yi's cell type labels")
ggsave("/mnt/data/projects/spatial-rads/analysis/figures/celltype_marker_validation.pdf",
       plot = p, width = 14, height = 8)
cat("Marker validation dot plot saved.\n")

# --- Characterize "a" population ---
epi_markers <- available_markers[available_markers %in% c("Epcam", "Krt8", "Krt18", "Krt14")]
cat(sprintf("\nEpithelial markers available: %s\n", paste(epi_markers, collapse = ", ")))

if (length(epi_markers) > 0) {
  p2 <- VlnPlot(obj, features = epi_markers,
                 group.by = "ImmuneAtlas_ImmGen_Main_cell_Types", pt.size = 0)
  ggsave("/mnt/data/projects/spatial-rads/analysis/figures/a_population_markers.pdf",
         plot = p2, width = 12, height = 6)
  cat("'a' population marker plot saved.\n")

  # Summary stats for "a" vs others (Seurat v5 compatible)
  a_cells <- which(obj$ImmuneAtlas_ImmGen_Main_cell_Types == "a")
  non_a <- which(obj$ImmuneAtlas_ImmGen_Main_cell_Types != "a")
  norm_data <- GetAssayData(obj, layer = "data")
  for (m in epi_markers) {
    a_mean <- mean(norm_data[m, a_cells])
    other_mean <- mean(norm_data[m, non_a])
    cat(sprintf("  %s: 'a' mean=%.3f, others mean=%.3f, ratio=%.1f\n",
                m, a_mean, other_mean, a_mean / max(other_mean, 0.001)))
  }
}
