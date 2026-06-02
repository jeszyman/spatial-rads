#!/usr/bin/env Rscript
# Master runner: executes all analysis scripts in order.
# Uses pre-computed seurat_clustered.rds, skipping QC/normalize/cluster steps.

cat("========== LOADING DATA ==========\n")
source("00_load_data.R")

# Scripts 01-03 are documentary — the cached object is already QC'd, normalized,
# clustered, and has validated cell types. Generate their plots from the cached object.

cat("\n========== QC SUMMARY (from cached object) ==========\n")
qc_summary <- data.frame(
  condition = names(table(obj$Condition)),
  post_filter = as.integer(table(obj$Condition))
)
print(qc_summary)
write_tsv(qc_summary, file.path(DATA_DIR, "qc_summary.tsv"))

cat("\nGenerating QC violin plots...\n")
p_qc <- VlnPlot(obj, features = c("nCount_RNA", "nFeature_RNA"),
                 group.by = "Condition", pt.size = 0, ncol = 2)
ggsave(file.path(PLOT_DIR, "qc_violins.png"), plot = p_qc, width = 16, height = 5, dpi = 150)

cat("\n========== UMAP LANDSCAPE ==========\n")
library(patchwork)
p1 <- DimPlot(obj, group.by = "seurat_clusters", raster = TRUE, label = TRUE) + labs(title = "Clusters")
p2 <- DimPlot(obj, group.by = "Condition", raster = TRUE) + labs(title = "Condition")
p3 <- DimPlot(obj, group.by = "Slide", raster = TRUE) + labs(title = "Slide (batch)")
p4 <- DimPlot(obj, group.by = "treatment", raster = TRUE) + labs(title = "Treatment")
p_umap <- (p1 | p2) / (p3 | p4)
ggsave(file.path(PLOT_DIR, "umap_landscape.png"), plot = p_umap, width = 16, height = 12, dpi = 150)

cat("\n========== CELL TYPE VALIDATION ==========\n")
marker_candidates <- c("Cd3e", "Cd3d", "Cd4", "Cd8a", "Cd8b1",
                        "Cd68", "Adgre1", "Csf1r", "Itgam",
                        "Epcam", "Krt8", "Krt18", "Krt14",
                        "Cd19", "Ms4a1", "Ncr1", "Klrb1c", "Pecam1", "Cdh5")
available_markers <- marker_candidates[marker_candidates %in% rownames(obj)]
cat(sprintf("Available markers: %d / %d\n", length(available_markers), length(marker_candidates)))

p_dot <- DotPlot(obj, features = available_markers,
                 group.by = "ImmuneAtlas_ImmGen_Main_cell_Types") +
  RotatedAxis() + labs(title = "Marker validation: Yi's cell type labels")
ggsave(file.path(PLOT_DIR, "celltype_marker_dotplot.png"), plot = p_dot, width = 14, height = 8, dpi = 150)

p_ct <- DimPlot(obj, group.by = "cell_type_validated", raster = TRUE, label = TRUE, repel = TRUE) +
  labs(title = "Cell types (validated)")
ggsave(file.path(PLOT_DIR, "umap_celltype_validated.png"), plot = p_ct, width = 10, height = 8, dpi = 150)

p_spatial <- obj@meta.data %>%
  as.data.frame() %>%
  filter(model == "flank") %>%
  ggplot(aes(x = x_slide_mm, y = y_slide_mm, color = cell_type_validated)) +
  geom_point(size = 0.1, alpha = 0.3) +
  facet_wrap(~Condition, scales = "free") +
  labs(title = "Spatial cell type distribution (flank)", color = "Cell type") +
  theme_bw() + theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))
ggsave(file.path(PLOT_DIR, "spatial_celltype_flank.png"), plot = p_spatial, width = 16, height = 12, dpi = 150)

run_script <- function(name, file) {
  cat(sprintf("\n========== %s ==========\n", name))
  tryCatch(source(file), error = function(e) {
    cat(sprintf("ERROR in %s: %s\n", file, e$message))
    cat("Continuing to next script...\n")
  })
}

run_script("SCRIPT 04: DEGs & KINETICS", "04_deg_kinetics.R")
run_script("SCRIPT 05: STRIPE DETECTION", "05_stripe_detection.R")
run_script("SCRIPT 06: H2AX VALIDATION", "06_h2ax_validation.R")
run_script("SCRIPT 07: PEAK/VALLEY CLASSIFICATION", "07_peak_valley_classify.R")
run_script("SCRIPT 08: PEAK VS VALLEY DEGs", "08_pv_degs.R")
run_script("SCRIPT 09: PEAK VS VALLEY PATHWAYS", "09_pv_pathways.R")
run_script("SCRIPT 10: PEAK VS VALLEY COMPOSITION", "10_pv_composition.R")
run_script("SCRIPT 11: SPATIAL NEIGHBORHOODS", "11_pv_spatial_neighborhoods.R")
run_script("SCRIPT 12: SBRT COMPARISON", "12_sbrt_comparison.R")
run_script("SCRIPT 13: SIGNATURE DEFINITION", "13_signatures_define.R")
run_script("SCRIPT 14: SIGNATURE SCORING", "14_signature_scoring.R")
run_script("SCRIPT 15: SIGNATURE KINETICS", "15_signature_kinetics.R")
run_script("SCRIPT 16: LAYER 4 SUMMARY", "16_summary.R")

cat("\n========== ALL SCRIPTS COMPLETE ==========\n")
cat(sprintf("Data files: %d\n", length(list.files(DATA_DIR, pattern = "\\.tsv$"))))
cat(sprintf("Plot files: %d\n", length(list.files(PLOT_DIR, pattern = "\\.png$"))))
