library(Seurat)
library(tidyverse)

obj <- readRDS("/mnt/data/projects/spatial-rads/analysis/objects/seurat_clustered.rds")

mbrt4h <- obj@meta.data %>%
  as.data.frame() %>%
  filter(Condition == "MBRT_4h") %>%
  select(cell_id, fov, x_slide_mm, y_slide_mm, cell_type_validated,
         TypeI_interferon, TypeII_interferon, STING, DNA_Damage_Repair)

# --- Clean FOV centroid map for comparison with H2AX image ---
fov_centroids <- mbrt4h %>%
  group_by(fov) %>%
  summarise(
    x = mean(x_slide_mm),
    y = mean(y_slide_mm),
    n = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(y))  # top to bottom like the PPTX image

# Assign FOV rows by y position
fov_centroids <- fov_centroids %>%
  mutate(row = as.integer(cut(y, breaks = seq(max(y) + 0.2, min(y) - 0.2, length.out = 10),
                               labels = FALSE)))

cat("=== FOV map (top to bottom, matching PPTX orientation) ===\n")
cat("FOV rows sorted high-y to low-y:\n\n")
for (r in sort(unique(fov_centroids$row))) {
  row_fovs <- fov_centroids %>% filter(row == r) %>% arrange(x)
  cat(sprintf("Row %d (y ~ %.1f mm): FOVs %s (%d FOVs)\n",
              r, mean(row_fovs$y), paste(row_fovs$fov, collapse = ", "), nrow(row_fovs)))
}

# --- Large, clean spatial map with FOV labels ---
p1 <- ggplot(fov_centroids, aes(x = x, y = y)) +
  geom_tile(width = 0.28, height = 0.28, fill = "gray90", color = "gray60") +
  geom_text(aes(label = fov), size = 2.5, fontface = "bold") +
  coord_fixed() +
  labs(title = "MBRT_4h FOV grid (spatial coordinates)",
       subtitle = "Compare with H2AX image to assign peak/valley per FOV row",
       x = "x_slide_mm", y = "y_slide_mm") +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank())
ggsave("/mnt/data/projects/spatial-rads/analysis/figures/mbrt4h_fov_grid_clean.pdf",
       plot = p1, width = 14, height = 10)

# --- Overlay: cell density background + FOV labels ---
p2 <- ggplot() +
  stat_density_2d_filled(data = mbrt4h, aes(x = x_slide_mm, y = y_slide_mm),
                          contour_var = "density", bins = 20) +
  scale_fill_viridis_d(option = "viridis", guide = "none") +
  geom_text(data = fov_centroids, aes(x = x, y = y, label = fov),
            size = 2, color = "white", fontface = "bold") +
  coord_fixed() +
  labs(title = "MBRT_4h: Cell density with FOV labels",
       x = "x_slide_mm", y = "y_slide_mm") +
  theme_bw(base_size = 14)
ggsave("/mnt/data/projects/spatial-rads/analysis/figures/mbrt4h_density_fov_overlay.pdf",
       plot = p2, width = 14, height = 10)

# --- Print FOV list in a format ready for manual peak/valley annotation ---
cat("\n\n=== FOV annotation template ===\n")
cat("Assign each FOV row based on H2AX image stripes:\n")
cat("P = peak (bright H2AX band, direct beam exposure)\n")
cat("V = valley (dark gap between beams)\n")
cat("M = mixed (FOV straddles peak/valley boundary)\n\n")

for (r in sort(unique(fov_centroids$row))) {
  row_fovs <- fov_centroids %>% filter(row == r) %>% arrange(x)
  cat(sprintf("Row %d (y=%.1f): [  ] FOVs %s\n",
              r, mean(row_fovs$y), paste(row_fovs$fov, collapse = ", ")))
}

cat("\n\nOnce you label each row P/V/M, I'll classify all 149K cells accordingly.\n")
