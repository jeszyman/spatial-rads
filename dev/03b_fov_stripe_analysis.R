library(Seurat)
library(tidyverse)

obj <- readRDS("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/objects/seurat_clustered.rds")

mbrt4h <- obj@meta.data %>%
  as.data.frame() %>%
  filter(Condition == "MBRT_4h") %>%
  select(cell_id, fov, x_slide_mm, y_slide_mm, cell_type_validated,
         TypeI_interferon, TypeII_interferon, STING, DNA_Damage_Repair)

# Get Cdkn1a (p21) expression per cell
p21_expr <- GetAssayData(obj, layer = "data")["Cdkn1a", rownames(obj@meta.data[obj@meta.data$Condition == "MBRT_4h", ])]
mbrt4h$p21 <- as.numeric(p21_expr[rownames(mbrt4h)])

# --- Per-FOV summary with spatial positions ---
fov_summary <- mbrt4h %>%
  group_by(fov) %>%
  summarise(
    n_cells = n(),
    mean_y = mean(y_slide_mm),
    mean_x = mean(x_slide_mm),
    min_y = min(y_slide_mm),
    max_y = max(y_slide_mm),
    mean_p21 = mean(p21, na.rm = TRUE),
    mean_ddr = mean(DNA_Damage_Repair, na.rm = TRUE),
    mean_sting = mean(STING, na.rm = TRUE),
    mean_ifn1 = mean(TypeI_interferon, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(mean_y)

cat("=== FOV summary sorted by y-position (top to bottom) ===\n")
print(fov_summary %>% select(fov, n_cells, mean_y, mean_p21, mean_ddr), n = 100)

# --- Identify FOV rows by y-position clustering ---
# FOVs in the same row should have similar y-positions
cat("\n=== FOV rows by y-position (0.3mm binning) ===\n")
fov_summary <- fov_summary %>%
  mutate(y_row = round(mean_y / 0.3) * 0.3)

row_summary <- fov_summary %>%
  group_by(y_row) %>%
  summarise(
    fovs = paste(fov, collapse = ","),
    n_fovs = n(),
    mean_p21 = mean(mean_p21),
    mean_ddr = mean(mean_ddr),
    mean_sting = mean(mean_sting),
    .groups = "drop"
  ) %>%
  arrange(y_row)

cat("\nRow-level p21 and DDR scores:\n")
print(row_summary, n = 30)

# --- Heatmap: FOV grid colored by p21 ---
p1 <- ggplot(fov_summary, aes(x = mean_x, y = mean_y, fill = mean_p21)) +
  geom_tile(width = 0.25, height = 0.25) +
  geom_text(aes(label = fov), size = 2, color = "white") +
  scale_fill_viridis_c(option = "inferno") +
  labs(title = "MBRT_4h: Mean p21 by FOV (spatial grid)",
       x = "x_slide_mm", y = "y_slide_mm", fill = "Mean p21") +
  coord_fixed() +
  theme_bw()
ggsave("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/figures/mbrt4h_fov_p21_grid.pdf",
       plot = p1, width = 10, height = 8)

# --- Heatmap: FOV grid colored by DDR ---
p2 <- ggplot(fov_summary, aes(x = mean_x, y = mean_y, fill = mean_ddr)) +
  geom_tile(width = 0.25, height = 0.25) +
  geom_text(aes(label = fov), size = 2, color = "white") +
  scale_fill_viridis_c(option = "inferno") +
  labs(title = "MBRT_4h: Mean DDR by FOV (spatial grid)",
       x = "x_slide_mm", y = "y_slide_mm", fill = "Mean DDR") +
  coord_fixed() +
  theme_bw()
ggsave("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/figures/mbrt4h_fov_ddr_grid.pdf",
       plot = p2, width = 10, height = 8)

# --- Line profile: mean p21 by FOV row ---
p3 <- ggplot(row_summary, aes(x = y_row, y = mean_p21)) +
  geom_line(linewidth = 1) +
  geom_point(aes(size = n_fovs), color = "red") +
  labs(x = "y position (mm)", y = "Mean p21",
       title = "p21 profile by FOV row — looking for stripe oscillation") +
  theme_bw()
ggsave("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/figures/mbrt4h_p21_by_fov_row.pdf",
       plot = p3, width = 10, height = 5)

# --- Also STING by FOV row ---
p4 <- ggplot(row_summary, aes(x = y_row, y = mean_sting)) +
  geom_line(linewidth = 1) +
  geom_point(aes(size = n_fovs), color = "blue") +
  labs(x = "y position (mm)", y = "Mean STING",
       title = "STING profile by FOV row") +
  theme_bw()
ggsave("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/figures/mbrt4h_sting_by_fov_row.pdf",
       plot = p4, width = 10, height = 5)

cat("\nAll FOV grid plots saved.\n")

# --- Now look at the H2AX image stripe geometry ---
# From image7.png, stripes appear roughly horizontal
# The MBRT beam spacing is typically 3.2mm center-to-center (400um peak width)
# Let's check if the y-range spans enough for ~2-3 stripes
cat(sprintf("\nMBRT_4h y-range: %.2f to %.2f mm (span: %.2f mm)\n",
            min(mbrt4h$y_slide_mm), max(mbrt4h$y_slide_mm),
            max(mbrt4h$y_slide_mm) - min(mbrt4h$y_slide_mm)))
cat(sprintf("MBRT_4h x-range: %.2f to %.2f mm (span: %.2f mm)\n",
            min(mbrt4h$x_slide_mm), max(mbrt4h$x_slide_mm),
            max(mbrt4h$x_slide_mm) - min(mbrt4h$x_slide_mm)))
cat("Expected ~3.2mm beam spacing -> ~1 stripe per 3.2mm of tissue\n")
