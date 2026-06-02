library(Seurat)
library(tidyverse)

obj <- readRDS("/mnt/data/projects/spatial-rads/analysis/objects/seurat_clustered.rds")

# --- Extract MBRT_4h spatial coordinates ---
mbrt4h <- obj@meta.data %>%
  as.data.frame() %>%
  filter(Condition == "MBRT_4h") %>%
  select(cell_id, fov, x_slide_mm, y_slide_mm, cell_type_validated,
         TypeI_interferon, TypeII_interferon, STING, DNA_Damage_Repair)

cat(sprintf("MBRT_4h: %d cells, %d FOVs\n", nrow(mbrt4h), length(unique(mbrt4h$fov))))
cat(sprintf("x range: %.2f - %.2f mm\n", min(mbrt4h$x_slide_mm), max(mbrt4h$x_slide_mm)))
cat(sprintf("y range: %.2f - %.2f mm\n", min(mbrt4h$y_slide_mm), max(mbrt4h$y_slide_mm)))

# --- Visualize spatial layout with DDR score as proxy ---
# The H2AX stripes should show up as bands of high DDR in the spatial plot
p1 <- ggplot(mbrt4h, aes(x = x_slide_mm, y = y_slide_mm, color = DNA_Damage_Repair)) +
  geom_point(size = 0.05, alpha = 0.3) +
  scale_color_viridis_c(option = "inferno") +
  labs(title = "MBRT_4h: DNA Damage Repair score (spatial)",
       color = "DDR score") +
  coord_fixed() +
  theme_bw()
ggsave("/mnt/data/projects/spatial-rads/analysis/figures/mbrt4h_ddr_spatial.pdf",
       plot = p1, width = 10, height = 6)
cat("DDR spatial plot saved.\n")

# --- Profile DDR along y-axis to find stripe pattern ---
# Bin cells by y position and compute mean DDR
y_profile <- mbrt4h %>%
  mutate(y_bin = round(y_slide_mm, 2)) %>%
  group_by(y_bin) %>%
  summarise(mean_ddr = mean(DNA_Damage_Repair, na.rm = TRUE),
            n_cells = n(), .groups = "drop")

p2 <- ggplot(y_profile, aes(x = y_bin, y = mean_ddr)) +
  geom_line() +
  geom_point(aes(size = n_cells), alpha = 0.5) +
  labs(x = "y_slide_mm", y = "Mean DDR score",
       title = "DDR profile along y-axis (MBRT_4h) -- looking for stripes") +
  theme_bw()
ggsave("/mnt/data/projects/spatial-rads/analysis/figures/mbrt4h_ddr_y_profile.pdf",
       plot = p2, width = 10, height = 5)
cat("Y-axis DDR profile saved.\n")

# --- Profile DDR along x-axis ---
x_profile <- mbrt4h %>%
  mutate(x_bin = round(x_slide_mm, 2)) %>%
  group_by(x_bin) %>%
  summarise(mean_ddr = mean(DNA_Damage_Repair, na.rm = TRUE),
            n_cells = n(), .groups = "drop")

p3 <- ggplot(x_profile, aes(x = x_bin, y = mean_ddr)) +
  geom_line() +
  geom_point(aes(size = n_cells), alpha = 0.5) +
  labs(x = "x_slide_mm", y = "Mean DDR score",
       title = "DDR profile along x-axis (MBRT_4h) -- looking for stripes") +
  theme_bw()
ggsave("/mnt/data/projects/spatial-rads/analysis/figures/mbrt4h_ddr_x_profile.pdf",
       plot = p3, width = 10, height = 5)
cat("X-axis DDR profile saved.\n")

# --- Also check STING along both axes (may show bystander valleys) ---
p4 <- ggplot(mbrt4h, aes(x = x_slide_mm, y = y_slide_mm, color = STING)) +
  geom_point(size = 0.05, alpha = 0.3) +
  scale_color_viridis_c(option = "plasma") +
  labs(title = "MBRT_4h: STING score (spatial)", color = "STING") +
  coord_fixed() +
  theme_bw()
ggsave("/mnt/data/projects/spatial-rads/analysis/figures/mbrt4h_sting_spatial.pdf",
       plot = p4, width = 10, height = 6)
cat("STING spatial plot saved.\n")

# --- Type I IFN spatial ---
p5 <- ggplot(mbrt4h, aes(x = x_slide_mm, y = y_slide_mm, color = TypeI_interferon)) +
  geom_point(size = 0.05, alpha = 0.3) +
  scale_color_viridis_c(option = "magma") +
  labs(title = "MBRT_4h: Type I IFN score (spatial)", color = "TypeI IFN") +
  coord_fixed() +
  theme_bw()
ggsave("/mnt/data/projects/spatial-rads/analysis/figures/mbrt4h_ifn1_spatial.pdf",
       plot = p5, width = 10, height = 6)
cat("Type I IFN spatial plot saved.\n")

cat("\nDone. Review DDR profiles to determine stripe orientation (x or y axis).\n")
cat("Then identify peak centers and classify cells.\n")
