library(Seurat)
library(tidyverse)

obj <- readRDS("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/objects/seurat_clustered.rds")

mbrt4h <- obj@meta.data %>%
  as.data.frame() %>%
  filter(Condition == "MBRT_4h") %>%
  select(cell_id, fov, x_slide_mm, y_slide_mm, cell_type_validated,
         TypeI_interferon, TypeII_interferon, STING, DNA_Damage_Repair)

p21_expr <- GetAssayData(obj, layer = "data")["Cdkn1a", rownames(obj@meta.data[obj@meta.data$Condition == "MBRT_4h", ])]
mbrt4h$p21 <- as.numeric(p21_expr[rownames(mbrt4h)])

# --- Fine-grained y-axis p21 profile (0.1mm bins) ---
y_profile <- mbrt4h %>%
  mutate(y_bin = round(y_slide_mm / 0.1) * 0.1) %>%
  group_by(y_bin) %>%
  summarise(mean_p21 = mean(p21, na.rm = TRUE),
            n_cells = n(), .groups = "drop") %>%
  filter(n_cells >= 50)  # drop sparse bins

cat("=== Y-axis p21 profile (0.1mm bins, n>=50 cells) ===\n")
print(y_profile, n = 50)

# --- Smooth the profile to find peaks and valleys ---
# Use a rolling mean with window of 5 bins (0.5mm)
library(zoo)
y_profile <- y_profile %>%
  arrange(y_bin) %>%
  mutate(p21_smooth = rollmean(mean_p21, k = 5, fill = NA, align = "center"))

p1 <- ggplot(y_profile, aes(x = y_bin)) +
  geom_point(aes(y = mean_p21, size = n_cells), alpha = 0.3) +
  geom_line(aes(y = p21_smooth), color = "red", linewidth = 1.2, na.rm = TRUE) +
  labs(x = "y position (mm)", y = "Mean p21 (Cdkn1a)",
       title = "MBRT_4h: p21 y-axis profile with smoothing (0.5mm window)") +
  theme_bw()
ggsave("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/figures/mbrt4h_p21_y_smooth.pdf",
       plot = p1, width = 12, height = 5)

# --- Also do per-FOV with explicit row assignment ---
fov_stats <- mbrt4h %>%
  group_by(fov) %>%
  summarise(
    mean_y = mean(y_slide_mm),
    mean_x = mean(x_slide_mm),
    mean_p21 = mean(p21, na.rm = TRUE),
    mean_ddr = mean(DNA_Damage_Repair, na.rm = TRUE),
    mean_sting = mean(STING, na.rm = TRUE),
    mean_ifn1 = mean(TypeI_interferon, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

# --- Spatial heatmap: cells colored by p21 with FOV boundaries ---
# Use 2D binning for cleaner view
p2 <- ggplot(mbrt4h, aes(x = x_slide_mm, y = y_slide_mm)) +
  stat_summary_2d(aes(z = p21), fun = mean, binwidth = c(0.1, 0.1)) +
  scale_fill_viridis_c(option = "inferno", name = "Mean p21") +
  coord_fixed() +
  labs(title = "MBRT_4h: p21 spatial heatmap (0.1mm bins)") +
  theme_bw()
ggsave("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/figures/mbrt4h_p21_heatmap_fine.pdf",
       plot = p2, width = 12, height = 8)

# --- Multiple pathway 2D heatmaps for comparison ---
# DDR
p3 <- ggplot(mbrt4h, aes(x = x_slide_mm, y = y_slide_mm)) +
  stat_summary_2d(aes(z = DNA_Damage_Repair), fun = mean, binwidth = c(0.1, 0.1)) +
  scale_fill_viridis_c(option = "inferno", name = "Mean DDR") +
  coord_fixed() +
  labs(title = "MBRT_4h: DDR spatial heatmap (0.1mm bins)") +
  theme_bw()

# STING
p4 <- ggplot(mbrt4h, aes(x = x_slide_mm, y = y_slide_mm)) +
  stat_summary_2d(aes(z = STING), fun = mean, binwidth = c(0.1, 0.1)) +
  scale_fill_viridis_c(option = "plasma", name = "Mean STING") +
  coord_fixed() +
  labs(title = "MBRT_4h: STING spatial heatmap (0.1mm bins)") +
  theme_bw()

# Type I IFN
p5 <- ggplot(mbrt4h, aes(x = x_slide_mm, y = y_slide_mm)) +
  stat_summary_2d(aes(z = TypeI_interferon), fun = mean, binwidth = c(0.1, 0.1)) +
  scale_fill_viridis_c(option = "magma", name = "Mean IFN-I") +
  coord_fixed() +
  labs(title = "MBRT_4h: Type I IFN spatial heatmap (0.1mm bins)") +
  theme_bw()

library(patchwork)
combo <- (p2 | p3) / (p4 | p5)
ggsave("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/figures/mbrt4h_4pathway_heatmaps.pdf",
       plot = combo, width = 20, height = 16)
cat("4-pathway spatial heatmap saved.\n")

# --- Attempt peak/valley classification ---
# From the per-FOV data, valleys appear at y ~ 0.5-0.7 and y ~ 1.1-1.3
# Let's use the smoothed profile to identify local minima
cat("\n=== Smoothed p21 profile ===\n")
y_profile %>%
  filter(!is.na(p21_smooth)) %>%
  select(y_bin, mean_p21, p21_smooth, n_cells) %>%
  print(n = 50)

# Find local minima in smoothed profile
smooth_data <- y_profile %>% filter(!is.na(p21_smooth))
diffs <- diff(smooth_data$p21_smooth)
# A local min is where the derivative goes from negative to positive
local_mins <- which(diffs[-length(diffs)] < 0 & diffs[-1] > 0) + 1
local_maxs <- which(diffs[-length(diffs)] > 0 & diffs[-1] < 0) + 1

cat("\nLocal minima (valleys) in smoothed p21:\n")
print(smooth_data[local_mins, c("y_bin", "p21_smooth", "n_cells")])
cat("\nLocal maxima (peaks) in smoothed p21:\n")
print(smooth_data[local_maxs, c("y_bin", "p21_smooth", "n_cells")])

# --- Annotate profile with detected peaks/valleys ---
p6 <- ggplot(y_profile, aes(x = y_bin)) +
  geom_point(aes(y = mean_p21, size = n_cells), alpha = 0.2) +
  geom_line(aes(y = p21_smooth), color = "red", linewidth = 1.2, na.rm = TRUE) +
  geom_vline(xintercept = smooth_data$y_bin[local_mins], color = "blue",
             linetype = "dashed", linewidth = 0.8) +
  geom_vline(xintercept = smooth_data$y_bin[local_maxs], color = "darkgreen",
             linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = smooth_data$y_bin[local_mins],
           y = max(y_profile$mean_p21, na.rm = TRUE) * 0.95,
           label = "V", color = "blue", size = 5) +
  annotate("text", x = smooth_data$y_bin[local_maxs],
           y = max(y_profile$mean_p21, na.rm = TRUE),
           label = "P", color = "darkgreen", size = 5) +
  labs(x = "y position (mm)", y = "Mean p21 (Cdkn1a)",
       title = "MBRT_4h: p21 profile with detected peaks (P) and valleys (V)") +
  theme_bw()
ggsave("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/figures/mbrt4h_p21_peaks_valleys.pdf",
       plot = p6, width = 14, height = 6)
cat("\nPeak/valley annotated profile saved.\n")
