library(Seurat)
library(tidyverse)
library(zoo)

obj <- readRDS("/mnt/data/projects/spatial-rads/analysis/objects/seurat_clustered.rds")

mbrt4h <- obj@meta.data %>%
  as.data.frame() %>%
  filter(Condition == "MBRT_4h") %>%
  select(cell_id, fov, x_slide_mm, y_slide_mm, cell_type_validated,
         TypeI_interferon, TypeII_interferon, STING, DNA_Damage_Repair)

p21_expr <- GetAssayData(obj, layer = "data")["Cdkn1a", rownames(obj@meta.data[obj@meta.data$Condition == "MBRT_4h", ])]
mbrt4h$p21 <- as.numeric(p21_expr[rownames(mbrt4h)])

# --- Tilt correction ---
# Stripes tilt downward ~15 degrees left to right
# A stripe at x=0,y=y0 is at y = y0 - x*tan(15) at position x
# Corrected coordinate: y_corr = y + x*tan(15)
# Cells at the same y_corr are on the same stripe

tilt_deg <- 15
tilt_rad <- tilt_deg * pi / 180
mbrt4h <- mbrt4h %>%
  mutate(y_corr = y_slide_mm + x_slide_mm * tan(tilt_rad))

cat(sprintf("Tilt correction: %.1f degrees, tan = %.4f\n", tilt_deg, tan(tilt_rad)))
cat(sprintf("y_corr range: %.2f to %.2f mm (span: %.2f)\n",
            min(mbrt4h$y_corr), max(mbrt4h$y_corr),
            max(mbrt4h$y_corr) - min(mbrt4h$y_corr)))

# --- Smoothed p21 profile along corrected axis ---
y_profile <- mbrt4h %>%
  mutate(y_bin = round(y_corr / 0.1) * 0.1) %>%
  group_by(y_bin) %>%
  summarise(mean_p21 = mean(p21, na.rm = TRUE),
            n_cells = n(), .groups = "drop") %>%
  filter(n_cells >= 50) %>%
  arrange(y_bin) %>%
  mutate(p21_smooth = rollmean(mean_p21, k = 5, fill = NA, align = "center"))

# Find peaks and valleys
smooth_data <- y_profile %>% filter(!is.na(p21_smooth))
diffs <- diff(smooth_data$p21_smooth)
local_mins <- which(diffs[-length(diffs)] < 0 & diffs[-1] > 0) + 1
local_maxs <- which(diffs[-length(diffs)] > 0 & diffs[-1] < 0) + 1

cat("\n=== Tilt-corrected p21 profile: local minima (valleys) ===\n")
print(smooth_data[local_mins, c("y_bin", "p21_smooth", "n_cells")])
cat("\n=== Local maxima (peaks) ===\n")
print(smooth_data[local_maxs, c("y_bin", "p21_smooth", "n_cells")])

p1 <- ggplot(y_profile, aes(x = y_bin)) +
  geom_point(aes(y = mean_p21, size = n_cells), alpha = 0.2) +
  geom_line(aes(y = p21_smooth), color = "red", linewidth = 1.2, na.rm = TRUE) +
  geom_vline(xintercept = smooth_data$y_bin[local_mins], color = "blue",
             linetype = "dashed", linewidth = 0.8) +
  geom_vline(xintercept = smooth_data$y_bin[local_maxs], color = "darkgreen",
             linetype = "dashed", linewidth = 0.8) +
  labs(x = "Tilt-corrected y (mm)", y = "Mean p21",
       title = sprintf("p21 profile along tilt-corrected axis (%d deg)", tilt_deg)) +
  theme_bw()
ggsave("/mnt/data/projects/spatial-rads/analysis/figures/mbrt4h_p21_tilted_profile.pdf",
       plot = p1, width = 14, height = 6)

# --- Now classify FOVs using corrected coordinate ---
fov_stats <- mbrt4h %>%
  group_by(fov) %>%
  summarise(
    y_corr_mean = mean(y_corr),
    x_mean = mean(x_slide_mm),
    y_mean = mean(y_slide_mm),
    mean_p21 = mean(p21, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  arrange(y_corr_mean)

# Bin FOVs by corrected y
fov_stats <- fov_stats %>%
  mutate(corr_row = as.integer(cut(y_corr_mean,
                                    breaks = seq(min(y_corr_mean) - 0.1,
                                                 max(y_corr_mean) + 0.5,
                                                 by = 0.5))))

cat("\n=== FOV rows by tilt-corrected y ===\n")
fov_stats %>%
  group_by(corr_row) %>%
  summarise(
    y_corr = mean(y_corr_mean),
    fovs = paste(fov, collapse = ","),
    n_fovs = n(),
    mean_p21 = mean(mean_p21),
    .groups = "drop"
  ) %>%
  print(n = 20)

# --- Grid search for best spacing using corrected coordinate ---
spacings <- seq(0.5, 1.5, by = 0.05)
offsets <- seq(min(mbrt4h$y_corr), min(mbrt4h$y_corr) + 1.5, by = 0.05)

results <- expand_grid(spacing = spacings, offset = offsets) %>%
  mutate(contrast = map2_dbl(spacing, offset, function(s, o) {
    beam_centers <- seq(o, max(mbrt4h$y_corr) + s, by = s)
    dist_to_peak <- sapply(mbrt4h$y_corr, function(yc) min(abs(yc - beam_centers)))
    quarter <- s / 4
    peak_p21 <- mean(mbrt4h$p21[dist_to_peak < quarter], na.rm = TRUE)
    valley_p21 <- mean(mbrt4h$p21[dist_to_peak > quarter * 2], na.rm = TRUE)
    n_peak <- sum(dist_to_peak < quarter)
    n_valley <- sum(dist_to_peak > quarter * 2)
    if (n_peak < 500 || n_valley < 500) return(0)
    peak_p21 - valley_p21
  }))

best <- results %>% arrange(desc(contrast)) %>% head(10)
cat("\n=== Top 10 stripe models (tilt-corrected, POSITIVE contrast = peaks have higher p21) ===\n")
print(best)

# Use best model
best_spacing <- best$spacing[1]
best_offset <- best$offset[1]
beam_centers <- seq(best_offset, max(mbrt4h$y_corr) + best_spacing, by = best_spacing)
cat(sprintf("\nBest: spacing=%.2fmm, offset=%.2f, contrast=%.4f\n",
            best_spacing, best_offset, best$contrast[1]))
cat(sprintf("Beam centers (corrected y): %s\n", paste(round(beam_centers, 2), collapse = ", ")))

# Classify
mbrt4h <- mbrt4h %>%
  mutate(dist_to_peak = sapply(y_corr, function(yc) min(abs(yc - beam_centers)))) %>%
  mutate(zone = case_when(
    dist_to_peak < best_spacing / 4 ~ "peak",
    dist_to_peak > best_spacing / 2 ~ "valley",
    TRUE ~ "boundary"
  ))

cat("\n=== Cell classification ===\n")
mbrt4h %>% count(zone) %>% mutate(pct = round(n / sum(n) * 100, 1)) %>% print()

cat("\n=== Pathway scores by zone ===\n")
mbrt4h %>%
  filter(zone != "boundary") %>%
  group_by(zone) %>%
  summarise(
    n = n(),
    mean_p21 = mean(p21, na.rm = TRUE),
    mean_DDR = mean(DNA_Damage_Repair, na.rm = TRUE),
    mean_STING = mean(STING, na.rm = TRUE),
    mean_IFN1 = mean(TypeI_interferon, na.rm = TRUE),
    mean_IFN2 = mean(TypeII_interferon, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print()

# --- Verification: spatial plot with tilted stripe lines ---
p2 <- ggplot(mbrt4h, aes(x = x_slide_mm, y = y_slide_mm, color = zone)) +
  geom_point(size = 0.02, alpha = 0.3) +
  scale_color_manual(values = c("peak" = "#E74C3C", "boundary" = "gray70", "valley" = "#3498DB")) +
  coord_fixed() +
  labs(title = sprintf("Peak/Valley with %d-deg tilt (spacing=%.2fmm)", tilt_deg, best_spacing),
       subtitle = "Compare tilted bands with H2AX image",
       color = "Zone") +
  theme_bw(base_size = 14)

# Add tilted beam center lines
for (bc in beam_centers) {
  x_range <- range(mbrt4h$x_slide_mm)
  # y_corr = y + x*tan(tilt) = bc => y = bc - x*tan(tilt)
  p2 <- p2 + geom_abline(intercept = bc, slope = -tan(tilt_rad),
                           color = "yellow", linewidth = 0.5, alpha = 0.7)
}
ggsave("/mnt/data/projects/spatial-rads/analysis/figures/mbrt4h_tilted_zones.pdf",
       plot = p2, width = 14, height = 10)

# --- With FOV labels ---
fov_centroids <- mbrt4h %>%
  group_by(fov, zone) %>%
  summarise(x = mean(x_slide_mm), y = mean(y_slide_mm), .groups = "drop") %>%
  group_by(fov) %>%
  slice_max(n = 1, order_by = zone) %>%  # take one entry per FOV
  ungroup()

p3 <- ggplot(mbrt4h, aes(x = x_slide_mm, y = y_slide_mm, color = zone)) +
  geom_point(size = 0.02, alpha = 0.2) +
  scale_color_manual(values = c("peak" = "#E74C3C", "boundary" = "gray70", "valley" = "#3498DB")) +
  geom_text(data = fov_centroids, aes(x = x, y = y, label = fov),
            color = "black", size = 2, fontface = "bold") +
  coord_fixed() +
  labs(title = "Tilted peak/valley zones with FOV labels",
       subtitle = "CHECK: do tilted bands match H2AX image?",
       color = "Zone") +
  theme_bw(base_size = 14)
for (bc in beam_centers) {
  p3 <- p3 + geom_abline(intercept = bc, slope = -tan(tilt_rad),
                           color = "yellow", linewidth = 0.5, alpha = 0.7)
}
ggsave("/mnt/data/projects/spatial-rads/analysis/figures/mbrt4h_tilted_zones_labels.pdf",
       plot = p3, width = 14, height = 10)

cat("\nTilted zone plots saved. Check mbrt4h_tilted_zones_labels.pdf against H2AX image.\n")
