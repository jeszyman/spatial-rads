library(Seurat)
library(tidyverse)
library(zoo)

obj <- readRDS("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/objects/seurat_clustered.rds")

mbrt4h <- obj@meta.data %>%
  as.data.frame() %>%
  filter(Condition == "MBRT_4h") %>%
  select(cell_id, fov, x_slide_mm, y_slide_mm, cell_type_validated,
         TypeI_interferon, TypeII_interferon, STING, DNA_Damage_Repair)

p21_expr <- GetAssayData(obj, layer = "data")["Cdkn1a", rownames(obj@meta.data[obj@meta.data$Condition == "MBRT_4h", ])]
mbrt4h$p21 <- as.numeric(p21_expr[rownames(mbrt4h)])

# --- Scan tilt angles 5-35 degrees ---
# For each angle, find the best spacing/offset that maximizes peak-valley p21 contrast
cat("Scanning tilt angles...\n")

angle_results <- map_dfr(seq(5, 35, by = 2), function(deg) {
  cat(sprintf("  %d deg...\n", deg))
  rad <- deg * pi / 180
  y_corr <- mbrt4h$y_slide_mm + mbrt4h$x_slide_mm * tan(rad)

  # Grid search spacing and offset
  best_contrast <- 0
  best_s <- NA
  best_o <- NA

  for (s in seq(0.6, 1.8, by = 0.1)) {
    for (o in seq(min(y_corr) + s/4, min(y_corr) + s * 1.5, by = 0.1)) {
      centers <- seq(o, max(y_corr) + s, by = s)
      dist <- sapply(y_corr, function(yc) min(abs(yc - centers)))
      q <- s / 4
      pk <- which(dist < q)
      vl <- which(dist > q * 2)
      if (length(pk) < 500 || length(vl) < 500) next
      contrast <- mean(mbrt4h$p21[pk], na.rm = TRUE) - mean(mbrt4h$p21[vl], na.rm = TRUE)
      if (contrast > best_contrast) {
        best_contrast <- contrast
        best_s <- s
        best_o <- o
      }
    }
  }

  tibble(angle = deg, best_spacing = best_s, best_offset = best_o, contrast = best_contrast)
})

cat("\n=== Tilt angle scan results ===\n")
print(angle_results, n = 20)

best_angle <- angle_results$angle[which.max(angle_results$contrast)]
best_spacing <- angle_results$best_spacing[which.max(angle_results$contrast)]
best_offset <- angle_results$best_offset[which.max(angle_results$contrast)]
cat(sprintf("\nOptimal: angle=%d deg, spacing=%.1fmm, contrast=%.4f\n",
            best_angle, best_spacing, max(angle_results$contrast)))

# --- Plot angle scan ---
p0 <- ggplot(angle_results, aes(x = angle, y = contrast)) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  geom_vline(xintercept = best_angle, color = "red", linetype = "dashed") +
  labs(x = "Tilt angle (degrees)", y = "Peak - Valley p21 contrast",
       title = "Optimal tilt angle for stripe detection") +
  theme_bw()
ggsave("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/figures/mbrt4h_tilt_angle_scan.pdf",
       plot = p0, width = 10, height = 6)

# --- Apply optimal model ---
rad <- best_angle * pi / 180
mbrt4h <- mbrt4h %>%
  mutate(y_corr = y_slide_mm + x_slide_mm * tan(rad))

beam_centers <- seq(best_offset, max(mbrt4h$y_corr) + best_spacing, by = best_spacing)
quarter <- best_spacing / 4

mbrt4h <- mbrt4h %>%
  mutate(dist_to_peak = sapply(y_corr, function(yc) min(abs(yc - beam_centers)))) %>%
  mutate(zone = case_when(
    dist_to_peak < quarter ~ "peak",
    dist_to_peak > quarter * 2 ~ "valley",
    TRUE ~ "boundary"
  ))

cat(sprintf("\nBeam centers (corrected y): %s\n", paste(round(beam_centers, 2), collapse = ", ")))
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

# --- Spatial plot ---
p1 <- ggplot(mbrt4h, aes(x = x_slide_mm, y = y_slide_mm, color = zone)) +
  geom_point(size = 0.02, alpha = 0.3) +
  scale_color_manual(values = c("peak" = "#E74C3C", "boundary" = "gray70", "valley" = "#3498DB")) +
  coord_fixed() +
  labs(title = sprintf("Peak/Valley: %d-deg tilt, %.1fmm spacing", best_angle, best_spacing),
       color = "Zone") +
  theme_bw(base_size = 14)
for (bc in beam_centers) {
  p1 <- p1 + geom_abline(intercept = bc, slope = -tan(rad),
                           color = "yellow", linewidth = 0.5, alpha = 0.7)
}
ggsave("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/figures/mbrt4h_optimal_zones.pdf",
       plot = p1, width = 14, height = 10)

# --- With FOV labels ---
fov_centroids <- mbrt4h %>%
  group_by(fov) %>%
  summarise(x = mean(x_slide_mm), y = mean(y_slide_mm),
            zone = names(which.max(table(zone))), .groups = "drop")

p2 <- ggplot(mbrt4h, aes(x = x_slide_mm, y = y_slide_mm, color = zone)) +
  geom_point(size = 0.02, alpha = 0.2) +
  scale_color_manual(values = c("peak" = "#E74C3C", "boundary" = "gray70", "valley" = "#3498DB")) +
  geom_text(data = fov_centroids, aes(x = x, y = y, label = fov),
            color = "black", size = 2, fontface = "bold") +
  coord_fixed() +
  labs(title = sprintf("Optimal zones (%d deg, %.1fmm) with FOV labels", best_angle, best_spacing),
       subtitle = "CHECK against H2AX image",
       color = "Zone") +
  theme_bw(base_size = 14)
for (bc in beam_centers) {
  p2 <- p2 + geom_abline(intercept = bc, slope = -tan(rad),
                           color = "yellow", linewidth = 0.5, alpha = 0.7)
}
ggsave("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/figures/mbrt4h_optimal_zones_labels.pdf",
       plot = p2, width = 14, height = 10)

cat("\nOptimal zone plots saved.\n")
