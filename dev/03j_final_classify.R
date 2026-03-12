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

# --- 15-degree tilt, grid search for best spacing ---
tilt_deg <- 15
rad <- tilt_deg * pi / 180
mbrt4h <- mbrt4h %>%
  mutate(y_corr = y_slide_mm + x_slide_mm * tan(rad))

# Fine grid search at fixed 15 degrees
results <- expand_grid(
  spacing = seq(0.6, 1.8, by = 0.05),
  offset = seq(min(mbrt4h$y_corr) + 0.1, min(mbrt4h$y_corr) + 2.0, by = 0.05)
) %>%
  mutate(contrast = map2_dbl(spacing, offset, function(s, o) {
    centers <- seq(o, max(mbrt4h$y_corr) + s, by = s)
    dist <- sapply(mbrt4h$y_corr, function(yc) min(abs(yc - centers)))
    # Simple 50/50 split: closer to center = peak, farther = valley
    half <- s / 4  # half of half-spacing
    pk <- which(dist < half)
    vl <- which(dist > half)
    if (length(pk) < 1000 || length(vl) < 1000) return(0)
    mean(mbrt4h$p21[pk], na.rm = TRUE) - mean(mbrt4h$p21[vl], na.rm = TRUE)
  }))

best <- results %>% arrange(desc(contrast)) %>% head(5)
cat("=== Top 5 at 15 degrees ===\n")
print(best)

best_spacing <- best$spacing[1]
best_offset <- best$offset[1]
beam_centers <- seq(best_offset, max(mbrt4h$y_corr) + best_spacing, by = best_spacing)

cat(sprintf("\nSpacing: %.2fmm, offset: %.2f\n", best_spacing, best_offset))
cat(sprintf("Beam centers: %s\n", paste(round(beam_centers, 2), collapse = ", ")))

# --- Classify: simple 50/50 split at half-spacing ---
mbrt4h <- mbrt4h %>%
  mutate(dist_to_peak = sapply(y_corr, function(yc) min(abs(yc - beam_centers)))) %>%
  mutate(zone = ifelse(dist_to_peak < best_spacing / 4, "peak", "valley"))

cat("\n=== Cell counts ===\n")
mbrt4h %>% count(zone) %>% mutate(pct = round(n / sum(n) * 100, 1)) %>% print()

cat("\n=== Pathway scores ===\n")
mbrt4h %>%
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

# --- Spatial plot with tilted lines ---
p1 <- ggplot(mbrt4h, aes(x = x_slide_mm, y = y_slide_mm, color = zone)) +
  geom_point(size = 0.02, alpha = 0.3) +
  scale_color_manual(values = c("peak" = "#E74C3C", "valley" = "#3498DB")) +
  coord_fixed() +
  labs(title = sprintf("Peak/Valley: 15-deg tilt, %.2fmm spacing", best_spacing),
       color = "Zone") +
  theme_bw(base_size = 14)
for (bc in beam_centers) {
  p1 <- p1 + geom_abline(intercept = bc, slope = -tan(rad),
                           color = "yellow", linewidth = 0.6, alpha = 0.8)
}
ggsave("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/figures/mbrt4h_final_zones.pdf",
       plot = p1, width = 14, height = 10)

# With FOV labels
fov_centroids <- mbrt4h %>%
  group_by(fov) %>%
  summarise(x = mean(x_slide_mm), y = mean(y_slide_mm),
            zone = names(which.max(table(zone))), .groups = "drop")

p2 <- p1 +
  geom_text(data = fov_centroids, aes(x = x, y = y, label = fov),
            color = "black", size = 2, fontface = "bold", inherit.aes = FALSE)
ggsave("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/figures/mbrt4h_final_zones_labels.pdf",
       plot = p2, width = 14, height = 10)

cat("\nPlots saved: mbrt4h_final_zones.pdf and mbrt4h_final_zones_labels.pdf\n")
