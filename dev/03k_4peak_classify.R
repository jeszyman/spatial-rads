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

# 15-degree tilt
rad <- 15 * pi / 180
mbrt4h <- mbrt4h %>%
  mutate(y_corr = y_slide_mm + x_slide_mm * tan(rad))

y_range <- range(mbrt4h$y_corr)
cat(sprintf("Corrected y range: %.2f to %.2f (span %.2f)\n", y_range[1], y_range[2], diff(y_range)))

# --- Force exactly 4 peaks: search spacing 0.8-1.3mm ---
results <- expand_grid(
  spacing = seq(0.8, 1.3, by = 0.02),
  offset = seq(y_range[1] + 0.1, y_range[1] + 1.3, by = 0.05)
) %>%
  mutate(n_peaks_in_range = map2_int(spacing, offset, function(s, o) {
    centers <- seq(o, y_range[2] + 0.1, by = s)
    # Count peaks that actually overlap with tissue cells
    sum(sapply(centers, function(c) sum(abs(mbrt4h$y_corr - c) < s/4) > 100))
  })) %>%
  filter(n_peaks_in_range == 4) %>%
  mutate(contrast = map2_dbl(spacing, offset, function(s, o) {
    centers <- seq(o, y_range[2] + 0.1, by = s)
    dist <- sapply(mbrt4h$y_corr, function(yc) min(abs(yc - centers)))
    pk <- which(dist < s / 4)
    vl <- which(dist > s / 4)
    if (length(pk) < 1000 || length(vl) < 1000) return(0)
    mean(mbrt4h$p21[pk], na.rm = TRUE) - mean(mbrt4h$p21[vl], na.rm = TRUE)
  }))

best <- results %>% arrange(desc(contrast)) %>% head(5)
cat("\n=== Top 5 models with exactly 4 peaks ===\n")
print(best)

best_spacing <- best$spacing[1]
best_offset <- best$offset[1]
beam_centers <- seq(best_offset, y_range[2] + 0.1, by = best_spacing)
# Trim to only centers with cells nearby
beam_centers <- beam_centers[sapply(beam_centers, function(c) sum(abs(mbrt4h$y_corr - c) < best_spacing/4) > 100)]

cat(sprintf("\nSpacing: %.2fmm, %d peaks\n", best_spacing, length(beam_centers)))
cat(sprintf("Beam centers (corrected y): %s\n", paste(round(beam_centers, 2), collapse = ", ")))

# Classify
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

# --- Spatial plot ---
p1 <- ggplot(mbrt4h, aes(x = x_slide_mm, y = y_slide_mm, color = zone)) +
  geom_point(size = 0.02, alpha = 0.3) +
  scale_color_manual(values = c("peak" = "#E74C3C", "valley" = "#3498DB")) +
  coord_fixed() +
  labs(title = sprintf("4 peaks, 15-deg tilt, %.2fmm spacing", best_spacing),
       color = "Zone") +
  theme_bw(base_size = 14)
for (bc in beam_centers) {
  p1 <- p1 + geom_abline(intercept = bc, slope = -tan(rad),
                           color = "yellow", linewidth = 0.6, alpha = 0.8)
}
ggsave("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/figures/mbrt4h_4peak_zones.pdf",
       plot = p1, width = 14, height = 10)

# With FOV labels
fov_centroids <- mbrt4h %>%
  group_by(fov) %>%
  summarise(x = mean(x_slide_mm), y = mean(y_slide_mm), .groups = "drop")

p2 <- p1 +
  geom_text(data = fov_centroids, aes(x = x, y = y, label = fov),
            color = "black", size = 2, fontface = "bold", inherit.aes = FALSE)
ggsave("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/figures/mbrt4h_4peak_zones_labels.pdf",
       plot = p2, width = 14, height = 10)

cat("\nCheck mbrt4h_4peak_zones_labels.pdf against H2AX image.\n")
