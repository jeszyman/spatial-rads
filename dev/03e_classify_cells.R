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

# --- From the H2AX image (image7.png), count ~5 bright horizontal bands ---
# H2AX image stripes are roughly evenly spaced across the tissue
# CosMx y-range: -0.18 to 3.55 mm (3.73mm span)
# ~5 stripes suggests center-to-center spacing ~0.75mm
# The tissue edge confounds — edge cells have higher p21 regardless of beam

# --- Strategy: use the FOV-level banding the user identified ---
# Two dark bands visible at y~0.6 and y~1.2 in per-FOV p21 heatmap
# These are valleys between beams

# Let's build a stripe model from the H2AX image:
# From image7.png, counting bright bands from bottom to top of tissue:
# Beam peaks roughly at y = 0.1, 0.85, 1.6, 2.35, 3.1 (5 beams, ~0.75mm spacing)
# Valleys between beams at: y = 0.47, 1.22, 1.97, 2.72

# But let's verify this by looking at the ACTUAL p21 data with multiple candidate spacings
spacings <- seq(0.5, 1.2, by = 0.05)
offsets <- seq(-0.2, 0.5, by = 0.05)

# For each spacing+offset, compute the contrast between peak and valley cells
results <- expand_grid(spacing = spacings, offset = offsets) %>%
  mutate(contrast = map2_dbl(spacing, offset, function(s, o) {
    # Define beam centers at offset, offset+s, offset+2s, ...
    beam_centers <- seq(o, max(mbrt4h$y_slide_mm) + s, by = s)
    # Distance to nearest beam center
    mbrt4h_dist <- mbrt4h %>%
      mutate(dist_to_peak = {
        d <- outer(y_slide_mm, beam_centers, function(y, c) abs(y - c))
        apply(d, 1, min)
      })
    # Classify: peak if within 1/4 of spacing, valley if in outer 1/4
    quarter <- s / 4
    peaks <- mbrt4h_dist %>% filter(dist_to_peak < quarter)
    valleys <- mbrt4h_dist %>% filter(dist_to_peak > quarter * 2)

    if (nrow(peaks) < 100 || nrow(valleys) < 100) return(0)
    mean(peaks$p21, na.rm = TRUE) - mean(valleys$p21, na.rm = TRUE)
  }))

best <- results %>% arrange(desc(abs(contrast))) %>% head(10)
cat("=== Top 10 spacing/offset combos by p21 peak-valley contrast ===\n")
print(best)

best_spacing <- best$spacing[1]
best_offset <- best$offset[1]
cat(sprintf("\nBest fit: spacing=%.2fmm, offset=%.2fmm, contrast=%.4f\n",
            best_spacing, best_offset, best$contrast[1]))

# --- Apply best-fit stripe model ---
beam_centers <- seq(best_offset, max(mbrt4h$y_slide_mm) + best_spacing, by = best_spacing)
cat(sprintf("\nBeam centers: %s\n", paste(round(beam_centers, 2), collapse = ", ")))

mbrt4h <- mbrt4h %>%
  mutate(dist_to_peak = {
    d <- outer(y_slide_mm, beam_centers, function(y, c) abs(y - c))
    apply(d, 1, min)
  })

quarter <- best_spacing / 4
mbrt4h <- mbrt4h %>%
  mutate(zone = case_when(
    dist_to_peak < quarter ~ "peak",
    dist_to_peak > quarter * 2 ~ "valley",
    TRUE ~ "boundary"
  ))

cat("\n=== Cell classification ===\n")
mbrt4h %>% count(zone) %>% mutate(pct = n / sum(n) * 100) %>% print()

cat("\n=== p21 by zone ===\n")
mbrt4h %>%
  group_by(zone) %>%
  summarise(mean_p21 = mean(p21, na.rm = TRUE),
            mean_ddr = mean(DNA_Damage_Repair, na.rm = TRUE),
            mean_sting = mean(STING, na.rm = TRUE),
            mean_ifn1 = mean(TypeI_interferon, na.rm = TRUE),
            n = n(), .groups = "drop") %>%
  print()

# --- Spatial plot with stripe zones ---
p1 <- ggplot(mbrt4h, aes(x = x_slide_mm, y = y_slide_mm, color = zone)) +
  geom_point(size = 0.02, alpha = 0.3) +
  scale_color_manual(values = c("peak" = "#E74C3C", "valley" = "#3498DB", "boundary" = "gray70")) +
  geom_hline(yintercept = beam_centers, color = "red", linetype = "dashed", alpha = 0.5) +
  coord_fixed() +
  labs(title = sprintf("MBRT_4h: Peak/Valley classification (spacing=%.2fmm)", best_spacing),
       color = "Zone") +
  theme_bw()
ggsave("/mnt/data/projects/spatial-rads/analysis/figures/mbrt4h_peak_valley_zones.pdf",
       plot = p1, width = 12, height = 8)

# --- p21 spatial with zone overlay ---
p2 <- ggplot(mbrt4h, aes(x = x_slide_mm, y = y_slide_mm)) +
  stat_summary_2d(aes(z = p21), fun = mean, binwidth = c(0.1, 0.1)) +
  scale_fill_viridis_c(option = "inferno", name = "Mean p21") +
  geom_hline(yintercept = beam_centers, color = "cyan", linetype = "dashed",
             linewidth = 0.5, alpha = 0.7) +
  coord_fixed() +
  labs(title = "p21 heatmap with beam centers (dashed lines)") +
  theme_bw()
ggsave("/mnt/data/projects/spatial-rads/analysis/figures/mbrt4h_p21_with_beams.pdf",
       plot = p2, width = 12, height = 8)

# --- Pathway comparison: peak vs valley ---
library(patchwork)

pathway_by_zone <- mbrt4h %>%
  pivot_longer(c(TypeI_interferon, TypeII_interferon, STING, DNA_Damage_Repair, p21),
               names_to = "pathway", values_to = "score") %>%
  group_by(zone, pathway) %>%
  summarise(mean_score = mean(score, na.rm = TRUE),
            se = sd(score, na.rm = TRUE) / sqrt(n()),
            n = n(), .groups = "drop")

p3 <- ggplot(pathway_by_zone %>% filter(zone != "boundary"),
             aes(x = pathway, y = mean_score, fill = zone)) +
  geom_col(position = "dodge") +
  geom_errorbar(aes(ymin = mean_score - se, ymax = mean_score + se),
                position = position_dodge(0.9), width = 0.2) +
  scale_fill_manual(values = c("peak" = "#E74C3C", "valley" = "#3498DB")) +
  labs(title = "Pathway scores: Peak vs Valley", y = "Mean score", x = NULL) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/mnt/data/projects/spatial-rads/analysis/figures/mbrt4h_peak_valley_pathways.pdf",
       plot = p3, width = 10, height = 6)

# --- Cell type composition by zone ---
comp <- mbrt4h %>%
  count(zone, cell_type_validated) %>%
  group_by(zone) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

p4 <- ggplot(comp %>% filter(zone != "boundary"),
             aes(x = zone, y = prop, fill = cell_type_validated)) +
  geom_col() +
  labs(title = "Cell type composition: Peak vs Valley", fill = "Cell type",
       y = "Proportion") +
  theme_bw()
ggsave("/mnt/data/projects/spatial-rads/analysis/figures/mbrt4h_peak_valley_composition.pdf",
       plot = p4, width = 10, height = 6)

# --- Save classification ---
write.csv(mbrt4h %>% select(cell_id, fov, x_slide_mm, y_slide_mm, zone, dist_to_peak,
                             cell_type_validated, p21, DNA_Damage_Repair, STING,
                             TypeI_interferon, TypeII_interferon),
          "/mnt/data/projects/spatial-rads/analysis/tables/mbrt4h_peak_valley.csv",
          row.names = FALSE)
cat("\nClassification saved to mbrt4h_peak_valley.csv\n")

# --- Print stripe geometry summary ---
cat(sprintf("\n=== STRIPE MODEL SUMMARY ===\n"))
cat(sprintf("Orientation: horizontal (y-axis variation)\n"))
cat(sprintf("Beam spacing: %.2f mm center-to-center\n", best_spacing))
cat(sprintf("Beam centers: %s mm\n", paste(round(beam_centers, 2), collapse = ", ")))
cat(sprintf("Peak zone: within %.3f mm of beam center\n", quarter))
cat(sprintf("Valley zone: > %.3f mm from beam center\n", quarter * 2))
cat(sprintf("Boundary: %.3f - %.3f mm from beam center\n", quarter, quarter * 2))
