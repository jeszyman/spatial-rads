library(Seurat)
library(tidyverse)

obj <- readRDS("/mnt/data/projects/spatial-rads/analysis/objects/seurat_clustered.rds")

mbrt4h <- obj@meta.data %>%
  as.data.frame() %>%
  filter(Condition == "MBRT_4h") %>%
  select(cell_id, fov, x_slide_mm, y_slide_mm, cell_type_validated,
         TypeI_interferon, TypeII_interferon, STING, DNA_Damage_Repair)

p21_expr <- GetAssayData(obj, layer = "data")["Cdkn1a", rownames(obj@meta.data[obj@meta.data$Condition == "MBRT_4h", ])]
mbrt4h$p21 <- as.numeric(p21_expr[rownames(mbrt4h)])

# --- Image-based FOV row assignment ---
# From H2AX image (image7.png): ~5 brown horizontal bands (peaks)
# separated by ~4 lighter gaps (valleys), top to bottom of tissue.
# FOV grid (image6.png) overlays on same tissue.
#
# Assignment (top of tissue = high y, bottom = low y):
#   Row 9 (y~3.2): FOVs 154-160         -> P  (top bright band)
#   Row 8 (y~2.8): FOVs 162-171         -> V  (gap)
#   Row 7 (y~2.3): FOVs 161,173-182     -> P  (bright band)
#   Row 6 (y~1.9): FOVs 184,172,187,189,183 -> V (gap)
#   Row 5 (y~1.7): FOVs 185-196         -> P  (central bright band)
#   Row 4 (y~1.3): FOVs 197-208         -> V  (gap)
#   Row 3 (y~0.8): FOVs 209-220         -> P  (bright band)
#   Row 2 (y~0.3): FOVs 221-231         -> V  (gap)
#   Row 1 (y~0.0): FOVs 232,233         -> P  (bottom edge, only 2 FOVs)

peak_fovs <- c(154:160,                           # Row 9
               161, 173:182,                       # Row 7
               185:196,                            # Row 5
               209:220,                            # Row 3
               232, 233)                           # Row 1

valley_fovs <- c(162:171,                          # Row 8
                 184, 172, 187, 189, 183,          # Row 6
                 197:208,                           # Row 4
                 221:231)                           # Row 2

# Verify all MBRT_4h FOVs are accounted for
all_fovs <- sort(unique(mbrt4h$fov))
assigned <- sort(unique(c(peak_fovs, valley_fovs)))
missing <- setdiff(all_fovs, assigned)
extra <- setdiff(assigned, all_fovs)
cat(sprintf("Total FOVs in data: %d\n", length(all_fovs)))
cat(sprintf("Assigned FOVs: %d\n", length(intersect(assigned, all_fovs))))
if (length(missing) > 0) cat(sprintf("MISSING FOVs: %s\n", paste(missing, collapse = ", ")))
if (length(extra) > 0) cat(sprintf("Extra FOVs (not in data): %s\n", paste(extra, collapse = ", ")))

# Classify cells
mbrt4h <- mbrt4h %>%
  mutate(zone = case_when(
    fov %in% peak_fovs ~ "peak",
    fov %in% valley_fovs ~ "valley",
    TRUE ~ "unassigned"
  ))

cat("\n=== Cell counts by zone ===\n")
mbrt4h %>% count(zone) %>% mutate(pct = round(n / sum(n) * 100, 1)) %>% print()

# --- VERIFICATION PLOT: spatial map colored by zone ---
p1 <- ggplot(mbrt4h, aes(x = x_slide_mm, y = y_slide_mm, color = zone)) +
  geom_point(size = 0.02, alpha = 0.4) +
  scale_color_manual(values = c("peak" = "#E74C3C", "valley" = "#3498DB")) +
  coord_fixed() +
  labs(title = "MBRT_4h: Image-based Peak/Valley FOV classification",
       subtitle = "Red = peak (H2AX+ band), Blue = valley (gap between beams)",
       color = "Zone") +
  theme_bw(base_size = 14)
ggsave("/mnt/data/projects/spatial-rads/analysis/figures/mbrt4h_image_based_zones.pdf",
       plot = p1, width = 14, height = 10)

# --- Same plot with FOV labels ---
fov_centroids <- mbrt4h %>%
  group_by(fov, zone) %>%
  summarise(x = mean(x_slide_mm), y = mean(y_slide_mm), .groups = "drop")

p2 <- ggplot(mbrt4h, aes(x = x_slide_mm, y = y_slide_mm, color = zone)) +
  geom_point(size = 0.02, alpha = 0.2) +
  scale_color_manual(values = c("peak" = "#E74C3C", "valley" = "#3498DB")) +
  geom_text(data = fov_centroids, aes(x = x, y = y, label = fov),
            color = "black", size = 2, fontface = "bold") +
  coord_fixed() +
  labs(title = "MBRT_4h: Peak/Valley zones with FOV labels",
       subtitle = "CHECK: do red/blue bands match H2AX image stripes?",
       color = "Zone") +
  theme_bw(base_size = 14)
ggsave("/mnt/data/projects/spatial-rads/analysis/figures/mbrt4h_zones_with_fov_labels.pdf",
       plot = p2, width = 14, height = 10)

cat("\nVerification plots saved. Compare mbrt4h_zones_with_fov_labels.pdf with H2AX image.\n")

# --- Quick pathway summary ---
cat("\n=== Pathway scores by zone ===\n")
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
