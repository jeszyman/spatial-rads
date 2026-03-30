if (!exists("PATHWAY_COLS")) source("00_load_data.R")

obj <- readRDS(file.path(ANALYSIS_DIR, "objects", "seurat_clustered.rds"))
stripe_model <- readRDS(file.path(DATA_DIR, "stripe_model.rds"))

mbrt4h <- obj@meta.data %>%
  as.data.frame() %>%
  filter(Condition == "MBRT_4h") %>%
  select(cell_id, fov, x_slide_mm, y_slide_mm, cell_type_validated,
         all_of(PATHWAY_COLS))

p21_expr <- GetAssayData(obj, layer = "data")["Cdkn1a",
  rownames(obj@meta.data[obj@meta.data$Condition == "MBRT_4h", ])]
mbrt4h$p21 <- as.numeric(p21_expr[rownames(mbrt4h)])

# --- Apply stripe model ---
rad <- stripe_model$tilt_deg * pi / 180
mbrt4h <- mbrt4h %>%
  mutate(
    y_corr = y_slide_mm + x_slide_mm * tan(rad),
    dist_to_peak = sapply(y_corr, function(yc) min(abs(yc - stripe_model$beam_centers))),
    zone = ifelse(dist_to_peak < stripe_model$spacing_mm / 4, "peak", "valley")
  )

cat("=== Peak/Valley Classification ===\n")
cat(sprintf("Model: %d-deg tilt, %.2fmm spacing, %d peaks\n",
            stripe_model$tilt_deg, stripe_model$spacing_mm, stripe_model$n_peaks))
cat(sprintf("Beam centers: %s\n", paste(round(stripe_model$beam_centers, 2), collapse = ", ")))

zone_counts <- mbrt4h %>% count(zone) %>% mutate(pct = round(n / sum(n) * 100, 1))
print(zone_counts)

cat("\n=== p21 validation ===\n")
mbrt4h %>%
  group_by(zone) %>%
  summarise(n = n(), mean_p21 = mean(p21, na.rm = TRUE),
            mean_DDR = mean(DNA_Damage_Repair, na.rm = TRUE), .groups = "drop") %>%
  print()

# --- Save classification ---
pv_class <- mbrt4h %>%
  select(cell_id, fov, x_slide_mm, y_slide_mm, y_corr, dist_to_peak, zone,
         cell_type_validated, p21, all_of(PATHWAY_COLS))
write_tsv(pv_class, file.path(DATA_DIR, "mbrt4h_peak_valley.tsv"))
cat(sprintf("Classification saved: %d cells\n", nrow(pv_class)))

# --- Spatial zone map ---
fov_centroids <- mbrt4h %>%
  group_by(fov) %>%
  summarise(x = mean(x_slide_mm), y = mean(y_slide_mm), .groups = "drop")

p <- ggplot(mbrt4h, aes(x = x_slide_mm, y = y_slide_mm, color = zone)) +
  geom_point(size = 0.02, alpha = 0.3) +
  scale_color_manual(values = ZONE_COLORS) +
  geom_text(data = fov_centroids, aes(x = x, y = y, label = fov),
            color = "black", size = 2, fontface = "bold", inherit.aes = FALSE) +
  coord_fixed() +
  labs(title = sprintf("Peak/Valley: %d-deg tilt, %.2fmm, %d peaks",
                       stripe_model$tilt_deg, stripe_model$spacing_mm, stripe_model$n_peaks),
       color = "Zone") +
  theme_bw(base_size = 14)
for (bc in stripe_model$beam_centers) {
  p <- p + geom_abline(intercept = bc, slope = -tan(rad), color = "yellow", linewidth = 0.6, alpha = 0.8)
}
ggsave(file.path(PLOT_DIR, "peak_valley_zones.png"), plot = p, width = 14, height = 10, dpi = 150)

# --- Add zone to Seurat object for downstream use ---
zone_lookup <- setNames(mbrt4h$zone, mbrt4h$cell_id)
obj$peak_valley_zone <- ifelse(colnames(obj) %in% names(zone_lookup),
                               zone_lookup[colnames(obj)],
                               NA_character_)
saveRDS(obj, file.path(ANALYSIS_DIR, "objects", "seurat_clustered.rds"))
cat("Zone added to Seurat metadata and saved.\n")
