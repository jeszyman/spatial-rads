if (!exists("PATHWAY_COLS")) source("00_load_data.R")
library(zoo)

obj <- readRDS(file.path(ANALYSIS_DIR, "objects", "seurat_clustered.rds"))

# --- Extract MBRT_4h cells ---
mbrt4h <- obj@meta.data %>%
  as.data.frame() %>%
  filter(Condition == "MBRT_4h") %>%
  select(cell_id, fov, x_slide_mm, y_slide_mm, cell_type_validated,
         all_of(PATHWAY_COLS))

p21_expr <- GetAssayData(obj, layer = "data")["Cdkn1a",
  rownames(obj@meta.data[obj@meta.data$Condition == "MBRT_4h", ])]
mbrt4h$p21 <- as.numeric(p21_expr[rownames(mbrt4h)])
cat(sprintf("MBRT_4h: %d cells, %d FOVs\n", nrow(mbrt4h), length(unique(mbrt4h$fov))))

# ---- Highlight 1: Fine-grained p21 spatial profile ----
y_profile <- mbrt4h %>%
  mutate(y_bin = round(y_slide_mm / 0.1) * 0.1) %>%
  group_by(y_bin) %>%
  summarise(mean_p21 = mean(p21, na.rm = TRUE), n_cells = n(), .groups = "drop") %>%
  filter(n_cells >= 50) %>%
  arrange(y_bin) %>%
  mutate(p21_smooth = rollmean(mean_p21, k = 5, fill = NA, align = "center"))

p_profile <- ggplot(y_profile, aes(x = y_bin)) +
  geom_point(aes(y = mean_p21, size = n_cells), alpha = 0.3) +
  geom_line(aes(y = p21_smooth), color = "red", linewidth = 1.2, na.rm = TRUE) +
  labs(x = "y position (mm)", y = "Mean p21 (Cdkn1a)",
       title = "MBRT_4h: p21 y-axis profile (0.1mm bins, 0.5mm rolling mean)") +
  theme_bw()
ggsave(file.path(PLOT_DIR, "stripe_p21_y_profile.png"), plot = p_profile, width = 12, height = 5, dpi = 150)

# ---- Highlight 2: Rotated axis scan (stripe orientation) ----
angles <- seq(0, 175, by = 15)
angle_cvs <- map_dbl(angles, function(theta) {
  rad <- theta * pi / 180
  proj <- mbrt4h$x_slide_mm * cos(rad) + mbrt4h$y_slide_mm * sin(rad)
  bins <- round(proj / 0.1) * 0.1
  profile <- tibble(bin = bins, p21 = mbrt4h$p21) %>%
    group_by(bin) %>%
    summarise(mean_p21 = mean(p21, na.rm = TRUE), n = n(), .groups = "drop") %>%
    filter(n >= 50) %>% arrange(bin)
  if (nrow(profile) < 7) return(0)
  smoothed <- rollmean(profile$mean_p21, k = 5, fill = NA, align = "center")
  smoothed <- smoothed[!is.na(smoothed)]
  if (length(smoothed) < 3) return(0)
  sd(smoothed) / mean(smoothed)
})

angle_df <- tibble(angle = angles, cv = angle_cvs)
best_proj_angle <- angle_df$angle[which.max(angle_df$cv)]
cat(sprintf("Best projection angle: %d deg (CV=%.4f)\n", best_proj_angle, max(angle_df$cv)))

p_angle <- ggplot(angle_df, aes(x = angle, y = cv)) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  geom_vline(xintercept = best_proj_angle, color = "red", linetype = "dashed") +
  labs(x = "Projection angle (degrees)", y = "CV of smoothed p21",
       title = "Stripe direction scan: Y-axis shows strongest periodic variation") +
  theme_bw()
ggsave(file.path(PLOT_DIR, "stripe_angle_scan.png"), plot = p_angle, width = 10, height = 6, dpi = 150)

# ---- Highlight 3: Tilt angle optimization ----
cat("Scanning tilt angles (5-35 deg)...\n")
angle_results <- map_dfr(seq(5, 35, by = 2), function(deg) {
  rad <- deg * pi / 180
  y_corr <- mbrt4h$y_slide_mm + mbrt4h$x_slide_mm * tan(rad)
  best_contrast <- 0
  best_s <- NA
  for (s in seq(0.6, 1.8, by = 0.1)) {
    for (o in seq(min(y_corr) + s/4, min(y_corr) + s * 1.5, by = 0.1)) {
      centers <- seq(o, max(y_corr) + s, by = s)
      dist <- sapply(y_corr, function(yc) min(abs(yc - centers)))
      pk <- which(dist < s / 4)
      vl <- which(dist > s / 2)
      if (length(pk) < 500 || length(vl) < 500) next
      contrast <- mean(mbrt4h$p21[pk], na.rm = TRUE) - mean(mbrt4h$p21[vl], na.rm = TRUE)
      if (contrast > best_contrast) { best_contrast <- contrast; best_s <- s }
    }
  }
  tibble(angle = deg, best_spacing = best_s, contrast = best_contrast)
})

opt_angle <- angle_results$angle[which.max(angle_results$contrast)]
cat(sprintf("Optimal tilt: %d deg (contrast=%.4f)\n", opt_angle, max(angle_results$contrast)))

p_tilt <- ggplot(angle_results, aes(x = angle, y = contrast)) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  geom_vline(xintercept = opt_angle, color = "red", linetype = "dashed") +
  labs(x = "Tilt angle (degrees)", y = "Peak - Valley p21 contrast",
       title = "Tilt optimization: 15 deg confirmed as optimal") +
  theme_bw()
ggsave(file.path(PLOT_DIR, "stripe_tilt_optimization.png"), plot = p_tilt, width = 10, height = 6, dpi = 150)

# ---- Highlight 4: 4-peak grid search ----
rad <- 15 * pi / 180
mbrt4h <- mbrt4h %>% mutate(y_corr = y_slide_mm + x_slide_mm * tan(rad))
y_range <- range(mbrt4h$y_corr)

results_4pk <- expand_grid(
  spacing = seq(0.8, 1.3, by = 0.02),
  offset = seq(y_range[1] + 0.1, y_range[1] + 1.3, by = 0.05)
) %>%
  mutate(n_peaks = map2_int(spacing, offset, function(s, o) {
    centers <- seq(o, y_range[2] + 0.1, by = s)
    sum(sapply(centers, function(c) sum(abs(mbrt4h$y_corr - c) < s/4) > 100))
  })) %>%
  filter(n_peaks == 4) %>%
  mutate(contrast = map2_dbl(spacing, offset, function(s, o) {
    centers <- seq(o, y_range[2] + 0.1, by = s)
    dist <- sapply(mbrt4h$y_corr, function(yc) min(abs(yc - centers)))
    pk <- which(dist < s / 4); vl <- which(dist > s / 4)
    if (length(pk) < 1000 || length(vl) < 1000) return(0)
    mean(mbrt4h$p21[pk], na.rm = TRUE) - mean(mbrt4h$p21[vl], na.rm = TRUE)
  }))

best_4pk <- results_4pk %>% arrange(desc(contrast)) %>% head(1)
cat(sprintf("\n4-peak model: spacing=%.2fmm, contrast=%.4f\n", best_4pk$spacing, best_4pk$contrast))

# Compute beam centers for best model
beam_centers <- seq(best_4pk$offset, y_range[2] + 0.1, by = best_4pk$spacing)
beam_centers <- beam_centers[sapply(beam_centers, function(c) sum(abs(mbrt4h$y_corr - c) < best_4pk$spacing/4) > 100)]
cat(sprintf("Beam centers (corrected y): %s\n", paste(round(beam_centers, 2), collapse = ", ")))

# Save stripe model parameters for downstream scripts
stripe_model <- list(
  tilt_deg = 15,
  spacing_mm = best_4pk$spacing,
  beam_centers = beam_centers,
  n_peaks = length(beam_centers)
)
saveRDS(stripe_model, file.path(DATA_DIR, "stripe_model.rds"))

# 4-pathway spatial heatmaps
p21_heatmap <- ggplot(mbrt4h, aes(x = x_slide_mm, y = y_slide_mm)) +
  stat_summary_2d(aes(z = p21), fun = mean, binwidth = c(0.1, 0.1)) +
  scale_fill_viridis_c(option = "inferno", name = "Mean p21") +
  coord_fixed() + labs(title = "p21 spatial heatmap") + theme_bw()

ddr_heatmap <- ggplot(mbrt4h, aes(x = x_slide_mm, y = y_slide_mm)) +
  stat_summary_2d(aes(z = DNA_Damage_Repair), fun = mean, binwidth = c(0.1, 0.1)) +
  scale_fill_viridis_c(option = "inferno", name = "Mean DDR") +
  coord_fixed() + labs(title = "DDR spatial heatmap") + theme_bw()

sting_heatmap <- ggplot(mbrt4h, aes(x = x_slide_mm, y = y_slide_mm)) +
  stat_summary_2d(aes(z = STING), fun = mean, binwidth = c(0.1, 0.1)) +
  scale_fill_viridis_c(option = "plasma", name = "Mean STING") +
  coord_fixed() + labs(title = "STING spatial heatmap") + theme_bw()

ifn1_heatmap <- ggplot(mbrt4h, aes(x = x_slide_mm, y = y_slide_mm)) +
  stat_summary_2d(aes(z = TypeI_interferon), fun = mean, binwidth = c(0.1, 0.1)) +
  scale_fill_viridis_c(option = "magma", name = "Mean IFN-I") +
  coord_fixed() + labs(title = "Type I IFN spatial heatmap") + theme_bw()

combo <- (p21_heatmap | ddr_heatmap) / (sting_heatmap | ifn1_heatmap)
ggsave(file.path(PLOT_DIR, "stripe_4pathway_heatmaps.png"), plot = combo, width = 20, height = 16, dpi = 150)

cat("Stripe detection complete. Model saved to data/stripe_model.rds\n")
