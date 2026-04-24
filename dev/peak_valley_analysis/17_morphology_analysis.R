if (!exists("PATHWAY_COLS")) source("00_load_data.R")

library(arrow)

# --- Load morphology data from parquet ---
meta <- read_parquet(
  file.path(INPUT_DIR, "Analysis_Mutter_01_CosMmR_Mutter_updated_metadata.parquet"),
  col_select = c("cell_id", "fov", "Slide", "Condition", "x_slide_mm", "y_slide_mm",
                  "Area.um2", "NucArea", "Circularity", "Solidity", "AspectRatio",
                  "qcFlagsCell", "ImmuneAtlas_ImmGen_Main_cell_Types")
) %>%
  as.data.frame() %>%
  filter(qcFlagsCell == "Pass") %>%
  mutate(cell_type = ifelse(ImmuneAtlas_ImmGen_Main_cell_Types == "a",
                            "tumor_epithelial",
                            ImmuneAtlas_ImmGen_Main_cell_Types))

# --- Peak/valley classification via FOV extrapolation ---
BLOCK <- "Block_21"; SLIDE_ID <- "20250529_214712_S4"
peak_fovs   <- c(224, 209, 201, 186, 172, 225, 229, 215, 175, 176, 206, 181, 168, 167)
valley_fovs <- c(199, 213, 227, 188, 204, 218, 178, 165, 170, 231, 189, 173, 154)

mbrt4h_all <- meta %>% filter(Condition == "MBRT_4h")

peak_cells <- mbrt4h_all %>% filter(fov %in% peak_fovs)
peak_fov_centroids <- peak_cells %>%
  group_by(fov) %>%
  summarise(cx = mean(x_slide_mm), cy = mean(y_slide_mm), .groups = "drop")

search_theta <- function(cx, cy, n_stripes = 4,
                         theta_grid = seq(-pi/2, pi/2, length.out = 181)) {
  best <- list(theta = NA, within_var = Inf)
  for (th in theta_grid) {
    d <- -sin(th) * cx + cos(th) * cy
    km <- kmeans(d, centers = n_stripes, nstart = 10)
    if (km$tot.withinss < best$within_var)
      best <- list(theta = th, within_var = km$tot.withinss,
                   centers = sort(km$centers[, 1]), d = d, cluster = km$cluster)
  }
  best
}

set.seed(1)
fit <- search_theta(peak_fov_centroids$cx, peak_fov_centroids$cy)
theta <- fit$theta
stripe_centers <- sort(fit$centers)
beam_spacing <- median(diff(stripe_centers))

peak_fov_centroids$d_perp <- -sin(theta) * peak_fov_centroids$cx + cos(theta) * peak_fov_centroids$cy
peak_fov_centroids$stripe <- fit$cluster
within_stripe_sd <- peak_fov_centroids %>%
  mutate(resid = d_perp - stripe_centers[stripe]) %>%
  group_by(stripe) %>% summarise(sd = sd(resid), .groups = "drop") %>%
  pull(sd) %>% mean()
peak_half <- within_stripe_sd + 0.15

d_all <- -sin(theta) * mbrt4h_all$x_slide_mm + cos(theta) * mbrt4h_all$y_slide_mm
d_range <- range(d_all)
all_centers <- seq(stripe_centers[1] - ceiling((stripe_centers[1] - d_range[1]) / beam_spacing) * beam_spacing,
                   stripe_centers[length(stripe_centers)] + ceiling((d_range[2] - stripe_centers[length(stripe_centers)]) / beam_spacing) * beam_spacing,
                   by = beam_spacing)

dist_to_nearest <- sapply(d_all, function(d) min(abs(d - all_centers)))
zone <- rep("transition", length(d_all))
zone[dist_to_nearest < peak_half] <- "peak"
zone[dist_to_nearest > (beam_spacing / 2 - peak_half)] <- "valley"

mbrt4h_all$zone <- zone
mbrt4h_all$d_to_peak <- dist_to_nearest

cat(sprintf("FOV extrapolation: theta=%.1f deg, spacing=%.3f mm, peak_half=%.3f mm\n",
            theta * 180 / pi, beam_spacing, peak_half))
cat("Zone counts:\n"); print(table(zone))

mbrt4h <- mbrt4h_all %>%
  filter(zone %in% c("peak", "valley")) %>%
  mutate(group = paste0("MBRT_", zone))

sbrt4h <- meta %>%
  filter(Condition == "SBRT_4h") %>%
  mutate(group = "SBRT_4h")

control <- meta %>%
  filter(Condition == "Control") %>%
  mutate(group = "Control")

morph <- bind_rows(control, sbrt4h, mbrt4h) %>%
  mutate(group = factor(group,
    levels = c("Control", "SBRT_4h", "MBRT_valley", "MBRT_peak")))

GROUP_COLORS <- c("Control" = "#7F8C8D", "SBRT_4h" = "#3498DB",
                  "MBRT_valley" = "#85C1E9", "MBRT_peak" = "#E74C3C")

cat(sprintf("Morphology dataset: %d cells across %d groups\n", nrow(morph),
            n_distinct(morph$group)))
morph %>% count(group) %>% print()

# --- Violin plots: 5 metrics x 4 groups ---
MORPH_METRICS <- c("Area.um2", "NucArea", "Circularity", "Solidity", "AspectRatio")
MORPH_LABELS <- c("Area.um2" = "Cell Area (µm²)", "NucArea" = "Nuclear Area (px)",
                  "Circularity" = "Circularity", "Solidity" = "Solidity",
                  "AspectRatio" = "Aspect Ratio")

morph_long <- morph %>%
  pivot_longer(cols = all_of(MORPH_METRICS), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = MORPH_METRICS)) %>%
  group_by(metric) %>%
  filter(value <= quantile(value, 0.99, na.rm = TRUE)) %>%
  ungroup()

p_violin <- ggplot(morph_long, aes(x = group, y = value, fill = group)) +
  geom_violin(scale = "width", alpha = 0.7, linewidth = 0.3) +
  geom_boxplot(width = 0.1, outlier.size = 0.2, alpha = 0.8) +
  facet_wrap(~metric, scales = "free_y",
             labeller = labeller(metric = MORPH_LABELS)) +
  scale_fill_manual(values = GROUP_COLORS) +
  labs(x = NULL, y = NULL, title = "Cell Morphology: MBRT Peak/Valley vs SBRT vs Control (4h)") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "none")

ggsave(file.path(PLOT_DIR, "morphology_violins.png"), plot = p_violin,
       width = 14, height = 10, dpi = 200)
cat("Saved: plots/morphology_violins.png\n")

# --- Cell-type-stratified violin for Area.um2 ---
top_ct <- morph %>%
  count(cell_type, sort = TRUE) %>%
  slice_head(n = 5) %>%
  pull(cell_type)

morph_ct <- morph %>% filter(cell_type %in% top_ct)

p_ct_area <- ggplot(morph_ct, aes(x = group, y = Area.um2, fill = group)) +
  geom_violin(scale = "width", alpha = 0.7, linewidth = 0.3) +
  geom_boxplot(width = 0.1, outlier.size = 0.2, alpha = 0.8) +
  facet_wrap(~cell_type, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = GROUP_COLORS) +
  labs(x = NULL, y = "Cell Area (µm²)",
       title = "Cell Area by Cell Type: MBRT Peak/Valley vs SBRT vs Control (4h)") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "none")

ggsave(file.path(PLOT_DIR, "morphology_area_by_celltype.png"), plot = p_ct_area,
       width = 16, height = 5, dpi = 200)
cat("Saved: plots/morphology_area_by_celltype.png\n")

# --- Cell type composition table ---
comp <- morph %>%
  count(group, cell_type) %>%
  group_by(group) %>%
  mutate(pct = round(n / sum(n) * 100, 2)) %>%
  ungroup() %>%
  arrange(group, desc(n))

write_tsv(comp, file.path(DATA_DIR, "morphology_celltype_composition.tsv"))
cat("Saved: data/morphology_celltype_composition.tsv\n")

cat("\n=== Cell type composition (%) ===\n")
comp %>%
  select(group, cell_type, pct) %>%
  pivot_wider(names_from = group, values_from = pct, values_fill = 0) %>%
  arrange(desc(Control)) %>%
  print(n = 20)

# --- Spatial scatter: MBRT_4h colored by Area.um2 ---
p_spatial <- ggplot(mbrt4h_all, aes(x = x_slide_mm, y = y_slide_mm, color = Area.um2)) +
  geom_point(size = 0.05, alpha = 0.4) +
  scale_color_viridis_c(option = "inferno", limits = c(0, quantile(mbrt4h_all$Area.um2, 0.99))) +
  coord_fixed() +
  labs(title = "MBRT 4h: Cell Area Spatial Distribution (FOV extrapolation)",
       color = "Area (µm²)") +
  theme_bw(base_size = 14)

for (ac in all_centers) {
  p_spatial <- p_spatial +
    geom_abline(intercept = ac / cos(theta), slope = sin(theta) / cos(theta),
                color = "white", linewidth = 0.5, alpha = 0.7, linetype = "dashed")
}

ggsave(file.path(PLOT_DIR, "morphology_spatial_area.png"), plot = p_spatial,
       width = 14, height = 10, dpi = 200)
cat("Saved: plots/morphology_spatial_area.png\n")

# --- Summary statistics ---
summary_stats <- morph %>%
  pivot_longer(cols = all_of(MORPH_METRICS), names_to = "metric", values_to = "value") %>%
  group_by(group, metric) %>%
  summarise(
    n = n(),
    mean = round(mean(value, na.rm = TRUE), 3),
    median = round(median(value, na.rm = TRUE), 3),
    sd = round(sd(value, na.rm = TRUE), 3),
    q25 = round(quantile(value, 0.25, na.rm = TRUE), 3),
    q75 = round(quantile(value, 0.75, na.rm = TRUE), 3),
    .groups = "drop"
  )

write_tsv(summary_stats, file.path(DATA_DIR, "morphology_summary.tsv"))
cat("Saved: data/morphology_summary.tsv\n")

cat("\n=== Morphology Summary ===\n")
summary_stats %>%
  filter(metric == "Area.um2") %>%
  print()

# --- Morphology kinetics: all timepoints ---
kinetics <- meta %>%
  filter(Condition %in% CONDITION_LEVELS) %>%
  mutate(
    treatment = case_when(
      Condition == "Control" ~ "NT",
      grepl("MBRT", Condition) ~ "MBRT",
      grepl("SBRT", Condition) ~ "SBRT"
    ),
    timepoint_h = case_when(
      Condition == "Control" ~ 0,
      grepl("1h", Condition) ~ 1,
      grepl("4h", Condition) ~ 4,
      grepl("day2", Condition) ~ 48,
      grepl("day6", Condition) ~ 144,
      grepl("day8", Condition) ~ 192,
      grepl("day10", Condition) ~ 240
    ),
    model = ifelse(grepl("^Tongue", Condition), "tongue", "flank")
  )

flank <- kinetics %>% filter(model == "flank")

flank_summary <- flank %>%
  group_by(treatment, timepoint_h) %>%
  summarise(
    n = n(),
    median_area = median(Area.um2, na.rm = TRUE),
    mean_area = mean(Area.um2, na.rm = TRUE),
    q25 = quantile(Area.um2, 0.25, na.rm = TRUE),
    q75 = quantile(Area.um2, 0.75, na.rm = TRUE),
    median_nuc = median(NucArea, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n=== Area Kinetics (flank) ===\n")
flank_summary %>% arrange(timepoint_h, treatment) %>% print(n = 30)

p_kinetics <- ggplot(flank_summary, aes(x = timepoint_h, y = median_area,
                                        color = treatment)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_ribbon(aes(ymin = q25, ymax = q75, fill = treatment),
              alpha = 0.15, color = NA) +
  scale_color_manual(values = TREATMENT_COLORS) +
  scale_fill_manual(values = TREATMENT_COLORS) +
  scale_x_continuous(breaks = c(0, 1, 4, 48, 144),
                     labels = c("0", "1h", "4h", "day2", "day6")) +
  labs(x = "Time post-irradiation", y = "Median Cell Area (µm²)",
       title = "Cell Area Kinetics: MBRT vs SBRT (flank)",
       subtitle = "Ribbon = IQR",
       color = "Treatment", fill = "Treatment") +
  theme_bw(base_size = 14)

ggsave(file.path(PLOT_DIR, "morphology_kinetics.png"), plot = p_kinetics,
       width = 10, height = 6, dpi = 200)
cat("Saved: plots/morphology_kinetics.png\n")

p_nuc_kinetics <- ggplot(flank_summary, aes(x = timepoint_h, y = median_nuc,
                                            color = treatment)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = TREATMENT_COLORS) +
  scale_x_continuous(breaks = c(0, 1, 4, 48, 144),
                     labels = c("0", "1h", "4h", "day2", "day6")) +
  labs(x = "Time post-irradiation", y = "Median Nuclear Area (px)",
       title = "Nuclear Area Kinetics: MBRT vs SBRT (flank)",
       color = "Treatment") +
  theme_bw(base_size = 14)

ggsave(file.path(PLOT_DIR, "morphology_nuc_kinetics.png"), plot = p_nuc_kinetics,
       width = 10, height = 6, dpi = 200)
cat("Saved: plots/morphology_nuc_kinetics.png\n")
