if (!exists("PATHWAY_COLS")) source("00_load_data.R")

library(arrow)
library(MASS)
library(patchwork)

INPUT_DIR <- "/mnt/data/projects/spatial-rads/inputs"

meta <- read_parquet(
  file.path(INPUT_DIR, "Analysis_Mutter_01_CosMmR_Mutter_updated_metadata.parquet"),
  col_select = c("cell_id", "fov", "Condition", "x_slide_mm", "y_slide_mm",
                  "qcFlagsCell", "ImmuneAtlas_ImmGen_Main_cell_Types")
) %>% as.data.frame() %>% filter(qcFlagsCell == "Pass")

expr <- read_parquet(
  file.path(INPUT_DIR, "projects_Mutter_01_CosMmR_app_data_raw_tx_counts_matrix.parquet"),
  col_select = c("cell_id", "Cdkn1a")
) %>% as.data.frame()

mbrt4h <- meta %>% filter(Condition == "MBRT_4h")
mbrt4h <- left_join(mbrt4h, expr, by = "cell_id")

peak_fovs <- c(224, 209, 201, 186, 172, 225, 229, 215, 175, 176, 206, 181, 168, 167)

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

peak_cells <- mbrt4h %>% filter(fov %in% peak_fovs)
peak_centroids <- peak_cells %>% group_by(fov) %>%
  summarise(cx = mean(x_slide_mm), cy = mean(y_slide_mm), .groups = "drop")

set.seed(1)
fit <- search_theta(peak_centroids$cx, peak_centroids$cy)
theta <- fit$theta
stripe_centers <- sort(fit$centers)
beam_spacing <- median(diff(stripe_centers))

peak_centroids$d_perp <- -sin(theta) * peak_centroids$cx + cos(theta) * peak_centroids$cy
peak_centroids$stripe <- fit$cluster
within_stripe_sd <- peak_centroids %>%
  mutate(resid = d_perp - stripe_centers[stripe]) %>%
  group_by(stripe) %>% summarise(sd = sd(resid), .groups = "drop") %>%
  pull(sd) %>% mean()
peak_half <- within_stripe_sd + 0.15

d_all <- -sin(theta) * mbrt4h$x_slide_mm + cos(theta) * mbrt4h$y_slide_mm
d_range <- range(d_all)
all_centers <- seq(
  stripe_centers[1] - ceiling((stripe_centers[1] - d_range[1]) / beam_spacing) * beam_spacing,
  stripe_centers[length(stripe_centers)] + ceiling((d_range[2] - stripe_centers[length(stripe_centers)]) / beam_spacing) * beam_spacing,
  by = beam_spacing)

dist_to_nearest <- sapply(d_all, function(d) min(abs(d - all_centers)))
zone <- rep("transition", length(d_all))
zone[dist_to_nearest < peak_half] <- "peak"
zone[dist_to_nearest > (beam_spacing / 2 - peak_half)] <- "valley"
mbrt4h$zone <- zone

bcells <- mbrt4h %>%
  filter(ImmuneAtlas_ImmGen_Main_cell_Types %in% c("B.cell", "Memory.B", "Plasmablast"))
bcells$p21_pos <- bcells$Cdkn1a > 0
cat(sprintf("B cells: %d total, %d Cdkn1a+ (%.1f%%)\n",
            nrow(bcells), sum(bcells$p21_pos), 100 * mean(bcells$p21_pos)))

# --- Panel 1: Bold threshold dots ---
p1 <- ggplot() +
  geom_point(data = mbrt4h, aes(x = x_slide_mm, y = y_slide_mm),
             color = "grey20", size = 0.01, alpha = 0.08) +
  geom_point(data = bcells %>% filter(p21_pos) %>% arrange(Cdkn1a),
             aes(x = x_slide_mm, y = y_slide_mm, color = Cdkn1a),
             size = 5, alpha = 0.95) +
  scale_color_gradient(low = "#00BFFF", high = "#FFD700", name = "Cdkn1a") +
  coord_fixed() +
  labs(title = "p21+ B cells (large dots) over all cells",
       subtitle = sprintf("%d Cdkn1a+ B cells / %d total", sum(bcells$p21_pos), nrow(bcells))) +
  theme_void(base_size = 14) +
  theme(plot.background = element_rect(fill = "grey5", color = NA),
        plot.title = element_text(color = "white", face = "bold"),
        plot.subtitle = element_text(color = "grey60"),
        legend.text = element_text(color = "white"),
        legend.title = element_text(color = "white"))

for (ac in all_centers) {
  p1 <- p1 + geom_abline(intercept = ac / cos(theta), slope = sin(theta) / cos(theta),
                          color = "red", linewidth = 0.6, alpha = 0.6, linetype = "dashed")
}

# --- Panel 2: Spatial bin p21+ fraction ---
bin_size <- 0.35
bcells$bx <- round(bcells$x_slide_mm / bin_size) * bin_size
bcells$by <- round(bcells$y_slide_mm / bin_size) * bin_size

hex_ratio <- bcells %>%
  group_by(bx, by) %>%
  summarise(n_total = n(), n_pos = sum(p21_pos), frac = n_pos / n_total, .groups = "drop") %>%
  filter(n_total >= 3)

p2 <- ggplot() +
  geom_tile(data = hex_ratio, aes(x = bx, y = by, fill = frac),
            width = bin_size, height = bin_size, alpha = 0.9) +
  scale_fill_gradient2(low = "#1a1a2e", mid = "#e94560", high = "#FFD700",
                       midpoint = median(hex_ratio$frac),
                       name = "Cdkn1a+ frac") +
  coord_fixed() +
  labs(title = "Fraction of Cdkn1a+ B cells per spatial bin",
       subtitle = sprintf("bin=%.0fµm, min 3 B cells/bin, %d bins", bin_size * 1000, nrow(hex_ratio))) +
  theme_void(base_size = 14) +
  theme(plot.background = element_rect(fill = "grey5", color = NA),
        plot.title = element_text(color = "white", face = "bold"),
        plot.subtitle = element_text(color = "grey60"),
        legend.text = element_text(color = "white"),
        legend.title = element_text(color = "white"))

for (ac in all_centers) {
  p2 <- p2 + geom_abline(intercept = ac / cos(theta), slope = sin(theta) / cos(theta),
                          color = "cyan", linewidth = 0.5, alpha = 0.5, linetype = "dashed")
}

# --- Panel 3: KDE enrichment ---
bw <- 0.4
xlim <- range(bcells$x_slide_mm) + c(-0.5, 0.5)
ylim <- range(bcells$y_slide_mm) + c(-0.5, 0.5)

kde_all <- kde2d(bcells$x_slide_mm, bcells$y_slide_mm, h = bw, n = 100,
                 lims = c(xlim, ylim))
kde_pos <- kde2d(bcells$x_slide_mm[bcells$p21_pos], bcells$y_slide_mm[bcells$p21_pos],
                 h = bw, n = 100, lims = c(xlim, ylim))

ratio_grid <- expand.grid(x = kde_all$x, y = kde_all$y)
ratio_grid$d_all <- as.vector(kde_all$z)
ratio_grid$d_pos <- as.vector(kde_pos$z)
ratio_grid$enrichment <- ratio_grid$d_pos / (ratio_grid$d_all + 1e-10)
ratio_grid <- ratio_grid %>% filter(d_all > quantile(d_all, 0.15))

p3 <- ggplot(ratio_grid, aes(x = x, y = y, fill = enrichment)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_gradient2(low = "#1a1a2e", mid = "#e94560", high = "#FFD700",
                       midpoint = median(ratio_grid$enrichment),
                       name = "p21+ enrichment") +
  coord_fixed() +
  labs(title = "KDE enrichment: p21+ B cell density / all B cell density",
       subtitle = sprintf("bandwidth=%.1fmm", bw)) +
  theme_void(base_size = 14) +
  theme(plot.background = element_rect(fill = "grey5", color = NA),
        plot.title = element_text(color = "white", face = "bold"),
        plot.subtitle = element_text(color = "grey60"),
        legend.text = element_text(color = "white"),
        legend.title = element_text(color = "white"))

for (ac in all_centers) {
  p3 <- p3 + geom_abline(intercept = ac / cos(theta), slope = sin(theta) / cos(theta),
                          color = "cyan", linewidth = 0.5, alpha = 0.6, linetype = "dashed")
}

combined <- p1 / p2 / p3 +
  plot_annotation(
    title = "MBRT 4h: B cell Cdkn1a spatial signal — three views",
    theme = theme(plot.title = element_text(color = "white", size = 16, face = "bold"),
                  plot.background = element_rect(fill = "grey5", color = NA)))

ggsave(file.path(PLOT_DIR, "bcell_p21_spatial_v2.png"), plot = combined,
       width = 14, height = 24, dpi = 200)
cat("Saved: plots/bcell_p21_spatial_v2.png\n")
