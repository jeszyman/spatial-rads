library(arrow)
library(tidyverse)
library(MASS)

INPUT_DIR <- "/mnt/data/projects/spatial-rads/inputs"
PLOT_DIR  <- "plots"

meta <- read_parquet(
  file.path(INPUT_DIR, "Analysis_Mutter_01_CosMmR_Mutter_updated_metadata.parquet"),
  col_select = c("cell_id", "fov", "Condition", "x_slide_mm", "y_slide_mm",
                  "qcFlagsCell", "ImmuneAtlas_ImmGen_Main_cell_Types")
) %>% as.data.frame() %>% filter(qcFlagsCell == "Pass", Condition == "MBRT_4h")

peak_genes <- c("Hsd17b2", "Clec7a", "Itgb1", "Tlr8", "Tlr2",
                "Tnfrsf10b", "Cxcr4", "Ccr7", "Aif1", "Hif1a",
                "Atf3", "Mki67", "Cdkn1a", "Gadd45b", "Bax")

expr <- read_parquet(
  file.path(INPUT_DIR, "projects_Mutter_01_CosMmR_app_data_raw_tx_counts_matrix.parquet"),
  col_select = c("cell_id", all_of(peak_genes))
) %>% as.data.frame()

mbrt4h <- left_join(meta, expr, by = "cell_id")

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

bcells <- mbrt4h %>%
  filter(ImmuneAtlas_ImmGen_Main_cell_Types %in% c("B.cell", "Memory.B", "Plasmablast"))
cat(sprintf("B cells: %d\n", nrow(bcells)))

for (g in peak_genes) {
  bcells[[g]][is.na(bcells[[g]])] <- 0
  mu <- mean(bcells[[g]])
  s <- sd(bcells[[g]])
  if (s > 0) bcells[[paste0(g, "_z")]] <- (bcells[[g]] - mu) / s
  else bcells[[paste0(g, "_z")]] <- 0
}
z_cols <- paste0(peak_genes, "_z")
bcells$composite <- rowMeans(bcells[, z_cols, drop = FALSE], na.rm = TRUE)

thresh <- quantile(bcells$composite, 0.75)
bcells$damage_high <- bcells$composite > thresh
cat(sprintf("Damage-high: %d cells\n", sum(bcells$damage_high)))

d_bcell <- sapply(-sin(theta) * bcells$x_slide_mm + cos(theta) * bcells$y_slide_mm,
                  function(d) min(abs(d - all_centers)))
rho <- cor(bcells$composite, d_bcell, method = "spearman")
cat(sprintf("Composite rho vs distance: %.3f\n", rho))

bw <- 0.25
xlim <- range(bcells$x_slide_mm) + c(-0.3, 0.3)
ylim <- range(bcells$y_slide_mm) + c(-0.3, 0.3)
ngrid <- 200

kde_all <- kde2d(bcells$x_slide_mm, bcells$y_slide_mm,
                 h = bw, n = ngrid, lims = c(xlim, ylim))
kde_hi <- kde2d(bcells$x_slide_mm[bcells$damage_high],
                bcells$y_slide_mm[bcells$damage_high],
                h = bw, n = ngrid, lims = c(xlim, ylim))
cat("KDE done\n")

grid <- expand.grid(x = kde_all$x, y = kde_all$y)
grid$d_all <- as.vector(kde_all$z)
grid$d_hi <- as.vector(kde_hi$z)
grid$enrichment <- grid$d_hi / (grid$d_all + 1e-12)
density_thresh <- quantile(grid$d_all, 0.25)
grid <- grid %>% filter(d_all > density_thresh)

xr <- range(mbrt4h$x_slide_mm) + c(-0.5, 0.5)
band_df <- do.call(rbind, lapply(all_centers, function(ac) {
  y_lo_l <- (ac - peak_half + sin(theta) * xr[1]) / cos(theta)
  y_lo_r <- (ac - peak_half + sin(theta) * xr[2]) / cos(theta)
  y_hi_l <- (ac + peak_half + sin(theta) * xr[1]) / cos(theta)
  y_hi_r <- (ac + peak_half + sin(theta) * xr[2]) / cos(theta)
  data.frame(x = c(xr[1], xr[2], xr[2], xr[1]),
             y = c(y_lo_l, y_lo_r, y_hi_r, y_hi_l),
             group = paste0("p_", round(ac, 3)))
}))

p <- ggplot() +
  geom_polygon(data = band_df, aes(x = x, y = y, group = group),
               fill = "#E74C3C", alpha = 0.12) +
  geom_raster(data = grid, aes(x = x, y = y, fill = enrichment), interpolate = TRUE) +
  scale_fill_viridis_c(option = "inferno", name = "Damage-high\nenrichment") +
  coord_fixed(xlim = range(mbrt4h$x_slide_mm) + c(-0.1, 0.1),
              ylim = range(mbrt4h$y_slide_mm) + c(-0.1, 0.1)) +
  labs(title = "MBRT 4h: B cell multi-gene damage score (smoothed enrichment)",
       subtitle = sprintf("15-gene composite | Top quartile enrichment | rho=%.3f | bw=250um", rho)) +
  theme_void(base_size = 15) +
  theme(plot.background = element_rect(fill = "grey8", color = NA),
        plot.title = element_text(color = "white", face = "bold", size = 16),
        plot.subtitle = element_text(color = "grey60", size = 10),
        legend.text = element_text(color = "white"),
        legend.title = element_text(color = "white"),
        plot.margin = margin(15, 15, 15, 15))

for (ac in all_centers) {
  p <- p + geom_abline(intercept = ac / cos(theta), slope = sin(theta) / cos(theta),
                        color = "white", linewidth = 0.3, alpha = 0.35, linetype = "dotted")
}

ggsave(file.path(PLOT_DIR, "bcell_multigene_damage_map.png"), plot = p,
       width = 14, height = 10, dpi = 250)
cat("Saved: plots/bcell_multigene_damage_map.png\n")
