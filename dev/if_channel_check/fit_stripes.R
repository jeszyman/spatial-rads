library(arrow)
library(data.table)
library(ggplot2)

peak_fovs  <- unique(c(224, 209, 201, 186, 172, 225, 229, 215, 175, 176, 206, 181, 168, 167))
valley_fovs <- c(199, 213, 227, 188, 204, 218, 178, 165, 170, 231, 189, 173, 154)

meta <- as.data.table(read_parquet("/tmp/mutter01_meta.parquet"))
s4 <- meta[Slide == "20250529_214712_S4" & Block == "Block_21"]

# FOV centroids
fov_ct <- s4[, .(cx = mean(x_slide_mm), cy = mean(y_slide_mm)), by = fov]
fov_ct[, label := fifelse(fov %in% peak_fovs, "peak",
                          fifelse(fov %in% valley_fovs, "valley", "unlabeled"))]

cat("Label counts:\n"); print(table(fov_ct$label))

# ---- Fit stripe model ----
# Model: lines at angle theta, spacing s, with 4 peak lines at positions y0 + k*s (k=0..3) along axis perpendicular to stripes.
# Project each labeled centroid onto the perpendicular axis: u = -sin(theta)*cx + cos(theta)*cy
# Peak: u closest to one of 4 peak positions. Valley: u closest to midpoint between peak positions.
# Parameters: theta, s, y0

project <- function(cx, cy, theta) -sin(theta) * cx + cos(theta) * cy

loss <- function(par, pk, vl) {
  theta <- par[1]; s <- par[2]; y0 <- par[3]
  u_pk <- project(pk$cx, pk$cy, theta)
  u_vl <- project(vl$cx, vl$cy, theta)
  peak_lines <- y0 + 0:3 * s
  # peak centroid should be close to a peak line
  d_pk <- sapply(u_pk, function(u) min(abs(u - peak_lines)))
  # valley should be close to midpoint between adjacent peak lines
  valley_lines <- y0 + (0:2 + 0.5) * s
  d_vl <- sapply(u_vl, function(u) min(abs(u - valley_lines)))
  sum(d_pk^2) + sum(d_vl^2)
}

# Grid search for good init, then optim
theta_grid <- seq(-pi/4, pi/4, length.out = 81)
s_init <- 1.02
pk <- fov_ct[label == "peak"]
vl <- fov_ct[label == "valley"]

# For a given theta, best y0 is a 1D search across a peak line that explains peaks
best <- list(val = Inf)
for (theta in theta_grid) {
  u_pk <- project(pk$cx, pk$cy, theta)
  u_vl <- project(vl$cx, vl$cy, theta)
  # init y0 as median of u_pk mod s
  y0_candidates <- seq(min(c(u_pk, u_vl)) - 0.5, min(c(u_pk, u_vl)) + 1.5, by = 0.05)
  for (y0 in y0_candidates) {
    v <- loss(c(theta, s_init, y0), pk, vl)
    if (v < best$val) best <- list(val = v, par = c(theta, s_init, y0))
  }
}
cat("Grid best loss:", round(best$val, 4), "init par:", round(best$par, 4), "\n")

opt <- optim(best$par, loss, pk = pk, vl = vl, method = "Nelder-Mead",
             control = list(reltol = 1e-10, maxit = 2000))
theta_fit <- opt$par[1]; s_fit <- opt$par[2]; y0_fit <- opt$par[3]
cat(sprintf("Fit: angle = %.2f deg, spacing = %.3f mm, y0 = %.3f mm\n",
            theta_fit * 180/pi, s_fit, y0_fit))
cat("Final loss (sum sq residual in mm):", round(opt$value, 5), "\n")

# Residuals for each labeled FOV
fov_ct[, u := project(cx, cy, theta_fit)]
peak_lines <- y0_fit + 0:3 * s_fit
valley_lines <- y0_fit + (0:2 + 0.5) * s_fit
fov_ct[, dist_to_nearest_peak := sapply(u, function(x) min(abs(x - peak_lines)))]
fov_ct[, dist_to_nearest_valley := sapply(u, function(x) min(abs(x - valley_lines)))]

cat("\nPeak-labeled residuals (um, to nearest peak line):\n")
print(summary(1000 * fov_ct[label == "peak", dist_to_nearest_peak]))
cat("\nValley-labeled residuals (um, to nearest valley line):\n")
print(summary(1000 * fov_ct[label == "valley", dist_to_nearest_valley]))

# ---- Per-cell classification ----
s4[, u := project(x_slide_mm, y_slide_mm, theta_fit)]
s4[, d_peak := sapply(u, function(x) min(abs(x - peak_lines)))]
s4[, d_valley := sapply(u, function(x) min(abs(x - valley_lines)))]
s4[, class := fifelse(d_peak < 0.130, "peak",
                      fifelse(d_valley < 0.130, "valley", "ambiguous"))]
# Buffer 130 um from each line; anything further in either direction is ambiguous

cat("\nPer-cell class counts:\n")
print(table(s4$class))
cat(sprintf("Peak fraction: %.1f%%\n", 100 * mean(s4$class == "peak")))
cat(sprintf("Valley fraction: %.1f%%\n", 100 * mean(s4$class == "valley")))
cat(sprintf("Ambiguous fraction: %.1f%%\n", 100 * mean(s4$class == "ambiguous")))

# Save per-cell labels
out <- s4[, .(cell_id, fov, x_slide_mm, y_slide_mm, u, d_peak, d_valley, class)]
fwrite(out, "/tmp/mutter01_block21_peak_valley_cells.tsv", sep = "\t")
cat("\nWrote /tmp/mutter01_block21_peak_valley_cells.tsv (", nrow(out), " cells)\n")

# ---- Diagnostic plot ----
# Build peak-line segments for drawing
xr <- range(s4$x_slide_mm); yr <- range(s4$y_slide_mm)
# For each line (peak or valley), find start/end by intersecting with tumor x range
line_coords <- function(u_val, xr, theta) {
  # y = (u_val + sin(theta) * x) / cos(theta)
  x0 <- xr[1]; x1 <- xr[2]
  y0 <- (u_val + sin(theta) * x0) / cos(theta)
  y1 <- (u_val + sin(theta) * x1) / cos(theta)
  data.table(x = c(x0, x1), y = c(y0, y1))
}
peak_lines_df <- rbindlist(lapply(seq_along(peak_lines), function(i) {
  cbind(line = paste0("peak_", i), line_coords(peak_lines[i], xr, theta_fit))
}))
valley_lines_df <- rbindlist(lapply(seq_along(valley_lines), function(i) {
  cbind(line = paste0("valley_", i), line_coords(valley_lines[i], xr, theta_fit))
}))

p <- ggplot() +
  geom_point(data = fov_ct, aes(cx, cy, color = label), size = 3) +
  geom_text(data = fov_ct, aes(cx, cy, label = fov), size = 2, vjust = -1.3) +
  geom_line(data = peak_lines_df, aes(x, y, group = line), color = "red", linewidth = 0.8) +
  geom_line(data = valley_lines_df, aes(x, y, group = line), color = "blue", linetype = "dashed", linewidth = 0.5) +
  scale_color_manual(values = c(peak = "darkred", valley = "darkblue", unlabeled = "gray70")) +
  coord_fixed() +
  labs(title = sprintf("Stripe fit: angle = %.2f deg, spacing = %.3f mm", theta_fit * 180/pi, s_fit),
       x = "x_slide (mm)", y = "y_slide (mm)",
       subtitle = "Red solid = peak lines | Blue dashed = valley midlines") +
  theme_minimal()

ggsave("/tmp/stripe_fit_diagnostic.png", p, width = 10, height = 7, dpi = 150)
cat("Wrote /tmp/stripe_fit_diagnostic.png\n")
