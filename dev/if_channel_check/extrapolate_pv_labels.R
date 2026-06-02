suppressPackageStartupMessages({
  library(arrow); library(data.table); library(Matrix); library(MERINGUE)
})

OUT <- "/home/jeszyman/repos/spatial-rads/dev/if_channel_check"
BLOCK <- "Block_21"; SLIDE <- "20250529_214712_S4"
N_SUB <- 8000  # more cells since we're labeling all of them
set.seed(1)

peak_fovs   <- c(224, 209, 201, 186, 172, 225, 229, 215, 175,
                 176, 206, 181, 168, 167)
valley_fovs <- c(199, 213, 227, 188, 204, 218, 178, 165, 170, 231, 189,
                 173, 154)

meta <- as.data.table(read_parquet("/tmp/mutter01_meta.parquet"))
counts_df <- as.data.table(read_parquet("/tmp/mutter01_counts.parquet"))

m <- meta[Slide == SLIDE & Block == BLOCK]
counts_df <- counts_df[cell_id %in% m$cell_id]
setkeyv(m, "cell_id"); setkeyv(counts_df, "cell_id")
m <- m[cell_id %in% counts_df$cell_id]
idx <- sample(nrow(m), min(N_SUB, nrow(m)))
m <- m[idx]; counts_df <- counts_df[cell_id %in% m$cell_id]
setorderv(m, "cell_id"); setorderv(counts_df, "cell_id")
cat("Cells:", nrow(m), "\n")

# ---- Step 1: fit stripe angle from peak FOV centroids ----
peak_cells <- m[fov %in% peak_fovs]
peak_fov_centroids <- peak_cells[, .(cx = mean(x_slide_mm),
                                     cy = mean(y_slide_mm)),
                                 by = fov]
search_theta <- function(cx, cy, n_stripes = 4,
                         theta_grid = seq(-pi/2, pi/2, length.out = 181)) {
  best <- list(theta = NA, within_var = Inf)
  for (th in theta_grid) {
    d <- -sin(th) * cx + cos(th) * cy
    km <- kmeans(d, centers = n_stripes, nstart = 10)
    if (km$tot.withinss < best$within_var) {
      best <- list(theta = th, within_var = km$tot.withinss,
                   centers = sort(km$centers[, 1]), d = d,
                   cluster = km$cluster)
    }
  }
  best
}
fit <- search_theta(peak_fov_centroids$cx, peak_fov_centroids$cy)
theta <- fit$theta
peak_fov_centroids[, d_perp := -sin(theta) * cx + cos(theta) * cy]

# Estimate beam half-width from spread of peak FOVs *within* each stripe
peak_fov_centroids[, stripe := fit$cluster]
stripe_centers <- sort(fit$centers)
beam_spacing <- median(diff(stripe_centers))
within_stripe_sd <- peak_fov_centroids[,
  .(sd = sd(d_perp - stripe_centers[stripe])), by = stripe][, mean(sd)]
# Peak half-width: 1 SD of peak-FOV centroids around stripe centerline
# (FOVs are ~0.5 mm, so this captures FOV size + beam width together)
peak_half <- within_stripe_sd + 0.15  # + a 150-um buffer
cat(sprintf("Stripe angle: %.1f deg | spacing: %.3f mm | beam half-width: %.3f mm\n",
            theta * 180 / pi, beam_spacing, peak_half))
cat("Stripe centerlines (perpendicular coord):",
    paste(round(stripe_centers, 3), collapse = ", "), "\n")

# ---- Step 2: extend centerlines beyond the 4 labeled stripes ----
# Extrapolate to cover full tumor extent
pos_all <- as.matrix(m[, .(x_slide_mm, y_slide_mm)])
d_all <- -sin(theta) * pos_all[, 1] + cos(theta) * pos_all[, 2]
d_range <- range(d_all)
all_centers <- seq(stripe_centers[1] - ceiling((stripe_centers[1] - d_range[1]) / beam_spacing) * beam_spacing,
                   stripe_centers[length(stripe_centers)] + ceiling((d_range[2] - stripe_centers[length(stripe_centers)]) / beam_spacing) * beam_spacing,
                   by = beam_spacing)
cat("Extrapolated stripe centers:", paste(round(all_centers, 3), collapse = ", "), "\n")

# ---- Step 3: label every cell ----
dist_to_nearest_peak <- sapply(d_all, function(d) min(abs(d - all_centers)))
label <- rep("transition", length(d_all))
label[dist_to_nearest_peak < peak_half] <- "peak"
# Valley = midway between peaks (>~0.6 from any peak for spacing 1.02)
valley_min <- (beam_spacing / 2) - peak_half  # distance from midline to valley core
label[dist_to_nearest_peak > (beam_spacing / 2 - peak_half)] <- "valley"
cat("Label counts (extrapolated):\n"); print(table(label))

m[, pv_extrap := label]
m[, d_to_peak := dist_to_nearest_peak]

# Compare to manual labels
m[, pv_manual := ifelse(fov %in% peak_fovs, "peak",
                         ifelse(fov %in% valley_fovs, "valley", NA))]
cat("\nManual vs extrapolated (only rows with manual label):\n")
print(m[!is.na(pv_manual), table(manual = pv_manual, extrap = pv_extrap)])
concordance <- m[!is.na(pv_manual),
                 sum(pv_manual == pv_extrap) / .N]
cat(sprintf("Manual-extrapolated concordance: %.1f%%\n", 100 * concordance))

# Save cell labels
fwrite(m[, .(cell_id, Slide, Block, fov, x_slide_mm, y_slide_mm,
             d_to_peak, pv_extrap, pv_manual)],
       "/tmp/block21_extrapolated_pv.tsv", sep = "\t")
cat("Saved /tmp/block21_extrapolated_pv.tsv\n")

# ---- Visualize ----
pos <- pos_all
x_range <- diff(range(pos[, 1])); y_range <- diff(range(pos[, 2]))
asp <- y_range / x_range
panel_w <- 900; panel_h <- round(panel_w * asp)
png(file.path(OUT, "block21_extrapolated_pv.png"),
    width = panel_w * 2 + 80, height = (panel_h + 60) * 2 + 40, res = 130)
par(mfrow = c(2, 2), mar = c(2, 2, 3, 1))

add_all_lines <- function() {
  for (c in all_centers) {
    # Line: -sin(theta)*x + cos(theta)*y = c -> y = (c + sin(theta)*x) / cos(theta)
    abline(a = c / cos(theta), b = sin(theta) / cos(theta),
           col = "black", lty = 2, lwd = 1.2)
  }
}

# Manual labels
col_manual <- ifelse(is.na(m$pv_manual), "grey85",
                     ifelse(m$pv_manual == "peak", "#e41a1c", "#377eb8"))
plot(pos, col = col_manual, pch = 19, cex = 0.3, asp = 1,
     main = "Manual peak/valley FOV labels", xlab = "x mm", ylab = "y mm")
add_all_lines()
legend("topright", legend = c("peak (manual)", "valley (manual)", "unlabeled"),
       col = c("#e41a1c", "#377eb8", "grey85"), pch = 19, bty = "n", cex = 0.8)

# Extrapolated labels
col_extrap <- c(peak = "#e41a1c", valley = "#377eb8",
                transition = "grey85")[m$pv_extrap]
plot(pos, col = col_extrap, pch = 19, cex = 0.3, asp = 1,
     main = sprintf("Extrapolated labels (peak half = %.2f mm)", peak_half),
     xlab = "x mm", ylab = "y mm")
add_all_lines()
legend("topright", legend = c("peak", "valley", "transition"),
       col = c("#e41a1c", "#377eb8", "grey85"), pch = 19, bty = "n", cex = 0.8)

# Perpendicular distance heatmap
dist_cols <- colorRampPalette(c("#e41a1c", "white", "#377eb8"))(100)
d_norm <- pmin(pmax((m$d_to_peak - 0) / (beam_spacing / 2), 0), 1)
col_dist <- dist_cols[round(d_norm * 99) + 1]
plot(pos, col = col_dist, pch = 19, cex = 0.3, asp = 1,
     main = "Perpendicular distance to nearest peak (mm)",
     xlab = "x mm", ylab = "y mm")
add_all_lines()

# Histogram of distance, colored by manual label for calibration
hist(m[pv_manual == "peak", d_to_peak], breaks = 30,
     col = adjustcolor("#e41a1c", 0.5), border = NA,
     xlab = "Distance to nearest peak centerline (mm)",
     main = "Manual peak (red) vs valley (blue) FOVs,\ndistance to extrapolated peak line",
     xlim = c(0, beam_spacing / 2 + 0.1))
hist(m[pv_manual == "valley", d_to_peak], breaks = 30,
     col = adjustcolor("#377eb8", 0.5), border = NA, add = TRUE)
abline(v = peak_half, lty = 2, col = "#e41a1c", lwd = 2)
abline(v = beam_spacing / 2 - peak_half, lty = 2, col = "#377eb8", lwd = 2)
legend("topright", legend = c("manual peak", "manual valley",
                              "peak cutoff", "valley cutoff"),
       fill = c(adjustcolor("#e41a1c", 0.5),
                adjustcolor("#377eb8", 0.5), NA, NA),
       border = c(NA, NA, NA, NA),
       lty = c(NA, NA, 2, 2),
       col = c(NA, NA, "#e41a1c", "#377eb8"), bty = "n", cex = 0.7)
dev.off()
cat("Saved block21_extrapolated_pv.png\n")
