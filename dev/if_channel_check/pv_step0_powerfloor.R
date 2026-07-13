# Step 0 of the rotation-null reconciliation: phase-null power floor + unit counts.
# Geometry-only, no counts needed. Decides whether a phase null can resolve anything
# on this block before any expensive confound-isolation work is run.
suppressPackageStartupMessages({ library(arrow); library(data.table) })

META <- "/mnt/data/projects/spatial-rads/inputs/mutter01/Analysis_Mutter_01_CosMmR_Mutter_updated_metadata.parquet"
SLIDE <- "20250529_214712_S4"
MBRT_BLOCK <- "Block_21"
set.seed(1)

peak_fovs <- c(224, 209, 201, 186, 172, 225, 229, 215, 175,
               176, 206, 181, 168, 167)

meta <- as.data.table(read_parquet(META))
setnames(meta, "ImmuneAtlas_ImmGen_Main_cell_Types", "main_type")
blk <- meta[Slide == SLIDE & Block == MBRT_BLOCK]
tum <- blk[main_type == "a"]                    # matches the contested test's cell set
cat(sprintf("Block_21: %d cells total, %d 'a'-bucket tumor cells\n", nrow(blk), nrow(tum)))

## ---- Fit real stripe geometry (identical to pv_rotation_null.R) ----
peak_centroids <- blk[fov %in% peak_fovs,
  .(cx = mean(x_slide_mm), cy = mean(y_slide_mm)), by = fov]
search_theta <- function(cx, cy, n = 4, grid = seq(-pi/2, pi/2, length.out = 181)) {
  best <- list(theta = NA, within = Inf)
  for (th in grid) {
    d <- -sin(th) * cx + cos(th) * cy
    km <- kmeans(d, centers = n, nstart = 10)
    if (km$tot.withinss < best$within)
      best <- list(theta = th, within = km$tot.withinss,
                   centers = sort(km$centers[,1]), cluster = km$cluster)
  }
  best
}
fit <- search_theta(peak_centroids$cx, peak_centroids$cy)
theta_real <- fit$theta
spacing <- median(diff(fit$centers))
peak_centroids[, d_perp := -sin(theta_real) * cx + cos(theta_real) * cy]
peak_centroids[, stripe := fit$cluster]
within_sd <- peak_centroids[, .(sd = sd(d_perp - fit$centers[stripe])), by = stripe][, mean(sd)]
peak_half <- within_sd + 0.15
cat(sprintf("Stripe fit: theta=%.1f deg, spacing=%.3f mm, within_sd=%.3f mm, peak_half=%.3f mm\n",
            theta_real * 180/pi, spacing, within_sd, peak_half))

## ---- Tumor extent along the beam-perpendicular axis -> number of periods ----
pos <- as.matrix(tum[, .(x_slide_mm, y_slide_mm)])
d_cells <- -sin(theta_real) * pos[,1] + cos(theta_real) * pos[,2]
d_span <- diff(range(d_cells))
n_periods <- d_span / spacing
cat(sprintf("Tumor perpendicular extent: %.3f mm = %.2f beam periods\n", d_span, n_periods))

## ---- label_cells identical to pv_rotation_null.R ----
label_cells <- function(d, phase_offset, spacing, half_w) {
  d_range <- range(d)
  centers_seq <- seq(d_range[1] - spacing, d_range[2] + spacing, by = spacing) + phase_offset
  dist_to_peak <- sapply(d, function(x) min(abs(x - centers_seq)))
  ifelse(dist_to_peak < half_w, "peak",
         ifelse(dist_to_peak > (spacing / 2 - half_w), "valley", "transition"))
}

## ---- Independent phase-offset count via peak-set decorrelation ----
# delta(phi) is a function of the peak/valley label sets, so delta decorrelates no
# slower than the peak SET does. The phi-lag at which the peak set becomes disjoint
# from phi=0 is thus an OPTIMISTIC (upper-bound) count of independent offsets ->
# an optimistic (lower-bound) p-floor. If even this floor is too high, the phase
# null is definitively futile.
phi_grid <- seq(0, spacing/2, length.out = 200)
peak0 <- which(label_cells(d_cells, 0, spacing, peak_half) == "peak")
jacc <- sapply(phi_grid, function(ph) {
  pk <- which(label_cells(d_cells, ph, spacing, peak_half) == "peak")
  length(intersect(pk, peak0)) / length(union(pk, peak0))
})
# decorrelation lag = first phi where Jaccard <= 1/e
lag_idx <- which(jacc <= exp(-1))[1]
decorr_lag <- if (is.na(lag_idx)) spacing/2 else phi_grid[lag_idx]
n_indep_decorr <- (spacing/2) / decorr_lag
cat(sprintf("Peak-set decorrelation lag (Jaccard<=1/e): %.3f mm -> %.2f independent offsets in [0, spacing/2)\n",
            decorr_lag, n_indep_decorr))

## ---- p-floor under the two-sided |delta| phase null ----
# smallest empirical p = 1 / (n independent offsets) when the real config is most extreme
p_floor_decorr <- 1 / n_indep_decorr
p_floor_periods <- 1 / n_periods
cat(sprintf("\n=== PHASE-NULL POWER FLOOR ===\n"))
cat(sprintf("  optimistic (decorrelation): p_floor ~ %.3f\n", p_floor_decorr))
cat(sprintf("  coarse (n periods):         p_floor ~ %.3f\n", p_floor_periods))

## ---- FOV pseudobulk-unit count (for the FOV-aggregated arm) ----
# Assign every block FOV a mean d, label by the real stripe grid, count peak vs valley FOVs.
fov_d <- blk[, .(cx = mean(x_slide_mm), cy = mean(y_slide_mm)), by = fov]
fov_d[, d := -sin(theta_real) * cx + cos(theta_real) * cy]
centers_seq <- seq(min(fov_d$d) - spacing, max(fov_d$d) + spacing, by = spacing) + fit$centers[1]
fov_d[, dist_peak := sapply(d, function(x) min(abs(x - centers_seq)))]
fov_d[, zone := ifelse(dist_peak < peak_half, "peak",
                ifelse(dist_peak > (spacing/2 - peak_half), "valley", "transition"))]
cat(sprintf("\nFOV pseudobulk units (real stripe grid): "))
print(table(fov_d$zone))
cat(sprintf("  (n peak FOVs used to FIT geometry: %d)\n", length(peak_fovs)))

## ---- Registration proxy ----
cat(sprintf("\nRegistration proxy: peak-FOV-centroid scatter about fitted stripe = %.3f mm (within_sd).\n", within_sd))
cat("  NOTE: true H2AX-overlay-to-centroid registration uncertainty is NOT measured here;\n")
cat("  it requires the H2AX IF overlay and must be measured before the phase-null band is set.\n")

## ---- Verdict line ----
plausible_effect_p <- 0.05  # a real dose effect ought to clear at least this under a working null
cat(sprintf("\n=== GATE ===\n"))
if (p_floor_decorr > plausible_effect_p) {
  cat(sprintf("  p_floor (%.3f) EXCEEDS 0.05: phase null CANNOT reach significance regardless of truth.\n", p_floor_decorr))
  cat("  -> SKIP the phase null (Steps 3-3.5); let effect-size stability + FOV-pseudobulk carry the verdict.\n")
} else {
  cat(sprintf("  p_floor (%.3f) is below 0.05: phase null can in principle resolve; proceed to Steps 1-3.\n", p_floor_decorr))
}
