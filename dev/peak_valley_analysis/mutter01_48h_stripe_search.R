library(arrow)
library(data.table)
library(ggplot2)
library(ggrepel)

PLOT_DIR <- "plots"
DATA_DIR <- "data"
INPUT_DIR <- "/mnt/data/projects/spatial-rads/inputs/mutter01"

# --- Load Mutter_01 metadata + counts for MBRT_day2 (48h) ---
cat("Loading Mutter_01 metadata...\n")
meta <- as.data.table(read_parquet(
  file.path(INPUT_DIR, "Analysis_Mutter_01_CosMmR_Mutter_updated_metadata.parquet"),
  col_select = c("cell_id", "fov", "Slide", "Condition", "Block",
                 "x_slide_mm", "y_slide_mm", "qcFlagsCell",
                 "ImmuneAtlas_ImmGen_Main_cell_Types")))
meta <- meta[qcFlagsCell == "Pass"]

cat("Loading counts...\n")
counts <- as.data.table(read_parquet(
  file.path(INPUT_DIR, "projects_Mutter_01_CosMmR_app_data_raw_tx_counts_matrix.parquet")))

gene_cols <- setdiff(names(counts), c("Slide", "fov", "cell_id"))

# --- Extract MBRT_day2 (48h, Block_36) ---
mbrt48 <- meta[Condition == "MBRT_day2"]
cat(sprintf("MBRT_day2: %d cells, %d FOVs, Block %s\n",
            nrow(mbrt48), uniqueN(mbrt48$fov), unique(mbrt48$Block)))
cat(sprintf("Spatial extent: x=[%.1f, %.1f], y=[%.1f, %.1f] mm\n",
            min(mbrt48$x_slide_mm), max(mbrt48$x_slide_mm),
            min(mbrt48$y_slide_mm), max(mbrt48$y_slide_mm)))

# Merge counts
counts48 <- counts[cell_id %in% mbrt48$cell_id]
setkeyv(counts48, "cell_id"); setkeyv(mbrt48, "cell_id")
counts48 <- counts48[mbrt48$cell_id]

# CP10K normalize target genes
sig_genes <- c("Cdkn1a", "Clu", "Ctsd")
lib_size <- rowSums(counts48[, ..gene_cols])
for (g in sig_genes) {
  mbrt48[, (g) := log2(counts48[[g]] / pmax(lib_size, 1) * 1e4 + 1)]
}
mbrt48[, sig_score := (Cdkn1a + Clu + Ctsd) / 3]

# --- Theta search: find stripe angle ---
# Known: spacing ~1.05mm. Search angle, score by autocorrelation at expected lag.
BIN_WIDTH <- 0.05  # mm
EXPECTED_SPACING <- 1.05  # mm
EXPECTED_LAG <- round(EXPECTED_SPACING / BIN_WIDTH)

theta_grid <- seq(-60, 60, by = 1) * pi / 180

score_periodicity <- function(x, y, expr, theta, bin_width, expected_lag) {
  d_perp <- -sin(theta) * x + cos(theta) * y
  bins <- cut(d_perp, breaks = seq(min(d_perp) - bin_width, max(d_perp) + bin_width, by = bin_width))
  bin_means <- tapply(expr, bins, mean, na.rm = TRUE)
  bin_means <- bin_means[!is.na(bin_means)]
  if (length(bin_means) < expected_lag + 5) return(NA_real_)
  ac <- acf(bin_means, lag.max = expected_lag + 5, plot = FALSE)$acf
  if (length(ac) > expected_lag) ac[expected_lag + 1] else NA_real_
}

cat("\nTheta search using p21 (Cdkn1a)...\n")
p21_scores <- sapply(theta_grid, function(th) {
  score_periodicity(mbrt48$x_slide_mm, mbrt48$y_slide_mm, mbrt48$Cdkn1a,
                    th, BIN_WIDTH, EXPECTED_LAG)
})

cat("Theta search using 3-gene signature...\n")
sig_scores <- sapply(theta_grid, function(th) {
  score_periodicity(mbrt48$x_slide_mm, mbrt48$y_slide_mm, mbrt48$sig_score,
                    th, BIN_WIDTH, EXPECTED_LAG)
})

theta_deg <- theta_grid * 180 / pi
best_p21_idx <- which.max(p21_scores)
best_sig_idx <- which.max(sig_scores)
cat(sprintf("Best p21 angle: %.0f deg (autocorr = %.4f)\n",
            theta_deg[best_p21_idx], p21_scores[best_p21_idx]))
cat(sprintf("Best signature angle: %.0f deg (autocorr = %.4f)\n",
            theta_deg[best_sig_idx], sig_scores[best_sig_idx]))

# Reference from 4h Block_21
cat(sprintf("Reference (4h Block_21): -25 deg\n"))

# --- Permutation test ---
N_PERM <- 500
cat(sprintf("\nPermutation test (%d permutations)...\n", N_PERM))

best_real_p21 <- max(p21_scores, na.rm = TRUE)
best_real_sig <- max(sig_scores, na.rm = TRUE)

set.seed(42)
perm_max_p21 <- numeric(N_PERM)
perm_max_sig <- numeric(N_PERM)
for (i in seq_len(N_PERM)) {
  shuf_idx <- sample(nrow(mbrt48))
  p21_shuf <- mbrt48$Cdkn1a[shuf_idx]
  sig_shuf <- mbrt48$sig_score[shuf_idx]
  perm_p21 <- sapply(theta_grid, function(th)
    score_periodicity(mbrt48$x_slide_mm, mbrt48$y_slide_mm, p21_shuf,
                      th, BIN_WIDTH, EXPECTED_LAG))
  perm_sig <- sapply(theta_grid, function(th)
    score_periodicity(mbrt48$x_slide_mm, mbrt48$y_slide_mm, sig_shuf,
                      th, BIN_WIDTH, EXPECTED_LAG))
  perm_max_p21[i] <- max(perm_p21, na.rm = TRUE)
  perm_max_sig[i] <- max(perm_sig, na.rm = TRUE)
  if (i %% 100 == 0) cat(sprintf("  permutation %d/%d\n", i, N_PERM))
}

perm_p_p21 <- mean(perm_max_p21 >= best_real_p21)
perm_p_sig <- mean(perm_max_sig >= best_real_sig)
cat(sprintf("\nPermutation p-values:\n"))
cat(sprintf("  p21: observed=%.4f, perm p=%.3f\n", best_real_p21, perm_p_p21))
cat(sprintf("  signature: observed=%.4f, perm p=%.3f\n", best_real_sig, perm_p_sig))

# --- Plot: theta scan ---
scan_df <- data.frame(
  theta = rep(theta_deg, 2),
  autocorr = c(p21_scores, sig_scores),
  signal = rep(c("Cdkn1a (p21)", "3-gene signature"), each = length(theta_deg))
)

p_scan <- ggplot(scan_df, aes(x = theta, y = autocorr, color = signal)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = -25, linetype = "dashed", color = "grey40") +
  annotate("text", x = -25, y = max(scan_df$autocorr, na.rm = TRUE) * 0.95,
           label = "4h reference (-25°)", size = 3, hjust = -0.05) +
  geom_point(data = scan_df[scan_df$theta == theta_deg[best_p21_idx] & scan_df$signal == "Cdkn1a (p21)", ],
             size = 3) +
  geom_point(data = scan_df[scan_df$theta == theta_deg[best_sig_idx] & scan_df$signal == "3-gene signature", ],
             size = 3) +
  labs(x = "Stripe angle (degrees)", y = "Autocorrelation at 1.05mm lag",
       title = "MBRT day 2 (48h) stripe angle search — Mutter_01 Block_36",
       subtitle = sprintf("p21 best: %.0f° (perm p=%.3f); signature best: %.0f° (perm p=%.3f)",
                          theta_deg[best_p21_idx], perm_p_p21,
                          theta_deg[best_sig_idx], perm_p_sig),
       color = NULL) +
  theme_bw(base_size = 12)
ggsave(file.path(PLOT_DIR, "mutter01_48h_theta_scan.png"), p_scan, width = 10, height = 4, dpi = 200)
cat("Saved: mutter01_48h_theta_scan.png\n")

# --- Use best angle to classify cells ---
use_sig <- if (perm_p_sig < perm_p_p21) "signature" else "p21"
best_theta <- if (use_sig == "signature") theta_grid[best_sig_idx] else theta_grid[best_p21_idx]
cat(sprintf("\nUsing %s angle (%.0f deg) for classification\n", use_sig, best_theta * 180 / pi))

d_perp <- -sin(best_theta) * mbrt48$x_slide_mm + cos(best_theta) * mbrt48$y_slide_mm

# Find beam centers via binned expression peaks
bins <- cut(d_perp, breaks = seq(min(d_perp) - BIN_WIDTH, max(d_perp) + BIN_WIDTH, by = BIN_WIDTH))
bin_centers <- tapply(d_perp, bins, mean, na.rm = TRUE)
bin_expr <- tapply(mbrt48$Cdkn1a, bins, mean, na.rm = TRUE)
valid <- !is.na(bin_expr) & !is.na(bin_centers)
bin_centers <- bin_centers[valid]
bin_expr <- bin_expr[valid]

# Seed beam centers at expected spacing from FFT/autocorrelation
d_range <- range(d_perp)
n_beams <- round(diff(d_range) / EXPECTED_SPACING)
beam_centers <- seq(d_range[1] + EXPECTED_SPACING / 2,
                    d_range[2] - EXPECTED_SPACING / 2,
                    by = EXPECTED_SPACING)

# Classify
dist_to_nearest <- sapply(d_perp, function(d) min(abs(d - beam_centers)))
peak_half <- EXPECTED_SPACING / 4
mbrt48[, pv_zone := fifelse(dist_to_nearest < peak_half, "peak",
                    fifelse(dist_to_nearest > EXPECTED_SPACING / 2 - peak_half, "valley", "transition"))]
mbrt48[, d_to_peak := dist_to_nearest]

cat("Zone counts:\n")
print(table(mbrt48$pv_zone))

# --- Per-cell spatial map ---
x_breaks <- seq(floor(min(mbrt48$x_slide_mm)), ceiling(max(mbrt48$x_slide_mm)), by = 1)
y_breaks <- seq(floor(min(mbrt48$y_slide_mm)), ceiling(max(mbrt48$y_slide_mm)), by = 1)

# Binary p21 expressing map
mbrt48[, p21_positive := Cdkn1a > 0]
p_spatial <- ggplot(mbrt48, aes(x = x_slide_mm, y = y_slide_mm, color = p21_positive)) +
  geom_point(size = 0.05, alpha = 0.3) +
  scale_color_manual(values = c("FALSE" = "grey80", "TRUE" = "#e41a1c"),
                     labels = c("Negative", "p21+")) +
  scale_x_continuous(breaks = x_breaks) +
  scale_y_continuous(breaks = y_breaks) +
  coord_fixed() +
  labs(x = "x (mm)", y = "y (mm)", color = NULL,
       title = "MBRT day 2 (48h) — per-cell Cdkn1a expression",
       subtitle = sprintf("Block_36; %d cells; beam angle = %.0f°",
                          nrow(mbrt48), round(best_theta * 180 / pi))) +
  theme_bw(base_size = 12) +
  theme(panel.background = element_rect(fill = "grey95"), panel.grid = element_blank())

# Add beam centerlines
x_lim <- range(mbrt48$x_slide_mm)
for (bc in beam_centers) {
  line_df <- data.frame(
    x = x_lim,
    y = (bc + sin(best_theta) * x_lim) / cos(best_theta)
  )
  p_spatial <- p_spatial +
    geom_line(data = line_df, aes(x = x, y = y), inherit.aes = FALSE,
              color = "cyan", linewidth = 0.5, alpha = 0.7, linetype = "dashed")
}
ggsave(file.path(PLOT_DIR, "mutter01_48h_p21_spatial.png"), p_spatial, width = 8, height = 8, dpi = 200)
cat("Saved: mutter01_48h_p21_spatial.png\n")

# --- Zone-colored spatial map ---
mbrt48_zoned <- mbrt48[pv_zone %in% c("peak", "valley")]
p_zones <- ggplot(mbrt48_zoned, aes(x = x_slide_mm, y = y_slide_mm, color = pv_zone)) +
  geom_point(size = 0.05, alpha = 0.3) +
  scale_color_manual(values = c("peak" = "#e41a1c", "valley" = "#377eb8")) +
  scale_x_continuous(breaks = x_breaks) +
  scale_y_continuous(breaks = y_breaks) +
  coord_fixed() +
  labs(x = "x (mm)", y = "y (mm)", color = "Zone",
       title = "MBRT day 2 (48h) — peak/valley classification",
       subtitle = sprintf("Block_36; angle=%.0f°, spacing=%.2fmm",
                          round(best_theta * 180 / pi), EXPECTED_SPACING)) +
  theme_bw(base_size = 12) +
  theme(panel.background = element_rect(fill = "grey95"), panel.grid = element_blank())
for (bc in beam_centers) {
  line_df <- data.frame(x = x_lim,
    y = (bc + sin(best_theta) * x_lim) / cos(best_theta))
  p_zones <- p_zones +
    geom_line(data = line_df, aes(x = x, y = y), inherit.aes = FALSE,
              color = "cyan", linewidth = 0.5, alpha = 0.7, linetype = "dashed")
}
ggsave(file.path(PLOT_DIR, "mutter01_48h_zones.png"), p_zones, width = 8, height = 8, dpi = 200)
cat("Saved: mutter01_48h_zones.png\n")

# --- Gradient plot ---
mbrt48_pv <- mbrt48[pv_zone %in% c("peak", "valley")]
bin_width_grad <- 0.02
mbrt48_pv[, dist_bin := round(d_to_peak / bin_width_grad) * bin_width_grad]
grad <- mbrt48_pv[, .(mean_p21 = mean(Cdkn1a), mean_sig = mean(sig_score),
                       pct_p21 = 100 * mean(Cdkn1a > 0), n = .N),
                   by = dist_bin][order(dist_bin)]

p_grad <- ggplot(grad, aes(x = dist_bin, y = mean_p21)) +
  geom_point(size = 0.8, alpha = 0.5) +
  geom_smooth(method = "loess", span = 0.5, color = "#e41a1c", se = TRUE) +
  geom_vline(xintercept = peak_half, linetype = "dashed", color = "#e41a1c", alpha = 0.6) +
  geom_vline(xintercept = EXPECTED_SPACING / 2 - peak_half, linetype = "dashed",
             color = "#377eb8", alpha = 0.6) +
  labs(x = "Distance from nearest beam centerline (mm)",
       y = "Mean Cdkn1a expression (log2 CP10K + 1)",
       title = "MBRT day 2 (48h) — p21 gradient from peak center",
       subtitle = sprintf("Block_36; angle=%.0f°; red/blue dashed = peak/valley boundaries",
                          round(best_theta * 180 / pi))) +
  theme_bw(base_size = 12)
ggsave(file.path(PLOT_DIR, "mutter01_48h_p21_gradient.png"), p_grad, width = 8, height = 4, dpi = 200)
cat("Saved: mutter01_48h_p21_gradient.png\n")

# --- Peak vs valley summary ---
cat("\n=== Peak vs Valley at 48h ===\n")
pv_summary <- mbrt48[pv_zone %in% c("peak", "valley"),
  .(n = .N,
    mean_p21 = round(mean(Cdkn1a), 4),
    pct_p21 = round(100 * mean(Cdkn1a > 0), 2),
    mean_clu = round(mean(Clu), 4),
    mean_ctsd = round(mean(Ctsd), 4),
    mean_sig = round(mean(sig_score), 4)),
  by = pv_zone]
print(pv_summary)

# FOV-level Wilcoxon
fov_p21 <- mbrt48[pv_zone %in% c("peak", "valley"),
  .(pct_p21 = 100 * mean(Cdkn1a > 0), n = .N), by = .(fov, pv_zone)][n >= 10]
pk_fovs <- fov_p21[pv_zone == "peak", pct_p21]
vl_fovs <- fov_p21[pv_zone == "valley", pct_p21]
if (length(pk_fovs) >= 3 && length(vl_fovs) >= 3) {
  wt <- wilcox.test(pk_fovs, vl_fovs)
  cat(sprintf("\nFOV-level Wilcoxon (p21 %% expressing): p = %.4f (n_pk=%d, n_vl=%d)\n",
              wt$p.value, length(pk_fovs), length(vl_fovs)))
}

cat("\nMutter_01 48h stripe search complete.\n")
