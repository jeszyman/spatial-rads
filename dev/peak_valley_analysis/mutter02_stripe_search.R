## Mutter_02 stripe detection: find MBRT beam angle and spacing at 48h
## Uses FFT-based periodicity search on Cdkn1a and a 3-gene signature
## (Cdkn1a + Clu + Gadd45b; Ctsd not in panel)

library(Seurat)
library(tidyverse)
library(zoo)
library(patchwork)

DATA_DIR  <- "data"
PLOT_DIR  <- "plots"
for (d in c(DATA_DIR, PLOT_DIR)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# Mutter_01 reference parameters
REF_THETA   <- -25
REF_SPACING <- 1.056

# Signature genes (Ctsd absent from panel; substitute Gadd45b)
SIG_GENES <- c("Cdkn1a", "Clu", "Gadd45b")

# ------------------------------------------------------------------
# Helper: FFT-based periodicity score for a given angle
# Projects cells onto perpendicular axis, bins, smooths, then FFT
# Returns: list(score, dominant_period, profile_df)
# ------------------------------------------------------------------
fft_periodicity <- function(x, y, expr, theta_deg, bin_width = 0.05,
                            min_cells = 30, smooth_k = 5) {
  rad <- theta_deg * pi / 180
  d_perp <- -sin(rad) * x + cos(rad) * y

  prof <- tibble(d = d_perp, val = expr) %>%
    mutate(bin = round(d / bin_width) * bin_width) %>%
    group_by(bin) %>%
    summarise(mean_val = mean(val, na.rm = TRUE), n = n(), .groups = "drop") %>%
    filter(n >= min_cells) %>%
    arrange(bin)

  if (nrow(prof) < 2 * smooth_k) {
    return(list(score = 0, period = NA_real_, profile = prof))
  }

  prof$smoothed <- rollmean(prof$mean_val, k = smooth_k, fill = NA, align = "center")
  vals <- prof$smoothed[!is.na(prof$smoothed)]

  if (length(vals) < 8) {
    return(list(score = 0, period = NA_real_, profile = prof))
  }

  # Detrend (remove linear trend)
  vals_detrend <- residuals(lm(vals ~ seq_along(vals)))

  # FFT
  n <- length(vals_detrend)
  ft <- Mod(fft(vals_detrend)[2:(n %/% 2 + 1)])
  freqs <- (1:(n %/% 2)) / (n * bin_width)  # cycles per mm

  # Restrict to plausible MBRT spacings: 0.5-2.0 mm -> freq 0.5-2.0 /mm
  valid <- freqs >= 0.5 & freqs <= 2.0
  if (!any(valid)) {
    return(list(score = 0, period = NA_real_, profile = prof))
  }

  peak_idx <- which(valid)[which.max(ft[valid])]
  score <- ft[peak_idx] / mean(ft)  # signal-to-noise ratio
  period <- 1 / freqs[peak_idx]

  list(score = score, period = period, profile = prof)
}

# ------------------------------------------------------------------
# Helper: classify cells into peak/valley given theta, spacing, beam centers
# ------------------------------------------------------------------
classify_pv <- function(x, y, theta_deg, beam_centers, spacing) {
  rad <- theta_deg * pi / 180
  d_perp <- -sin(rad) * x + cos(rad) * y
  dist_to_peak <- sapply(d_perp, function(d) min(abs(d - beam_centers)))
  zone <- ifelse(dist_to_peak < spacing / 4, "peak",
           ifelse(dist_to_peak > spacing / 2 - spacing / 4, "valley", "transition"))
  list(d_perp = d_perp, dist_to_peak = dist_to_peak, zone = zone)
}

# ------------------------------------------------------------------
# Helper: find beam centers via contrast optimization
# ------------------------------------------------------------------
find_beam_centers <- function(d_perp, expr, spacing, n_offsets = 200) {
  d_range <- range(d_perp)
  offsets <- seq(d_range[1], d_range[1] + spacing, length.out = n_offsets)

  best_contrast <- -Inf
  best_offset   <- NA
  best_centers  <- NULL

  for (o in offsets) {
    centers <- seq(o, d_range[2] + spacing, by = spacing)
    dist <- sapply(d_perp, function(d) min(abs(d - centers)))
    pk <- which(dist < spacing / 4)
    vl <- which(dist > spacing / 2 - spacing / 4)
    if (length(pk) < 200 || length(vl) < 200) next
    contrast <- mean(expr[pk], na.rm = TRUE) - mean(expr[vl], na.rm = TRUE)
    if (contrast > best_contrast) {
      best_contrast <- contrast
      best_offset   <- o
      best_centers  <- centers
    }
  }

  # Trim centers to those with enough cells
  if (!is.null(best_centers)) {
    keep <- sapply(best_centers, function(c) sum(abs(d_perp - c) < spacing / 4) > 50)
    best_centers <- best_centers[keep]
  }

  list(contrast = best_contrast, offset = best_offset, centers = best_centers)
}

# ==================================================================
# MAIN: process each slide
# ==================================================================
slide_results <- list()

for (slide_num in 1:4) {
  cat(sprintf("\n========== SLIDE %d ==========\n", slide_num))
  rds_path <- file.path(DATA_DIR, sprintf("seurat_mutter02_slide%d_qc.rds", slide_num))
  obj <- readRDS(rds_path)

  # Extract MBRT cells
  mbrt <- subset(obj, treatment == "MBRT")
  n_mbrt <- ncol(mbrt)
  cat(sprintf("MBRT cells: %d\n", n_mbrt))

  # Extract expression
  expr_mat <- GetAssayData(mbrt, layer = "data")
  p21 <- as.numeric(expr_mat["Cdkn1a", ])

  # 3-gene signature: mean of available genes
  sig_genes_present <- SIG_GENES[SIG_GENES %in% rownames(expr_mat)]
  cat(sprintf("Signature genes present: %s\n", paste(sig_genes_present, collapse = ", ")))
  sig_score <- colMeans(expr_mat[sig_genes_present, , drop = FALSE])

  x <- mbrt$x_slide_mm
  y <- mbrt$y_slide_mm

  # ------------------------------------------------------------------
  # Step 1: Theta scan — FFT periodicity for both markers
  # ------------------------------------------------------------------
  cat("Scanning angles...\n")
  angles <- seq(-45, 45, by = 1)

  theta_results <- map_dfr(angles, function(theta) {
    res_p21 <- fft_periodicity(x, y, p21, theta)
    res_sig <- fft_periodicity(x, y, sig_score, theta)
    tibble(
      theta       = theta,
      score_p21   = res_p21$score,
      period_p21  = res_p21$period,
      score_sig   = res_sig$score,
      period_sig  = res_sig$period
    )
  })

  # Best theta from p21
  best_p21_idx <- which.max(theta_results$score_p21)
  best_theta_p21 <- theta_results$theta[best_p21_idx]
  best_period_p21 <- theta_results$period_p21[best_p21_idx]

  # Best theta from signature
  best_sig_idx <- which.max(theta_results$score_sig)
  best_theta_sig <- theta_results$theta[best_sig_idx]
  best_period_sig <- theta_results$period_sig[best_sig_idx]

  # Use whichever has stronger signal
  if (theta_results$score_sig[best_sig_idx] > theta_results$score_p21[best_p21_idx]) {
    best_theta   <- best_theta_sig
    best_spacing <- best_period_sig
    best_marker  <- "signature"
  } else {
    best_theta   <- best_theta_p21
    best_spacing <- best_period_p21
    best_marker  <- "p21"
  }

  cat(sprintf("Best theta (p21):  %d deg, period=%.3f mm, score=%.2f\n",
              best_theta_p21, best_period_p21, theta_results$score_p21[best_p21_idx]))
  cat(sprintf("Best theta (sig):  %d deg, period=%.3f mm, score=%.2f\n",
              best_theta_sig, best_period_sig, theta_results$score_sig[best_sig_idx]))
  cat(sprintf("Using: %s (theta=%d, spacing=%.3f mm)\n", best_marker, best_theta, best_spacing))

  # ------------------------------------------------------------------
  # Step 2: Fine-tune spacing via contrast optimization
  # ------------------------------------------------------------------
  cat("Optimizing beam centers...\n")
  rad <- best_theta * pi / 180
  d_perp <- -sin(rad) * x + cos(rad) * y

  # Try spacing grid near FFT estimate
  spacing_grid <- seq(max(0.6, best_spacing - 0.3),
                      min(2.0, best_spacing + 0.3), by = 0.02)

  spacing_results <- map_dfr(spacing_grid, function(s) {
    res <- find_beam_centers(d_perp, p21, s)
    tibble(spacing = s, contrast = res$contrast,
           n_peaks = length(res$centers))
  })

  opt_row <- spacing_results %>% arrange(desc(contrast)) %>% slice(1)
  opt_spacing <- opt_row$spacing
  cat(sprintf("Optimized spacing: %.3f mm (contrast=%.4f, n_peaks=%d)\n",
              opt_spacing, opt_row$contrast, opt_row$n_peaks))

  # Final beam centers with optimized spacing
  final_bc <- find_beam_centers(d_perp, p21, opt_spacing)
  beam_centers <- final_bc$centers
  cat(sprintf("Beam centers (d_perp): %s\n",
              paste(round(beam_centers, 3), collapse = ", ")))

  # ------------------------------------------------------------------
  # Step 3: Classify peak/valley
  # ------------------------------------------------------------------
  pv <- classify_pv(x, y, best_theta, beam_centers, opt_spacing)
  zone_tab <- table(pv$zone)
  cat("Zone counts:\n"); print(zone_tab)

  mean_p21_peak   <- mean(p21[pv$zone == "peak"], na.rm = TRUE)
  mean_p21_valley <- mean(p21[pv$zone == "valley"], na.rm = TRUE)
  log2fc <- log2((mean_p21_peak + 0.01) / (mean_p21_valley + 0.01))
  cat(sprintf("p21: peak=%.4f, valley=%.4f, log2FC=%.3f\n",
              mean_p21_peak, mean_p21_valley, log2fc))

  mean_sig_peak   <- mean(sig_score[pv$zone == "peak"], na.rm = TRUE)
  mean_sig_valley <- mean(sig_score[pv$zone == "valley"], na.rm = TRUE)
  sig_log2fc <- log2((mean_sig_peak + 0.01) / (mean_sig_valley + 0.01))
  cat(sprintf("sig: peak=%.4f, valley=%.4f, log2FC=%.3f\n",
              mean_sig_peak, mean_sig_valley, sig_log2fc))

  # Store results
  slide_results[[slide_num]] <- list(
    slide         = slide_num,
    n_cells       = n_mbrt,
    best_theta    = best_theta,
    opt_spacing   = opt_spacing,
    n_peaks       = length(beam_centers),
    beam_centers  = beam_centers,
    contrast      = final_bc$contrast,
    p21_peak      = mean_p21_peak,
    p21_valley    = mean_p21_valley,
    p21_log2fc    = log2fc,
    sig_peak      = mean_sig_peak,
    sig_valley    = mean_sig_valley,
    sig_log2fc    = sig_log2fc,
    best_marker   = best_marker,
    theta_results = theta_results,
    pv_zone       = pv$zone,
    d_perp        = pv$d_perp,
    dist_to_peak  = pv$dist_to_peak,
    x             = x,
    y             = y,
    p21           = p21,
    sig_score     = sig_score
  )

  # ------------------------------------------------------------------
  # Plot 1: Theta scan
  # ------------------------------------------------------------------
  p_theta <- ggplot(theta_results, aes(x = theta)) +
    geom_line(aes(y = score_p21, color = "Cdkn1a"), linewidth = 1) +
    geom_line(aes(y = score_sig, color = "3-gene sig"), linewidth = 1) +
    geom_vline(xintercept = best_theta, color = "red", linetype = "dashed") +
    geom_vline(xintercept = REF_THETA, color = "grey50", linetype = "dotted") +
    annotate("text", x = REF_THETA - 2, y = max(theta_results$score_p21) * 0.95,
             label = sprintf("Mutter_01 ref: %d°", REF_THETA),
             hjust = 1, size = 3, color = "grey40") +
    annotate("text", x = best_theta + 2, y = max(theta_results$score_p21) * 0.90,
             label = sprintf("Best: %d° (%.2f mm)", best_theta, opt_spacing),
             hjust = 0, size = 3, color = "red") +
    scale_color_manual(values = c("Cdkn1a" = "steelblue", "3-gene sig" = "darkgreen")) +
    labs(x = "Beam angle (degrees)", y = "FFT periodicity score (SNR)",
         title = sprintf("Slide %d — Stripe angle scan (48h MBRT)", slide_num),
         color = "Marker") +
    theme_bw(base_size = 12)
  ggsave(file.path(PLOT_DIR, sprintf("mutter02_slide%d_theta_scan.png", slide_num)),
         plot = p_theta, width = 10, height = 6, dpi = 150)

  # ------------------------------------------------------------------
  # Plot 2: Stripe map with beam centerlines
  # ------------------------------------------------------------------
  plot_df <- tibble(x = x, y = y, zone = pv$zone) %>%
    mutate(zone = factor(zone, levels = c("peak", "transition", "valley")))

  p_map <- ggplot(plot_df, aes(x = x, y = y, color = zone)) +
    geom_point(size = 0.01, alpha = 0.2) +
    scale_color_manual(values = c(peak = "#d73027", transition = "grey70",
                                  valley = "#4575b4"),
                       drop = FALSE) +
    coord_fixed() +
    labs(title = sprintf("Slide %d — Peak/valley zones (theta=%d°, spacing=%.2f mm, %d peaks)",
                         slide_num, best_theta, opt_spacing, length(beam_centers)),
         x = "x (mm)", y = "y (mm)", color = "Zone") +
    theme_bw(base_size = 12) +
    theme(legend.position = "right")

  # Add beam centerlines
  for (bc in beam_centers) {
    # d_perp = -sin(theta)*x + cos(theta)*y = bc
    # => y = (bc + sin(theta)*x) / cos(theta)
    slope <- sin(rad) / cos(rad)
    intercept <- bc / cos(rad)
    p_map <- p_map + geom_abline(intercept = intercept, slope = slope,
                                  color = "yellow", linewidth = 0.5, alpha = 0.8)
  }
  ggsave(file.path(PLOT_DIR, sprintf("mutter02_slide%d_stripe_map.png", slide_num)),
         plot = p_map, width = 12, height = 10, dpi = 150)

  # ------------------------------------------------------------------
  # Plot 3: p21 gradient vs distance from beam center
  # ------------------------------------------------------------------
  grad_df <- tibble(dist = pv$dist_to_peak, p21 = p21, sig = sig_score) %>%
    filter(dist <= opt_spacing / 2) %>%
    mutate(dist_bin = round(dist / 0.02) * 0.02) %>%
    group_by(dist_bin) %>%
    summarise(mean_p21 = mean(p21, na.rm = TRUE),
              mean_sig = mean(sig, na.rm = TRUE),
              n = n(), .groups = "drop") %>%
    filter(n >= 20)

  p_grad <- ggplot(grad_df, aes(x = dist_bin)) +
    geom_point(aes(y = mean_p21, color = "Cdkn1a"), size = 2) +
    geom_smooth(aes(y = mean_p21, color = "Cdkn1a"),
                method = "loess", se = FALSE, linewidth = 1.2) +
    geom_point(aes(y = mean_sig, color = "3-gene sig"), size = 2, shape = 17) +
    geom_smooth(aes(y = mean_sig, color = "3-gene sig"),
                method = "loess", se = FALSE, linewidth = 1.2, linetype = "dashed") +
    geom_vline(xintercept = opt_spacing / 4, linetype = "dotted", color = "grey50") +
    annotate("text", x = opt_spacing / 4, y = max(grad_df$mean_p21, na.rm = TRUE),
             label = "peak|transition boundary", angle = 90, vjust = -0.5,
             size = 3, color = "grey40") +
    scale_color_manual(values = c("Cdkn1a" = "steelblue", "3-gene sig" = "darkgreen")) +
    labs(x = "Distance from beam center (mm)", y = "Mean expression",
         title = sprintf("Slide %d — Expression gradient (48h MBRT)", slide_num),
         color = "Marker") +
    theme_bw(base_size = 12)
  ggsave(file.path(PLOT_DIR, sprintf("mutter02_slide%d_p21_gradient.png", slide_num)),
         plot = p_grad, width = 8, height = 6, dpi = 150)

  rm(obj, mbrt, expr_mat)
  gc(verbose = FALSE)
  cat(sprintf("Slide %d complete.\n", slide_num))
}

# ==================================================================
# Summary across all slides
# ==================================================================
cat("\n\n========== CROSS-SLIDE SUMMARY ==========\n")
summary_df <- map_dfr(slide_results, function(r) {
  tibble(
    slide       = r$slide,
    n_cells     = r$n_cells,
    theta       = r$best_theta,
    spacing_mm  = r$opt_spacing,
    n_peaks     = r$n_peaks,
    p21_peak    = r$p21_peak,
    p21_valley  = r$p21_valley,
    p21_log2fc  = r$p21_log2fc,
    sig_peak    = r$sig_peak,
    sig_valley  = r$sig_valley,
    sig_log2fc  = r$sig_log2fc,
    contrast    = r$contrast,
    best_marker = r$best_marker
  )
})

cat("\nPer-slide results:\n")
print(as.data.frame(summary_df), digits = 3)

cat(sprintf("\nMutter_01 reference: theta=%d deg, spacing=%.3f mm\n", REF_THETA, REF_SPACING))
cat(sprintf("Mutter_02 mean theta: %.1f deg (sd=%.1f)\n",
            mean(summary_df$theta), sd(summary_df$theta)))
cat(sprintf("Mutter_02 mean spacing: %.3f mm (sd=%.3f)\n",
            mean(summary_df$spacing_mm), sd(summary_df$spacing_mm)))
cat(sprintf("Mutter_02 mean p21 log2FC (peak/valley): %.3f (sd=%.3f)\n",
            mean(summary_df$p21_log2fc), sd(summary_df$p21_log2fc)))

# Write summary TSV
write_tsv(summary_df, file.path(DATA_DIR, "mutter02_stripe_summary.tsv"))
cat("\nWrote: data/mutter02_stripe_summary.tsv\n")

# ------------------------------------------------------------------
# Summary figure: 4-panel comparison
# ------------------------------------------------------------------

# Panel A: theta scan overlay
theta_all <- map_dfr(slide_results, function(r) {
  r$theta_results %>% mutate(slide = factor(r$slide))
})

pA <- ggplot(theta_all, aes(x = theta, y = score_p21, color = slide)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = REF_THETA, linetype = "dotted", color = "grey50") +
  annotate("text", x = REF_THETA, y = max(theta_all$score_p21, na.rm = TRUE) * 0.95,
           label = "Mutter_01", hjust = 1.1, size = 3, color = "grey40") +
  labs(x = "Beam angle (deg)", y = "FFT score (p21)",
       title = "A) Theta scan", color = "Slide") +
  theme_bw(base_size = 11)

# Panel B: spacing comparison
pB <- ggplot(summary_df, aes(x = factor(slide), y = spacing_mm)) +
  geom_col(fill = "steelblue", width = 0.6) +
  geom_hline(yintercept = REF_SPACING, linetype = "dashed", color = "red") +
  annotate("text", x = 0.5, y = REF_SPACING + 0.02,
           label = sprintf("Mutter_01: %.3f mm", REF_SPACING),
           hjust = 0, size = 3, color = "red") +
  labs(x = "Slide", y = "Spacing (mm)", title = "B) Beam spacing") +
  theme_bw(base_size = 11)

# Panel C: p21 peak vs valley
pv_df <- summary_df %>%
  select(slide, p21_peak, p21_valley) %>%
  pivot_longer(cols = c(p21_peak, p21_valley), names_to = "zone",
               values_to = "mean_p21") %>%
  mutate(zone = gsub("p21_", "", zone))

pC <- ggplot(pv_df, aes(x = factor(slide), y = mean_p21, fill = zone)) +
  geom_col(position = "dodge", width = 0.6) +
  scale_fill_manual(values = c(peak = "#d73027", valley = "#4575b4")) +
  labs(x = "Slide", y = "Mean Cdkn1a", title = "C) Peak vs valley p21",
       fill = "Zone") +
  theme_bw(base_size = 11)

# Panel D: log2FC effect sizes
pD <- ggplot(summary_df, aes(x = factor(slide))) +
  geom_col(aes(y = p21_log2fc, fill = "Cdkn1a"), width = 0.4,
           position = position_nudge(x = -0.2)) +
  geom_col(aes(y = sig_log2fc, fill = "3-gene sig"), width = 0.4,
           position = position_nudge(x = 0.2)) +
  geom_hline(yintercept = 0, color = "grey30") +
  scale_fill_manual(values = c("Cdkn1a" = "steelblue", "3-gene sig" = "darkgreen")) +
  labs(x = "Slide", y = "log2FC (peak/valley)",
       title = "D) Effect sizes", fill = "Marker") +
  theme_bw(base_size = 11)

combo <- (pA | pB) / (pC | pD) +
  plot_annotation(
    title = "Mutter_02: MBRT stripe detection at 48h",
    subtitle = sprintf("Reference (Mutter_01 4h): theta=%d°, spacing=%.3f mm",
                       REF_THETA, REF_SPACING),
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )
ggsave(file.path(PLOT_DIR, "mutter02_stripe_summary.png"),
       plot = combo, width = 14, height = 10, dpi = 150)

cat("\nAll plots saved. Done.\n")
