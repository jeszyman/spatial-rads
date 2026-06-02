suppressPackageStartupMessages({
  library(arrow); library(data.table); library(Matrix); library(ggplot2)
  library(ggrepel)
})

OUT <- "/home/jeszyman/repos/spatial-rads/dev/if_channel_check"
SLIDE <- "20250529_214712_S4"
MBRT_BLOCK <- "Block_21"   # MBRT 4h
SBRT_BLOCK <- "Block_17"   # SBRT 4h
set.seed(1)

peak_fovs <- c(224, 209, 201, 186, 172, 225, 229, 215, 175,
               176, 206, 181, 168, 167)

meta <- as.data.table(read_parquet("/tmp/mutter01_meta.parquet"))
counts_df <- as.data.table(read_parquet("/tmp/mutter01_counts.parquet"))
setnames(meta, "ImmuneAtlas_ImmGen_Main_cell_Types", "main_type")

# ---- Fit stripe geometry from MBRT Block_21 peak FOV centroids ----
mbrt_meta <- meta[Slide == SLIDE & Block == MBRT_BLOCK]
peak_centroids <- mbrt_meta[fov %in% peak_fovs,
  .(cx = mean(x_slide_mm), cy = mean(y_slide_mm)), by = fov]
search_theta <- function(cx, cy, n = 4,
                         grid = seq(-pi/2, pi/2, length.out = 181)) {
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
theta <- fit$theta
beam_spacing <- median(diff(fit$centers))
peak_centroids[, d_perp := -sin(theta) * cx + cos(theta) * cy]
peak_centroids[, stripe := fit$cluster]
within_sd <- peak_centroids[, .(sd = sd(d_perp - fit$centers[stripe])),
                            by = stripe][, mean(sd)]
peak_half <- within_sd + 0.15
cat(sprintf("MBRT stripe: theta=%.1f deg, spacing=%.3f, peak_half=%.3f\n",
            theta * 180/pi, beam_spacing, peak_half))

label_cells <- function(pos, fit_centers, theta, beam_spacing, peak_half) {
  d <- -sin(theta) * pos[,1] + cos(theta) * pos[,2]
  d_range <- range(d)
  all_centers <- seq(fit_centers[1] -
                       ceiling((fit_centers[1] - d_range[1]) / beam_spacing) * beam_spacing,
                     fit_centers[length(fit_centers)] +
                       ceiling((d_range[2] - fit_centers[length(fit_centers)]) / beam_spacing) * beam_spacing,
                     by = beam_spacing)
  dist_to_peak <- sapply(d, function(x) min(abs(x - all_centers)))
  ifelse(dist_to_peak < peak_half, "peak",
         ifelse(dist_to_peak > (beam_spacing/2 - peak_half), "valley", "transition"))
}

# ---- For each block: label tumor cells, compute per-gene peak/valley mean ----
compute_pv <- function(block_label, block_name) {
  m <- meta[Slide == SLIDE & Block == block_name & main_type == "a"]
  c_df <- counts_df[cell_id %in% m$cell_id]
  setkeyv(m, "cell_id"); setkeyv(c_df, "cell_id")
  m <- m[cell_id %in% c_df$cell_id]; setorderv(m, "cell_id"); setorderv(c_df, "cell_id")
  cat(block_label, block_name, "tumor cells:", nrow(m), "\n")

  pos <- as.matrix(m[, .(x_slide_mm, y_slide_mm)])
  # Note: SBRT uses MBRT-fitted stripe geometry (null null; arbitrary pattern)
  pv <- label_cells(pos, fit$centers, theta, beam_spacing, peak_half)
  m[, pv := pv]

  gene_cols <- setdiff(colnames(c_df), c("Slide","fov","cell_id"))
  mat <- t(as.matrix(c_df[, ..gene_cols]))
  lib <- colSums(mat); lib[lib == 0] <- 1
  normed <- log2(t(t(mat) / lib) * 1e4 + 1)

  peak_idx <- which(pv == "peak"); valley_idx <- which(pv == "valley")
  cat("  peak:", length(peak_idx), "valley:", length(valley_idx), "\n")
  peak_mean <- rowMeans(normed[, peak_idx, drop = FALSE])
  valley_mean <- rowMeans(normed[, valley_idx, drop = FALSE])
  pct_expr <- rowSums(mat > 0) / ncol(mat)
  # Wilcoxon per gene (sampled for speed)
  n_sample <- min(length(peak_idx), length(valley_idx), 3000)
  p_sample <- sample(peak_idx, n_sample); v_sample <- sample(valley_idx, n_sample)
  wilcox_p <- sapply(seq_len(nrow(normed)), function(g) {
    if (pct_expr[g] < 0.02) return(NA_real_)
    wilcox.test(normed[g, p_sample], normed[g, v_sample])$p.value
  })
  data.table(gene = rownames(normed),
             peak_mean = peak_mean,
             valley_mean = valley_mean,
             peak_minus_valley = peak_mean - valley_mean,
             pct_expr = pct_expr,
             wilcox_p = wilcox_p,
             arm = block_label)
}

mbrt <- compute_pv("MBRT", MBRT_BLOCK)
sbrt <- compute_pv("SBRT", SBRT_BLOCK)
merged <- merge(mbrt, sbrt, by = "gene", suffixes = c("_mbrt", "_sbrt"))
merged[, delta := peak_minus_valley_mbrt - peak_minus_valley_sbrt]
setorder(merged, -delta)
cat("\nTop 15 peak-enriched MBRT-specific (MBRT − SBRT peak-valley delta):\n")
print(merged[pct_expr_mbrt > 0.05 & pct_expr_sbrt > 0.05,
             .(gene, MBRT_PV = round(peak_minus_valley_mbrt, 3),
               SBRT_PV = round(peak_minus_valley_sbrt, 3),
               delta = round(delta, 3),
               mbrt_p = sprintf("%.1e", wilcox_p_mbrt))][1:15])
cat("\nTop 15 valley-enriched MBRT-specific:\n")
print(merged[pct_expr_mbrt > 0.05 & pct_expr_sbrt > 0.05][
  order(delta)][1:15,
    .(gene, MBRT_PV = round(peak_minus_valley_mbrt, 3),
      SBRT_PV = round(peak_minus_valley_sbrt, 3),
      delta = round(delta, 3),
      mbrt_p = sprintf("%.1e", wilcox_p_mbrt))])

fwrite(merged, "/tmp/tumor_pv_mbrt_vs_sbrt_null.tsv", sep = "\t")
cat("Saved /tmp/tumor_pv_mbrt_vs_sbrt_null.tsv\n")

# ---- Plot ----
keep <- merged[pct_expr_mbrt > 0.05 & pct_expr_sbrt > 0.05]
# Null distribution: SBRT peak-valley
sd_null <- sd(keep$peak_minus_valley_sbrt)
keep[, mbrt_z := peak_minus_valley_mbrt / sd_null]
keep[, is_hit := abs(delta) > 2 * sd_null & wilcox_p_mbrt < 1e-4]
label_genes <- keep[is_hit == TRUE, ][order(-abs(delta))][1:15, gene]
keep[, label := ifelse(gene %in% label_genes, gene, NA)]

p1 <- ggplot(keep, aes(x = peak_minus_valley_sbrt,
                       y = peak_minus_valley_mbrt,
                       color = is_hit)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey70") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = c(-2, 2) * sd_null,
             linetype = "dotted", color = "grey50") +
  geom_point(alpha = 0.6, size = 0.8) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 30,
                  color = "black") +
  scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "grey60")) +
  labs(x = "SBRT tumor peak − valley (mean log-norm)",
       y = "MBRT tumor peak − valley (mean log-norm)",
       title = "MBRT 4h vs SBRT 4h tumor peak/valley delta (Block_21 vs Block_17)",
       subtitle = sprintf("Same S4 slide, same extrapolated stripe geometry applied. SBRT serves as geometry null (sd=%.3f). Hits: |delta| > 2sd AND wilcox p<1e-4.", sd_null)) +
  theme_minimal() +
  theme(legend.position = "none", plot.title = element_text(size = 11),
        plot.subtitle = element_text(size = 9))

ggsave(file.path(OUT, "tumor_pv_mbrt_vs_sbrt_null.png"), p1,
       width = 10, height = 8, dpi = 130)
cat("Saved tumor_pv_mbrt_vs_sbrt_null.png\n")
