suppressPackageStartupMessages({
  library(arrow); library(data.table); library(Matrix); library(ggplot2)
  library(ggrepel); library(patchwork)
})

OUT <- "/home/jeszyman/repos/spatial-rads/dev/if_channel_check"
SLIDE <- "20250529_214712_S4"
MBRT_BLOCK <- "Block_21"
set.seed(1)

peak_fovs <- c(224, 209, 201, 186, 172, 225, 229, 215, 175,
               176, 206, 181, 168, 167)

meta <- as.data.table(read_parquet("/tmp/mutter01_meta.parquet"))
counts_df <- as.data.table(read_parquet("/tmp/mutter01_counts.parquet"))
setnames(meta, "ImmuneAtlas_ImmGen_Main_cell_Types", "main_type")

# ---- Load MBRT Block_21 tumor cells ----
m <- meta[Slide == SLIDE & Block == MBRT_BLOCK & main_type == "a"]
c_df <- counts_df[cell_id %in% m$cell_id]
setkeyv(m, "cell_id"); setkeyv(c_df, "cell_id")
m <- m[cell_id %in% c_df$cell_id]
setorderv(m, "cell_id"); setorderv(c_df, "cell_id")
cat("MBRT Block_21 tumor cells:", nrow(m), "\n")

gene_cols <- setdiff(colnames(c_df), c("Slide","fov","cell_id"))
mat <- t(as.matrix(c_df[, ..gene_cols]))
lib <- colSums(mat); lib[lib == 0] <- 1
normed <- log2(t(t(mat) / lib) * 1e4 + 1)
pct_expr <- rowSums(mat > 0) / ncol(mat)
keep_genes <- pct_expr > 0.05
normed <- normed[keep_genes, ]
cat("Genes kept:", nrow(normed), "\n")

pos <- as.matrix(m[, .(x_slide_mm, y_slide_mm)])

# ---- Fit real stripe geometry ----
peak_centroids <- meta[Slide == SLIDE & Block == MBRT_BLOCK & fov %in% peak_fovs,
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
theta_real <- fit$theta
beam_spacing <- median(diff(fit$centers))
peak_centroids[, d_perp := -sin(theta_real) * cx + cos(theta_real) * cy]
peak_centroids[, stripe := fit$cluster]
within_sd <- peak_centroids[, .(sd = sd(d_perp - fit$centers[stripe])),
                            by = stripe][, mean(sd)]
peak_half <- within_sd + 0.15
cat(sprintf("Real stripe: theta=%.1f deg, spacing=%.3f, peak_half=%.3f\n",
            theta_real * 180/pi, beam_spacing, peak_half))

# ---- Function: given theta, phase offset, compute peak/valley labels on pos ----
label_cells <- function(pos, theta, phase_offset, spacing, half_w) {
  d <- -sin(theta) * pos[, 1] + cos(theta) * pos[, 2]
  # Pick stripe centers densely across d range
  d_range <- range(d)
  centers_seq <- seq(d_range[1] - spacing, d_range[2] + spacing, by = spacing)
  # Shift by phase_offset (fraction of spacing)
  centers_seq <- centers_seq + phase_offset
  dist_to_peak <- sapply(d, function(x) min(abs(x - centers_seq)))
  ifelse(dist_to_peak < half_w, "peak",
         ifelse(dist_to_peak > (spacing / 2 - half_w), "valley", "transition"))
}

# ---- Function: compute per-gene peak-valley mean delta given labels ----
compute_delta <- function(labels, normed) {
  peak_idx <- which(labels == "peak")
  valley_idx <- which(labels == "valley")
  if (length(peak_idx) < 50 || length(valley_idx) < 50) return(NULL)
  rowMeans(normed[, peak_idx, drop = FALSE]) -
    rowMeans(normed[, valley_idx, drop = FALSE])
}

# ---- Real MBRT labels ----
# Anchor phase so that fit$centers[1] lands on a peak
phase_real <- fit$centers[1]  # first stripe center = phase reference
d_cells_real <- -sin(theta_real) * pos[,1] + cos(theta_real) * pos[,2]
d_range_real <- range(d_cells_real)
# Extrapolate centers through full range
centers_real <- seq(fit$centers[1] -
                      ceiling((fit$centers[1] - d_range_real[1]) / beam_spacing) * beam_spacing,
                    fit$centers[length(fit$centers)] +
                      ceiling((d_range_real[2] - fit$centers[length(fit$centers)]) / beam_spacing) * beam_spacing,
                    by = beam_spacing)
real_labels <- ifelse(
  sapply(d_cells_real, function(x) min(abs(x - centers_real))) < peak_half, "peak",
  ifelse(sapply(d_cells_real, function(x) min(abs(x - centers_real))) >
           (beam_spacing / 2 - peak_half), "valley", "transition")
)
cat("Real label counts:\n"); print(table(real_labels))

real_delta <- compute_delta(real_labels, normed)

# ---- Rotation null: sample many angles, keep spacing and half-width identical ----
n_null <- 200
# Exclude angles within ±30 deg of real (stripes still partially overlap within that range).
# Generate more candidates so enough pass the exclusion.
theta_raw <- runif(n_null * 3, -pi/2, pi/2)
min_sep <- pi/6   # 30 deg minimum separation from real angle
theta_nulls <- theta_raw[
  pmin(abs(theta_raw - theta_real),
       abs(theta_raw - (theta_real + pi)),
       abs(theta_raw - (theta_real - pi))) > min_sep
]
theta_nulls <- head(theta_nulls, n_null)
# Diagnostic: report concordance of null labels with real labels per rotation
concord <- numeric(length(theta_nulls))
for (i in seq_along(theta_nulls)) {
  phase <- runif(1, 0, beam_spacing)
  labs <- label_cells(pos, theta_nulls[i], phase, beam_spacing, peak_half)
  # Concordance ignoring transition cells in either set
  both <- !grepl("transition", labs) & !grepl("transition", real_labels)
  concord[i] <- mean(labs[both] == real_labels[both])
}
cat(sprintf("Null rotations: %d iterations, min sep %.0f deg from real. Label concordance range: %.2f-%.2f (mean %.2f)\n",
            length(theta_nulls), min_sep * 180/pi,
            min(concord), max(concord), mean(concord)))

null_deltas <- matrix(NA_real_, nrow = length(real_delta), ncol = length(theta_nulls))
rownames(null_deltas) <- names(real_delta)
for (i in seq_along(theta_nulls)) {
  # Random phase
  phase <- runif(1, 0, beam_spacing)
  labels <- label_cells(pos, theta_nulls[i], phase, beam_spacing, peak_half)
  d <- compute_delta(labels, normed)
  if (!is.null(d)) null_deltas[, i] <- d
  if (i %% 25 == 0) cat("  null iter", i, "\n")
}

# Per-gene empirical p-value: fraction of null deltas with |delta| >= |real|
emp_p <- sapply(seq_along(real_delta), function(g) {
  nd <- null_deltas[g, ]; nd <- nd[!is.na(nd)]
  if (length(nd) < 10) return(NA_real_)
  sum(abs(nd) >= abs(real_delta[g])) / length(nd)
})
names(emp_p) <- names(real_delta)

# Also pull rotation_null sd per gene and overall
null_sd <- apply(null_deltas, 1, sd, na.rm = TRUE)
overall_null_sd <- sd(as.vector(null_deltas), na.rm = TRUE)
cat(sprintf("Rotation-null overall sd: %.3f (vs SBRT null sd: ~0.016)\n",
            overall_null_sd))

# ---- Join with SBRT null from prior analysis ----
sbrt <- fread("/tmp/tumor_pv_mbrt_vs_sbrt_null.tsv")
combined <- data.table(gene = names(real_delta),
                       mbrt_delta = real_delta,
                       rot_null_sd = null_sd,
                       rot_null_p = emp_p)
combined <- merge(combined, sbrt[, .(gene, sbrt_delta = peak_minus_valley_sbrt,
                                      pct_expr_mbrt, pct_expr_sbrt)],
                  by = "gene", all.x = TRUE)
combined[, survives_rot := rot_null_p < 0.05 & abs(mbrt_delta) > 2 * overall_null_sd]
setorder(combined, -mbrt_delta)
fwrite(combined, "/tmp/tumor_pv_rotation_null.tsv", sep = "\t")

cat("\n=== Top 20 peak-enriched with rotation-null validation ===\n")
print(combined[pct_expr_mbrt > 0.05][1:20,
  .(gene, mbrt_delta = round(mbrt_delta, 3),
    sbrt_delta = round(sbrt_delta, 3),
    rot_null_p = sprintf("%.3f", rot_null_p),
    survives = survives_rot)])
cat(sprintf("\nTotal hits surviving rotation null (p<0.05 AND |delta|>2*null sd): %d\n",
            sum(combined$survives_rot, na.rm = TRUE)))

# ---- Plot: 3-curve density ----
densdata <- rbindlist(list(
  data.table(arm = "MBRT real (θ=−25°)", val = real_delta),
  data.table(arm = "MBRT rotation null (200 random θ)",
             val = as.vector(null_deltas)),
  data.table(arm = "SBRT null (different tumor)",
             val = combined$sbrt_delta[!is.na(combined$sbrt_delta)])
))
densdata[, arm := factor(arm, levels = c(
  "SBRT null (different tumor)",
  "MBRT rotation null (200 random θ)",
  "MBRT real (θ=−25°)"))]

p_dens <- ggplot(densdata, aes(x = val, fill = arm, color = arm)) +
  geom_density(alpha = 0.3, linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  scale_fill_manual(values = c(
    "SBRT null (different tumor)" = "grey50",
    "MBRT rotation null (200 random θ)" = "darkgreen",
    "MBRT real (θ=−25°)" = "#e41a1c")) +
  scale_color_manual(values = c(
    "SBRT null (different tumor)" = "grey30",
    "MBRT rotation null (200 random θ)" = "darkgreen",
    "MBRT real (θ=−25°)" = "#b30000")) +
  labs(x = "Peak − Valley (log-norm expression)", y = "Density of genes",
       fill = NULL, color = NULL,
       title = "Gene-level peak−valley distribution under three null models",
       subtitle = paste0("Rotation null is the right control: MBRT tumor, real stripe spacing + half-width, but random stripe angle.\n",
                         "If MBRT peak-valley is real biology (dose-driven), it should exceed the rotation null.")) +
  theme_minimal() +
  theme(legend.position = "bottom", plot.title = element_text(size = 11),
        plot.subtitle = element_text(size = 9))

# Scatter: real MBRT delta vs rotation-null sd per gene
scat <- combined[!is.na(rot_null_sd) & pct_expr_mbrt > 0.05]
top_peak <- scat[order(-mbrt_delta)][1:12, gene]
top_valley <- scat[order(mbrt_delta)][1:5, gene]
scat[, label := ifelse(gene %in% c(top_peak, top_valley), gene, NA_character_)]
scat[, category := fifelse(survives_rot == TRUE,
                           fifelse(mbrt_delta > 0, "survives (peak)", "survives (valley)"),
                           "does not survive")]

# Volcano-style: real MBRT delta on x, -log10(empirical rotation-null p) on y
scat[, neg_log_p := -log10(pmax(rot_null_p, 1/200))]
p_scat <- ggplot(scat, aes(x = mbrt_delta, y = neg_log_p,
                           color = category)) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  annotate("text", x = max(scat$mbrt_delta) * 0.95,
           y = -log10(0.05) + 0.08,
           label = "rotation-null p = 0.05", hjust = 1, size = 3, color = "grey40") +
  geom_point(alpha = 0.75, size = 1.1) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 50,
                  color = "black", min.segment.length = 0,
                  box.padding = 0.3, seed = 1) +
  scale_color_manual(values = c("survives (peak)" = "#e41a1c",
                                "survives (valley)" = "#377eb8",
                                "does not survive" = "grey60")) +
  labs(x = "Real MBRT peak − valley (log-norm) in tumor cells at 4h",
       y = "−log10 (rotation-null p-value)   higher = more surprising",
       color = NULL,
       title = "Volcano: real MBRT peak/valley effect vs rotation-null significance",
       subtitle = paste0(
         "x = observed peak-valley difference under real MBRT stripes. Right = peak-enriched, left = valley-enriched.\n",
         "y = how unlikely this effect is under the rotation null (higher = more significant). Dashed = p=0.05 cutoff.\n",
         "Red dots (upper-right) = real peak enrichment. Grey dots near 0 = no signal or architecture artifact.")) +
  theme_minimal() +
  theme(legend.position = "bottom", plot.title = element_text(size = 11),
        plot.subtitle = element_text(size = 9))

combined_plot <- p_dens / p_scat + plot_layout(heights = c(1, 1.4))
ggsave(file.path(OUT, "tumor_pv_rotation_null.png"),
       combined_plot, width = 12, height = 13, dpi = 130)
cat("Saved tumor_pv_rotation_null.png\n")
