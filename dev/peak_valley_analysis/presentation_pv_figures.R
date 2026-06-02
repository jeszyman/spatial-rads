suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(ggplot2)
  library(dplyr)
})

INPUT_DIR  <- "/mnt/data/projects/spatial-rads/inputs/mutter01"
PLOT_DIR   <- file.path(dirname(normalizePath(".")), "peak_valley_analysis", "plots")
if (!dir.exists(PLOT_DIR)) dir.create(PLOT_DIR, recursive = TRUE)

BLOCK <- "Block_21"
SLIDE <- "20250529_214712_S4"
set.seed(1)

peak_fovs   <- c(224, 209, 201, 186, 172, 225, 229, 215, 175,
                  176, 206, 181, 168, 167)
valley_fovs <- c(199, 213, 227, 188, 204, 218, 178, 165, 170,
                  231, 189, 173, 154)

# --- Load data ---
meta <- as.data.table(read_parquet(
  file.path(INPUT_DIR, "Analysis_Mutter_01_CosMmR_Mutter_updated_metadata.parquet")))
counts_full <- as.data.table(read_parquet(
  file.path(INPUT_DIR, "projects_Mutter_01_CosMmR_app_data_raw_tx_counts_matrix.parquet")))

m <- meta[Slide == SLIDE & Block == BLOCK & qcFlagsCell == "Pass"]
counts_b21 <- counts_full[cell_id %in% m$cell_id]
setkeyv(m, "cell_id"); setkeyv(counts_b21, "cell_id")
m <- m[cell_id %in% counts_b21$cell_id]
cat("Block 21 cells:", nrow(m), "\n")

# ===========================================================
# STRIPE MODEL FIT (from manually called peak FOVs)
# ===========================================================
peak_cells <- m[fov %in% peak_fovs]
peak_fov_centroids <- peak_cells[, .(cx = mean(x_slide_mm),
                                     cy = mean(y_slide_mm)), by = fov]

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

fit <- search_theta(peak_fov_centroids$cx, peak_fov_centroids$cy)
theta <- fit$theta
stripe_centers <- sort(fit$centers)
beam_spacing <- median(diff(stripe_centers))

peak_fov_centroids[, d_perp := -sin(theta) * cx + cos(theta) * cy]
peak_fov_centroids[, stripe := fit$cluster]
within_stripe_sd <- peak_fov_centroids[,
  .(sd = sd(d_perp - stripe_centers[stripe])), by = stripe][, mean(sd)]
peak_half <- within_stripe_sd + 0.15

d_all <- -sin(theta) * m$x_slide_mm + cos(theta) * m$y_slide_mm
d_range <- range(d_all)
all_centers <- seq(
  stripe_centers[1] - ceiling((stripe_centers[1] - d_range[1]) / beam_spacing) * beam_spacing,
  stripe_centers[length(stripe_centers)] + ceiling((d_range[2] - stripe_centers[length(stripe_centers)]) / beam_spacing) * beam_spacing,
  by = beam_spacing)

dist_to_nearest <- sapply(d_all, function(d) min(abs(d - all_centers)))
label <- rep("transition", length(d_all))
label[dist_to_nearest < peak_half] <- "peak"
label[dist_to_nearest > (beam_spacing / 2 - peak_half)] <- "valley"

m[, pv_extrap := label]
m[, d_to_peak := dist_to_nearest]
m[, pv_manual := fifelse(fov %in% peak_fovs, "peak",
                         fifelse(fov %in% valley_fovs, "valley", NA_character_))]

cat(sprintf("Stripe angle: %.1f deg | spacing: %.3f mm | peak half: %.3f mm\n",
            theta * 180 / pi, beam_spacing, peak_half))
cat("Zone counts:\n"); print(table(label))

# --- Beam centerlines for overlay ---
x_lim <- range(m$x_slide_mm)
beam_lines <- do.call(rbind, lapply(all_centers, function(c) {
  data.frame(
    x = x_lim,
    y = (c + sin(theta) * x_lim) / cos(theta),
    center = c)
}))

zone_colors <- c("peak" = "#e41a1c", "valley" = "#377eb8", "transition" = "grey80")

# --- Consistent axis breaks for spatial plots ---
x_range <- range(m$x_slide_mm)
y_range <- range(m$y_slide_mm)
tick_interval <- 1  # 1 mm ticks
x_breaks <- seq(floor(x_range[1]), ceiling(x_range[2]), by = tick_interval)
y_breaks <- seq(floor(y_range[1]), ceiling(y_range[2]), by = tick_interval)

# ===========================================================
# FIGURE 1: Manual FOV labels
# ===========================================================
m[, manual_col := fifelse(is.na(pv_manual), "unlabeled", pv_manual)]
manual_colors <- c("peak" = "#e41a1c", "valley" = "#377eb8", "unlabeled" = "grey85")

p1 <- ggplot(m, aes(x = x_slide_mm, y = y_slide_mm, color = manual_col)) +
  geom_point(size = 0.05, alpha = 0.3) +
  geom_line(data = beam_lines, aes(x = x, y = y, group = center),
            color = "black", linetype = "dashed", linewidth = 0.5,
            inherit.aes = FALSE) +
  scale_color_manual(values = manual_colors,
                     labels = c("Peak (manual)", "Valley (manual)", "Unlabeled"),
                     name = NULL) +
  scale_x_continuous(breaks = x_breaks) +
  scale_y_continuous(breaks = y_breaks) +
  coord_fixed() +
  labs(x = "x (mm)", y = "y (mm)",
       title = "Manual peak/valley FOV classification from gamma-H2AX IHC") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

ggsave(file.path(PLOT_DIR, "pres_manual_fov_labels.png"), p1,
       width = 8, height = 10, dpi = 200)
cat("Saved: pres_manual_fov_labels.png\n")

# ===========================================================
# FIGURE 2: Extrapolated labels
# ===========================================================
m[, pv_extrap_f := factor(pv_extrap, levels = c("peak", "valley", "transition"))]

p2 <- ggplot(m, aes(x = x_slide_mm, y = y_slide_mm, color = pv_extrap_f)) +
  geom_point(size = 0.05, alpha = 0.3) +
  geom_line(data = beam_lines, aes(x = x, y = y, group = center),
            color = "black", linetype = "dashed", linewidth = 0.5,
            inherit.aes = FALSE) +
  scale_color_manual(values = zone_colors, name = NULL) +
  scale_x_continuous(breaks = x_breaks) +
  scale_y_continuous(breaks = y_breaks) +
  coord_fixed() +
  labs(x = "x (mm)", y = "y (mm)",
       title = sprintf("Extrapolated peak/valley (%.0f-deg tilt, %.2f mm spacing)",
                       theta * 180 / pi, beam_spacing)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

ggsave(file.path(PLOT_DIR, "pres_extrapolated_zones.png"), p2,
       width = 8, height = 10, dpi = 200)
cat("Saved: pres_extrapolated_zones.png\n")

# ===========================================================
# FIGURE 3: Validation histogram
# ===========================================================
hist_data <- m[!is.na(pv_manual), .(d_to_peak, pv_manual)]

p3 <- ggplot(hist_data, aes(x = d_to_peak, fill = pv_manual)) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
  geom_vline(xintercept = peak_half, linetype = "dashed",
             color = "#e41a1c", linewidth = 0.8) +
  geom_vline(xintercept = beam_spacing / 2 - peak_half, linetype = "dashed",
             color = "#377eb8", linewidth = 0.8) +
  scale_fill_manual(values = c("peak" = "#e41a1c", "valley" = "#377eb8"),
                    labels = c("Manual peak FOVs", "Manual valley FOVs"),
                    name = NULL) +
  labs(x = "Distance to nearest beam centerline (mm)",
       y = "Cell count",
       title = "Stripe model validation",
       subtitle = "Manually called peak/valley FOVs separate cleanly by modeled beam distance") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(PLOT_DIR, "pres_validation_histogram.png"), p3,
       width = 8, height = 5, dpi = 200)
cat("Saved: pres_validation_histogram.png\n")

# ===========================================================
# FIGURE 4: p21 % expressing kinetics (all conditions + PV overlay)
# ===========================================================
design <- as.data.table(readxl::read_excel(
  "/home/jeszyman/repos/spatial-rads/data/metadata.xlsx"))
design <- design[, .(slide, block_id, condition, model, treatment, timepoint_h)]

meta_all <- merge(meta[qcFlagsCell == "Pass"], design,
                  by.x = c("Slide", "Block"), by.y = c("slide", "block_id"))
meta_all <- meta_all[is.na(model) | model == "flank"]

counts_all <- counts_full[cell_id %in% meta_all$cell_id]
setkeyv(meta_all, "cell_id"); setkeyv(counts_all, "cell_id")
meta_all <- meta_all[cell_id %in% counts_all$cell_id]

lib <- rowSums(counts_all[, !c("Slide", "fov", "cell_id"), with = FALSE])
meta_all[, Cdkn1a := log2((counts_all$Cdkn1a / pmax(lib, 1)) * 1e4 + 1)]

classify_cells <- function(x) {
  out <- rep(NA_character_, length(x))
  out[x == "a"] <- "Tumor (4T1)"
  out[x %in% c("B.cell", "Memory.B", "Plasma", "Plasmablast",
               "GC_centroblasts", "GC_centrocyes", "Spleen.CD19")] <- "B cells"
  out[x %in% c("CD8.T.cell", "Spleen.Naive.CD4", "Spleen.Naive.CD8",
               "Spleen.CD4Act.48hrs", "Spleen.LN.Naive.CD4",
               "Spleen.Treg", "Colon.Treg.Nrplo")] <- "T cells"
  out[x %in% c("Macrophage", "PerC.macrophage", "spleen.red.pulp.macs",
               "small.peritoneal.macs", "Microglia")] <- "Macrophages"
  out[x %in% c("Ly6Clo.blood.monocytes", "Ly6Chi.blood.monocytes",
               "Dendritic")] <- "Monocytes/DC"
  out[x %in% c("Fibroblastic.reticular", "Pericyte")] <- "Stroma"
  out
}

meta_all[, cell_class := classify_cells(ImmuneAtlas_ImmGen_Main_cell_Types)]
meta_all <- meta_all[!is.na(cell_class)]
meta_all[, cell_class := factor(cell_class,
  levels = c("B cells", "T cells", "Monocytes/DC", "Macrophages", "Tumor (4T1)", "Stroma"))]

meta_all[, treat := fifelse(is.na(treatment) | treatment == "NT", "Control", treatment)]
meta_all[condition == "Control", treat := "Control"]

# Kinetics summary
summ <- meta_all[, .(
  n_cells = .N,
  pct_expr = 100 * mean(Cdkn1a > 0, na.rm = TRUE)
), by = .(cell_class, treat, timepoint_h, condition)]

# Peak/valley overlay at 4h using EXTRAPOLATED labels
mbrt4h_all <- meta_all[condition == "MBRT_4h"]
d_4h <- -sin(theta) * mbrt4h_all$x_slide_mm + cos(theta) * mbrt4h_all$y_slide_mm
dist_4h <- sapply(d_4h, function(d) min(abs(d - all_centers)))
pv_label_4h <- rep(NA_character_, nrow(mbrt4h_all))
pv_label_4h[dist_4h < peak_half] <- "MBRT peak"
pv_label_4h[dist_4h > (beam_spacing / 2 - peak_half)] <- "MBRT valley"
mbrt4h_all[, pv_label := pv_label_4h]

pv_summ <- mbrt4h_all[!is.na(pv_label), .(
  n_cells = .N,
  pct_expr = 100 * mean(Cdkn1a > 0, na.rm = TRUE),
  treat = pv_label,
  timepoint_h = 4,
  condition = pv_label
), by = .(cell_class, pv_label)]
pv_summ[, pv_label := NULL]
summ <- rbind(summ, pv_summ, fill = TRUE)

tp_map <- c("0" = "0 Gy", "1" = "1h", "4" = "4h", "48" = "2d", "144" = "6d")
tp_levels <- c("0 Gy", "1h", "4h", "2d", "6d")
summ[, tp := factor(tp_map[as.character(timepoint_h)], levels = tp_levels)]

treat_levels <- c("Control", "MBRT", "SBRT", "MBRT peak", "MBRT valley")
summ[, treat := factor(treat, levels = treat_levels)]
treat_colors <- c("Control" = "gray40", "MBRT" = "#e41a1c",
                  "SBRT" = "#377eb8", "MBRT peak" = "#b30000",
                  "MBRT valley" = "#fc9272")
treat_shapes <- c("Control" = 16, "MBRT" = 17, "SBRT" = 15,
                  "MBRT peak" = 8, "MBRT valley" = 4)

line_series <- summ[treat %in% c("Control", "MBRT", "SBRT")]
pv_points   <- summ[treat %in% c("MBRT peak", "MBRT valley")]

p4 <- ggplot() +
  geom_line(data = line_series,
            aes(x = tp, y = pct_expr, color = treat, group = treat),
            linewidth = 0.6) +
  geom_point(data = line_series,
             aes(x = tp, y = pct_expr, color = treat, shape = treat),
             size = 3) +
  geom_point(data = pv_points,
             aes(x = tp, y = pct_expr, color = treat, shape = treat),
             size = 4, stroke = 1.2) +
  facet_wrap(~ cell_class, scales = "free_y") +
  scale_color_manual(values = treat_colors) +
  scale_shape_manual(values = treat_shapes) +
  labs(x = "Timepoint",
       y = "% cells expressing Cdkn1a (p21)",
       title = "p21 expression kinetics by cell type",
       subtitle = "Mutter_01 flank tumors (n = 1 per condition)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 11, face = "bold"))

ggsave(file.path(PLOT_DIR, "pres_p21_pctexpr_kinetics.png"), p4,
       width = 11, height = 7, dpi = 200)
cat("Saved: pres_p21_pctexpr_kinetics.png\n")

# ===========================================================
# FIGURE 5: p21 gradient vs absolute distance from beam center
# ===========================================================
# Use MBRT_4h Block_21 cells with extrapolated distance
mbrt4h_grad <- meta[Slide == SLIDE & Block == BLOCK & qcFlagsCell == "Pass"]
mbrt4h_grad <- mbrt4h_grad[cell_id %in% counts_full$cell_id]
setkeyv(mbrt4h_grad, "cell_id")
counts_grad <- counts_full[cell_id %in% mbrt4h_grad$cell_id]
setkeyv(counts_grad, "cell_id")
counts_grad <- counts_grad[mbrt4h_grad$cell_id]

lib_grad <- rowSums(counts_grad[, !c("Slide", "fov", "cell_id"), with = FALSE])
mbrt4h_grad[, Cdkn1a := log2((counts_grad$Cdkn1a / pmax(lib_grad, 1)) * 1e4 + 1)]
mbrt4h_grad[, Gadd45b := log2((counts_grad$Gadd45b / pmax(lib_grad, 1)) * 1e4 + 1)]
mbrt4h_grad[, Bax := log2((counts_grad$Bax / pmax(lib_grad, 1)) * 1e4 + 1)]

d_grad <- -sin(theta) * mbrt4h_grad$x_slide_mm + cos(theta) * mbrt4h_grad$y_slide_mm
mbrt4h_grad[, d_to_peak := sapply(d_grad, function(d) min(abs(d - all_centers)))]

mbrt4h_grad[, cell_class := classify_cells(ImmuneAtlas_ImmGen_Main_cell_Types)]

# Radiosensitivity ordering (most sensitive first)
rs_levels <- c("All cells", "B cells", "T cells", "Monocytes/DC",
               "Macrophages", "Stroma", "Tumor (4T1)")

grad_typed <- mbrt4h_grad[!is.na(cell_class)]
grad_all <- copy(grad_typed)
grad_all[, cell_class := "All cells"]
grad_combined <- rbind(grad_typed, grad_all)
grad_combined[, cell_class := factor(cell_class, levels = rs_levels)]

MIN_CELLS_PER_FOV <- 10
MIN_CELLS_PER_BIN <- 20
bin_width <- 0.02

mbrt4h_grad[, pv_zone := fifelse(d_to_peak < peak_half, "peak",
  fifelse(d_to_peak > (beam_spacing / 2 - peak_half), "valley", NA_character_))]

# --- Helper: gradient + FOV Wilcoxon for any gene ---
make_gene_plots <- function(dt, gene_col, gene_name, gene_label, plot_dir) {
  gc <- copy(dt[!is.na(cell_class)])
  gc[, expr := get(gene_col)]
  gc_all <- copy(gc)
  gc_all[, cell_class := "All cells"]
  gc_comb <- rbind(gc, gc_all)
  gc_comb[, cell_class := factor(cell_class, levels = rs_levels)]

  gc_comb[, d_bin := floor(d_to_peak / bin_width) * bin_width + bin_width / 2]
  binned <- gc_comb[, .(mean_expr = mean(expr), n_cells = .N),
                    by = .(cell_class, d_bin)]
  binned <- binned[n_cells >= MIN_CELLS_PER_BIN]

  # Drop facets with <5 bins (loess needs enough points)
  bin_counts <- binned[, .N, by = cell_class]
  keep_classes <- bin_counts[N >= 5, cell_class]
  binned <- binned[cell_class %in% keep_classes]

  # Gradient plot
  pg <- ggplot(binned, aes(x = d_bin, y = mean_expr)) +
    geom_point(size = 0.6, alpha = 0.5, color = "grey40") +
    geom_smooth(method = "loess", span = 0.5, se = TRUE,
                color = "#e41a1c", fill = "#e41a1c", alpha = 0.2) +
    geom_vline(xintercept = peak_half, linetype = "dashed",
               color = "#e41a1c", linewidth = 0.4, alpha = 0.6) +
    geom_vline(xintercept = beam_spacing / 2 - peak_half, linetype = "dashed",
               color = "#377eb8", linewidth = 0.4, alpha = 0.6) +
    facet_wrap(~ cell_class, scales = "free_y") +
    labs(x = "Distance from nearest beam centerline (mm)",
         y = sprintf("Mean %s expression (log2 CP10K + 1)", gene_name),
         title = sprintf("%s expression gradient from peak center", gene_label),
         subtitle = "MBRT 4h, Block 21 — ordered by radiosensitivity; red/blue dashed = peak/valley boundaries") +
    theme_bw(base_size = 12) +
    theme(strip.text = element_text(size = 11, face = "bold"))

  grad_file <- file.path(plot_dir, sprintf("pres_%s_gradient.png", tolower(gene_col)))
  ggsave(grad_file, pg, width = 12, height = 8, dpi = 200)
  cat(sprintf("Saved: %s\n", basename(grad_file)))

  # FOV-level % expressing
  fov_typed <- dt[!is.na(cell_class) & !is.na(pv_zone), .(
    pct_expr = 100 * mean(get(gene_col) > 0), n_cells = .N
  ), by = .(cell_class, fov, pv_zone)]
  fov_typed <- fov_typed[n_cells >= MIN_CELLS_PER_FOV]

  fov_allcells <- dt[!is.na(cell_class) & !is.na(pv_zone), .(
    pct_expr = 100 * mean(get(gene_col) > 0), n_cells = .N
  ), by = .(fov, pv_zone)]
  fov_allcells <- fov_allcells[n_cells >= MIN_CELLS_PER_FOV]
  fov_allcells[, cell_class := "All cells"]

  fov_all <- rbind(fov_typed, fov_allcells, use.names = TRUE, fill = TRUE)
  fov_all[, cell_class := factor(cell_class, levels = rs_levels)]
  fov_all[, pv_zone := factor(pv_zone, levels = c("peak", "valley"))]

  # Wilcoxon per cell class
  wres <- fov_all[, {
    pk <- pct_expr[pv_zone == "peak"]
    vl <- pct_expr[pv_zone == "valley"]
    if (length(pk) >= 3 && length(vl) >= 3) {
      wt <- suppressWarnings(wilcox.test(pk, vl, alternative = "greater"))
      .(p = wt$p.value, n_pk = length(pk), n_vl = length(vl))
    } else {
      .(p = NA_real_, n_pk = length(pk), n_vl = length(vl))
    }
  }, by = cell_class]
  wres[, p_label := fifelse(is.na(p), "",
    fifelse(p < 0.001, "p < 0.001",
    fifelse(p < 0.01, sprintf("p = %.3f", p),
    fifelse(p < 0.05, sprintf("p = %.2f", p),
    sprintf("p = %.2f", p)))))]

  cat(sprintf("\n%s FOV-level Wilcoxon (peak > valley):\n", gene_label))
  print(wres)

  annot <- merge(wres[, .(cell_class, p_label)],
                 fov_all[, .(ymax = max(pct_expr)), by = cell_class],
                 by = "cell_class")

  pf <- ggplot(fov_all, aes(x = pv_zone, y = pct_expr, fill = pv_zone)) +
    geom_boxplot(outlier.size = 0.8, width = 0.6, alpha = 0.7) +
    geom_jitter(width = 0.15, size = 0.8, alpha = 0.5, color = "grey30") +
    geom_text(data = annot, aes(x = 1.5, y = ymax * 1.08, label = p_label),
              size = 3.5, inherit.aes = FALSE) +
    scale_fill_manual(values = c("peak" = "#e41a1c", "valley" = "#377eb8")) +
    facet_wrap(~ cell_class, scales = "free_y", nrow = 1) +
    labs(x = NULL, y = sprintf("%% cells expressing %s per FOV", gene_name),
         title = sprintf("Peak vs valley %s expression (FOV-level)", gene_label),
         subtitle = sprintf("Each point = one FOV (>=%d cells of type); Wilcoxon rank-sum (peak > valley)",
                            MIN_CELLS_PER_FOV)) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none",
          strip.text = element_text(size = 9, face = "bold"),
          axis.text.x = element_text(size = 10))

  fov_file <- file.path(plot_dir, sprintf("pres_%s_fov_wilcoxon.png", tolower(gene_col)))
  ggsave(fov_file, pf, width = 14, height = 4.5, dpi = 200)
  cat(sprintf("Saved: %s\n", basename(fov_file)))
}

# Generate for all three genes
make_gene_plots(mbrt4h_grad, "Cdkn1a", "Cdkn1a", "p21 (Cdkn1a)", PLOT_DIR)
make_gene_plots(mbrt4h_grad, "Gadd45b", "Gadd45b", "Gadd45b", PLOT_DIR)
make_gene_plots(mbrt4h_grad, "Bax", "Bax", "Bax", PLOT_DIR)

# ===========================================================
# FIGURE 7: p21 spatial heatmap (KNN-smoothed)
# ===========================================================
library(RANN)

hm <- copy(mbrt4h_grad[!is.na(cell_class)])
pos_mat <- as.matrix(hm[, .(x_slide_mm, y_slide_mm)])

K_SMOOTH <- 150
cat(sprintf("KNN smoothing p21 over k=%d neighbors for %d cells...\n", K_SMOOTH, nrow(hm)))
nn <- nn2(pos_mat, pos_mat, k = K_SMOOTH)
p21_vec <- hm$Cdkn1a
hm[, p21_smooth := rowMeans(matrix(p21_vec[nn$nn.idx], ncol = K_SMOOTH))]

# Subsample for plotting (full 149K cells is heavy for ggplot)
set.seed(42)
hm_sub <- hm[sample(.N, min(30000, .N))]

x_lim_hm <- range(hm$x_slide_mm)
beam_lines_hm <- do.call(rbind, lapply(all_centers, function(c) {
  data.frame(x = x_lim_hm,
             y = (c + sin(theta) * x_lim_hm) / cos(theta),
             center = c)
}))

p7 <- ggplot(hm_sub, aes(x = x_slide_mm, y = y_slide_mm, color = p21_smooth)) +
  geom_point(size = 0.15, alpha = 0.8) +
  geom_line(data = beam_lines_hm, aes(x = x, y = y, group = center),
            color = "cyan", linewidth = 0.6, linetype = "solid",
            inherit.aes = FALSE, alpha = 0.8) +
  scale_color_gradientn(
    colors = c("grey10", "#4a1486", "#e41a1c", "#FFD700"),
    name = "Smoothed\nCdkn1a\n(k=150 NN)") +
  scale_x_continuous(breaks = x_breaks) +
  scale_y_continuous(breaks = y_breaks) +
  coord_fixed() +
  labs(x = "x (mm)", y = "y (mm)",
       title = "Cdkn1a (p21) spatial expression — MBRT 4h",
       subtitle = "Each cell colored by mean p21 of 150 nearest neighbors; cyan = beam centerlines") +
  theme_bw(base_size = 12) +
  theme(panel.background = element_rect(fill = "grey10"),
        panel.grid = element_blank())

ggsave(file.path(PLOT_DIR, "pres_p21_spatial_heatmap.png"), p7,
       width = 8, height = 10, dpi = 250)
cat("Saved: pres_p21_spatial_heatmap.png\n")

# ===========================================================
# FIGURE 8: All-cells peak vs valley volcano (FOV-level p-values)
# ===========================================================
library(ggrepel)

cat("\n=== Peak/Valley DGE: all-cells volcano ===\n")
gene_cols <- setdiff(names(counts_grad), c("Slide", "fov", "cell_id"))
peak_cells <- mbrt4h_grad[pv_zone == "peak"]
valley_cells <- mbrt4h_grad[pv_zone == "valley"]

counts_peak <- counts_full[cell_id %in% peak_cells$cell_id]
counts_valley <- counts_full[cell_id %in% valley_cells$cell_id]
setkeyv(counts_peak, "cell_id"); setkeyv(counts_valley, "cell_id")
counts_peak <- counts_peak[peak_cells$cell_id]
counts_valley <- counts_valley[valley_cells$cell_id]

lib_peak <- rowSums(counts_peak[, ..gene_cols])
lib_valley <- rowSums(counts_valley[, ..gene_cols])
norm_peak <- log2(sweep(as.matrix(counts_peak[, ..gene_cols]), 1, pmax(lib_peak, 1), "/") * 1e4 + 1)
norm_valley <- log2(sweep(as.matrix(counts_valley[, ..gene_cols]), 1, pmax(lib_valley, 1), "/") * 1e4 + 1)

# Per-gene summary
dge <- data.table(
  gene = gene_cols,
  mean_peak = colMeans(norm_peak),
  mean_valley = colMeans(norm_valley),
  pct_peak = colMeans(norm_peak > 0) * 100,
  pct_valley = colMeans(norm_valley > 0) * 100
)
dge[, log2FC := mean_peak - mean_valley]
dge[, delta_pct := pct_peak - pct_valley]

# FOV-level Wilcoxon for each gene (honest p-values)
cat("Running FOV-level Wilcoxon for", length(gene_cols), "genes...\n")
pk_fov_vec <- peak_cells$fov
vl_fov_vec <- valley_cells$fov

fov_gene_pvals <- sapply(gene_cols, function(g) {
  pk_fov_pct <- data.table(fov = pk_fov_vec, expr = counts_peak[[g]] > 0)[, .(pct = 100 * mean(expr)), by = fov]
  vl_fov_pct <- data.table(fov = vl_fov_vec, expr = counts_valley[[g]] > 0)[, .(pct = 100 * mean(expr)), by = fov]
  if (nrow(pk_fov_pct) < 3 || nrow(vl_fov_pct) < 3) return(NA_real_)
  suppressWarnings(wilcox.test(pk_fov_pct$pct, vl_fov_pct$pct)$p.value)
})

dge[, fov_p := fov_gene_pvals]
dge[, fov_padj := p.adjust(fov_p, method = "BH")]
dge[, neg_log_p := -log10(pmax(fov_p, 1e-20))]

# Classify using raw (nominal) p-values
dge[, category := fifelse(is.na(fov_p), "n.s.",
  fifelse(fov_p < 0.05 & log2FC > 0, "Peak-enriched (p < 0.05)",
  fifelse(fov_p < 0.05 & log2FC < 0, "Valley-enriched (p < 0.05)", "n.s.")))]

# Label significant hits + large-FC non-significant genes
n_label <- 15
top_peak <- dge[category == "Peak-enriched (p < 0.05)"][order(-log2FC)][1:min(n_label, .N), gene]
top_valley <- dge[category == "Valley-enriched (p < 0.05)"][order(log2FC)][1:min(n_label, .N), gene]
large_fc <- dge[category == "n.s." & abs(log2FC) > 0.2, gene]
dge[, label := fifelse(gene %in% c(top_peak, top_valley, large_fc), gene, NA_character_)]

cat(sprintf("Nominal p < 0.05: %d peak-enriched, %d valley-enriched (none survive BH correction)\n",
    sum(dge$category == "Peak-enriched (p < 0.05)", na.rm = TRUE),
    sum(dge$category == "Valley-enriched (p < 0.05)", na.rm = TRUE)))

p8 <- ggplot(dge, aes(x = log2FC, y = neg_log_p, color = category)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_point(alpha = 0.7, size = 1.5) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 30,
                  color = "black", min.segment.length = 0,
                  box.padding = 0.3, seed = 1) +
  scale_color_manual(values = c("Peak-enriched (p < 0.05)" = "#e41a1c",
                                "Valley-enriched (p < 0.05)" = "#377eb8",
                                "n.s." = "grey60")) +
  labs(x = "log2 fold change (peak / valley)",
       y = "-log10(FOV-level Wilcoxon p-value)",
       color = NULL,
       title = "Peak vs valley DGE — all cells",
       subtitle = "MBRT 4h Block 21; nominal p-values (not FDR-corrected); dashed line = p = 0.05") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(PLOT_DIR, "pres_pv_volcano.png"), p8,
       width = 10, height = 4, dpi = 200)
cat("Saved: pres_pv_volcano.png\n")

# Save DGE table
fwrite(dge[order(fov_p)], file.path(PLOT_DIR, "..", "data", "pres_pv_dge_allcells.tsv"), sep = "\t")
cat("Saved: data/pres_pv_dge_allcells.tsv\n")

# ===========================================================
# FIGURE 9: UpSet plots — cell type decomposition
# ===========================================================
library(UpSetR)

cat("\n=== Cell-type decomposition ===\n")
ct_levels <- c("B cells", "T cells", "Monocytes/DC", "Macrophages", "Stroma", "Tumor (4T1)")

# Per-cell-type DGE: log2FC + FOV-level Wilcoxon
ct_dge_list <- list()
for (ct in ct_levels) {
  pk <- mbrt4h_grad[cell_class == ct & pv_zone == "peak"]
  vl <- mbrt4h_grad[cell_class == ct & pv_zone == "valley"]
  if (nrow(pk) < 50 || nrow(vl) < 50) next

  cpk <- counts_full[cell_id %in% pk$cell_id]
  cvl <- counts_full[cell_id %in% vl$cell_id]
  setkeyv(cpk, "cell_id"); setkeyv(cvl, "cell_id")
  cpk <- cpk[pk$cell_id]; cvl <- cvl[vl$cell_id]

  lpk <- rowSums(cpk[, ..gene_cols]); lvl <- rowSums(cvl[, ..gene_cols])
  npk <- log2(sweep(as.matrix(cpk[, ..gene_cols]), 1, pmax(lpk, 1), "/") * 1e4 + 1)
  nvl <- log2(sweep(as.matrix(cvl[, ..gene_cols]), 1, pmax(lvl, 1), "/") * 1e4 + 1)

  ct_res <- data.table(
    gene = gene_cols,
    cell_type = ct,
    log2FC = colMeans(npk) - colMeans(nvl),
    pct_peak = colMeans(npk > 0) * 100,
    pct_valley = colMeans(nvl > 0) * 100
  )
  ct_res[, delta_pct := pct_peak - pct_valley]

  # FOV-level p-values per cell type
  ct_fov_p <- sapply(gene_cols, function(g) {
    pk_fov <- data.table(fov = pk$fov, expr = cpk[[g]] > 0)[, .(pct = 100 * mean(expr)), by = fov]
    vl_fov <- data.table(fov = vl$fov, expr = cvl[[g]] > 0)[, .(pct = 100 * mean(expr)), by = fov]
    pk_fov <- pk_fov[, .SD[.N >= 1], by = fov]  # keep all FOVs
    vl_fov <- vl_fov[, .SD[.N >= 1], by = fov]
    if (nrow(pk_fov) < 3 || nrow(vl_fov) < 3) return(NA_real_)
    suppressWarnings(wilcox.test(pk_fov$pct, vl_fov$pct)$p.value)
  })
  ct_res[, fov_p := ct_fov_p]
  ct_res[, fov_padj := p.adjust(fov_p, method = "BH")]

  ct_dge_list[[ct]] <- ct_res
  cat(sprintf("  %s: %d peak-up, %d valley-up (nominal p < 0.05)\n", ct,
      sum(ct_res$fov_p < 0.05 & ct_res$log2FC > 0, na.rm = TRUE),
      sum(ct_res$fov_p < 0.05 & ct_res$log2FC < 0, na.rm = TRUE)))
}
ct_dge <- rbindlist(ct_dge_list)
fwrite(ct_dge[order(fov_p)], file.path(PLOT_DIR, "..", "data", "pres_pv_dge_by_celltype.tsv"), sep = "\t")

# Build UpSet input: for each direction, which cell types call each gene significant?
make_upset_data <- function(ct_dge, direction = "peak") {
  if (direction == "peak") {
    sig <- ct_dge[fov_p < 0.05 & log2FC > 0]
  } else {
    sig <- ct_dge[fov_p < 0.05 & log2FC < 0]
  }
  if (nrow(sig) == 0) return(NULL)

  all_genes <- unique(sig$gene)
  all_cts <- unique(sig$cell_type)
  mat <- matrix(0L, nrow = length(all_genes), ncol = length(all_cts),
                dimnames = list(all_genes, all_cts))
  for (i in seq_len(nrow(sig))) {
    mat[sig$gene[i], sig$cell_type[i]] <- 1L
  }
  as.data.frame(mat)
}

# Peak-enriched UpSet
peak_upset_data <- make_upset_data(ct_dge, "peak")
if (!is.null(peak_upset_data) && ncol(peak_upset_data) >= 2 && nrow(peak_upset_data) >= 2) {
  png(file.path(PLOT_DIR, "pres_pv_upset_peak.png"), width = 10, height = 6, units = "in", res = 200)
  print(upset(peak_upset_data, sets = names(peak_upset_data),
              order.by = "freq", mb.ratio = c(0.6, 0.4),
              main.bar.color = "#e41a1c", sets.bar.color = "#e41a1c",
              text.scale = 1.3,
              mainbar.y.label = "Shared peak-enriched genes",
              sets.x.label = "Genes per cell type"))
  dev.off()
  cat("Saved: pres_pv_upset_peak.png\n")
} else {
  cat("Not enough peak-enriched genes across cell types for UpSet plot\n")
}

# Valley-enriched UpSet
valley_upset_data <- make_upset_data(ct_dge, "valley")
if (!is.null(valley_upset_data) && ncol(valley_upset_data) >= 2 && nrow(valley_upset_data) >= 2) {
  png(file.path(PLOT_DIR, "pres_pv_upset_valley.png"), width = 10, height = 6, units = "in", res = 200)
  print(upset(valley_upset_data, sets = names(valley_upset_data),
              order.by = "freq", mb.ratio = c(0.6, 0.4),
              main.bar.color = "#377eb8", sets.bar.color = "#377eb8",
              text.scale = 1.3,
              mainbar.y.label = "Shared valley-enriched genes",
              sets.x.label = "Genes per cell type"))
  dev.off()
  cat("Saved: pres_pv_upset_valley.png\n")
} else {
  cat("Not enough valley-enriched genes across cell types for UpSet plot\n")
}

cat("\nAll presentation figures generated.\n")
