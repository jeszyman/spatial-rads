suppressPackageStartupMessages({
  library(arrow); library(data.table); library(Matrix); library(MERINGUE)
})

OUT <- "/home/jeszyman/repos/spatial-rads/dev/if_channel_check"
SLIDE <- "20250529_214712_S4"; BLOCK <- "Block_21"
set.seed(1)

peak_fovs   <- c(224, 209, 201, 186, 172, 225, 229, 215, 175,
                 176, 206, 181, 168, 167)
valley_fovs <- c(199, 213, 227, 188, 204, 218, 178, 165, 170, 231, 189,
                 173, 154)

meta <- as.data.table(read_parquet("/tmp/mutter01_meta.parquet"))
counts_df <- as.data.table(read_parquet("/tmp/mutter01_counts.parquet"))
setnames(meta, "ImmuneAtlas_ImmGen_Main_cell_Types", "main_type")

mac_types <- c("Macrophage","PerC.macrophage","spleen.red.pulp.macs",
               "small.peritoneal.macs","Microglia")
m <- meta[Slide == SLIDE & Block == BLOCK & main_type %in% mac_types]
counts_df <- counts_df[cell_id %in% m$cell_id]
setkeyv(m, "cell_id"); setkeyv(counts_df, "cell_id")
m <- m[cell_id %in% counts_df$cell_id]
idx <- sample(nrow(m), min(3500, nrow(m)))
m <- m[idx]; counts_df <- counts_df[cell_id %in% m$cell_id]
setorderv(m, "cell_id"); setorderv(counts_df, "cell_id")

gene_cols <- setdiff(colnames(counts_df), c("Slide","fov","cell_id"))
mat <- t(as.matrix(counts_df[, ..gene_cols]))
colnames(mat) <- counts_df$cell_id
lib <- colSums(mat); lib[lib == 0] <- 1
normed <- log2(t(t(mat) / lib) * 1e4 + 1)

pos <- as.matrix(m[, .(x_slide_mm, y_slide_mm)])
rownames(pos) <- m$cell_id

# ---- Fit stripe geometry from peak FOV centroids (full Block_21 meta, not just macs) ----
all_block_cells <- meta[Slide == SLIDE & Block == BLOCK]
peak_centroids <- all_block_cells[fov %in% peak_fovs,
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
cat(sprintf("Stripe: theta=%.1f deg, spacing=%.3f mm, peak_half=%.3f mm\n",
            theta * 180/pi, beam_spacing, peak_half))

# Extrapolate labels to macrophage cells
d_mac <- -sin(theta) * pos[,1] + cos(theta) * pos[,2]
d_range <- range(d_mac)
all_centers <- seq(fit$centers[1] -
                     ceiling((fit$centers[1] - d_range[1]) / beam_spacing) * beam_spacing,
                   fit$centers[length(fit$centers)] +
                     ceiling((d_range[2] - fit$centers[length(fit$centers)]) / beam_spacing) * beam_spacing,
                   by = beam_spacing)
dist_to_peak <- sapply(d_mac, function(d) min(abs(d - all_centers)))
pv_extrap <- ifelse(dist_to_peak < peak_half, "peak",
                    ifelse(dist_to_peak > (beam_spacing/2 - peak_half),
                           "valley", "transition"))
cat("Macrophage label counts:\n"); print(table(pv_extrap))

# ---- Smoothing and pattern scoring with extrapolated labels ----
w <- getSpatialNeighbors(pos, filterDist = 0.05)
w_smooth <- getSpatialNeighbors(pos, filterDist = 0.15)
w_smooth <- w_smooth + Matrix::Diagonal(nrow(w_smooth))
rs <- Matrix::rowSums(w_smooth); rs[rs == 0] <- 1
w_smooth <- w_smooth / rs
smooth_expr <- function(v) as.numeric(w_smooth %*% v)

grp <- fread("/tmp/meringue_macrophage_gene_pattern.tsv")
patterns <- sort(unique(grp$pattern))

pattern_stats <- rbindlist(lapply(patterns, function(pp) {
  genes <- grp[pattern == pp, gene]
  genes <- intersect(genes, rownames(normed))
  expr <- colMeans(normed[genes, , drop = FALSE])
  # Extrapolated peak vs valley
  pk <- mean(expr[pv_extrap == "peak"])
  vl <- mean(expr[pv_extrap == "valley"])
  # Wilcoxon peak vs valley, both extrapolated
  wt <- wilcox.test(expr[pv_extrap == "peak"], expr[pv_extrap == "valley"])
  # Moran's I on mean
  mt <- moranTest(expr, weight = w)
  data.table(pattern = pp, n_genes = length(genes),
             morans_I = round(mt["observed"], 3),
             morans_p = mt["p.value"],
             peak = round(pk, 3), valley = round(vl, 3),
             peak_minus_valley = round(pk - vl, 3),
             wilcox_p = wt$p.value,
             top_genes = paste(head(genes, 10), collapse = ", "))
}))
fwrite(pattern_stats, "/tmp/meringue_macrophage_pattern_stats_extrap.tsv", sep = "\t")
cat("\n=== Pattern stats with extrapolated labels ===\n")
print(pattern_stats[, .(pattern, n_genes, morans_I,
                        peak, valley, peak_minus_valley,
                        wilcox_p = sprintf("%.1e", wilcox_p))])

# ---- Plot: correct aspect, with peak stripe lines ----
x_rng <- diff(range(pos[,1])); y_rng <- diff(range(pos[,2]))
asp <- y_rng / x_rng
panel_w <- 800
panel_h <- round(panel_w * asp)
ncol_plot <- 2
nrow_plot <- ceiling(length(patterns) / ncol_plot)

png(file.path(OUT, "meringue_macrophage_refined.png"),
    width = panel_w * ncol_plot + 80,
    height = (panel_h + 80) * nrow_plot + 40, res = 130)
par(mfrow = c(nrow_plot, ncol_plot), mar = c(2, 2, 3.5, 1), pty = "m")

add_lines <- function() {
  for (c in all_centers) {
    abline(a = c / cos(theta), b = sin(theta) / cos(theta),
           col = "black", lty = 2, lwd = 1.1)
  }
}

for (pp in patterns) {
  genes <- grp[pattern == pp, gene]
  genes <- intersect(genes, rownames(normed))
  expr <- smooth_expr(colMeans(normed[genes, , drop = FALSE]))
  names(expr) <- colnames(normed)
  ps <- pattern_stats[pattern == pp]
  title <- sprintf("Pattern %d (%d genes, I=%.2f, P−V=%+.2f, wilcox p=%.1e)\n%s",
                   pp, length(genes), ps$morans_I,
                   ps$peak_minus_valley, ps$wilcox_p,
                   substr(ps$top_genes, 1, 80))
  interpolate(pos, expr, scale = TRUE, fill = TRUE,
              zlim = c(-1.3, 1.3), binSize = 50,
              main = title, plot = TRUE)
  add_lines()
}
dev.off()
cat("Saved meringue_macrophage_refined.png\n")
