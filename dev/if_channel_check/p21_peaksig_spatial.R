suppressPackageStartupMessages({
  library(arrow); library(data.table); library(Matrix); library(MERINGUE)
})

OUT <- "/home/jeszyman/repos/spatial-rads/dev/if_channel_check"
BLOCK <- "Block_21"; SLIDE <- "20250529_214712_S4"
N_SUB <- 5000
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

gene_cols <- setdiff(colnames(counts_df), c("Slide", "fov", "cell_id"))
mat <- t(as.matrix(counts_df[, ..gene_cols]))
colnames(mat) <- counts_df$cell_id
lib <- colSums(mat); lib[lib == 0] <- 1
normed <- log2(t(t(mat) / lib) * 1e4 + 1)

pos <- as.matrix(m[, .(x_slide_mm, y_slide_mm)])
rownames(pos) <- m$cell_id

# Fit stripe axis by searching over stripe angles and finding the angle
# that best clusters peak-FOV centroids into 4 discrete stripes.
peak_cells <- m[fov %in% peak_fovs]
peak_fov_centroids <- peak_cells[, .(cx = mean(x_slide_mm),
                                     cy = mean(y_slide_mm)),
                                 by = fov]
# theta = stripe-axis angle from horizontal (radians).
# Projection onto perpendicular direction: d = -sin(theta)*x + cos(theta)*y
search_theta <- function(cx, cy, n_stripes = 4, theta_grid = seq(-pi/2, pi/2, length.out = 181)) {
  best <- list(theta = NA, within_var = Inf, d = NULL, km = NULL)
  for (th in theta_grid) {
    d <- -sin(th) * cx + cos(th) * cy
    km <- kmeans(d, centers = n_stripes, nstart = 10)
    if (km$tot.withinss < best$within_var) {
      best <- list(theta = th, within_var = km$tot.withinss, d = d, km = km)
    }
  }
  best
}
fit <- search_theta(peak_fov_centroids$cx, peak_fov_centroids$cy)
peak_slope <- tan(fit$theta)  # slope of stripes themselves
cat(sprintf("Fitted stripe angle: %.1f deg from horizontal | slope = %.3f | tilt from vertical = %.1f deg\n",
            fit$theta * 180/pi, peak_slope, 90 - fit$theta * 180/pi))
peak_fov_centroids[, stripe := fit$km$cluster]
peak_fov_centroids[, d_perp := fit$d]
peak_line_params <- peak_fov_centroids[, .(
  cx_mean = mean(cx),
  cy_mean = mean(cy)
), by = stripe]
peak_line_params[, slope := peak_slope]
peak_line_params[, intercept := cy_mean - slope * cx_mean]
cat("Stripe intercepts (sorted by perpendicular offset):\n")
print(peak_line_params[order(intercept)])

add_peak_lines <- function() {
  ylim <- par("usr")[3:4]
  for (i in seq_len(nrow(peak_line_params))) {
    abline(a = peak_line_params$intercept[i],
           b = peak_line_params$slope[i],
           col = "black", lty = 2, lwd = 1.5)
  }
}

# Spatial neighbor smoothing
w_smooth <- getSpatialNeighbors(pos, filterDist = 0.15)
w_smooth <- w_smooth + Matrix::Diagonal(nrow(w_smooth))
rs <- Matrix::rowSums(w_smooth); rs[rs == 0] <- 1
w_smooth <- w_smooth / rs
smooth_expr <- function(v) as.numeric(w_smooth %*% v)

# p21 = Cdkn1a
if (!("Cdkn1a" %in% rownames(normed))) stop("Cdkn1a not in matrix")
p21 <- normed["Cdkn1a", ]
p21_sm <- smooth_expr(p21)
names(p21_sm) <- colnames(normed)

# Peak signature: top 10 tumor peak-UP from DGE
dge <- fread("/tmp/dge_pv_by_celltype_fixed.tsv")
tumor_peak <- dge[cell_class == "Tumor cells (a)" & avg_log2FC > 0.1
                 ][order(-avg_log2FC)][1:15, gene]
tumor_peak <- intersect(tumor_peak, rownames(normed))
cat("Tumor peak-sig genes used:", length(tumor_peak), "->",
    paste(tumor_peak, collapse = ", "), "\n")
ps <- colMeans(normed[tumor_peak, , drop = FALSE])
ps_sm <- smooth_expr(ps); names(ps_sm) <- colnames(normed)

# Valley sig
tumor_val <- dge[cell_class == "Tumor cells (a)" & avg_log2FC < -0.1
                ][order(avg_log2FC)][1:15, gene]
tumor_val <- intersect(tumor_val, rownames(normed))
vs <- colMeans(normed[tumor_val, , drop = FALSE])
vs_sm <- smooth_expr(vs); names(vs_sm) <- colnames(normed)

# Plot
# Match panel aspect to data aspect ratio
x_range <- diff(range(pos[, 1])); y_range <- diff(range(pos[, 2]))
asp <- y_range / x_range  # ~0.58 for this tumor
panel_w <- 900; panel_h <- round(panel_w * asp)
png(file.path(OUT, "p21_peaksig_spatial.png"),
    width = panel_w * 2 + 80, height = (panel_h + 60) * 2 + 40, res = 130)
par(mfrow = c(2, 2), mar = c(2, 2, 3, 1), pty = "m")

# Peak/valley FOV overlay
pv <- ifelse(m$fov %in% peak_fovs, "peak",
             ifelse(m$fov %in% valley_fovs, "valley", "other"))
col <- c(peak = "#e41a1c", valley = "#377eb8", other = "grey85")[pv]
plot(pos, col = col, pch = 19, cex = 0.3, asp = 1,
     main = "Peak (red) / Valley (blue) FOV assignment",
     xlab = "x mm", ylab = "y mm")
add_peak_lines()
legend("topright", legend = c("peak", "valley", "other"),
       col = c("#e41a1c", "#377eb8", "grey85"), pch = 19, bty = "n", cex = 0.8)

ZL <- c(-1.3, 1.3)  # tighter saturation for better contrast
interpolate(pos, p21_sm, scale = TRUE, fill = TRUE,
            zlim = ZL, binSize = 60,
            main = "Cdkn1a (p21) expression", plot = TRUE)
add_peak_lines()

interpolate(pos, ps_sm, scale = TRUE, fill = TRUE,
            zlim = ZL, binSize = 60,
            main = sprintf("Tumor peak signature (%d genes)", length(tumor_peak)),
            plot = TRUE)
add_peak_lines()

interpolate(pos, vs_sm, scale = TRUE, fill = TRUE,
            zlim = ZL, binSize = 60,
            main = sprintf("Tumor valley signature (%d genes)", length(tumor_val)),
            plot = TRUE)
add_peak_lines()

dev.off()
cat("Saved p21_peaksig_spatial.png\n")
