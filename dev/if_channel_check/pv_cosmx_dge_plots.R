suppressPackageStartupMessages({
  library(arrow); library(data.table); library(Matrix); library(Seurat)
  library(ggplot2); library(patchwork)
})

OUT <- "/home/jeszyman/repos/spatial-rads/dev/if_channel_check"
SLIDE <- "20250529_214712_S4"; BLOCK <- "Block_21"
set.seed(1)

peak_fovs <- c(224, 209, 201, 186, 172, 225, 229, 215, 175,
               176, 206, 181, 168, 167)

meta <- as.data.table(read_parquet("/tmp/mutter01_meta.parquet"))
counts_df <- as.data.table(read_parquet("/tmp/mutter01_counts.parquet"))
setnames(meta, "ImmuneAtlas_ImmGen_Main_cell_Types", "main_type")

m <- meta[Slide == SLIDE & Block == BLOCK & main_type == "a"]
c_df <- counts_df[cell_id %in% m$cell_id]
setkeyv(m, "cell_id"); setkeyv(c_df, "cell_id")
m <- m[cell_id %in% c_df$cell_id]
setorderv(m, "cell_id"); setorderv(c_df, "cell_id")

gene_cols <- setdiff(colnames(c_df), c("Slide","fov","cell_id"))
mat <- t(as.matrix(c_df[, ..gene_cols]))
colnames(mat) <- c_df$cell_id
lib <- colSums(mat); lib[lib == 0] <- 1
normed <- log2(t(t(mat) / lib) * 1e4 + 1)

# 4T1 tumor filter
epith <- intersect(c("Krt8","Krt18","Epcam","Cdh1"), rownames(normed))
endo  <- intersect(c("Pecam1","Cdh5","Vwf"), rownames(normed))
ep_s <- colMeans(normed[epith, , drop = FALSE])
en_s <- if (length(endo) > 0) colMeans(normed[endo, , drop = FALSE]) else 0
keep4t1 <- ep_s > 0.3 & en_s < 0.1
m <- m[keep4t1]; mat <- mat[, keep4t1]; normed <- normed[, keep4t1]
cat("4T1 tumor cells:", ncol(mat), "\n")

# Stripe geometry + peak/valley extrapolation
peak_centroids <- meta[Slide == SLIDE & Block == BLOCK & fov %in% peak_fovs,
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
spacing <- median(diff(fit$centers))
peak_centroids[, d_perp := -sin(theta) * cx + cos(theta) * cy]
peak_centroids[, stripe := fit$cluster]
within_sd <- peak_centroids[, .(sd = sd(d_perp - fit$centers[stripe])),
                            by = stripe][, mean(sd)]
peak_half <- within_sd + 0.15

pos <- as.matrix(m[, .(x_slide_mm, y_slide_mm)])
d <- -sin(theta) * pos[,1] + cos(theta) * pos[,2]
d_range <- range(d)
all_ctr <- seq(fit$centers[1] - ceiling((fit$centers[1] - d_range[1]) / spacing) * spacing,
               fit$centers[length(fit$centers)] +
                 ceiling((d_range[2] - fit$centers[length(fit$centers)]) / spacing) * spacing,
               by = spacing)
dist_to_peak <- sapply(d, function(x) min(abs(x - all_ctr)))
label <- ifelse(dist_to_peak < peak_half, "peak",
                ifelse(dist_to_peak > (spacing/2 - peak_half), "valley", "transition"))
keep <- label %in% c("peak","valley")
mat <- mat[, keep]; m <- m[keep]; normed <- normed[, keep]; label <- label[keep]
pos <- pos[keep, ]
cat("Final cells — peak:", sum(label == "peak"), "valley:", sum(label == "valley"), "\n")

# Seurat workflow
obj <- CreateSeuratObject(counts = mat)
obj <- NormalizeData(obj, verbose = FALSE)
obj$pv <- factor(label, levels = c("valley","peak"))
Idents(obj) <- obj$pv

# FindMarkers, rank by effect size (not p)
de <- FindMarkers(obj, ident.1 = "peak", ident.2 = "valley",
                  test.use = "wilcox", logfc.threshold = 0,
                  min.pct = 0.05, verbose = FALSE)
de$gene <- rownames(de)
de <- as.data.table(de)
setorder(de, -avg_log2FC)

top_peak   <- de[avg_log2FC > 0][1:15, gene]
top_valley <- de[order(avg_log2FC)][1:10, gene]
top_genes  <- c(top_peak, top_valley)

# ---- Dotplot ----
p_dot <- DotPlot(obj, features = unique(top_genes)) +
  scale_color_gradient2(low = "#377eb8", mid = "white", high = "#e41a1c",
                        midpoint = 0) +
  coord_flip() +
  labs(x = "Gene", y = "Region",
       title = "Top 15 peak-UP + 10 valley-UP genes (by avg log2FC)",
       subtitle = "4T1 tumor cells, Block_21 MBRT 4h. Dot size = % expressing, color = scaled expression.") +
  theme(axis.text.x = element_text(angle = 0),
        plot.title = element_text(size = 11),
        plot.subtitle = element_text(size = 9))

# ---- Heatmap of top 25 genes across cells ordered by peak/valley ----
obj_sub <- subset(obj, downsample = 500)  # 500 per group for readable heatmap
obj_sub <- ScaleData(obj_sub, features = unique(top_genes), verbose = FALSE)
p_heat <- DoHeatmap(obj_sub, features = unique(top_genes),
                    group.by = "pv", raster = TRUE) +
  labs(title = "Heatmap: top peak+valley markers, 500 cells subsample per group") +
  theme(axis.text.y = element_text(size = 7),
        plot.title = element_text(size = 11))

# ---- Spatial plots of top 4 genes (peak panel) + top 2 valley genes ----
x_rng <- diff(range(pos[,1])); y_rng <- diff(range(pos[,2]))
asp <- y_rng / x_rng

# Peak stripe lines for overlay
add_peak_lines <- function(p) {
  for (c in all_ctr) {
    p <- p + geom_abline(intercept = c / cos(theta),
                         slope = sin(theta) / cos(theta),
                         linetype = "dashed", color = "black",
                         linewidth = 0.3, alpha = 0.5)
  }
  p
}

plot_gene_spatial <- function(gene) {
  if (!(gene %in% rownames(normed))) return(ggplot() + theme_void() +
                                             labs(title = paste(gene, "not found")))
  expr <- normed[gene, ]
  cap_hi <- quantile(expr, 0.98); cap_lo <- quantile(expr, 0.02)
  df <- data.table(x = pos[,1], y = pos[,2],
                   expr = pmin(pmax(expr, cap_lo), cap_hi))
  p <- ggplot(df, aes(x, y, color = expr)) +
    geom_point(size = 0.4, alpha = 0.7) +
    scale_color_viridis_c(option = "magma", name = "log-norm") +
    coord_fixed() + theme_void() +
    theme(plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
          legend.position = "right") +
    labs(title = gene)
  add_peak_lines(p)
}

top4_peak   <- de[avg_log2FC > 0][1:4, gene]
top2_valley <- de[order(avg_log2FC)][1:2, gene]
spatial_plots <- lapply(c(top4_peak, top2_valley), plot_gene_spatial)
p_spatial <- wrap_plots(spatial_plots, ncol = 3) +
  plot_annotation(title = "Spatial expression of top peak+valley genes",
                  subtitle = "Dashed lines = fitted peak centerlines")

# Combine and save
combined_top <- p_dot / p_heat + plot_layout(heights = c(1, 1.3))
ggsave(file.path(OUT, "tumor_pv_dotplot_heatmap.png"),
       combined_top, width = 12, height = 14, dpi = 130)
ggsave(file.path(OUT, "tumor_pv_spatial.png"),
       p_spatial, width = 16, height = 10, dpi = 130)
cat("Saved tumor_pv_dotplot_heatmap.png and tumor_pv_spatial.png\n")
