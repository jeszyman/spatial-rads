library(arrow)
library(data.table)
library(ggplot2)
library(patchwork)

meta <- as.data.table(read_parquet("/tmp/mutter01_meta.parquet"))

# Pick a dense FOV with good cell type diversity (S3 fov 1: 1208 cells)
fov_meta <- meta[Slide == "20250529_214712_S3" & fov == 1]
cat("FOV cells:", nrow(fov_meta), "\n")

if_cols <- c("Mean.DAPI", "Mean.PanCK", "Mean.CD45", "Mean.CD298.B2M", "Mean.G")

# Clip to 99th percentile per channel to prevent outliers dominating color scale
plot_if <- function(col) {
  vals <- fov_meta[[col]]
  cap <- quantile(vals, 0.99, na.rm = TRUE)
  ggplot(fov_meta, aes(x = x_FOV_px, y = y_FOV_px, color = pmin(.data[[col]], cap))) +
    geom_point(size = 0.3, alpha = 0.8) +
    scale_color_viridis_c(name = col, option = "magma") +
    coord_fixed() +
    theme_void() +
    theme(legend.position = "bottom", legend.key.width = unit(1, "cm")) +
    labs(title = col)
}

panels <- lapply(if_cols, plot_if)

# Add a 6th panel: cell type
fov_meta[, ct := ImmuneAtlas_ImmGen_cellFamily_Main_cell_Types]
ct_panel <- ggplot(fov_meta, aes(x = x_FOV_px, y = y_FOV_px, color = ct)) +
  geom_point(size = 0.3, alpha = 0.8) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "bottom", legend.text = element_text(size = 6)) +
  guides(color = guide_legend(override.aes = list(size = 2), nrow = 3)) +
  labs(title = "Cell family (insitutype)")

combined <- wrap_plots(c(panels, list(ct_panel)), ncol = 3)
ggsave("/home/jeszyman/repos/spatial-rads/dev/if_channel_check/if_channels_S3_fov1.png",
       combined, width = 15, height = 11, dpi = 150)

# ---- Quantitative: which cell family shows highest Mean.G? ----
cat("\n=== Mean.G by cell family (top 15) ===\n")
g_by_ct <- meta[!is.na(ImmuneAtlas_ImmGen_cellFamily_Main_cell_Types),
                .(mean_G = mean(Mean.G, na.rm = TRUE),
                  median_G = median(Mean.G, na.rm = TRUE),
                  n = .N),
                by = ImmuneAtlas_ImmGen_cellFamily_Main_cell_Types][order(-mean_G)]
print(head(g_by_ct, 15))

cat("\n=== Correlation between IF channels (Pearson) ===\n")
cor_mat <- cor(meta[, .(Mean.DAPI, Mean.PanCK, Mean.CD45, Mean.CD298.B2M, Mean.G)],
               use = "pairwise.complete.obs")
print(round(cor_mat, 3))

cat("\n=== Per-main-cell-type: G-high fraction ===\n")
g_thresh <- quantile(meta$Mean.G, 0.95, na.rm = TRUE)
cat("G-high threshold (95th pct):", round(g_thresh, 1), "\n")
g_high_by_ct <- meta[, .(g_high_frac = mean(Mean.G > g_thresh, na.rm = TRUE), n = .N),
                     by = ImmuneAtlas_ImmGen_Main_cell_Types][order(-g_high_frac)]
print(head(g_high_by_ct, 15))
