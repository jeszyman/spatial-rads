suppressPackageStartupMessages({
  library(arrow); library(data.table); library(Matrix)
  library(ggplot2); library(patchwork)
})
OUT <- "/home/jeszyman/repos/spatial-rads/dev/if_channel_check"

grp <- fread("/tmp/meringue_spatial_patterns_block21.tsv")
cat("Pattern table:\n"); print(grp[, .N, by = pattern])

meta <- as.data.table(read_parquet("/tmp/mutter01_meta.parquet"))
counts <- as.data.table(read_parquet("/tmp/mutter01_counts.parquet"))

m <- meta[Slide == "20250529_214712_S4" & Block == "Block_21"]
c_df <- counts[cell_id %in% m$cell_id]
setkeyv(m, "cell_id"); setkeyv(c_df, "cell_id")
m <- m[cell_id %in% c_df$cell_id]

set.seed(1)
idx <- sample(nrow(m), 15000)
m <- m[idx]; c_df <- c_df[cell_id %in% m$cell_id]
setorderv(m, "cell_id"); setorderv(c_df, "cell_id")

gene_cols <- setdiff(colnames(c_df), c("Slide", "fov", "cell_id"))
mat <- t(as.matrix(c_df[, ..gene_cols]))
colnames(mat) <- c_df$cell_id; rownames(mat) <- gene_cols
lib <- colSums(mat); lib[lib == 0] <- 1
normed <- log2(t(t(mat) / lib) * 1e4 + 1)

pos <- as.matrix(m[, .(x_slide_mm, y_slide_mm)])

plot_pattern <- function(pattern_id, genes) {
  gens <- intersect(genes, rownames(normed))
  expr <- if (length(gens) > 1) colMeans(normed[gens, ]) else normed[gens, ]
  cap_hi <- quantile(expr, 0.98)
  cap_lo <- quantile(expr, 0.02)
  df <- data.table(x = pos[,1], y = pos[,2],
                   expr = pmin(pmax(expr, cap_lo), cap_hi))
  ggplot(df, aes(x, y, color = expr)) +
    geom_point(size = 1.0, alpha = 0.85) +
    scale_color_viridis_c(option = "magma", name = "mean") +
    coord_fixed() + theme_void() +
    theme(legend.position = "right", legend.key.height = unit(0.4, "cm"),
          plot.title = element_text(size = 11, hjust = 0.5)) +
    labs(title = sprintf("Pattern %d (%d genes): %s",
                         pattern_id, length(gens),
                         paste(head(gens, 6), collapse = ", ")))
}

pats <- sort(unique(grp$pattern))
plots <- lapply(pats, function(p) plot_pattern(p, grp[pattern == p, gene]))
combined <- wrap_plots(plots, ncol = 2)
ggsave(file.path(OUT, "meringue_patterns_block21.png"),
       combined, width = 14, height = 4 * ceiling(length(pats) / 2), dpi = 130)
cat("Saved meringue_patterns_block21.png\n")
