suppressPackageStartupMessages({
  library(arrow); library(data.table); library(Matrix); library(MERINGUE)
  library(ggplot2); library(patchwork)
})

OUT <- "/home/jeszyman/repos/spatial-rads/dev/if_channel_check"
BLOCK <- "Block_21"  # MBRT 4h on S4
SLIDE <- "20250529_214712_S4"
N_SUB <- 2500
set.seed(1)

cat("Loading data...\n")
meta <- as.data.table(read_parquet("/tmp/mutter01_meta.parquet"))
counts_df <- as.data.table(read_parquet("/tmp/mutter01_counts.parquet"))

m <- meta[Slide == SLIDE & Block == BLOCK]
counts_df <- counts_df[cell_id %in% m$cell_id]
setkeyv(m, "cell_id"); setkeyv(counts_df, "cell_id")
m <- m[cell_id %in% counts_df$cell_id]

idx <- sample(nrow(m), N_SUB)
m <- m[idx]; counts_df <- counts_df[cell_id %in% m$cell_id]
setorderv(m, "cell_id"); setorderv(counts_df, "cell_id")
cat("Cells:", nrow(m), "\n")

gene_cols <- setdiff(colnames(counts_df), c("Slide", "fov", "cell_id"))
mat <- t(as.matrix(counts_df[, ..gene_cols]))
colnames(mat) <- counts_df$cell_id; rownames(mat) <- gene_cols
stopifnot(all(m$cell_id == colnames(mat)))

# Normalize
lib <- colSums(mat); lib[lib == 0] <- 1
normed <- log2(t(t(mat) / lib) * 1e4 + 1)
pct_exp <- rowSums(mat > 0) / ncol(mat)
normed <- normed[pct_exp > 0.05, ]
cat("Genes after filter:", nrow(normed), "\n")

pos <- as.matrix(m[, .(x_slide_mm, y_slide_mm)])
rownames(pos) <- m$cell_id

# ---- Step 1: spatial neighbors (using filterDist in mm) ----
cat("Building spatial neighbor matrix (dist <= 0.05 mm = 50 um)...\n")
w <- getSpatialNeighbors(pos, filterDist = 0.05)
cat("Neighbors per cell (mean):", round(mean(rowSums(w > 0)), 1), "\n")

# ---- Step 2: spatial autocorrelation (Moran's I with significance) ----
cat("Computing Moran's I with significance...\n")
I <- getSpatialPatterns(normed, w)
cat("Top 10 by Moran's I:\n")
print(head(I[order(-I$observed), ], 10))

# ---- Step 3: filter to significant spatial genes ----
cat("\nFiltering to significant spatial genes...\n")
results.filter <- filterSpatialPatterns(mat = normed, I = I, w = w,
                                        adjustPv = TRUE, alpha = 0.05,
                                        minPercentCells = 0.05,
                                        verbose = TRUE)
cat("Significant spatial genes:", length(results.filter), "\n")
if (length(results.filter) < 10) stop("Too few significant genes; check params")

# ---- Step 4: spatial cross-correlation matrix ----
cat("Computing spatial cross-correlation matrix...\n")
scc <- spatialCrossCorMatrix(mat = as.matrix(normed[results.filter, ]), weight = w)
cat("SCC matrix:", dim(scc), "\n")

# ---- Step 5: group into spatial patterns ----
cat("Grouping genes into spatial patterns...\n")
ggroup <- groupSigSpatialPatterns(pos = pos,
                                  mat = as.matrix(normed[results.filter, ]),
                                  scc = scc,
                                  power = 1,
                                  hclustMethod = "ward.D",
                                  deepSplit = 2,
                                  zlim = c(-1.5, 1.5),
                                  plot = FALSE)
cat("Pattern group counts:\n"); print(table(ggroup$groups))

# Save group assignments
grp_df <- data.table(gene = names(ggroup$groups), pattern = ggroup$groups)
setorder(grp_df, pattern, gene)
fwrite(grp_df, "/tmp/meringue_spatial_patterns_block21.tsv", sep = "\t")
cat("Saved /tmp/meringue_spatial_patterns_block21.tsv\n")

# Print genes per pattern (truncated)
for (g in sort(unique(grp_df$pattern))) {
  genes <- grp_df[pattern == g, gene]
  cat(sprintf("\n=== Pattern %d (%d genes) ===\n", g, length(genes)))
  cat(head(genes, 25), sep = ", "); cat("\n")
}

# ---- Step 6: visualize each pattern as a mean-expression spatial map ----
plot_pattern <- function(pattern_id, genes) {
  if (length(genes) < 1) return(ggplot() + theme_void())
  pattern_expr <- colMeans(normed[genes, , drop = FALSE])
  cap_hi <- quantile(pattern_expr, 0.99)
  cap_lo <- quantile(pattern_expr, 0.01)
  df <- data.table(x = pos[, 1], y = pos[, 2],
                   expr = pmin(pmax(pattern_expr, cap_lo), cap_hi))
  ggplot(df, aes(x, y, color = expr)) +
    geom_point(size = 1.2, alpha = 0.9) +
    scale_color_viridis_c(option = "magma", name = "mean\nexpr") +
    coord_fixed() + theme_void() +
    theme(legend.position = "right", legend.key.height = unit(0.4, "cm"),
          plot.title = element_text(size = 11, hjust = 0.5)) +
    labs(title = sprintf("Pattern %d — %d genes (e.g. %s)",
                         pattern_id, length(genes),
                         paste(head(genes, 4), collapse = ", ")))
}

patterns <- sort(unique(grp_df$pattern))
plots <- lapply(patterns, function(p) plot_pattern(p, grp_df[pattern == p, gene]))
combined <- wrap_plots(plots, ncol = 2)
ggsave(file.path(OUT, "meringue_patterns_block21.png"), combined,
       width = 14, height = 4 * ceiling(length(patterns) / 2), dpi = 130)
cat("\nSaved", file.path(OUT, "meringue_patterns_block21.png"), "\n")

# Save scc heatmap
png(file.path(OUT, "meringue_scc_heatmap_block21.png"), width = 1200, height = 1200)
heatmap(scc, symm = TRUE, scale = "none", col = colorRampPalette(c("blue", "white", "red"))(50),
        main = "Spatial cross-correlation (Block_21 MBRT 4h)")
dev.off()
cat("Saved SCC heatmap.\n")
