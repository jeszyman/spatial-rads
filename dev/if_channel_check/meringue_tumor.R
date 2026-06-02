suppressPackageStartupMessages({
  library(arrow); library(data.table); library(Matrix); library(MERINGUE)
  library(ggplot2); library(patchwork)
})

OUT <- "/home/jeszyman/repos/spatial-rads/dev/if_channel_check"
BLOCK <- "Block_21"
SLIDE <- "20250529_214712_S4"
N_SUB <- 2500
set.seed(1)

peak_fovs   <- unique(c(224, 209, 201, 186, 172, 225, 229, 215, 175,
                        176, 206, 181, 168, 167))
valley_fovs <- c(199, 213, 227, 188, 204, 218, 178, 165, 170, 231, 189,
                 173, 154)

cat("Loading data...\n")
meta <- as.data.table(read_parquet("/tmp/mutter01_meta.parquet"))
counts_df <- as.data.table(read_parquet("/tmp/mutter01_counts.parquet"))

# Tumor-only in Block_21
m <- meta[Slide == SLIDE & Block == BLOCK &
          ImmuneAtlas_ImmGen_Main_cell_Types == "a"]
cat("Tumor cells in Block_21:", nrow(m), "\n")
counts_df <- counts_df[cell_id %in% m$cell_id]
setkeyv(m, "cell_id"); setkeyv(counts_df, "cell_id")
m <- m[cell_id %in% counts_df$cell_id]

idx <- sample(nrow(m), min(N_SUB, nrow(m)))
m <- m[idx]; counts_df <- counts_df[cell_id %in% m$cell_id]
setorderv(m, "cell_id"); setorderv(counts_df, "cell_id")
cat("Subsampled to:", nrow(m), "\n")

gene_cols <- setdiff(colnames(counts_df), c("Slide", "fov", "cell_id"))
mat <- t(as.matrix(counts_df[, ..gene_cols]))
colnames(mat) <- counts_df$cell_id; rownames(mat) <- gene_cols
stopifnot(all(m$cell_id == colnames(mat)))

lib <- colSums(mat); lib[lib == 0] <- 1
normed <- log2(t(t(mat) / lib) * 1e4 + 1)
pct_exp <- rowSums(mat > 0) / ncol(mat)
normed <- normed[pct_exp > 0.05, ]
cat("Genes after filter:", nrow(normed), "\n")

pos <- as.matrix(m[, .(x_slide_mm, y_slide_mm)])
rownames(pos) <- m$cell_id

cat("Building spatial neighbors (filterDist = 0.05 mm)...\n")
w <- getSpatialNeighbors(pos, filterDist = 0.05)
cat("Mean neighbors:", round(mean(rowSums(w > 0)), 1), "\n")

cat("Moran's I with significance...\n")
I <- getSpatialPatterns(normed, w)
cat("Top 10 by Moran's I:\n")
print(head(I[order(-I$observed), ], 10))

results.filter <- filterSpatialPatterns(mat = normed, I = I, w = w,
                                        adjustPv = TRUE, alpha = 0.05,
                                        minPercentCells = 0.05,
                                        verbose = TRUE)
cat("Significant spatial genes:", length(results.filter), "\n")
if (length(results.filter) < 5) stop("Too few significant genes")

cat("Spatial cross-correlation matrix...\n")
scc <- spatialCrossCorMatrix(mat = as.matrix(normed[results.filter, ]), weight = w)

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

grp_df <- data.table(gene = names(ggroup$groups),
                     pattern = as.integer(as.character(ggroup$groups)))
setorder(grp_df, pattern, gene)
fwrite(grp_df, "/tmp/meringue_spatial_patterns_tumor.tsv", sep = "\t")
cat("Saved /tmp/meringue_spatial_patterns_tumor.tsv\n")

for (g in sort(unique(grp_df$pattern))) {
  genes <- grp_df[pattern == g, gene]
  cat(sprintf("\n=== Pattern %d (%d genes) ===\n", g, length(genes)))
  cat(head(genes, 25), sep = ", "); cat("\n")
}

# Pattern spatial maps, with peak/valley FOV annotation
m[, pv := fifelse(fov %in% peak_fovs, "peak",
                  fifelse(fov %in% valley_fovs, "valley", "other"))]
pv_vec <- m$pv

plot_pattern <- function(pattern_id, genes) {
  if (length(genes) < 1) return(ggplot() + theme_void())
  expr <- colMeans(normed[genes, , drop = FALSE])
  cap_hi <- quantile(expr, 0.99); cap_lo <- quantile(expr, 0.01)
  df <- data.table(x = pos[, 1], y = pos[, 2],
                   expr = pmin(pmax(expr, cap_lo), cap_hi),
                   pv = pv_vec)
  ggplot(df, aes(x, y, color = expr)) +
    geom_point(size = 1.2, alpha = 0.9) +
    scale_color_viridis_c(option = "magma", name = "mean") +
    coord_fixed() + theme_void() +
    theme(legend.position = "right", legend.key.height = unit(0.4, "cm"),
          plot.title = element_text(size = 10, hjust = 0.5)) +
    labs(title = sprintf("Pattern %d (%d genes): %s",
                         pattern_id, length(genes),
                         paste(head(genes, 5), collapse = ", ")))
}

# Peak/valley cell overlay (separate small panel)
pv_plot <- ggplot(data.table(x = pos[, 1], y = pos[, 2], pv = pv_vec),
                  aes(x, y, color = pv)) +
  geom_point(size = 0.8, alpha = 0.8) +
  scale_color_manual(values = c(peak = "#e41a1c", valley = "#377eb8",
                                other = "grey85")) +
  coord_fixed() + theme_void() +
  theme(legend.position = "right",
        plot.title = element_text(size = 10, hjust = 0.5)) +
  labs(title = "Peak (red) / Valley (blue) FOV assignment")

patterns <- sort(unique(grp_df$pattern))
plots <- lapply(patterns, function(p) plot_pattern(p, grp_df[pattern == p, gene]))
plots <- c(list(pv_plot), plots)
combined <- wrap_plots(plots, ncol = 2)
ggsave(file.path(OUT, "meringue_tumor_patterns.png"), combined,
       width = 14, height = 4 * ceiling(length(plots) / 2), dpi = 130)
cat("\nSaved meringue_tumor_patterns.png\n")

# Mean pattern score per peak/valley
score_by_pv <- rbindlist(lapply(patterns, function(pp) {
  genes <- grp_df[pattern == pp, gene]
  expr <- colMeans(normed[genes, , drop = FALSE])
  data.table(pattern = pp,
             pv = pv_vec,
             score = expr)[, .(mean = mean(score), n = .N),
                           by = .(pattern, pv)]
}))
fwrite(score_by_pv, "/tmp/meringue_tumor_pattern_by_pv.tsv", sep = "\t")
cat("Pattern mean score by peak/valley:\n")
print(dcast(score_by_pv, pattern ~ pv, value.var = "mean"))
