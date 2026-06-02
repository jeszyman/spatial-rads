suppressPackageStartupMessages({
  library(arrow); library(data.table); library(Matrix); library(MERINGUE)
})

OUT <- "/home/jeszyman/repos/spatial-rads/dev/if_channel_check"
BLOCK <- "Block_21"; SLIDE <- "20250529_214712_S4"; N_SUB <- 2500
set.seed(1)

# Re-load cells and gene patterns (deterministic via seed)
meta <- as.data.table(read_parquet("/tmp/mutter01_meta.parquet"))
counts_df <- as.data.table(read_parquet("/tmp/mutter01_counts.parquet"))
m <- meta[Slide == SLIDE & Block == BLOCK &
          ImmuneAtlas_ImmGen_Main_cell_Types == "a"]
counts_df <- counts_df[cell_id %in% m$cell_id]
setkeyv(m, "cell_id"); setkeyv(counts_df, "cell_id")
m <- m[cell_id %in% counts_df$cell_id]
idx <- sample(nrow(m), min(N_SUB, nrow(m)))
m <- m[idx]; counts_df <- counts_df[cell_id %in% m$cell_id]
setorderv(m, "cell_id"); setorderv(counts_df, "cell_id")

gene_cols <- setdiff(colnames(counts_df), c("Slide", "fov", "cell_id"))
mat <- t(as.matrix(counts_df[, ..gene_cols]))
colnames(mat) <- counts_df$cell_id
lib <- colSums(mat); lib[lib == 0] <- 1
normed <- log2(t(t(mat) / lib) * 1e4 + 1)
pct_exp <- rowSums(mat > 0) / ncol(mat)
normed <- normed[pct_exp > 0.05, ]

pos <- as.matrix(m[, .(x_slide_mm, y_slide_mm)])
rownames(pos) <- m$cell_id
w <- getSpatialNeighbors(pos, filterDist = 0.05)

# Smoothing operator: wider radius averages each cell with ~30 spatial neighbors
w_smooth <- getSpatialNeighbors(pos, filterDist = 0.15)
w_smooth <- w_smooth + Matrix::Diagonal(nrow(w_smooth))  # include self
row_sums <- Matrix::rowSums(w_smooth); row_sums[row_sums == 0] <- 1
w_smooth <- w_smooth / row_sums
smooth_expr <- function(vec) as.numeric(w_smooth %*% vec)

grp <- fread("/tmp/meringue_spatial_patterns_tumor.tsv")
patterns <- sort(unique(grp$pattern))

# ---- Interpolated smooth pattern heatmaps (tutorial style) ----
png(file.path(OUT, "meringue_tumor_patterns_smooth.png"),
    width = 1400, height = 900, res = 130)
par(mfrow = c(2, ceiling(length(patterns) / 2)),
    mar = c(2, 2, 3, 1))
for (pp in patterns) {
  genes <- grp[pattern == pp, gene]
  expr <- colMeans(normed[genes, , drop = FALSE])
  expr <- smooth_expr(expr)
  names(expr) <- colnames(normed)
  interpolate(pos, expr, scale = TRUE, fill = TRUE,
              zlim = c(-2, 2), binSize = 60,
              main = sprintf("Pattern %d (%d genes): %s",
                             pp, length(genes),
                             paste(head(genes, 4), collapse = ", ")),
              plot = TRUE)
}
dev.off()
cat("Saved meringue_tumor_patterns_smooth.png\n")

# ---- Spatial neighbor network ----
png(file.path(OUT, "meringue_tumor_network.png"),
    width = 900, height = 900, res = 130)
par(mar = c(2, 2, 3, 1))
plotNetwork(pos, w, main = "Spatial neighbor network (filterDist = 50 um)",
            line.col = "grey70")
dev.off()
cat("Saved meringue_tumor_network.png\n")

# ---- Peak/valley FOV overlay using base plot ----
peak_fovs   <- c(224, 209, 201, 186, 172, 225, 229, 215, 175,
                 176, 206, 181, 168, 167)
valley_fovs <- c(199, 213, 227, 188, 204, 218, 178, 165, 170, 231, 189,
                 173, 154)
pv <- ifelse(m$fov %in% peak_fovs, "peak",
             ifelse(m$fov %in% valley_fovs, "valley", "other"))
col <- c(peak = "#e41a1c", valley = "#377eb8", other = "grey85")[pv]
png(file.path(OUT, "meringue_tumor_pv_overlay.png"),
    width = 900, height = 900, res = 130)
par(mar = c(3, 3, 3, 1))
plot(pos, col = col, pch = 19, cex = 0.4, asp = 1,
     xlab = "x (mm)", ylab = "y (mm)",
     main = "Tumor cells: peak (red) vs valley (blue) FOVs")
legend("topright", legend = c("peak", "valley", "other"),
       col = c("#e41a1c", "#377eb8", "grey85"), pch = 19, bty = "n")
dev.off()
cat("Saved meringue_tumor_pv_overlay.png\n")
