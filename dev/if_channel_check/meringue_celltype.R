suppressPackageStartupMessages({
  library(arrow); library(data.table); library(Matrix); library(MERINGUE)
})

OUT <- "/home/jeszyman/repos/spatial-rads/dev/if_channel_check"
BLOCK <- "Block_21"; SLIDE <- "20250529_214712_S4"
N_SUB <- 3500
set.seed(1)

args <- commandArgs(trailingOnly = TRUE)
CT <- if (length(args) >= 1) args[1] else "tumor"

classify <- function(x) {
  out <- rep(NA_character_, length(x))
  out[x == "a"] <- "tumor"
  out[x %in% c("Macrophage","PerC.macrophage","spleen.red.pulp.macs",
               "small.peritoneal.macs","Microglia")] <- "macrophage"
  out[x %in% c("B.cell","Memory.B","Plasma","Plasmablast",
               "GC_centroblasts","GC_centrocyes","Spleen.CD19")] <- "bcell"
  out[x %in% c("CD8.T.cell","Spleen.Naive.CD4","Spleen.Naive.CD8",
               "Spleen.CD4Act.48hrs","Spleen.LN.Naive.CD4",
               "Spleen.Treg","Colon.Treg.Nrplo")] <- "tcell"
  out
}

peak_fovs   <- c(224, 209, 201, 186, 172, 225, 229, 215, 175,
                 176, 206, 181, 168, 167)
valley_fovs <- c(199, 213, 227, 188, 204, 218, 178, 165, 170, 231, 189,
                 173, 154)

meta <- as.data.table(read_parquet("/tmp/mutter01_meta.parquet"))
counts_df <- as.data.table(read_parquet("/tmp/mutter01_counts.parquet"))
setnames(meta, "ImmuneAtlas_ImmGen_Main_cell_Types", "main_type")
meta[, ct := classify(main_type)]

m <- meta[Slide == SLIDE & Block == BLOCK & ct == CT]
cat(CT, "cells in Block_21:", nrow(m), "\n")
if (nrow(m) < 300) stop("Too few cells")
counts_df <- counts_df[cell_id %in% m$cell_id]
setkeyv(m, "cell_id"); setkeyv(counts_df, "cell_id")
m <- m[cell_id %in% counts_df$cell_id]

idx <- sample(nrow(m), min(N_SUB, nrow(m)))
m <- m[idx]; counts_df <- counts_df[cell_id %in% m$cell_id]
setorderv(m, "cell_id"); setorderv(counts_df, "cell_id")
cat("Subsampled:", nrow(m), "\n")

gene_cols <- setdiff(colnames(counts_df), c("Slide", "fov", "cell_id"))
mat <- t(as.matrix(counts_df[, ..gene_cols]))
colnames(mat) <- counts_df$cell_id

lib <- colSums(mat); lib[lib == 0] <- 1
normed <- log2(t(t(mat) / lib) * 1e4 + 1)
pct_exp <- rowSums(mat > 0) / ncol(mat)
normed <- normed[pct_exp > 0.02, ]
cat("Genes after filter:", nrow(normed), "\n")

pos <- as.matrix(m[, .(x_slide_mm, y_slide_mm)])
rownames(pos) <- m$cell_id
w <- getSpatialNeighbors(pos, filterDist = 0.05)
cat("Mean neighbors:", round(mean(rowSums(w > 0)), 1), "\n")

# Moran's I
I <- getSpatialPatterns(normed, w)

# Relaxed filter: raw p-value, looser alpha
sig_genes <- rownames(I)[!is.na(I$p.value) & I$p.value < 0.01]
# keep only genes with minPercentCells
keep <- rownames(normed)[rowSums(normed > 0) / ncol(normed) > 0.02]
sig_genes <- intersect(sig_genes, keep)
cat("Significant spatial genes (raw p<0.01, pct>2%):", length(sig_genes), "\n")
if (length(sig_genes) < 10) stop("Too few genes")

# Spatial cross-correlation and grouping
scc <- spatialCrossCorMatrix(mat = as.matrix(normed[sig_genes, ]), weight = w)
ggroup <- groupSigSpatialPatterns(pos = pos,
                                  mat = as.matrix(normed[sig_genes, ]),
                                  scc = scc,
                                  power = 1,
                                  hclustMethod = "ward.D",
                                  deepSplit = 2,
                                  zlim = c(-1.5, 1.5),
                                  plot = FALSE)
grp <- data.table(gene = names(ggroup$groups),
                  pattern = as.integer(as.character(ggroup$groups)))
setorder(grp, pattern, gene)
cat("Pattern sizes:\n"); print(table(grp$pattern))

# Smoothing operator for viz
w_smooth <- getSpatialNeighbors(pos, filterDist = 0.15)
w_smooth <- w_smooth + Matrix::Diagonal(nrow(w_smooth))
rs <- Matrix::rowSums(w_smooth); rs[rs == 0] <- 1
w_smooth <- w_smooth / rs
smooth_expr <- function(vec) as.numeric(w_smooth %*% vec)

# Per-pattern Moran's I on mean score
pattern_summary <- rbindlist(lapply(sort(unique(grp$pattern)), function(pp) {
  genes <- grp[pattern == pp, gene]
  expr <- colMeans(normed[genes, , drop = FALSE])
  mt <- moranTest(expr, weight = w)
  # Peak/valley scores
  pv <- ifelse(m$fov %in% peak_fovs, "peak",
               ifelse(m$fov %in% valley_fovs, "valley", "other"))
  peak_mean <- mean(expr[pv == "peak"])
  valley_mean <- mean(expr[pv == "valley"])
  data.table(pattern = pp, n_genes = length(genes),
             morans_I = mt["observed"], p_value = mt["p.value"],
             peak = peak_mean, valley = valley_mean,
             peak_minus_valley = peak_mean - valley_mean,
             top_genes = paste(head(genes, 8), collapse = ", "))
}))
fwrite(pattern_summary, sprintf("/tmp/meringue_%s_patterns.tsv", CT), sep = "\t")
fwrite(grp, sprintf("/tmp/meringue_%s_gene_pattern.tsv", CT), sep = "\t")
cat("\n=== Pattern summary (", CT, ") ===\n", sep = "")
print(pattern_summary[, .(pattern, n_genes, morans_I = round(morans_I, 3),
                          p_value, peak = round(peak, 3),
                          valley = round(valley, 3),
                          peak_minus_valley = round(peak_minus_valley, 3))])

# Plot
patterns <- sort(unique(grp$pattern))
nrow_plot <- ceiling(length(patterns) / 3)
png(sprintf("%s/meringue_%s_smooth.png", OUT, CT),
    width = 1500, height = 500 * nrow_plot, res = 130)
par(mfrow = c(nrow_plot, 3), mar = c(2, 2, 3.2, 1))
for (pp in patterns) {
  genes <- grp[pattern == pp, gene]
  expr <- smooth_expr(colMeans(normed[genes, , drop = FALSE]))
  names(expr) <- colnames(normed)
  mt <- pattern_summary[pattern == pp]
  title <- sprintf("Pattern %d (%d genes, I=%.2f, p=%.1e)\n%s",
                   pp, length(genes), mt$morans_I, mt$p_value,
                   substr(mt$top_genes, 1, 70))
  interpolate(pos, expr, scale = TRUE, fill = TRUE,
              zlim = c(-2, 2), binSize = 60,
              main = title, plot = TRUE)
}
dev.off()
cat("Saved", sprintf("meringue_%s_smooth.png", CT), "\n")

# Peak/valley overlay
pv <- ifelse(m$fov %in% peak_fovs, "peak",
             ifelse(m$fov %in% valley_fovs, "valley", "other"))
col <- c(peak = "#e41a1c", valley = "#377eb8", other = "grey85")[pv]
png(sprintf("%s/meringue_%s_pv_overlay.png", OUT, CT),
    width = 900, height = 900, res = 130)
par(mar = c(3, 3, 3, 1))
plot(pos, col = col, pch = 19, cex = 0.4, asp = 1,
     xlab = "x (mm)", ylab = "y (mm)",
     main = sprintf("%s cells: peak (red) / valley (blue)", CT))
legend("topright", legend = c("peak", "valley", "other"),
       col = c("#e41a1c", "#377eb8", "grey85"), pch = 19, bty = "n")
dev.off()
cat("Saved", sprintf("meringue_%s_pv_overlay.png", CT), "\n")
