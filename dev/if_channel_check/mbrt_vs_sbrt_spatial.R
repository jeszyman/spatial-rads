suppressPackageStartupMessages({
  library(arrow); library(data.table); library(Matrix)
  library(FNN); library(ggplot2); library(patchwork)
})

OUT <- "/home/jeszyman/repos/spatial-rads/dev/if_channel_check"

meta <- as.data.table(read_parquet("/tmp/mutter01_meta.parquet"))
counts_df <- as.data.table(read_parquet("/tmp/mutter01_counts.parquet"))

# Helper: compute Moran's I on one block
analyze_block <- function(block_id, label) {
  m <- meta[Slide == "20250529_214712_S4" & Block == block_id]
  c_df <- counts_df[cell_id %in% m$cell_id]
  setkeyv(m, "cell_id"); setkeyv(c_df, "cell_id")
  m <- m[cell_id %in% c_df$cell_id]
  cat(sprintf("\n%s (%s): %d cells, %d FOVs\n", label, block_id, nrow(m),
              length(unique(m$fov))))

  gene_cols <- setdiff(colnames(c_df), c("Slide", "fov", "cell_id"))
  mat <- t(as.matrix(c_df[, ..gene_cols]))
  colnames(mat) <- c_df$cell_id; rownames(mat) <- gene_cols
  stopifnot(all(m$cell_id == colnames(mat)))

  lib <- colSums(mat); lib[lib == 0] <- 1
  normed <- log2(t(t(mat) / lib) * 1e4 + 1)
  pct_exp <- rowSums(mat > 0) / ncol(mat)
  normed <- normed[pct_exp > 0.03, ]

  pos <- as.matrix(m[, .(x_slide_mm, y_slide_mm)])
  rownames(pos) <- m$cell_id

  list(meta = m, normed = normed, pos = pos, label = label)
}

morans_I_quick <- function(normed, pos, k = 30) {
  nn <- get.knn(pos, k = k)
  n <- nrow(pos)
  ii <- rep(seq_len(n), each = k); jj <- as.vector(t(nn$nn.index))
  W <- sparseMatrix(i = ii, j = jj, x = 1, dims = c(n, n))
  W <- (W + t(W)) > 0; W <- as(W, "dMatrix") * 1
  row_sums <- rowSums(W); row_sums[row_sums == 0] <- 1
  Wn <- W / row_sums
  center <- rowMeans(normed)
  X <- normed - center
  Wx <- X %*% Wn
  numer <- rowSums(X * Wx)
  denom <- rowSums(X * X); denom[denom == 0] <- NA
  data.table(gene = rownames(normed), morans_I = numer / denom)[order(-morans_I)]
}

# ---- Subsample for speed (large enough to render clearly) ----
set.seed(1)
subsample <- function(obj, n_sub = 15000) {
  if (ncol(obj$normed) <= n_sub) return(obj)
  idx <- sample(ncol(obj$normed), n_sub)
  obj$normed <- obj$normed[, idx]
  obj$pos <- obj$pos[idx, ]
  obj$meta <- obj$meta[idx]
  obj
}

mbrt <- subsample(analyze_block("Block_21", "MBRT 4h"))
sbrt <- subsample(analyze_block("Block_17", "SBRT 4h"))

cat("\nComputing Moran's I (k=30)...\n")
mbrt$I <- morans_I_quick(mbrt$normed, mbrt$pos, k = 30)
sbrt$I <- morans_I_quick(sbrt$normed, sbrt$pos, k = 30)

cat("\n=== Top 15 Moran's I: MBRT 4h ===\n"); print(head(mbrt$I, 15))
cat("\n=== Top 15 Moran's I: SBRT 4h ===\n"); print(head(sbrt$I, 15))

# Compare matched genes
common <- intersect(mbrt$I$gene, sbrt$I$gene)
cmp <- merge(mbrt$I[gene %in% common], sbrt$I[gene %in% common],
             by = "gene", suffixes = c("_mbrt", "_sbrt"))
cmp[, delta := morans_I_mbrt - morans_I_sbrt]
setorder(cmp, -delta)
fwrite(cmp, "/tmp/moran_mbrt_vs_sbrt_4h.tsv", sep = "\t")
cat("\n=== Genes with greatest MBRT-SBRT Moran's I difference (top 15) ===\n")
print(head(cmp, 15))

# ---- Spatial plot: one gene per row, MBRT | SBRT side by side ----
plot_one <- function(obj, gene, point_size = 1.5, alpha = 0.85) {
  if (!gene %in% rownames(obj$normed)) return(ggplot() + theme_void())
  expr <- obj$normed[gene, ]
  cap <- quantile(expr, 0.99)
  df <- data.table(x = obj$pos[, 1], y = obj$pos[, 2], expr = pmin(expr, cap))
  I_val <- obj$I[gene == gene, morans_I][1]
  ggplot(df, aes(x, y, color = expr)) +
    geom_point(size = point_size, alpha = alpha) +
    scale_color_viridis_c(option = "magma", name = NULL) +
    coord_fixed() + theme_void() +
    theme(legend.position = "right", legend.key.height = unit(0.4, "cm"),
          plot.title = element_text(size = 11, hjust = 0.5)) +
    labs(title = sprintf("%s — %s (I=%.3f)", gene, obj$label, I_val))
}

genes_to_show <- c("S100a6", "Lgals1", "Vim", "Cdkn1a", "Cxcl10", "Arg1")
rows <- lapply(genes_to_show, function(g) {
  plot_one(mbrt, g) + plot_one(sbrt, g)
})
big <- wrap_plots(rows, ncol = 1)
ggsave(file.path(OUT, "mbrt_vs_sbrt_spatial.png"), big,
       width = 14, height = 4 * length(genes_to_show), dpi = 120)
cat("\nSaved", file.path(OUT, "mbrt_vs_sbrt_spatial.png"), "\n")

# ---- Single big plot: S100a6 MBRT vs SBRT with much bigger markers ----
g <- "S100a6"
p1 <- plot_one(mbrt, g, point_size = 2, alpha = 0.9)
p2 <- plot_one(sbrt, g, point_size = 2, alpha = 0.9)
ggsave(file.path(OUT, "S100a6_mbrt_vs_sbrt.png"), p1 + p2,
       width = 16, height = 7, dpi = 150)
cat("Saved S100a6 side-by-side.\n")
