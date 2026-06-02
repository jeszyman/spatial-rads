suppressPackageStartupMessages({
  library(arrow); library(data.table); library(Matrix)
  library(FNN); library(ggplot2); library(patchwork)
})

OUT <- "/home/jeszyman/repos/spatial-rads/dev/if_channel_check"
cat("Loading Block_21...\n")

meta <- as.data.table(read_parquet("/tmp/mutter01_meta.parquet"))
s4 <- meta[Slide == "20250529_214712_S4" & Block == "Block_21"]

counts_df <- as.data.table(read_parquet("/tmp/mutter01_counts.parquet"))
counts_df <- counts_df[cell_id %in% s4$cell_id]
setkey(counts_df, cell_id)
setkey(s4, cell_id)
s4 <- s4[cell_id %in% counts_df$cell_id]

# Subsample to 8000 cells
set.seed(1)
n_sub <- 8000
idx <- sample(nrow(s4), n_sub)
s4 <- s4[idx]
counts_df <- counts_df[cell_id %in% s4$cell_id]
setorderv(s4, "cell_id"); setorderv(counts_df, "cell_id")

gene_cols <- setdiff(colnames(counts_df), c("Slide", "fov", "cell_id"))
mat <- t(as.matrix(counts_df[, ..gene_cols]))
colnames(mat) <- counts_df$cell_id
rownames(mat) <- gene_cols

lib <- colSums(mat); lib[lib == 0] <- 1
normed <- log2(t(t(mat) / lib) * 1e4 + 1)
pct_exp <- rowSums(mat > 0) / ncol(mat)
normed <- normed[pct_exp > 0.03, ]
cat(sprintf("Cells: %d, genes: %d\n", ncol(normed), nrow(normed)))

# ---- Build sparse kNN weight matrix ----
cat("Building kNN graph...\n")
pos <- as.matrix(s4[, .(x_slide_mm, y_slide_mm)])
k <- 10
nn <- get.knn(pos, k = k)
n <- nrow(pos)
ii <- rep(seq_len(n), each = k)
jj <- as.vector(t(nn$nn.index))
W <- sparseMatrix(i = ii, j = jj, x = 1, dims = c(n, n))
W <- (W + t(W)) > 0  # symmetrize
W <- as(W, "dMatrix") * 1
row_sums <- rowSums(W); row_sums[row_sums == 0] <- 1
Wn <- W / row_sums  # row-normalized
cat("Weight matrix:", n, "x", n, "nnz:", length(W@x), "\n")

# ---- Moran's I per gene: I = (x' Wn x) / (x' x), x centered ----
cat("Computing Moran's I per gene...\n")
center <- rowMeans(normed)
X <- normed - center  # gene x cell
# For each gene, I = sum_i x_i * sum_j W_ij x_j / sum_i x_i^2
# Matrix form: Wx = X %*% t(Wn)  (rows are genes, cols are cells)
# But Wn is cell x cell; we want X %*% Wn^T for each gene
Wx <- X %*% Wn  # genes x cells
numer <- rowSums(X * Wx)
denom <- rowSums(X * X); denom[denom == 0] <- NA
morans_I <- numer / denom
I_df <- data.table(gene = rownames(normed), morans_I = morans_I)
setorder(I_df, -morans_I)

cat("Top 30 by Moran's I:\n")
print(head(I_df, 30))
fwrite(I_df, "/tmp/meringue_morans_I_block21.tsv", sep = "\t")

# ---- Spotlight DGE hits ----
dge_peak_up   <- c("Cdkn1a", "Ifitm3", "S100a6", "Lgals1", "Arg1", "Ctsd", "Col6a1", "B2m", "Slpi", "Junb", "Il24")
dge_valley_up <- c("Cxcl10", "Oasl1", "Igfbp6", "Ndrg1", "Il1a", "Cst7", "Apoe", "Ifnl2/3", "H2-Q10")
spotlight <- intersect(c(dge_peak_up, dge_valley_up), I_df$gene)
cat("\nDGE hits Moran's I:\n")
print(I_df[gene %in% spotlight][order(-morans_I)])

# ---- Spatial plots ----
plot_gene <- function(g, label) {
  if (!g %in% rownames(normed)) return(ggplot() + theme_void() + labs(title = paste(g, "NA")))
  expr <- normed[g, ]
  cap <- quantile(expr, 0.98)
  df <- data.table(x = pos[, 1], y = pos[, 2], expr = pmin(expr, cap))
  ggplot(df, aes(x, y, color = expr)) +
    geom_point(size = 0.3, alpha = 0.7) +
    scale_color_viridis_c(option = "magma") +
    coord_fixed() + theme_void() +
    theme(legend.position = "bottom", legend.key.width = unit(0.6, "cm"),
          plot.title = element_text(size = 9, hjust = 0.5)) +
    labs(title = sprintf("%s (%s, I=%.3f)", g, label, I_df[gene == g, morans_I]))
}

peak_plots   <- lapply(intersect(dge_peak_up, rownames(normed)), plot_gene, label = "peak UP")
valley_plots <- lapply(intersect(dge_valley_up, rownames(normed)), plot_gene, label = "valley UP")
g_all <- wrap_plots(c(peak_plots, valley_plots), ncol = 4)
ggsave(file.path(OUT, "spatial_dge_hits_block21.png"), g_all,
       width = 18, height = 3 * ceiling(length(c(peak_plots, valley_plots)) / 4), dpi = 150)
cat("Saved spatial_dge_hits_block21.png\n")

top10 <- head(I_df, 10)$gene
tp <- lapply(top10, plot_gene, label = "top Moran's I")
g_top <- wrap_plots(tp, ncol = 4)
ggsave(file.path(OUT, "spatial_top10_morans_block21.png"), g_top, width = 18, height = 12, dpi = 150)
cat("Saved spatial_top10_morans_block21.png\n")

cat("\nDone.\n")
