suppressPackageStartupMessages({
  library(arrow); library(data.table); library(Matrix); library(Seurat)
  library(ggplot2); library(ggrepel)
})

OUT <- "/home/jeszyman/repos/spatial-rads/dev/if_channel_check"
SLIDE <- "20250529_214712_S4"; BLOCK <- "Block_21"
set.seed(1)

peak_fovs <- c(224, 209, 201, 186, 172, 225, 229, 215, 175,
               176, 206, 181, 168, 167)

meta <- as.data.table(read_parquet("/tmp/mutter01_meta.parquet"))
counts_df <- as.data.table(read_parquet("/tmp/mutter01_counts.parquet"))
setnames(meta, "ImmuneAtlas_ImmGen_Main_cell_Types", "main_type")

# Block_21 "a" cells
m <- meta[Slide == SLIDE & Block == BLOCK & main_type == "a"]
c_df <- counts_df[cell_id %in% m$cell_id]
setkeyv(m, "cell_id"); setkeyv(c_df, "cell_id")
m <- m[cell_id %in% c_df$cell_id]
setorderv(m, "cell_id"); setorderv(c_df, "cell_id")
cat("a-bucket cells:", nrow(m), "\n")

gene_cols <- setdiff(colnames(c_df), c("Slide","fov","cell_id"))
mat <- t(as.matrix(c_df[, ..gene_cols]))
colnames(mat) <- c_df$cell_id
lib <- colSums(mat); lib[lib == 0] <- 1
normed <- log2(t(t(mat) / lib) * 1e4 + 1)

# 4T1 tumor filter
epith <- intersect(c("Krt8","Krt18","Epcam","Cdh1"), rownames(normed))
endo  <- intersect(c("Pecam1","Cdh5","Vwf"), rownames(normed))
ep_s  <- colMeans(normed[epith, , drop = FALSE])
en_s  <- if (length(endo) > 0) colMeans(normed[endo, , drop = FALSE]) else rep(0, ncol(normed))
keep4t1 <- ep_s > 0.3 & en_s < 0.1
m <- m[keep4t1]; mat <- mat[, keep4t1]
cat("4T1 cells:", ncol(mat), "\n")

pos <- as.matrix(m[, .(x_slide_mm, y_slide_mm)])

# Stripe geometry fit + cell label extrapolation
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
cat("Peak:", sum(label == "peak"), " Valley:", sum(label == "valley"),
    " Dropped transition:", sum(label == "transition"), "\n")

mat <- mat[, keep]; m <- m[keep]; label <- label[keep]

# Seurat FindMarkers, Wilcoxon, peak vs valley
obj <- CreateSeuratObject(counts = mat)
obj <- NormalizeData(obj, verbose = FALSE)
obj$pv <- factor(label, levels = c("valley","peak"))
Idents(obj) <- obj$pv
cat("Running FindMarkers (Wilcoxon, peak vs valley)...\n")
de <- FindMarkers(obj, ident.1 = "peak", ident.2 = "valley",
                  test.use = "wilcox", logfc.threshold = 0,
                  min.pct = 0.05, verbose = FALSE)
de$gene <- rownames(de)
de <- as.data.table(de)[, .(gene, avg_log2FC, pct.1, pct.2, p_val, p_val_adj)]
setorder(de, -avg_log2FC)
fwrite(de, "/tmp/tumor_pv_simple_dge.tsv", sep = "\t")

cat("\n=== Top 25 peak-UP (avg_log2FC > 0) ===\n")
print(de[avg_log2FC > 0][1:25])
cat("\n=== Top 15 valley-UP (avg_log2FC < 0) ===\n")
print(de[order(avg_log2FC)][1:15])
cat(sprintf("\nTotal genes tested: %d\n", nrow(de)))
cat(sprintf("FDR < 0.05 & |log2FC| > 0.1:  peak-UP = %d, valley-UP = %d\n",
            sum(de$p_val_adj < 0.05 & de$avg_log2FC >  0.1),
            sum(de$p_val_adj < 0.05 & de$avg_log2FC < -0.1)))

# Volcano
de[, neg_log_p := -log10(pmax(p_val_adj, 1e-300))]
de[, sig := fifelse(p_val_adj < 0.05 & abs(avg_log2FC) > 0.1,
                    fifelse(avg_log2FC > 0, "peak-UP", "valley-UP"),
                    "n.s.")]
top_up   <- de[avg_log2FC > 0.1 & p_val_adj < 0.05][order(-avg_log2FC)][1:20, gene]
top_down <- de[avg_log2FC < -0.1 & p_val_adj < 0.05][order(avg_log2FC)][1:10, gene]
de[, label := ifelse(gene %in% c(top_up, top_down), gene, NA_character_)]

p <- ggplot(de, aes(x = avg_log2FC, y = neg_log_p, color = sig)) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dotted", color = "grey60") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.8, size = 1.1) +
  geom_text_repel(aes(label = label), size = 3.2, max.overlaps = 40,
                  color = "black", min.segment.length = 0, box.padding = 0.3,
                  seed = 1) +
  scale_color_manual(values = c("peak-UP" = "#e41a1c",
                                "valley-UP" = "#377eb8",
                                "n.s." = "grey70")) +
  labs(x = "avg log2FC (peak vs valley)",
       y = "−log10 BH-adjusted p",
       color = NULL,
       title = "4T1 tumor cells: peak vs valley at 4h MBRT (Block_21)",
       subtitle = paste0(
         "Simple Wilcoxon DGE. n_peak=", sum(label == "peak"),
         ", n_valley=", sum(label == "valley"),
         ". Cutoffs: BH p_adj<0.05 AND |log2FC|>0.1.")) +
  theme_minimal() +
  theme(legend.position = "bottom", plot.title = element_text(size = 12))

ggsave(file.path(OUT, "tumor_pv_simple_dge.png"), p, width = 11, height = 8, dpi = 130)
cat("Saved tumor_pv_simple_dge.png\n")
