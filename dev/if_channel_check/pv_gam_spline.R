suppressPackageStartupMessages({
  library(arrow); library(data.table); library(Matrix); library(mgcv)
  library(ggplot2); library(ggrepel); library(patchwork)
})

OUT <- "/home/jeszyman/repos/spatial-rads/dev/if_channel_check"
SLIDE <- "20250529_214712_S4"; MBRT_BLOCK <- "Block_21"
N_CELLS <- 15000
set.seed(1)

peak_fovs <- c(224, 209, 201, 186, 172, 225, 229, 215, 175,
               176, 206, 181, 168, 167)

meta <- as.data.table(read_parquet("/tmp/mutter01_meta.parquet"))
counts_df <- as.data.table(read_parquet("/tmp/mutter01_counts.parquet"))
setnames(meta, "ImmuneAtlas_ImmGen_Main_cell_Types", "main_type")

m <- meta[Slide == SLIDE & Block == MBRT_BLOCK & main_type == "a"]
c_df <- counts_df[cell_id %in% m$cell_id]
setkeyv(m, "cell_id"); setkeyv(c_df, "cell_id")
m <- m[cell_id %in% c_df$cell_id]
setorderv(m, "cell_id"); setorderv(c_df, "cell_id")
cat("a-bucket cells:", nrow(m), "\n")

gene_cols <- setdiff(colnames(c_df), c("Slide","fov","cell_id"))
mat_full <- t(as.matrix(c_df[, ..gene_cols]))
lib_full <- colSums(mat_full); lib_full[lib_full == 0] <- 1
normed_full <- log2(t(t(mat_full) / lib_full) * 1e4 + 1)

# ---- 4T1 tumor filter: epithelial-positive AND endothelial-negative ----
epith_markers <- intersect(c("Krt8","Krt18","Epcam","Cdh1"), rownames(normed_full))
endo_markers  <- intersect(c("Pecam1","Cdh5","Vwf"), rownames(normed_full))
cat("Epith markers present:", paste(epith_markers, collapse=","), "\n")
cat("Endo markers present:", paste(endo_markers, collapse=","), "\n")
epith_score <- colMeans(normed_full[epith_markers, , drop = FALSE])
endo_score  <- if (length(endo_markers) > 0) colMeans(normed_full[endo_markers, , drop = FALSE]) else rep(0, ncol(normed_full))
is_4t1 <- epith_score > 0.3 & endo_score < 0.1
cat(sprintf("4T1-positive cells: %d / %d (%.0f%%)\n",
            sum(is_4t1), length(is_4t1), 100 * mean(is_4t1)))

m <- m[is_4t1]; c_df <- c_df[is_4t1]
normed_full <- normed_full[, is_4t1]

if (nrow(m) > N_CELLS) {
  idx <- sample(nrow(m), N_CELLS)
  m <- m[idx]; c_df <- c_df[idx]; normed_full <- normed_full[, idx]
}
setorderv(m, "cell_id"); setorderv(c_df, "cell_id")
cat("Cells after subsample:", nrow(m), "\n")

# Filter genes by pct_expr on subsetted, normalized matrix
mat_sub <- t(as.matrix(c_df[, ..gene_cols]))
pct_expr <- rowSums(mat_sub > 0) / ncol(mat_sub)
normed <- normed_full[pct_expr > 0.05, ]
cat("Genes:", nrow(normed), "\n")

pos <- as.matrix(m[, .(x_slide_mm, y_slide_mm)])

# Fit stripe geometry, label cells
peak_centroids <- meta[Slide == SLIDE & Block == MBRT_BLOCK & fov %in% peak_fovs,
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
beam_spacing <- median(diff(fit$centers))
peak_centroids[, d_perp := -sin(theta) * cx + cos(theta) * cy]
peak_centroids[, stripe := fit$cluster]
within_sd <- peak_centroids[, .(sd = sd(d_perp - fit$centers[stripe])),
                            by = stripe][, mean(sd)]
peak_half <- within_sd + 0.15

d_cells <- -sin(theta) * pos[,1] + cos(theta) * pos[,2]
d_range <- range(d_cells)
centers_real <- seq(fit$centers[1] -
                      ceiling((fit$centers[1] - d_range[1]) / beam_spacing) * beam_spacing,
                    fit$centers[length(fit$centers)] +
                      ceiling((d_range[2] - fit$centers[length(fit$centers)]) / beam_spacing) * beam_spacing,
                    by = beam_spacing)
dist_to_peak <- sapply(d_cells, function(x) min(abs(x - centers_real)))
label <- ifelse(dist_to_peak < peak_half, "peak",
                ifelse(dist_to_peak > (beam_spacing/2 - peak_half),
                       "valley", "transition"))
keep <- label %in% c("peak", "valley")
cat("Cells kept (peak+valley):", sum(keep), "\n")

df_base <- data.table(x = pos[keep, 1], y = pos[keep, 2],
                      pv = factor(label[keep], levels = c("valley", "peak")))
normed_keep <- normed[, keep]

# ---- Per-gene 2D GAM ----
# For each gene, fit: expr ~ s(x, y, k=50) + pv
# Record pv coefficient, its SE, and p-value.
# This asks: after smoothing out any tumor-wide spatial trend, does peak vs valley
# still add explanatory power?
K_BASIS <- 50  # spatial smoother complexity
n_genes <- nrow(normed_keep)
cat(sprintf("Fitting %d GAMs with s(x, y, k=%d)...\n", n_genes, K_BASIS))

res <- data.table(gene = rownames(normed_keep),
                  raw_delta = NA_real_,
                  adj_delta = NA_real_,
                  se = NA_real_,
                  p = NA_real_,
                  edf_spatial = NA_real_)

t0 <- Sys.time()
for (i in seq_len(n_genes)) {
  df <- copy(df_base)
  df$expr <- normed_keep[i, ]
  res$raw_delta[i] <- mean(df$expr[df$pv == "peak"]) -
                     mean(df$expr[df$pv == "valley"])
  # bam() is fast for n > 1000
  fit_gam <- tryCatch(
    bam(expr ~ s(x, y, k = K_BASIS) + pv, data = df, discrete = TRUE),
    error = function(e) NULL)
  if (is.null(fit_gam)) next
  s <- summary(fit_gam)
  # Coefficient on pvpeak (the peak-vs-valley effect after spline)
  pv_row <- rownames(s$p.table)[grepl("pv", rownames(s$p.table))][1]
  if (is.na(pv_row)) next
  res$adj_delta[i] <- s$p.table[pv_row, "Estimate"]
  res$se[i] <- s$p.table[pv_row, "Std. Error"]
  res$p[i] <- s$p.table[pv_row, "Pr(>|t|)"]
  res$edf_spatial[i] <- s$s.table[1, "edf"]
  if (i %% 50 == 0) {
    dt_elapsed <- as.numeric(Sys.time() - t0, units = "secs")
    eta <- dt_elapsed / i * (n_genes - i)
    cat(sprintf("  %d/%d  |  %.1fs elapsed  |  ETA %.0fs\n",
                i, n_genes, dt_elapsed, eta))
  }
}
res[, p_adj := p.adjust(p, method = "BH")]
setorder(res, p)
fwrite(res, "/tmp/tumor_pv_gam.tsv", sep = "\t")

cat("\n=== Top 25 genes (smallest p after 2D-spline adjustment) ===\n")
print(res[1:25, .(gene, raw_delta = round(raw_delta, 3),
                  adj_delta = round(adj_delta, 3),
                  se = round(se, 3),
                  p = sprintf("%.1e", p),
                  p_adj = sprintf("%.1e", p_adj),
                  edf = round(edf_spatial, 1))])

hits <- res[p_adj < 0.05]
cat(sprintf("\nHits (BH-adjusted p < 0.05): %d\n", nrow(hits)))
cat(sprintf("Peak-enriched hits: %d | Valley-enriched hits: %d\n",
            sum(hits$adj_delta > 0), sum(hits$adj_delta < 0)))

# Volcano
res[, neg_log_p := -log10(pmax(p, 1e-300))]
res[, category := ifelse(p_adj < 0.05,
                         ifelse(adj_delta > 0, "peak (BH p<0.05)", "valley (BH p<0.05)"),
                         "n.s.")]
top_labels <- c(res[adj_delta > 0][order(-adj_delta)][1:15, gene],
                res[adj_delta < 0][order(adj_delta)][1:5, gene])
res[, label := ifelse(gene %in% top_labels, gene, NA_character_)]

p_volc <- ggplot(res, aes(x = adj_delta, y = neg_log_p, color = category)) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.75, size = 1.1) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 50,
                  color = "black", min.segment.length = 0,
                  box.padding = 0.3, seed = 1) +
  scale_color_manual(values = c("peak (BH p<0.05)" = "#e41a1c",
                                "valley (BH p<0.05)" = "#377eb8",
                                "n.s." = "grey60")) +
  labs(x = "Architecture-adjusted MBRT peak − valley effect (from GAM with 2D spline)",
       y = "−log10(raw GAM p-value)",
       color = NULL,
       title = "Volcano: tumor peak/valley effect after adjusting for 2D spatial architecture",
       subtitle = paste0(
         "Each point is a gene. Model: expr ~ s(x, y, k=50) + pv. ",
         "The 2D spline absorbs tumor architecture; remaining pv effect = architecture-adjusted peak/valley.\n",
         "Dashed line = unadjusted p=0.05. Colored = BH-adjusted p<0.05.")) +
  theme_minimal() +
  theme(legend.position = "bottom", plot.title = element_text(size = 11),
        plot.subtitle = element_text(size = 9))

# Raw-vs-adjusted scatter
p_comp <- ggplot(res, aes(x = raw_delta, y = adj_delta, color = category)) +
  geom_abline(slope = 1, linetype = "dotted", color = "grey50") +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_point(alpha = 0.7, size = 0.9) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 50,
                  color = "black", min.segment.length = 0, seed = 1) +
  scale_color_manual(values = c("peak (BH p<0.05)" = "#e41a1c",
                                "valley (BH p<0.05)" = "#377eb8",
                                "n.s." = "grey60")) +
  labs(x = "Raw peak − valley delta (naive)",
       y = "Architecture-adjusted delta (from GAM)",
       color = NULL,
       title = "Raw vs architecture-adjusted effect size",
       subtitle = "Points below the dotted line lost magnitude after adjustment = partly architecture.") +
  theme_minimal() +
  theme(legend.position = "bottom", plot.title = element_text(size = 11),
        plot.subtitle = element_text(size = 9))

combined <- p_volc / p_comp + plot_layout(heights = c(1.3, 1))
ggsave(file.path(OUT, "tumor_pv_gam_4t1.png"), combined,
       width = 12, height = 14, dpi = 130)
cat("Saved tumor_pv_gam_4t1.png\n")
fwrite(res, "/tmp/tumor_pv_gam_4t1.tsv", sep = "\t")
