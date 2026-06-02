suppressPackageStartupMessages({
  library(arrow); library(data.table); library(Matrix); library(ggplot2)
  library(ggrepel); library(patchwork); library(RANN)
})

OUT <- "/home/jeszyman/repos/spatial-rads/dev/if_channel_check"
SLIDE <- "20250529_214712_S4"; MBRT_BLOCK <- "Block_21"
set.seed(1)

peak_fovs <- c(224, 209, 201, 186, 172, 225, 229, 215, 175,
               176, 206, 181, 168, 167)

meta <- as.data.table(read_parquet("/tmp/mutter01_meta.parquet"))
counts_df <- as.data.table(read_parquet("/tmp/mutter01_counts.parquet"))
setnames(meta, "ImmuneAtlas_ImmGen_Main_cell_Types", "main_type")

# ---- Load MBRT Block_21 tumor cells ----
m <- meta[Slide == SLIDE & Block == MBRT_BLOCK & main_type == "a"]
c_df <- counts_df[cell_id %in% m$cell_id]
setkeyv(m, "cell_id"); setkeyv(c_df, "cell_id")
m <- m[cell_id %in% c_df$cell_id]
setorderv(m, "cell_id"); setorderv(c_df, "cell_id")

# Subsample to speed up (null requires re-computing labels 200 times)
N_CELLS <- 15000
if (nrow(m) > N_CELLS) {
  idx <- sample(nrow(m), N_CELLS)
  m <- m[idx]; c_df <- c_df[cell_id %in% m$cell_id]
  setorderv(m, "cell_id"); setorderv(c_df, "cell_id")
}
cat("Tumor cells:", nrow(m), "\n")

gene_cols <- setdiff(colnames(c_df), c("Slide","fov","cell_id"))
mat <- t(as.matrix(c_df[, ..gene_cols]))
lib <- colSums(mat); lib[lib == 0] <- 1
normed <- log2(t(t(mat) / lib) * 1e4 + 1)
pct_expr <- rowSums(mat > 0) / ncol(mat)
normed <- normed[pct_expr > 0.05, ]
cat("Genes kept:", nrow(normed), "\n")

pos <- as.matrix(m[, .(x_slide_mm, y_slide_mm)])

# ---- Real peak/valley labels (extrapolated stripe geometry) ----
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
real_label <- ifelse(dist_to_peak < peak_half, "peak",
                     ifelse(dist_to_peak > (beam_spacing/2 - peak_half),
                            "valley", "transition"))
n_peak <- sum(real_label == "peak")
cat("Real peak cells:", n_peak, " valley:", sum(real_label == "valley"), "\n")

compute_delta <- function(is_peak, normed) {
  if (sum(is_peak) < 100 || sum(!is_peak) < 100) return(NULL)
  rowMeans(normed[, is_peak, drop = FALSE]) -
    rowMeans(normed[, !is_peak, drop = FALSE])
}

real_delta <- compute_delta(real_label == "peak", normed)

# ---- Build spatial kNN for contiguous region growing ----
cat("Building kNN...\n")
k_nn <- 30
nn <- nn2(pos, pos, k = k_nn + 1)$nn.idx[, -1]

grow_random_region <- function(target_size, nn, n_seeds = 4) {
  # Start from n_seeds random cells, grow region breadth-first to target_size
  seeds <- sample(nrow(nn), n_seeds)
  in_region <- rep(FALSE, nrow(nn))
  in_region[seeds] <- TRUE
  frontier <- seeds
  while (sum(in_region) < target_size && length(frontier) > 0) {
    # Expand frontier: add neighbors not yet in region
    candidates <- unique(as.vector(nn[frontier, ]))
    candidates <- candidates[!in_region[candidates]]
    if (length(candidates) == 0) {
      # Pick a new random seed
      free <- which(!in_region)
      if (length(free) == 0) break
      new_seed <- sample(free, 1)
      in_region[new_seed] <- TRUE
      frontier <- new_seed
      next
    }
    # Add candidates up to target
    need <- target_size - sum(in_region)
    to_add <- head(candidates, need)
    in_region[to_add] <- TRUE
    frontier <- to_add
  }
  in_region
}

# ---- Null: random contiguous regions of same size as peak ----
n_null <- 200
cat("Generating", n_null, "random-region nulls...\n")
null_deltas <- matrix(NA_real_, nrow = nrow(normed), ncol = n_null)
rownames(null_deltas) <- rownames(normed)
for (i in seq_len(n_null)) {
  is_fake_peak <- grow_random_region(n_peak, nn, n_seeds = 4)
  d <- compute_delta(is_fake_peak, normed)
  if (!is.null(d)) null_deltas[, i] <- d
  if (i %% 25 == 0) cat("  iter", i, "\n")
}

# Empirical p-value per gene
emp_p <- sapply(seq_along(real_delta), function(g) {
  nd <- null_deltas[g, ]; nd <- nd[!is.na(nd)]
  if (length(nd) < 10) return(NA_real_)
  sum(abs(nd) >= abs(real_delta[g])) / length(nd)
})
names(emp_p) <- names(real_delta)
overall_null_sd <- sd(as.vector(null_deltas), na.rm = TRUE)
cat(sprintf("Random-region null overall sd: %.3f\n", overall_null_sd))

# Prior rotation null for comparison
rot <- fread("/tmp/tumor_pv_rotation_null.tsv")
out <- data.table(gene = names(real_delta),
                  mbrt_delta = real_delta,
                  rand_region_p = emp_p)
out <- merge(out, rot[, .(gene, rot_null_p)], by = "gene", all.x = TRUE)
out[, survives_rand := rand_region_p < 0.05 & abs(mbrt_delta) > 2 * overall_null_sd]
setorder(out, -mbrt_delta)
fwrite(out, "/tmp/tumor_pv_random_region_null.tsv", sep = "\t")

cat("\n=== Top 20 peak-enriched with random-region null ===\n")
print(out[1:20, .(gene, mbrt_delta = round(mbrt_delta, 3),
                  rand_region_p = sprintf("%.3f", rand_region_p),
                  rot_null_p = sprintf("%.3f", rot_null_p),
                  survives_rand)])
cat(sprintf("\nTotal hits (random-region p<0.05 AND |delta|>2*null sd): %d\n",
            sum(out$survives_rand, na.rm = TRUE)))

# Plot: volcano using random-region null
out[, neg_log_p := -log10(pmax(rand_region_p, 1/n_null))]
top_peak <- out[order(-mbrt_delta)][1:15, gene]
top_valley <- out[order(mbrt_delta)][1:5, gene]
out[, label := ifelse(gene %in% c(top_peak, top_valley), gene, NA_character_)]
out[, category := fifelse(survives_rand == TRUE,
                          fifelse(mbrt_delta > 0, "survives (peak)", "survives (valley)"),
                          "does not survive")]

# 3-curve density: real MBRT, random-region null, rotation null
rot_deltas <- rot$mbrt_delta  # all genes' real MBRT (same gene, same data)
densdata <- rbindlist(list(
  data.table(arm = sprintf("Random-region null (%d iters)", n_null),
             val = as.vector(null_deltas)),
  data.table(arm = "Real MBRT (θ=−25° stripes)",
             val = out$mbrt_delta)
))
p_dens <- ggplot(densdata, aes(x = val, fill = arm, color = arm)) +
  geom_density(alpha = 0.35, linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  scale_fill_manual(values = setNames(c("steelblue", "#e41a1c"), unique(densdata$arm))) +
  scale_color_manual(values = setNames(c("steelblue4", "#b30000"), unique(densdata$arm))) +
  labs(x = "Peak − Valley (log-norm)", y = "Density of genes",
       fill = NULL, color = NULL,
       title = "Random-region null vs real MBRT peak/valley deltas",
       subtitle = paste0(
         "Null = fake 'peak' blobs grown from random seeds to cover the same cell count as real peaks. ",
         "Keeps contiguous regions (not scattered cells), like stripes but arbitrary location.\n",
         "If real MBRT effects exceed this null, they're not explained by 'any random patch of tumor'.")) +
  theme_minimal() +
  theme(legend.position = "bottom", plot.title = element_text(size = 11),
        plot.subtitle = element_text(size = 9))

p_volc <- ggplot(out, aes(x = mbrt_delta, y = neg_log_p, color = category)) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  annotate("text", x = max(out$mbrt_delta) * 0.95, y = -log10(0.05) + 0.08,
           label = "p = 0.05", hjust = 1, size = 3, color = "grey40") +
  geom_point(alpha = 0.75, size = 1.1) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 50,
                  color = "black", min.segment.length = 0,
                  box.padding = 0.3, seed = 1) +
  scale_color_manual(values = c("survives (peak)" = "#e41a1c",
                                "survives (valley)" = "#377eb8",
                                "does not survive" = "grey60")) +
  labs(x = "Real MBRT peak − valley (log-norm)",
       y = "−log10(random-region null p)",
       color = NULL,
       title = "Volcano: MBRT peak-enrichment vs random-region null",
       subtitle = "Top-right = peak-enriched AND unlikely under random blobs of same size.") +
  theme_minimal() +
  theme(legend.position = "bottom", plot.title = element_text(size = 11),
        plot.subtitle = element_text(size = 9))

combined <- p_dens / p_volc + plot_layout(heights = c(1, 1.3))
ggsave(file.path(OUT, "tumor_pv_random_region_null.png"),
       combined, width = 12, height = 13, dpi = 130)
cat("Saved tumor_pv_random_region_null.png\n")
