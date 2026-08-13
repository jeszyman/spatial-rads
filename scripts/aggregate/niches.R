#!/usr/bin/env Rscript
# aggregate.smk spatial niches. Each cell is described by the cell_subtype
# composition of its k=20 nearest neighbours (within its sample, the contiguous
# tissue unit); k-means (K=6) on those vectors defines data-driven niches N1..N6.
# Method ported from dev/mbrt_vs_sbrt/07_niche_clustering.R, re-run on the unified
# cross-dataset labels + new coords over all 20 flank samples (prior dev niches are
# re-derived, not imported). M02 day-2 niche frequency gets the same propeller test
# as composition.R: getTransformedProps(logit) -> lmFit(~0+condition+slide_id) ->
# 3 contrasts -> eBayes(robust) -> topTable with global BH.
# Args: <full_labels.parquet> <coords.parquet> <obs.parquet> <out_per_cell.parquet>
#       <out_centroids.tsv> <out_freq.tsv> <out_test.tsv> <plot_heatmap> <plot_freq>
suppressPackageStartupMessages({
  library(data.table); library(arrow); library(RANN)
  library(ggplot2)
})
a <- commandArgs(trailingOnly = TRUE)
labels_path <- a[1]; coords_path <- a[2]; obs_path <- a[3]
out_pc <- a[4]; out_cent <- a[5]; out_freq <- a[6]; out_lm_input <- a[7]
plot_heat <- a[8]; plot_freq <- a[9]
K <- 6L; SEED <- 42L; KNN <- 20L

dir.create(dirname(out_cent), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(plot_heat), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_pc), recursive = TRUE, showWarnings = FALSE)

lab <- as.data.table(read_parquet(labels_path))[, .(cell, cell_subtype)]
co  <- as.data.table(read_parquet(coords_path))[, .(cell, sample_id, x_slide_mm, y_slide_mm)]
ob  <- as.data.table(read_parquet(obs_path))[, .(cell, condition, slide_id, timepoint_h, dataset)]
d <- lab[co, on = "cell"]; d <- ob[d, on = "cell"]
d <- d[!is.na(cell_subtype)]
subt <- sort(unique(d$cell_subtype))

build_comp <- function(ds) {                       # k=20 NN cell_subtype fractions
  xy <- as.matrix(ds[, .(x_slide_mm, y_slide_mm)])
  k_use <- min(KNN + 1L, nrow(ds))
  nn <- RANN::nn2(xy, k = k_use)$nn.idx[, -1, drop = FALSE]
  codes <- match(ds$cell_subtype, subt)
  neigh <- matrix(codes[nn], nrow = nrow(nn))
  comp <- vapply(seq_along(subt),
                 function(ci) rowMeans(neigh == ci, na.rm = TRUE),
                 numeric(nrow(neigh)))
  colnames(comp) <- subt; comp
}
parts <- lapply(split(seq_len(nrow(d)), d$sample_id), function(idx) {
  list(cells = d$cell[idx], comp = build_comp(d[idx]))
})
all_comp  <- do.call(rbind, lapply(parts, `[[`, "comp"))
all_cells <- unlist(lapply(parts, `[[`, "cells"), use.names = FALSE)

set.seed(SEED)
km <- kmeans(all_comp, centers = K, nstart = 10, iter.max = 50)
d <- data.table(cell = all_cells, niche = paste0("N", km$cluster))[d, on = "cell"]

write_parquet(d[, .(cell, sample_id, niche, cell_subtype, x_slide_mm, y_slide_mm)], out_pc)

cent <- as.data.table(km$centers)[, niche := paste0("N", .I)]
cent_long <- melt(cent, id.vars = "niche", variable.name = "cell_subtype",
                  value.name = "mean_frac")
fwrite(cent_long, out_cent, sep = "\t")

freq <- d[, .(n = .N), by = .(sample_id, niche, condition, timepoint_h, dataset, slide_id)]
freq[, frac := n / sum(n), by = sample_id]
setorder(freq, dataset, sample_id, niche)
fwrite(freq, out_freq, sep = "\t")

# --- engine input: per-cell niche labels (cell, sample_id, label, condition, slide_id), all
# flank samples. The lm_engine applies cohort_samples.tsv whitelist filtering. ---
lm_input <- d[!grepl("^Tongue", condition), .(cell, sample_id, label = niche, condition, slide_id)]
fwrite(lm_input, out_lm_input, sep = "\t")

# --- plots ----------------------------------------------------------------------
p_heat <- ggplot(cent_long, aes(cell_subtype, niche, fill = mean_frac)) +
  geom_tile() +
  scale_fill_viridis_c(name = "mean\nNN frac") +
  labs(x = NULL, y = NULL, title = "Niche centroid composition (k=20 NN cell_subtype fractions)") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot_heat, p_heat, width = 9, height = 4.5, dpi = 150)

m2f <- freq[dataset == "Mutter_02"]
p_freq <- ggplot(m2f, aes(niche, frac, fill = condition)) +
  geom_boxplot(outlier.size = 0.5, position = position_dodge(width = 0.8)) +
  geom_point(size = 0.7, position = position_dodge(width = 0.8)) +
  scale_fill_brewer(palette = "Set1") +
  labs(x = NULL, y = "niche fraction of sample",
       title = "M02 niche frequency by arm") +
  theme_bw(base_size = 10)
ggsave(plot_freq, p_freq, width = 8, height = 4.5, dpi = 150)

cat(sprintf("niches: %d cells, K=%d | cluster sizes %s | lm input %d flank cells\n",
            nrow(d), K, paste(km$size, collapse = "/"), nrow(lm_input)))
