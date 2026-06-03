#!/usr/bin/env Rscript
# aggregate.smk tumor-immune spatial mixing. Per sample (contiguous tissue unit),
# k=20 NN; per cell count immune- and tumor-compartment neighbours. Two readouts:
# (1) immune-neighbour fraction per cell (dev 05_spatial_nn.R), and (2) the Keren
# et al. 2018 mixing score per sample = tumor->immune neighbour edges / immune->
# immune edges (dev 12_mixing_score.R). Necrosis exclusion: cells in necrosis_zone
# are dropped before aggregating at the day-2/day-6 timepoints (the dev convention),
# the NN graph itself is built on all cells. Compartment (tumor/stroma/immune) comes
# from the unified labels, so prior dev findings are re-tested, not imported. M02
# day-2 gets a limma test (~0+condition+slide_id, 3 contrasts, 95% CIs, BH) on the
# per-sample metric matrix; with few outcome rows eBayes moderation is negligible.
# Args: <full_labels.parquet> <coords.parquet> <obs.parquet> <out_per_sample.tsv>
#       <out_test.tsv> <out_per_cell.parquet> <plot_mixing>
suppressPackageStartupMessages({
  library(data.table); library(arrow); library(RANN); library(limma); library(ggplot2)
})
a <- commandArgs(trailingOnly = TRUE)
labels_path <- a[1]; coords_path <- a[2]; obs_path <- a[3]
out_ps <- a[4]; out_test <- a[5]; out_pc <- a[6]; plot_mix <- a[7]
KNN <- 20L

dir.create(dirname(out_ps), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(plot_mix), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_pc), recursive = TRUE, showWarnings = FALSE)

lab <- as.data.table(read_parquet(labels_path))[, .(cell, compartment)]
co  <- as.data.table(read_parquet(coords_path))[, .(cell, sample_id, x_slide_mm, y_slide_mm, necrosis_zone)]
ob  <- as.data.table(read_parquet(obs_path))[, .(cell, condition, slide_id, timepoint_h, dataset)]
d <- lab[co, on = "cell"]; d <- ob[d, on = "cell"]
d <- d[!is.na(compartment)]

# per-cell immune/tumor neighbour counts via per-sample k=20 NN (graph on all cells)
nn_counts <- function(ds) {
  xy <- as.matrix(ds[, .(x_slide_mm, y_slide_mm)])
  k_use <- min(KNN + 1L, nrow(ds))
  nn <- RANN::nn2(xy, k = k_use)$nn.idx[, -1, drop = FALSE]
  comp_neigh <- matrix(ds$compartment[nn], nrow = nrow(nn))
  data.table(cell = ds$cell,
             immune_neighbor_n = rowSums(comp_neigh == "immune"),
             tumor_neighbor_n  = rowSums(comp_neigh == "tumor"),
             total_neighbors   = ncol(comp_neigh))
}
cnt <- rbindlist(lapply(split(seq_len(nrow(d)), d$sample_id),
                        function(idx) nn_counts(d[idx])))
d <- cnt[d, on = "cell"]
d[, immune_frac := immune_neighbor_n / total_neighbors]

write_parquet(d[, .(cell, sample_id, compartment, immune_neighbor_n,
                    tumor_neighbor_n, total_neighbors, immune_frac, necrosis_zone)], out_pc)

# necrosis exclusion at day2/day6, then per-sample aggregation
dm <- d[!(timepoint_h %in% c(48L, 144L) & necrosis_zone == TRUE)]
per_sample <- dm[, {
  n_tot <- .N; n_tum <- sum(compartment == "tumor"); n_imm <- sum(compartment == "immune")
  ti <- sum(immune_neighbor_n[compartment == "tumor"])
  ii <- sum(immune_neighbor_n[compartment == "immune"])
  gif <- n_imm / n_tot
  .(n_cells = n_tot, n_tumor = n_tum, n_immune = n_imm,
    tumor_immune_edges = ti, immune_immune_edges = ii,
    mixing_score = ti / max(ii, 1),
    mean_immune_frac = mean(immune_frac),
    immune_per_tumor = ti / max(n_tum, 1),
    enrichment_over_random = (ti / max(n_tum, 1)) / max(KNN * gif, 1e-3))
}, by = .(sample_id, condition, slide_id, timepoint_h, dataset)]
setorder(per_sample, dataset, sample_id)
fwrite(per_sample, out_ps, sep = "\t")

# --- M02 day2 limma test on the per-sample metric matrix ------------------------
m2 <- per_sample[dataset == "Mutter_02" & timepoint_h == 48L]
m2[, condition := factor(condition, levels = c("Control", "MBRT_day2", "SBRT_day2"))]
setkey(m2, sample_id)
metrics <- c("mixing_score", "mean_immune_frac", "immune_per_tumor", "enrichment_over_random")
mat <- t(as.matrix(m2[, ..metrics])); colnames(mat) <- m2$sample_id
design <- model.matrix(~ 0 + condition + slide_id, data = m2)
colnames(design) <- make.names(colnames(design))
fit <- lmFit(mat, design)
cm  <- makeContrasts(
  MBRT_vs_Ctrl = conditionMBRT_day2 - conditionControl,
  SBRT_vs_Ctrl = conditionSBRT_day2 - conditionControl,
  MBRT_vs_SBRT = conditionMBRT_day2 - conditionSBRT_day2, levels = design)
fit2 <- eBayes(contrasts.fit(fit, cm))
test <- rbindlist(lapply(colnames(cm), function(cn) {
  tt <- topTable(fit2, coef = cn, number = Inf, sort.by = "none", confint = TRUE)
  data.table(metric = rownames(tt), contrast = cn, estimate = tt$logFC,
             ci_low = tt$CI.L, ci_high = tt$CI.R, t_stat = tt$t, pvalue = tt$P.Value)
}))
test[, padj := p.adjust(pvalue, method = "BH")]
test[, `:=`(method = "limma", n_samples_per_group = 4L, dataset = "Mutter_02")]
fwrite(test, out_test, sep = "\t")

# --- plot: mixing score + immune-neighbour fraction by arm ----------------------
pd <- melt(m2[, c("sample_id", "condition", "mixing_score", "mean_immune_frac")],
           id.vars = c("sample_id", "condition"), variable.name = "metric", value.name = "value")
p <- ggplot(pd, aes(condition, value, fill = condition)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_point(size = 1, position = position_jitter(width = 0.12)) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_fill_brewer(palette = "Set1") +
  labs(x = NULL, y = NULL, title = "M02 day2 tumor-immune mixing by arm (n=4/arm)") +
  theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(plot_mix, p, width = 8, height = 4.5, dpi = 150)

cat(sprintf("mixing: %d samples | M02 day2 test %d rows, %d padj<0.05 | mixing_score range %.3f-%.3f\n",
            nrow(per_sample), nrow(test), test[padj < 0.05, .N],
            min(m2$mixing_score), max(m2$mixing_score)))
