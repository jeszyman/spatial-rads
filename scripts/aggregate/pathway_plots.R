#!/usr/bin/env Rscript
# aggregate.smk pathway track -- PLOT half. Reads the cached score/test TSVs from
# pathway_scores.R and renders the three figures. Split from the ~5h compute so a
# plot bug never forces a recompute. Args: <summary.tsv> <test.tsv>
#   <plot_heatmap> <plot_timecourse> <plot_scatter>
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
a <- commandArgs(trailingOnly = TRUE)
summary_out <- fread(a[1]); test_out <- fread(a[2])
plot_heat <- a[3]; plot_tc <- a[4]; plot_scatter <- a[5]
dir.create(dirname(plot_heat), recursive = TRUE, showWarnings = FALSE)
cm_names <- c("MBRT_vs_Ctrl", "SBRT_vs_Ctrl", "MBRT_vs_SBRT")

# --- plot 1: M02 primary-pathway test heatmap (UCell), * = padj<0.05 ------------
ph <- test_out[score_type == "UCell" & tier == "primary"]
ph[, sig := !is.na(padj_bh) & padj_bh < 0.05]
ph[, contrast := factor(contrast, levels = cm_names)]
ph[, pathway_name := factor(pathway_name, levels = rev(sort(unique(pathway_name))))]
ct_ord <- ph[contrast == "SBRT_vs_Ctrl", .(m = mean(estimate, na.rm = TRUE)),
             by = cell_type][order(m), cell_type]
ph[, cell_type := factor(cell_type, levels = ct_ord)]
p1 <- ggplot(ph, aes(cell_type, pathway_name, fill = estimate)) +
  geom_tile(colour = "grey92") +
  geom_text(data = ph[sig == TRUE], aes(label = "*"), size = 3, vjust = 0.75) +
  facet_wrap(~ contrast, ncol = 1) +
  scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick",
                       midpoint = 0, name = "mean-score\nestimate") +
  labs(x = NULL, y = NULL,
       title = "M02 day2 primary-pathway scores (UCell, limma, global BH; * padj<0.05)") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot_heat, p1, width = 11, height = 6.5, dpi = 150)

# --- plot 2: M01 primary-pathway timecourse (UCell, cell-count-weighted) --------
m1  <- summary_out[dataset == "Mutter_01" & tier == "primary" & score_type == "UCell"]
m1w <- m1[, .(wmean = sum(mean * n_cells) / sum(n_cells)),
          by = .(sample_id, pathway_name, timepoint_h, treatment)]
p2 <- ggplot(m1w, aes(timepoint_h, wmean, colour = treatment)) +
  geom_line(aes(group = treatment), linewidth = 0.4) +
  geom_point(size = 1.1) +
  facet_wrap(~ pathway_name, scales = "free_y") +
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "timepoint (h)", y = "mean UCell score (cell-count weighted)",
       title = "M01 primary-pathway score over time (n=1, descriptive)") +
  theme_bw(base_size = 9)
ggsave(plot_tc, p2, width = 8, height = 6, dpi = 150)

# --- plot 3: UCell vs AMS per-sample-mean concordance ---------------------------
wide <- dcast(summary_out, dataset + cell_type + pathway_name + sample_id ~ score_type,
              value.var = "mean")
tier_lkp <- unique(summary_out[, .(pathway_name, tier)])
wide[, tier := tier_lkp$tier[match(pathway_name, tier_lkp$pathway_name)]]
r_all <- wide[is.finite(UCell) & is.finite(AMS), round(cor(UCell, AMS), 3)]
p3 <- ggplot(wide, aes(UCell, AMS, colour = tier)) +
  geom_point(size = 0.4, alpha = 0.15) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  scale_colour_manual(values = c(primary = "firebrick", exploratory = "grey50")) +
  labs(x = "per-sample mean UCell", y = "per-sample mean AddModuleScore",
       title = sprintf("UCell vs AddModuleScore concordance (overall Pearson r = %.3f)", r_all)) +
  theme_bw(base_size = 9)
ggsave(plot_scatter, p3, width = 7, height = 6, dpi = 150)
cat(sprintf("pathway_plots: heatmap %d rows | timecourse %d rows | scatter %d points\n",
            nrow(ph), nrow(m1w), nrow(wide)))
