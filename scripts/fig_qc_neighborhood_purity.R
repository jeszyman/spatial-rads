#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# fig_qc_neighborhood_purity.R
# Transcriptional neighborhood purity by final cell subtype (Plummer et al., Nat
# Biotechnol 2025, Fig 4d): a typing-QC view of whether the locked cross-dataset
# labels are transcriptionally coherent. Per subtype, the distribution of per-cell
# purity (fraction of k=30 scVI-latent neighbors sharing the label), one box per
# subtype, faceted by compartment, ordered by median. Reads the per-subtype box
# quantiles from celltype_neighborhood_purity.py.
# Panel A: per-subtype purity distribution. Panel B: median purity vs subtype abundance,
# showing purity tracks positively with cell count (rarer subtypes form looser neighborhoods),
# a trend with clear exceptions (SmoothMuscle, NK) so it reads as association, not a rule.
# Args (all optional, canonical defaults): <celltype_neighborhood_purity.tsv> <out.png>
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({library(tidyverse); library(patchwork)})
source("/home/jeszyman/repos/science/R/figure_schema.R")

a   <- commandArgs(trailingOnly = TRUE)
inp <- if (length(a) >= 1) a[1] else "results/aggregate/celltype_neighborhood_purity.tsv"
out <- if (length(a) >= 2) sub("\\.png$", "", a[2]) else "results/aggregate/plots/qc_neighborhood_purity"

# facet_grid(space = "free_x") shows only each compartment's own subtypes per panel, so a single
# factor ordered by (compartment, descending median) renders correctly ordered within every facet
# without a tidytext reorder_within dependency.
d <- read_tsv(inp, show_col_types = FALSE) %>%
  mutate(compartment = factor(compartment, levels = c("tumor", "stroma", "immune"))) %>%
  arrange(compartment, desc(middle)) %>%
  mutate(cell_subtype = fct_inorder(cell_subtype))

rho_k <- signif(d$k_robust_spearman[1], 3)

comp_cols <- c(tumor = "#cc0000", stroma = "#4a86e8", immune = "#e69138")

# purity-vs-abundance association across all subtypes (panel B annotation)
rho_n <- suppressWarnings(cor(log10(d$n), d$middle, method = "spearman"))

legend_text <- str_c(
  "Transcriptional neighborhood purity (Plummer et al., Nat Biotechnol 2025): for each cell, ",
  "the fraction of its 30 nearest neighbors in the scVI integrated latent that share its final ",
  "cell-subtype label; a coherence check on the locked cross-dataset typing, complementary to ",
  "marker recall. (A) Per-subtype purity, one box per subtype (5/25/50/75/95th percentile over ",
  "its cells), the same latent the labels were built in; coarse-compartment purity is near-ceiling ",
  "by construction. (B) Purity tracks positively with subtype abundance (Spearman ",
  sprintf("rho = %.2f", rho_n),
  " over per-subtype medians): rarer subtypes, several of them marker-rescued, form looser ",
  "neighborhoods. The association is not a rule, with smooth muscle and NK cells sitting purer ",
  "than more abundant subtypes. Neighborhood size is not load-bearing: per-subtype medians ",
  sprintf("rank-correlate rho = %s between k = 15 and k = 30.", rho_k))

pA <- ggplot(d, aes(cell_subtype, middle, fill = compartment)) +
  geom_boxplot(aes(ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax),
               stat = "identity", width = 0.7, linewidth = 0.3, alpha = 0.9,
               outlier.shape = NA) +
  facet_grid(~ compartment, scales = "free_x", space = "free_x") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_fill_manual(values = comp_cols, guide = "none") +
  labs(x = NULL, y = "neighborhood purity (k = 30)", tag = "A") +
  theme_scifig(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"))

pB <- ggplot(d, aes(n, middle, color = compartment)) +
  geom_point(size = 2.6, alpha = 0.9) +
  ggrepel::geom_text_repel(aes(label = cell_subtype), size = 3, show.legend = FALSE,
                           max.overlaps = 20, seed = 1) +
  scale_x_log10() +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_color_manual(values = comp_cols, name = NULL) +
  labs(x = "cells assigned (log scale)", y = "median neighborhood purity", tag = "B") +
  theme_scifig(base_size = 12)

p <- (pA | pB) +
  plot_layout(widths = c(1.5, 1)) +
  plot_annotation(caption = str_wrap(legend_text, width = round(11 * 15)),
                  theme = theme(plot.caption = element_text(size = 10, hjust = 0, color = "gray30",
                                                            lineheight = 1.3, margin = margin(t = 18)),
                                plot.caption.position = "plot"))

save_plot(p, out, w = 13, h = 6)
