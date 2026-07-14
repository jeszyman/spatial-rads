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

comp_cols <- c(tumor = "#cc0000", stroma = "#4a86e8", immune = "#e69138")

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
  plot_layout(widths = c(1.5, 1))

save_plot(p, out, w = 13, h = 6)
