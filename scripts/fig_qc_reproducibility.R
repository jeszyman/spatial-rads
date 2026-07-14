#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# fig_qc_reproducibility.R
# Replicate reproducibility (M02 day-2, n=4/arm): (left) PCA of the scaled per-sample
# technical metrics, (right) each sample's mean pseudobulk Spearman concordance to its
# arm-mates. Outlier check on the composition inference -- a sample sitting apart in
# PCA or poorly concordant with its arm-mates could drive an arm effect. Reads the
# per-sample table from qc_reproducibility.R.
# Args (all optional, canonical defaults): <qc_reproducibility.tsv> <out.png>
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({library(tidyverse); library(patchwork)})
source("/home/jeszyman/repos/science/R/figure_schema.R")

a   <- commandArgs(trailingOnly = TRUE)
inp <- if (length(a) >= 1) a[1] else "results/aggregate/qc_reproducibility.tsv"
out <- if (length(a) >= 2) sub("\\.png$", "", a[2]) else "results/aggregate/plots/qc_reproducibility"

d <- read_tsv(inp, show_col_types = FALSE) %>%
  mutate(arm = factor(arm, levels = c("Control", "MBRT", "SBRT")))

# PC variance-explained travels in the table; fall back to bare labels if absent.
pc1_lab <- if ("pc1_ve" %in% names(d)) sprintf("PC1 (%.0f%%)", d$pc1_ve[1]) else "PC1"
pc2_lab <- if ("pc2_ve" %in% names(d)) sprintf("PC2 (%.0f%%)", d$pc2_ve[1]) else "PC2"

arm_cols <- c(Control = "#4a86e8", MBRT = "#e69138", SBRT = "#cc0000")

p1 <- ggplot(d, aes(PC1, PC2, color = arm)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(aes(label = slide_id), size = 3, show.legend = FALSE, max.overlaps = 20) +
  scale_color_manual(values = arm_cols) +
  labs(x = pc1_lab, y = pc2_lab, subtitle = "Technical-metric PCA", color = NULL) +
  theme_scifig(base_size = 12)

p2 <- ggplot(d, aes(arm, mean_r_armmates, color = arm)) +
  geom_jitter(width = 0.12, height = 0, size = 3) +
  scale_color_manual(values = arm_cols, guide = "none") +
  labs(x = NULL, y = "mean Spearman r to arm-mates",
       subtitle = "Pseudobulk replicate concordance") +
  theme_scifig(base_size = 12)

p <- (p1 | p2)

save_plot(p, out, w = 11, h = 5.2)
