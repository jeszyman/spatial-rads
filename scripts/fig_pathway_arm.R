#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# fig_pathway_arm.R
# Arm-level pathway response (Mutter_02 day-2): UCell-score arm effect for the
# primary curated gene sets, cell type x pathway, one facet per contrast. Tile
# fill = effect (limma estimate on per-sample UCell means); star = padj_bh<0.05.
# The pathway view of the arm contrasts, on the locked labels. Reads the
# pathway test table (primary tier, UCell score type).
# Args (optional): <pathway_test_m02day2.tsv> <out.png>
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({library(tidyverse)})
source("/home/jeszyman/repos/science/R/figure_schema.R")

a   <- commandArgs(trailingOnly = TRUE)
tst <- if (length(a) >= 1) a[1] else "results/aggregate/pathway_test_m02day2.tsv"
out <- if (length(a) >= 2) sub("\\.png$", "", a[2]) else "results/aggregate/plots/pathway_arm"

clv <- c("MBRT_vs_Ctrl", "SBRT_vs_Ctrl", "MBRT_vs_SBRT")
# STING is demoted (2/3 on-panel genes shared with the Type-I/II IFN sets, so
# not an independent readout); the table still carries it, the figure excludes it.
d <- read_tsv(tst, show_col_types = FALSE) %>%
  filter(tier == "primary", score_type == "UCell", pathway_name != "STING") %>%
  mutate(contrast = factor(contrast, levels = clv),
         pathway_name = str_remove(pathway_name, "^HALLMARK_"),
         pathway_lab = str_wrap(str_replace_all(pathway_name, "_", " "), 12),
         sig = padj_bh < 0.05)

# Cap fill at a symmetric quantile so one outlier tile doesn't wash the scale.
lim <- quantile(abs(d$estimate), 0.98, na.rm = TRUE)

p <- ggplot(d, aes(contrast, cell_type, fill = estimate)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(sig, "*", "")), size = 4, vjust = 0.78) +
  facet_wrap(~ pathway_lab, nrow = 1) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
                       limits = c(-lim, lim), oob = scales::squish, name = "UCell\narm effect") +
  labs(x = NULL, y = NULL) +
  theme_scifig(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        strip.text = element_text(size = 8))

save_plot(p, out, w = 11, h = 5.5)
