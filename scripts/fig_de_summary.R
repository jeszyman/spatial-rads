#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# fig_de_summary.R
# Arm-level pseudobulk DE summary (Mutter_02 day-2), the modernized replacement
# for the retired deg_pseudobulk figures. MDE-aware and call-class-aware, not a
# volcano (invalid on sparse CosMx). (A) confirmatory-hit count per cell type x
# contrast that clears its n=4 MDE, so underpowered non-results don't read as
# nulls. (B) effect vs MDE for the confirmatory family, colored by call_class
# (composition fraction_shift vs per-cell regulation vs ambiguous). Reads the
# single canonical results_master.tsv (DE rows).
# Args (optional): <results_master.tsv> <out.png>
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({library(tidyverse); library(patchwork)})
source("/home/jeszyman/repos/science/R/figure_schema.R")

a   <- commandArgs(trailingOnly = TRUE)
mas <- if (length(a) >= 1) a[1] else "results/aggregate/results_master.tsv"
out <- if (length(a) >= 2) sub("\\.png$", "", a[2]) else "results/aggregate/plots/de_summary"

clv <- c("MBRT_vs_Ctrl", "SBRT_vs_Ctrl", "MBRT_vs_SBRT")
de  <- read_tsv(mas, show_col_types = FALSE) %>%
  filter(readout_class == "DE") %>%
  mutate(contrast = factor(contrast, levels = clv))

# Panel A: confirmatory hits (padj_confirmatory<0.05 AND clears its MDE) per
# cell type x contrast. Cell types ordered by total hits.
hits <- de %>% filter(tier == "confirmatory") %>%
  group_by(unit, contrast) %>%
  summarise(n_hit = sum(padj_confirmatory < 0.05 & clears_mde, na.rm = TRUE), .groups = "drop")
ord  <- hits %>% group_by(unit) %>% summarise(tot = sum(n_hit)) %>% arrange(tot) %>% pull(unit)
hits <- mutate(hits, unit = factor(unit, levels = ord))

pA <- ggplot(hits, aes(contrast, unit, fill = n_hit)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(n_hit > 0, n_hit, "")), size = 3.2) +
  scale_fill_gradient(low = "grey92", high = "#cc0000", name = "MDE-clearing\nconfirmatory hits") +
  labs(x = NULL, y = NULL, tag = "A") +
  theme_scifig(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

# Panel B: |effect| vs MDE, confirmatory family, colored by call_class. Points
# above the diagonal clear the MDE; color separates composition vs regulation.
cc_cols <- c(fraction_shift = "#4a86e8", regulation = "#cc0000",
             ambiguous = "grey65", detection_only = "#e69138")
pb <- de %>% filter(tier == "confirmatory") %>%
  mutate(abs_eff = abs(effect),
         call_class = factor(call_class, levels = names(cc_cols)))

pB <- ggplot(pb, aes(mde, abs_eff, color = call_class)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
  geom_point(size = 1.6, alpha = 0.8) +
  facet_wrap(~ contrast, nrow = 1) +
  scale_color_manual(values = cc_cols, name = "call class", drop = FALSE) +
  labs(x = "minimum detectable effect (n=4)", y = "|log2FC|", tag = "B") +
  theme_scifig(base_size = 11)

p <- pA / pB + plot_layout(heights = c(1.1, 1))
save_plot(p, out, w = 10, h = 9)
