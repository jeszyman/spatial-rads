#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# fig_marker_coverage.R
# Cell-lineage / substate marker panel coverage: fraction of each set's markers
# present on the common panel, coloured by tier, thin sets flagged. Reads the
# marker panel-coverage table. A lineage with off-panel markers is silently
# weakened; this figure makes that visible.
# Args: <marker_panel_coverage.tsv> <out.png>
# -----------------------------------------------------------------------------
suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
cov <- args[1]; out <- args[2]

d <- read_tsv(cov, show_col_types = FALSE) %>%
  arrange(tier, frac_panel)
d$set <- factor(d$set, levels = d$set)

p <- ggplot(d, aes(frac_panel, set)) +
  geom_segment(aes(x = 0, xend = frac_panel, y = set, yend = set), color = "grey70") +
  geom_point(aes(color = tier, shape = thin), size = 4) +
  geom_text(aes(label = sprintf("%d/%d", n_panel, n_total)), hjust = -0.35, size = 3.2) +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1.12),
                     breaks = seq(0, 1, 0.25)) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 4),
                     labels = c("FALSE" = "usable", "TRUE" = "thin"), name = NULL) +
  labs(x = "Fraction of set markers on the common 950-gene panel", y = NULL,
       title = "Cell-marker panel coverage",
       subtitle = "Off-panel markers silently weaken a lineage; x = thin set") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

ggsave(out, p, width = 8, height = 5, dpi = 300)
cat("wrote", out, "\n")
