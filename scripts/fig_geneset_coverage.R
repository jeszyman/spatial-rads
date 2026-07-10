#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# fig_geneset_coverage.R
# Pathway-set panel coverage: fraction of each curated confirmatory set's genes
# present on the common panel, flagged when below the usability floor. Reads the
# gene-set panel-coverage table. Interpretation aid: low coverage means a null in
# that set is a panel blind spot, not necessarily biology.
# Args: <gene_set_panel_coverage.tsv> <out.png>
# -----------------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
cov <- args[1]
out <- args[2]

d <- read_tsv(cov, show_col_types = FALSE, comment = "#") %>%   # drop msigdbr_version line
  filter(tier == "primary") %>%
  mutate(coverage = n_panel / n_total) %>%
  arrange(coverage)
d$set <- factor(d$set, levels = d$set)

p <- ggplot(d, aes(coverage, set)) +
  geom_segment(aes(x = 0, xend = coverage, y = set, yend = set), color = "grey70") +
  geom_point(aes(color = thin), size = 4.5) +
  geom_text(aes(label = sprintf("%d/%d", n_panel, n_total)), hjust = -0.35, size = 3.4) +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1.12),
                     breaks = seq(0, 1, 0.25)) +
  scale_color_manual(values = c("FALSE" = "#2a9d8f", "TRUE" = "#e76f51"),
                     labels = c("FALSE" = "usable", "TRUE" = "thin (below floor)"),
                     name = NULL) +
  labs(x = "Fraction of set genes on the common 950-gene panel", y = NULL,
       title = "Pathway-set panel coverage (curated confirmatory sets)",
       subtitle = "Low coverage = a null there is a panel blind spot, not necessarily biology") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

ggsave(out, p, width = 8, height = 4, dpi = 300)
cat("wrote", out, "\n")
