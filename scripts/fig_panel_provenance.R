#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# fig_panel_provenance.R
# Panel provenance: per-dataset decomposition of the delivered gene panel into
# shared-stock / dataset-only-stock / custom, against the UCC stock reference.
# Reads the data-derived panel_provenance membership table.
# Args: <panel_provenance.tsv> <out.png>
# -----------------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
pv  <- args[1]
out <- args[2]

d <- read_tsv(pv, show_col_types = FALSE)

shared     <- sum(d$mutter_01 & d$mutter_02)          # 950 shared stock
m01_only   <- sum(d$mutter_01 & !d$mutter_02)          # stock genes only in M01
m02_custom <- sum(d$mutter_02 & d$custom)              # M02 non-stock custom
ucc        <- sum(d$ucc_standard)                      # UCC stock size

df <- tribble(
  ~dataset,    ~category,          ~n,
  "Mutter_01", "Shared stock",     shared,
  "Mutter_01", "M01-only stock",   m01_only,
  "Mutter_02", "Shared stock",     shared,
  "Mutter_02", "M02 custom",       m02_custom) %>%
  mutate(category = factor(category,
                           levels = c("Shared stock", "M01-only stock", "M02 custom")))

p <- ggplot(df, aes(n, fct_rev(dataset), fill = category)) +
  geom_col(width = 0.6, color = "white") +
  geom_vline(xintercept = ucc, linetype = "dashed", color = "grey30") +
  annotate("text", x = ucc, y = 1.5, label = paste0("UCC stock = ", ucc),
           hjust = 1.1, size = 3.4, color = "grey30") +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5),
            color = "white", fontface = "bold", size = 4) +
  scale_fill_manual(values = c("Shared stock" = "#3c78d8",
                               "M01-only stock" = "#8e8e8e",
                               "M02 custom" = "#e06666")) +
  labs(x = "Genes", y = NULL, fill = NULL,
       title = "Panel provenance: stock vs custom composition",
       subtitle = "M01 = full 1000-gene UCC stock; M02 = 950 shared stock + 21 custom (~50 stock dropped)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.major.y = element_blank(), legend.position = "bottom")

ggsave(out, p, width = 8, height = 3.3, dpi = 300)
cat("wrote", out, "\n")
