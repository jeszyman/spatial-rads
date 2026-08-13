#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# fig_panel_provenance.R
# Panel provenance: per-dataset decomposition of the delivered gene panel into
# shared-stock / dataset-only-stock / custom, against the UCC stock reference.
# Reads the data-derived panel_provenance membership table.
# Args (all optional, canonical defaults): <panel_provenance.tsv> <out.png>
# -----------------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))
source("/home/jeszyman/repos/science/R/figure_schema.R")

a   <- commandArgs(trailingOnly = TRUE)
pv  <- if (length(a) >= 1) a[1] else "results/data_model/panel_provenance.tsv"
out <- if (length(a) >= 2) sub("\\.png$", "", a[2]) else "results/data_model/plots/panel_provenance"

d <- read_tsv(pv, show_col_types = FALSE)

shared     <- sum(d$mutter_01 & d$mutter_02)          # 950 shared stock
m01_only   <- sum(d$mutter_01 & !d$mutter_02)          # stock genes only in M01
m02_custom <- sum(d$mutter_02 & d$custom)              # M02 non-stock custom
ucc        <- sum(d$ucc_standard)                      # UCC stock size

df <- tribble(
  ~dataset,    ~category,          ~n,
  "Mutter_01", "Shared Stock",     shared,
  "Mutter_01", "M01-Only Stock",   m01_only,
  "Mutter_02", "Shared Stock",     shared,
  "Mutter_02", "M02 Custom",       m02_custom) %>%
  mutate(category = factor(category,
                           levels = c("Shared Stock", "M01-Only Stock", "M02 Custom")))

legend_text <- str_c(
  "The two cohorts share a 950-gene common panel, so cross-dataset differential analysis ",
  "rests on a matched feature space. Mutter_01 carries the full 1000-gene UCC stock panel; ",
  "Mutter_02 drops ~50 stock genes and adds 21 custom genes (Fields add-on). Bars decompose ",
  "each cohort's delivered panel into shared stock, cohort-only stock, and custom; the dashed ",
  "line marks the 1000-gene UCC stock reference.")

# shared-stock (950) label centred in its wide segment; the thin M01-only / M02-custom
# slivers cannot hold a centred label, so those counts read just past the sliver edge.
shared_lab <- df %>% filter(category == "Shared Stock") %>%
  mutate(lab_x = ifelse(dataset == "Mutter_01", m01_only, m02_custom) + shared / 2)
thin_lab <- df %>% filter(category != "Shared Stock") %>% mutate(lab_x = n + 10)

p <- ggplot(df, aes(n, fct_rev(dataset), fill = category)) +
  geom_col(width = 0.6, color = "white") +
  geom_vline(xintercept = ucc, linetype = "dashed", color = "grey30") +
  geom_text(data = shared_lab, aes(x = lab_x, label = n),
            color = "white", fontface = "bold", size = 4) +
  geom_text(data = thin_lab, aes(x = lab_x, label = n), hjust = 0,
            color = "white", fontface = "bold", size = 4) +
  scale_x_continuous(breaks = seq(0, 1000, 250),
                     sec.axis = dup_axis(breaks = ucc, labels = "UCC stock", name = NULL)) +
  scale_fill_manual(values = c("Shared Stock" = "#3c78d8",
                               "M01-Only Stock" = "#8e8e8e",
                               "M02 Custom" = "#e06666")) +
  labs(x = "Genes (n)", y = NULL, fill = NULL,
       caption = str_wrap(legend_text, width = round(8 * 15))) +
  theme_scifig(base_size = 12) +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.4, "cm"),
        axis.title.x = element_text(margin = margin(t = 8)),
        plot.caption = element_text(size = 10, hjust = 0, color = "gray30",
                                    lineheight = 1.3, margin = margin(t = 20)),
        plot.caption.position = "plot")

save_plot(p, out, w = 8, h = 4.2)
