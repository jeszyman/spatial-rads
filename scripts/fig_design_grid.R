#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# fig_design_grid.R
# Study-design grid: samples per treatment arm x timepoint, faceted by dataset x
# model. Empty cells (grey) expose the design gaps; the fill/label is n. Reads the
# master sample sheet only.
# Args (all optional, canonical defaults): <samples.tsv> <out.png>
# -----------------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))
source("/home/jeszyman/repos/science/R/figure_schema.R")

a   <- commandArgs(trailingOnly = TRUE)
ss  <- if (length(a) >= 1) a[1] else "results/data_model/samples.tsv"
out <- if (length(a) >= 2) sub("\\.png$", "", a[2]) else "results/data_model/plots/design_grid"

tp_lab <- c("0" = "0h", "1" = "1h", "4" = "4h", "48" = "2d",
            "144" = "6d", "192" = "8d", "240" = "10d")
tp_lvl <- c("0h", "1h", "4h", "2d", "6d", "8d", "10d")

d <- read_tsv(ss, show_col_types = FALSE) %>%
  mutate(
    treatment = recode(treatment, NT = "Control"),
    treatment = factor(treatment, levels = c("Control", "MBRT", "SBRT")),
    tp        = factor(tp_lab[as.character(timepoint_h)], levels = tp_lvl),
    panel     = paste0(name, " (", model, ")"))

grid <- d %>%
  count(panel, tp, treatment, name = "n") %>%
  group_by(panel) %>%
  complete(tp, treatment) %>%    # keep observed timepoints x all arms -> gaps show as NA
  ungroup() %>%
  mutate(nlab = ifelse(is.na(n), "", as.character(n)))

legend_text <- str_c(
  "Study design: Mutter_02 is a 2x2 (4h + day-2, n = 2 per arm at each timepoint); ",
  "Mutter_01 is a single-tumour kinetic series (n = 1 per arm x timepoint). ",
  "Cell fill and label give the number of samples; grey cells are ",
  "design gaps (no sample). Faceted by dataset and tumour model; arms are ",
  "Control (0 Gy), MBRT, and SBRT.")

p <- ggplot(grid, aes(treatment, fct_rev(tp), fill = n)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = nlab), color = "white", fontface = "bold", size = 4.2) +
  facet_wrap(~ panel, scales = "free_y") +
  scale_fill_viridis_c(option = "D", begin = 0.25, end = 0.9,
                       na.value = "grey92", breaks = 1:4) +
  labs(x = "Treatment Arm", y = "Timepoint", fill = "Samples (n)",
       caption = str_wrap(legend_text, width = round(9 * 15))) +
  theme_scifig(base_size = 12) +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "bottom",
        legend.key.height = unit(0.4, "cm"),
        axis.title.x = element_text(margin = margin(t = 8)),
        axis.title.y = element_text(margin = margin(r = 8)),
        plot.caption = element_text(size = 10, hjust = 0, color = "gray30",
                                    lineheight = 1.3, margin = margin(t = 20)),
        plot.caption.position = "plot")

save_plot(p, out, w = 9, h = 5.5)
