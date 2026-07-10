#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# fig_design_grid.R
# Study-design grid: samples per treatment arm x timepoint, faceted by dataset x
# model. Empty cells (grey) expose the design gaps; the fill/label is n. Reads the
# master sample sheet only.
# Args: <samples.tsv> <out.png>
# -----------------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
ss  <- args[1]
out <- args[2]

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

p <- ggplot(grid, aes(treatment, fct_rev(tp), fill = n)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = nlab), color = "white", fontface = "bold", size = 4.2) +
  facet_wrap(~ panel, scales = "free_y") +
  scale_fill_viridis_c(option = "D", begin = 0.25, end = 0.9,
                       na.value = "grey92", breaks = 1:4) +
  labs(x = "Treatment arm", y = "Timepoint", fill = "n samples",
       title = "Experimental design: samples per arm x timepoint",
       subtitle = "Mutter_01 flank kinetic series (n=1) + tongue; Mutter_02 day-2 cohort (n=4/arm). Grey = no sample.") +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"))

ggsave(out, p, width = 9, height = 4.5, dpi = 300)
cat("wrote", out, "\n")
