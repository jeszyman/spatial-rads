#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# fig_composition_arm_scatter.R
# Arm-vs-arm compositional scatter (Mutter_02 day-2): per-cell-type pooled
# proportion, treated arm (y) vs Control (x), one panel per irradiated arm.
# Points on the diagonal = no compositional shift; off-diagonal = a cell type
# enriched/depleted under that arm. The compositional analog of the cross-cohort
# census. Reads the per-sample composition table (locked labels).
# Args (optional): <composition_by_sample.tsv> <out.png>
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({library(tidyverse); library(patchwork)})
source("/home/jeszyman/repos/science/R/figure_schema.R")

a    <- commandArgs(trailingOnly = TRUE)
comp <- if (length(a) >= 1) a[1] else "results/aggregate/composition_by_sample.tsv"
out  <- if (length(a) >= 2) sub("\\.png$", "", a[2]) else "results/aggregate/plots/composition_arm_scatter"

comp_map <- c(
  Tumor = "tumor",
  Fibroblast = "stroma", SmoothMuscle = "stroma", Endothelial = "stroma",
  Adipocyte = "stroma", unassigned = "stroma",
  Macrophages = "immune", `T cells` = "immune", `NK cells` = "immune", ILC = "immune",
  DC = "immune", `Plasma cells` = "immune", `Mast cells` = "immune",
  Neutrophils = "immune")
comp_cols <- c(tumor = "#cc0000", stroma = "#4a86e8", immune = "#e69138")

# Pooled per-arm proportion (cells of a type / arm total).
frac <- read_tsv(comp, show_col_types = FALSE) %>%
  filter(dataset == "Mutter_02") %>%
  mutate(arm = recode(condition, MBRT_day2 = "MBRT", SBRT_day2 = "SBRT")) %>%
  group_by(arm, cell_type) %>% summarise(n = sum(n_cells), .groups = "drop") %>%
  group_by(arm) %>% mutate(frac = n / sum(n)) %>% ungroup() %>%
  select(arm, cell_type, frac) %>%
  pivot_wider(names_from = arm, values_from = frac, values_fill = 0) %>%
  mutate(compartment = factor(comp_map[cell_type], levels = c("tumor", "stroma", "immune")))

flr  <- min(frac %>% select(Control, MBRT, SBRT) %>% as.matrix() %>% .[. > 0]) / 2
lims <- range(frac %>% select(Control, MBRT, SBRT) %>% as.matrix() %>% pmax(flr))

panel <- function(treated, tag) {
  d <- frac %>% transmute(cell_type, compartment,
                          x = pmax(Control, flr), y = pmax(.data[[treated]], flr))
  ggplot(d, aes(x, y, color = compartment)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(aes(label = cell_type), size = 3, show.legend = FALSE,
                             max.overlaps = 20, seed = 1) +
    scale_x_log10(limits = lims) + scale_y_log10(limits = lims) +
    scale_color_manual(values = comp_cols, name = NULL) +
    coord_equal() +
    labs(x = "Control proportion", y = sprintf("%s proportion", treated), tag = tag) +
    theme_scifig(base_size = 12)
}

p <- (panel("MBRT", "A") | panel("SBRT", "B")) + plot_layout(guides = "collect")
save_plot(p, out, w = 12, h = 6)
