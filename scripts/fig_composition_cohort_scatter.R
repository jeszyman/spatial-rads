#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# fig_composition_cohort_scatter.R
# Cross-cohort cell-type census on the locked unified labels: M01 vs M02 pooled
# proportion of each atlas cell type, one point per type, colored by compartment.
# Agreement sits on the diagonal; a type off-diagonal is enriched in one cohort.
# Two panels: (A) all M01 timepoints pooled vs M02; (B) day-2 only, the sole
# timepoint M02 samples (so the M01 axis is its two day-2 slides). Reads the
# per-sample composition table (final cell_subtype labels + dataset + timepoint).
# Args (both optional, canonical defaults): <composition_by_sample.tsv> <out.png>
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({library(tidyverse); library(patchwork)})
source("/home/jeszyman/repos/science/R/figure_schema.R")

a   <- commandArgs(trailingOnly = TRUE)
inp <- if (length(a) >= 1) a[1] else "results/aggregate/composition_by_sample.tsv"
out <- if (length(a) >= 2) sub("\\.png$", "", a[2]) else "results/aggregate/plots/composition_cohort_scatter"

comp_map <- c(
  Tumor = "tumor",
  Fibroblast = "stroma", SmoothMuscle = "stroma", Endothelial = "stroma",
  Adipocyte = "stroma", unassigned = "stroma",
  Macrophages = "immune", `T cells` = "immune", `NK cells` = "immune", ILC = "immune",
  DC = "immune", `Plasma cells` = "immune", `Mast cells` = "immune",
  Neutrophils = "immune", `Epithelial cells` = "immune")
comp_cols <- c(tumor = "#cc0000", stroma = "#4a86e8", immune = "#e69138")

raw <- read_tsv(inp, show_col_types = FALSE)

# Pooled per-cohort proportion (total cells of a type / cohort total) over a
# sample subset, returned wide as m01/m02 with the compartment label.
cohort_frac <- function(df) {
  df %>%
    group_by(dataset, cell_type) %>% summarise(n = sum(n_cells), .groups = "drop") %>%
    group_by(dataset) %>% mutate(frac = n / sum(n)) %>% ungroup() %>%
    select(dataset, cell_type, frac) %>%
    pivot_wider(names_from = dataset, values_from = frac, values_fill = 0) %>%
    rename(m01 = Mutter_01, m02 = Mutter_02) %>%
    mutate(compartment = factor(comp_map[cell_type], levels = c("tumor", "stroma", "immune")))
}

d_all <- cohort_frac(raw)
d_d2  <- cohort_frac(raw %>% filter(timepoint_h == 48))

# Shared log-axis floor and limits across both panels so the diagonal is comparable.
flr  <- min(c(d_all$m01[d_all$m01 > 0], d_all$m02[d_all$m02 > 0],
              d_d2$m01[d_d2$m01 > 0],  d_d2$m02[d_d2$m02 > 0])) / 2
pin  <- function(d) mutate(d, m01p = pmax(m01, flr), m02p = pmax(m02, flr))
d_all <- pin(d_all); d_d2 <- pin(d_d2)
lims  <- range(c(d_all$m01p, d_all$m02p, d_d2$m01p, d_d2$m02p))

for (d in list(d_all, d_d2)) {
  z <- d %>% filter(m01 == 0 | m02 == 0)
  if (nrow(z)) cat("Cohort-absent (pinned to floor):", paste(z$cell_type, collapse = ", "), "\n")
}

panel <- function(d, subtitle, tag) {
  ggplot(d, aes(m01p, m02p, color = compartment)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(aes(label = cell_type), size = 3, show.legend = FALSE,
                             max.overlaps = 20, seed = 1) +
    scale_x_log10(limits = lims) + scale_y_log10(limits = lims) +
    scale_color_manual(values = comp_cols, name = NULL) +
    coord_equal() +
    labs(x = "Mutter_01 proportion", y = "Mutter_02 proportion",
         subtitle = subtitle, tag = tag) +
    theme_scifig(base_size = 12)
}

pA <- panel(d_all, "All timepoints (M01 n=8)", "A")
pB <- panel(d_d2,  "Day 2 only (M01 n=2)",     "B")

p <- (pA | pB) + plot_layout(guides = "collect")

save_plot(p, out, w = 12, h = 6)
