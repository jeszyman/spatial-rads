#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# fig_composition_arm.R
# Arm-level compositional descriptive (Mutter_02 day-2, n=4/arm): (A) stacked
# proportion of each atlas cell type per treatment arm; (B) per-cell-type
# propeller effect (logit log2FC) with 95% CI, faceted by contrast. The
# descriptive backdrop to the composition inference. Reads the per-sample
# composition table + the propeller test table (both on the locked labels).
# Args (optional): <composition_by_sample.tsv> <composition_test.tsv> <out.png>
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({library(tidyverse); library(patchwork)})
source("/home/jeszyman/repos/science/R/figure_schema.R")

a    <- commandArgs(trailingOnly = TRUE)
comp <- if (length(a) >= 1) a[1] else "results/aggregate/composition_by_sample.tsv"
test <- if (length(a) >= 2) a[2] else "results/aggregate/composition_test_m02day2.tsv"
out  <- if (length(a) >= 3) sub("\\.png$", "", a[3]) else "results/aggregate/plots/composition_arm"

comp_map <- c(
  Tumor = "tumor",
  Fibroblast = "stroma", SmoothMuscle = "stroma", Endothelial = "stroma",
  Adipocyte = "stroma", unassigned = "stroma",
  Macrophages = "immune", `T cells` = "immune", `NK cells` = "immune", ILC = "immune",
  DC = "immune", `Plasma cells` = "immune", `Mast cells` = "immune",
  Neutrophils = "immune", `Epithelial cells` = "immune")
arm_cols <- c(Control = "#4a86e8", MBRT = "#e69138", SBRT = "#cc0000")

# Order cell types by compartment then overall abundance for a stable stack/axis.
d0 <- read_tsv(comp, show_col_types = FALSE) %>%
  filter(dataset == "Mutter_02") %>%
  mutate(arm = recode(condition, MBRT_day2 = "MBRT", SBRT_day2 = "SBRT"),
         arm = factor(arm, levels = c("Control", "MBRT", "SBRT")),
         compartment = factor(comp_map[cell_type], levels = c("tumor", "stroma", "immune")))
ord <- d0 %>% group_by(cell_type, compartment) %>% summarise(n = sum(n_cells), .groups = "drop") %>%
  arrange(compartment, desc(n)) %>% pull(cell_type)
d0  <- mutate(d0, cell_type = factor(cell_type, levels = ord))

# Panel A: mean per-arm proportion (pooled cells within arm), stacked.
barsd <- d0 %>% group_by(arm, cell_type) %>% summarise(n = sum(n_cells), .groups = "drop") %>%
  group_by(arm) %>% mutate(frac = n / sum(n)) %>% ungroup()
ct_cols <- setNames(scales::hue_pal()(length(ord)), ord)

pA <- ggplot(barsd, aes(arm, frac, fill = cell_type)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.15) +
  scale_fill_manual(values = ct_cols, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "proportion of cells", tag = "A") +
  guides(fill = guide_legend(ncol = 1, keyheight = unit(9, "pt"))) +
  theme_scifig(base_size = 11) +
  theme(legend.text = element_text(size = 8))

# Panel B: propeller effect + 95% CI, one facet per contrast.
td <- read_tsv(test, show_col_types = FALSE) %>%
  mutate(contrast = factor(contrast, levels = c("MBRT_vs_Ctrl", "SBRT_vs_Ctrl", "MBRT_vs_SBRT")),
         cell_type = factor(cell_type, levels = rev(ord)),
         sig = padj < 0.05)

pB <- ggplot(td, aes(log2FC_logit, cell_type, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
  geom_errorbarh(aes(xmin = ci_low_log2, xmax = ci_high_log2), height = 0, linewidth = 0.4) +
  geom_point(size = 1.8) +
  facet_wrap(~ contrast, nrow = 1) +
  scale_color_manual(values = c(`FALSE` = "gray55", `TRUE` = "#cc0000"),
                     labels = c("n.s.", "padj < 0.05"), name = NULL) +
  labs(x = "composition effect (logit log2FC)", y = NULL, tag = "B") +
  theme_scifig(base_size = 11)

p <- pA / pB + plot_layout(heights = c(1, 1.4))
save_plot(p, out, w = 10, h = 9)
