## 13_necrosis_best_figure.R
## Best single figure for the central-necrosis story:
## Centerpiece: MBRT_day2 vs SBRT_day2 spatial maps (LARGE) showing fragmented vs consolidated necrosis.
## Supporting: Control reference + day6 convergence + kinetic + fragmentation index plots.

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
})

DATA_DIR <- "data"
PLOT_DIR <- "plots"
OBJECTS_DIR <- "/mnt/data/projects/spatial-rads/analysis/objects/mbrt_vs_sbrt"

obj <- readRDS(file.path(OBJECTS_DIR, "seurat_filtered.rds"))
cat(sprintf("Loaded: %d cells\n", ncol(obj)))

necr_q <- read_tsv(file.path(DATA_DIR, "necrosis_quantification.tsv"),
                   show_col_types = FALSE)
necr_cc <- read_tsv(file.path(DATA_DIR, "necrosis_connected_components.tsv"),
                    show_col_types = FALSE)

# Add fragmentation: % in largest cluster
necr_cc <- necr_cc %>%
  mutate(pct_in_largest = 100 * largest_region_bins / pmax(necrotic_bins, 1))

make_cell_df <- function(cond) {
  as_tibble(obj@meta.data) %>%
    filter(Condition == cond, !is.na(necrosis_zone)) %>%
    mutate(zone_label = ifelse(necrosis_zone, "necrosis", "viable"))
}

make_spatial_panel <- function(cond, title_prefix = "") {
  m <- make_cell_df(cond)
  if (nrow(m) == 0) return(ggplot() + theme_void())

  q_row <- necr_q %>% filter(Condition == cond)
  cc_row <- necr_cc %>% filter(Condition == cond)
  pct_calib  <- if (nrow(q_row) > 0) round(q_row$necr_pct_calib[1], 1) else NA
  largest_mm2 <- if (nrow(cc_row) > 0) round(cc_row$largest_region_mm2[1], 2) else NA
  pct_largest <- if (nrow(cc_row) > 0) round(cc_row$pct_in_largest[1], 0) else NA

  annotation <- sprintf("%.1f%% necrosis · largest %.2f mm² · %d%% in single cluster",
                       pct_calib, largest_mm2, pct_largest)

  ggplot(m, aes(x = x_slide_mm, y = y_slide_mm, color = zone_label)) +
    geom_point(size = 0.25, alpha = 0.6) +
    scale_color_manual(values = c("viable" = "grey85", "necrosis" = "#E74C3C"),
                       guide = "none") +
    labs(title = paste0(title_prefix, cond), subtitle = annotation) +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(size = 10, color = "grey20"),
          aspect.ratio = 1,
          plot.background = element_rect(fill = "white", color = NA))
}

# Day 2 as centerpiece (biggest panels)
p_m2 <- make_spatial_panel("MBRT_day2")
p_s2 <- make_spatial_panel("SBRT_day2")
# Day 6 supporting
p_m6 <- make_spatial_panel("MBRT_day6")
p_s6 <- make_spatial_panel("SBRT_day6")
# Control reference
p_ctrl <- make_spatial_panel("Control")

# Kinetics
TREATMENT_COLORS <- c("MBRT" = "#E74C3C", "SBRT" = "#3498DB", "NT" = "#7F8C8D")

# Filter MBRT_1h artifact for kinetics
necr_q_clean <- necr_q %>%
  filter(Condition != "MBRT_1h",
         Condition %in% c("Control","MBRT_4h","MBRT_day2","MBRT_day6",
                          "SBRT_4h","SBRT_day2","SBRT_day6"))
necr_cc_clean <- necr_cc %>%
  filter(Condition != "MBRT_1h",
         Condition %in% c("Control","MBRT_4h","MBRT_day2","MBRT_day6",
                          "SBRT_4h","SBRT_day2","SBRT_day6")) %>%
  mutate(treatment = case_when(Condition == "Control" ~ "NT",
                               grepl("MBRT", Condition) ~ "MBRT",
                               grepl("SBRT", Condition) ~ "SBRT"))

p_pct <- necr_q_clean %>%
  ggplot(aes(x = factor(timepoint_h), y = necr_pct_calib,
             color = treatment, group = treatment)) +
  geom_line(linewidth = 1.2) + geom_point(size = 4) +
  scale_color_manual(values = TREATMENT_COLORS, name = NULL) +
  labs(title = "Total necrosis fraction (Control-calibrated)",
       x = "Hours post-RT", y = "% cells in necrosis-zone") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", size = 12))

p_frag <- necr_cc_clean %>%
  ggplot(aes(x = factor(timepoint_h), y = pct_in_largest,
             color = treatment, group = treatment)) +
  geom_line(linewidth = 1.2) + geom_point(size = 4) +
  scale_color_manual(values = TREATMENT_COLORS, name = NULL) +
  labs(title = "Necrosis consolidation",
       subtitle = "% of necrotic bins in single largest cluster (high = consolidated central; low = fragmented)",
       x = "Hours post-RT", y = "% in largest cluster") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", size = 12))

# Build layout. Top: BIG day2 spatial maps side-by-side. Middle: smaller row of Control/day6.
# Bottom: kinetics.
top <- p_m2 | p_s2
middle <- p_ctrl | p_m6 | p_s6
bottom <- p_pct | p_frag

composite <- top / middle / bottom +
  plot_layout(heights = c(2.0, 1.0, 1.0)) +
  plot_annotation(
    title = "MBRT spares central tumor architecture vs SBRT at day 2",
    subtitle = "Day 2 SBRT: necrosis consolidated in single 4.46 mm² central region (84% in one cluster). MBRT: less necrosis (4.5% vs 10%) AND fragmented across small regions (33% in largest, 1.04 mm²). By day 6, SBRT and MBRT necrosis converge.",
    caption = "Necrosis-zone defined per cell as k=20 NN-distance > 95th percentile of Control. Connected components on 0.1mm grid (≥50% necrosis bins per bin, ≥3 cells per bin). n=1 per condition; descriptive only.",
    theme = theme(plot.title = element_text(face = "bold", size = 16),
                  plot.subtitle = element_text(size = 12),
                  plot.caption = element_text(size = 9, hjust = 0))
  )

ggsave(file.path(PLOT_DIR, "necrosis_best_figure.png"),
       plot = composite, width = 17, height = 18, dpi = 150)

cat("\n=== 13_necrosis_best_figure.R complete ===\n")
cat(sprintf("Saved: %s\n", file.path(PLOT_DIR, "necrosis_best_figure.png")))
