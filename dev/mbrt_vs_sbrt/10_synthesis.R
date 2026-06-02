## 10_synthesis.R
## Pull together headline MBRT-vs-SBRT findings into a single composite figure.
## Reads from existing TSVs; lightweight.

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

DATA_DIR <- "data"
PLOT_DIR <- "plots"

TREATMENT_COLORS <- c("MBRT" = "#E74C3C", "SBRT" = "#3498DB", "NT" = "#7F8C8D")
flank_conds <- c("Control", "MBRT_1h", "MBRT_4h", "MBRT_day2", "MBRT_day6",
                 "SBRT_4h", "SBRT_day2", "SBRT_day6")

# === Panel A: Myeloid composition trajectory ===
comp <- read_tsv(file.path(DATA_DIR, "composition_trajectories.tsv"),
                 show_col_types = FALSE) %>%
  filter(universe == "major9_only", Condition %in% flank_conds,
         treatment %in% c("MBRT","SBRT","NT"))

p_mye <- comp %>%
  filter(cell_type_major == "myeloid") %>%
  ggplot(aes(x = factor(timepoint_h), y = prop, color = treatment, group = treatment)) +
  geom_line(linewidth = 1.0) + geom_point(size = 3) +
  scale_color_manual(values = TREATMENT_COLORS) +
  labs(title = "A. Myeloid recruitment kinetics",
       x = "Hours post-RT", y = "Myeloid proportion",
       caption = "MBRT recruits myeloid earlier and stronger") +
  theme_bw(base_size = 11)

# === Panel B: Tumor immune-neighbor fraction ===
nn <- read_tsv(file.path(DATA_DIR, "spatial_nn_kinetics.tsv"),
               show_col_types = FALSE) %>%
  filter(Condition %in% flank_conds, treatment %in% c("MBRT","SBRT","NT"))

p_nn <- nn %>%
  filter(cell_type_major == "tumor_epithelial") %>%
  ggplot(aes(x = factor(timepoint_h), y = mean_immune_frac,
             color = treatment, group = treatment)) +
  geom_line(linewidth = 1.0) + geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_immune_frac - sd_immune_frac,
                    ymax = mean_immune_frac + sd_immune_frac),
                width = 0.2, alpha = 0.4) +
  scale_color_manual(values = TREATMENT_COLORS) +
  labs(title = "B. Tumor cells: spatial immune neighborhood (k=20)",
       x = "Hours post-RT", y = "Mean fraction immune neighbors",
       caption = "MBRT tumor cells are surrounded by immune cells earlier") +
  theme_bw(base_size = 11)

# === Panel C: Necrosis kinetics (Control-calibrated) ===
necr <- read_tsv(file.path(DATA_DIR, "necrosis_quantification.tsv"),
                 show_col_types = FALSE) %>%
  filter(Condition %in% flank_conds, treatment %in% c("MBRT","SBRT","NT"))

p_necr <- necr %>%
  ggplot(aes(x = factor(timepoint_h), y = necr_pct_calib,
             color = treatment, group = treatment)) +
  geom_line(linewidth = 1.0) + geom_point(size = 3) +
  scale_color_manual(values = TREATMENT_COLORS) +
  labs(title = "C. Necrosis fraction (Control-calibrated)",
       x = "Hours post-RT", y = "% cells > 95th pct of Control NN-distance",
       caption = "Day2 SBRT 2.2x MBRT; convergence by day6") +
  theme_bw(base_size = 11)

# === Panel D: Connected components — largest contiguous necrotic region ===
cc <- read_tsv(file.path(DATA_DIR, "necrosis_connected_components.tsv"),
               show_col_types = FALSE) %>%
  filter(Condition %in% flank_conds)

p_cc <- cc %>%
  ggplot(aes(x = factor(timepoint_h), y = largest_region_mm2,
             color = treatment, group = treatment)) +
  geom_line(linewidth = 1.0) + geom_point(size = 3) +
  scale_color_manual(values = TREATMENT_COLORS) +
  labs(title = "D. Largest contiguous necrotic region",
       x = "Hours post-RT", y = "Largest region size (mm²)",
       caption = "Day2 SBRT 4.46 mm² (single central region) vs MBRT 1.04 mm² (fragmented)") +
  theme_bw(base_size = 11)

# === Composite ===
composite <- (p_mye + p_nn) / (p_necr + p_cc) +
  plot_annotation(
    title = "MBRT vs SBRT: comprehensive comparison (4T1 flank, Mutter_01)",
    subtitle = "Earlier immune engagement, less central necrosis at day2",
    caption = "n=1 per condition; descriptive only. See spatial-rads.org Layer 2 Extended for full results.",
    theme = theme(plot.title = element_text(size = 14, face = "bold"),
                  plot.subtitle = element_text(size = 11)))

ggsave(file.path(PLOT_DIR, "synthesis_executive_summary.png"),
       plot = composite, width = 14, height = 9, dpi = 150)

# === Synthesis summary table ===
summary_df <- tribble(
  ~axis, ~metric, ~MBRT_value, ~SBRT_value, ~delta_or_ratio, ~timepoint, ~interpretation,
  "Composition", "myeloid_proportion", 0.42, 0.31, "+0.11 ΔMBRT−SBRT", "day2", "MBRT recruits more myeloid",
  "Composition", "myeloid_proportion", 0.20, 0.12, "+0.08 ΔMBRT−SBRT", "4h", "Earlier MBRT myeloid recruitment",
  "Composition", "tumor_proportion", 0.41, 0.51, "−0.10 ΔMBRT−SBRT", "4h", "SBRT preserves more tumor early",
  "Composition", "B_cell_proportion", 0.085, 0.020, "+0.07 ΔMBRT−SBRT", "day6", "MBRT-driven late B cell expansion",
  "Composition", "plasma_proportion", 0.045, 0.000, "+0.045 ΔMBRT−SBRT", "day6", "MBRT-driven plasmablast emergence",
  "Spatial NN", "tumor_immune_frac", 0.375, 0.308, "+0.067 ΔMBRT−SBRT", "4h", "MBRT tumor mixed with immune earlier",
  "Spatial NN", "tumor_immune_frac", 0.548, 0.457, "+0.091 ΔMBRT−SBRT", "day2", "Gap widens at day2",
  "Spatial NN", "tumor_immune_frac", 0.719, 0.728, "−0.009 ΔMBRT−SBRT", "day6", "Convergence by day6",
  "Necrosis", "calibrated_pct", 4.5, 10.0, "0.45x SBRT", "day2", "MBRT day2 has 2.2x less central necrosis",
  "Necrosis", "calibrated_pct", 12.1, 8.7, "1.39x SBRT", "day6", "MBRT necrosis catches up by day6",
  "Necrosis", "largest_region_mm2", 1.04, 4.46, "0.23x SBRT", "day2", "SBRT day2 consolidated central; MBRT fragmented",
  "Necrosis", "largest_region_mm2", 4.93, 3.57, "1.38x SBRT", "day6", "MBRT day6 region size now exceeds SBRT",
  "DEG (tumor)", "Cd74_log2FC_4h", 0.68, NA, "MBRT-up", "4h", "Antigen presentation invariant chain",
  "DEG (tumor)", "B2m_log2FC_4h", 0.54, NA, "MBRT-up", "4h", "MHC class I light chain",
  "DEG (tumor)", "Spp1_log2FC_4h", -0.86, NA, "SBRT-up", "4h", "Hypoxia/scar marker",
  "DEG (myeloid)", "Cd163_log2FC_day2", 2.51, NA, "MBRT-up", "day2", "M2 macrophage scavenger receptor",
  "DEG (myeloid)", "Mrc1_log2FC_day2", 2.12, NA, "MBRT-up", "day2", "M2 macrophage CD206",
  "DEG (myeloid)", "Ccl8_log2FC_day2", 2.50, NA, "MBRT-up", "day2", "T cell chemokine",
  "DEG (B.cell)", "Ighg1_log2FC_day6", 4.40, NA, "MBRT-up", "day6", "Class-switched IgG heavy chain",
  "DEG (B.cell)", "Igkc_log2FC_day6", 3.30, NA, "MBRT-up", "day6", "Antibody kappa light chain",
  "DEG (B.cell)", "Xbp1_log2FC_day6", 1.06, NA, "MBRT-up", "day6", "Plasma cell differentiation TF"
)
write_tsv(summary_df, file.path(DATA_DIR, "synthesis_summary.tsv"))

cat("\n=== Synthesis composite saved ===\n")
cat(sprintf("File: %s\n", file.path(PLOT_DIR, "synthesis_executive_summary.png")))
cat(sprintf("Summary table: %s\n", file.path(DATA_DIR, "synthesis_summary.tsv")))
cat("\n=== 10_synthesis.R complete ===\n")
