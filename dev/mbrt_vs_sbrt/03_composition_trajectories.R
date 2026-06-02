## 03_composition_trajectories.R
## Cell type proportions per (Condition, timepoint, treatment), trajectories,
## and MBRT vs SBRT delta heatmap. Composition uses ALL cells (no necrosis exclusion)
## to avoid biasing the universe.

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
})

DATA_DIR <- "data"
PLOT_DIR <- "plots"
OBJECTS_DIR <- "/mnt/data/projects/spatial-rads/analysis/objects/mbrt_vs_sbrt"

obj <- readRDS(file.path(OBJECTS_DIR, "seurat_filtered.rds"))
cat(sprintf("Loaded: %d cells\n", ncol(obj)))

# Universe (post-Layer-1) for the unfiltered denominator
counts_universe <- read_tsv(file.path(DATA_DIR, "cell_counts_universe.tsv"),
                            show_col_types = FALSE)

# 9-major-only view from filtered object
counts_major <- as_tibble(obj@meta.data) %>%
  count(Condition, treatment, timepoint_h, cell_type_major)

# Compute proportions for both universes
prop_major <- counts_major %>%
  group_by(Condition) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(universe = "major9_only")

prop_universe <- counts_universe %>%
  group_by(Condition) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(treatment = case_when(Condition == "Control" ~ "NT",
                               grepl("MBRT", Condition) ~ "MBRT",
                               grepl("SBRT", Condition) ~ "SBRT",
                               TRUE ~ "other"),
         timepoint_h = case_when(Condition == "Control" ~ 0,
                                 grepl("1h", Condition) ~ 1,
                                 grepl("4h", Condition) ~ 4,
                                 grepl("day2", Condition) ~ 48,
                                 grepl("day6", Condition) ~ 144,
                                 grepl("day8", Condition) ~ 192,
                                 grepl("day10", Condition) ~ 240,
                                 TRUE ~ NA_real_),
         universe = "post_layer1") %>%
  filter(!grepl("^Tongue", Condition))  # flank-only — drops anything that snuck in

composition_long <- bind_rows(prop_major, prop_universe)
write_tsv(composition_long, file.path(DATA_DIR, "composition_trajectories.tsv"))

# Restrict to flank conditions for plotting
flank_conditions <- c("Control",
                      "MBRT_1h", "MBRT_4h", "MBRT_day2", "MBRT_day6",
                      "SBRT_4h", "SBRT_day2", "SBRT_day6")

# === Plot 1: stacked bar (major9 only) ===
celltype_palette <- c(
  "tumor_epithelial" = "#E41A1C",
  "endothelial"      = "#984EA3",
  "fibroblast"       = "#FF7F00",
  "smooth_muscle"    = "#A65628",
  "T.cell"           = "#377EB8",
  "B.cell"           = "#4DAF4A",
  "myeloid"          = "#FFFF33",
  "NK"               = "#F781BF",
  "plasma"           = "#999999"
)

p_stacked <- prop_major %>%
  filter(Condition %in% flank_conditions) %>%
  mutate(Condition = factor(Condition, levels = flank_conditions),
         cell_type_major = factor(cell_type_major,
                                  levels = names(celltype_palette))) %>%
  ggplot(aes(x = Condition, y = prop, fill = cell_type_major)) +
  geom_col() +
  scale_fill_manual(values = celltype_palette, name = "Cell type") +
  labs(title = "Cell type composition (4T1 flank, 9 major types)",
       x = NULL, y = "Proportion") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(PLOT_DIR, "composition_stacked.png"),
       plot = p_stacked, width = 11, height = 5, dpi = 150)

# === Plot 2: trajectories — proportion over time, faceted by cell type ===
TREATMENT_COLORS <- c("MBRT" = "#E74C3C", "SBRT" = "#3498DB", "NT" = "#7F8C8D")

p_traj <- prop_major %>%
  filter(Condition %in% flank_conditions, treatment %in% c("MBRT","SBRT","NT")) %>%
  ggplot(aes(x = factor(timepoint_h), y = prop, color = treatment, group = treatment)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  facet_wrap(~cell_type_major, scales = "free_y", ncol = 3) +
  scale_color_manual(values = TREATMENT_COLORS) +
  labs(title = "Cell type proportions over time (flank)",
       x = "Hours post-RT", y = "Proportion of major-9 universe") +
  theme_bw()
ggsave(file.path(PLOT_DIR, "composition_trajectories.png"),
       plot = p_traj, width = 12, height = 8, dpi = 150)

# === Plot 3: delta MBRT - SBRT heatmap ===
delta_df <- prop_major %>%
  filter(treatment %in% c("MBRT", "SBRT"), Condition %in% flank_conditions) %>%
  select(cell_type_major, timepoint_h, treatment, prop) %>%
  pivot_wider(names_from = treatment, values_from = prop) %>%
  mutate(delta_MBRT_SBRT = MBRT - SBRT) %>%
  filter(!is.na(delta_MBRT_SBRT))

p_delta <- delta_df %>%
  ggplot(aes(x = factor(timepoint_h), y = cell_type_major, fill = delta_MBRT_SBRT)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%+.3f", delta_MBRT_SBRT)),
            size = 3, color = "black") +
  scale_fill_gradient2(low = "#3498DB", mid = "white", high = "#E74C3C",
                       midpoint = 0, name = "MBRT - SBRT\n(prop)") +
  labs(title = "Composition delta: MBRT - SBRT proportion (4T1 flank)",
       x = "Hours post-RT", y = NULL) +
  theme_bw() + theme(panel.grid = element_blank())
ggsave(file.path(PLOT_DIR, "composition_delta_mbrt_vs_sbrt.png"),
       plot = p_delta, width = 8, height = 6, dpi = 150)

write_tsv(delta_df, file.path(DATA_DIR, "composition_delta_mbrt_vs_sbrt.tsv"))

cat("\n=== Composition snapshot ===\n")
print(prop_major %>%
        filter(Condition %in% flank_conditions) %>%
        select(Condition, cell_type_major, n, prop) %>%
        pivot_wider(names_from = cell_type_major, values_from = prop) %>%
        select(Condition, tumor_epithelial, T.cell, myeloid, B.cell, NK,
               endothelial, fibroblast, smooth_muscle, plasma))

cat("\n=== 03_composition_trajectories.R complete ===\n")
