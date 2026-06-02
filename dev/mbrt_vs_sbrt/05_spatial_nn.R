## 05_spatial_nn.R
## k=20 nearest-neighbor immune fraction per cell, per (slide x condition).
## Aggregate per (cell_type_major x timepoint x condition).
## Necrosis_zone cells excluded at day2/day6 (default behavior).

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(RANN)
  library(future)
  library(future.apply)
})

DATA_DIR <- "data"
PLOT_DIR <- "plots"
OBJECTS_DIR <- "/mnt/data/projects/spatial-rads/analysis/objects/mbrt_vs_sbrt"
NN_CACHE_DIR <- file.path(OBJECTS_DIR, "nn_cache")  # NN indices on data disk
if (!dir.exists(NN_CACHE_DIR)) dir.create(NN_CACHE_DIR, recursive = TRUE)

obj <- readRDS(file.path(OBJECTS_DIR, "seurat_filtered.rds"))
cat(sprintf("Loaded: %d cells\n", ncol(obj)))

# Use full meta (do not subset) — we want NN within all cells of a slide.
meta <- as_tibble(obj@meta.data) %>%
  mutate(cell_id = colnames(obj),
         x_slide_mm = obj$x_slide_mm,
         y_slide_mm = obj$y_slide_mm)

slide_col <- if ("Slide" %in% colnames(meta)) "Slide" else "slide_ID_numeric"
slides <- unique(meta[[slide_col]])
cat(sprintf("Slides: %d\n", length(slides)))

immune_types <- c("T.cell", "B.cell", "myeloid", "NK", "plasma")

# Compute NN per slide
nn_results <- list()
t0 <- Sys.time()
for (s in slides) {
  m_s <- meta[meta[[slide_col]] == s, ]
  if (nrow(m_s) < 50) next
  cache_file <- file.path(NN_CACHE_DIR, sprintf("nn_slide_%s.rds", s))
  if (file.exists(cache_file)) {
    nn <- readRDS(cache_file)
  } else {
    coords <- as.matrix(m_s[, c("x_slide_mm", "y_slide_mm")])
    k_use <- min(21, nrow(m_s))
    nn <- RANN::nn2(coords, k = k_use)
    saveRDS(nn, cache_file)
  }
  nn_idx <- nn$nn.idx[, -1, drop = FALSE]  # drop self
  ct_neigh_mat <- matrix(m_s$cell_type_major[nn_idx], nrow = nrow(nn_idx), ncol = ncol(nn_idx))
  immune_logical <- matrix(ct_neigh_mat %in% immune_types,
                           nrow = nrow(ct_neigh_mat), ncol = ncol(ct_neigh_mat))
  tumor_logical  <- matrix(ct_neigh_mat == "tumor_epithelial",
                           nrow = nrow(ct_neigh_mat), ncol = ncol(ct_neigh_mat))
  immune_frac <- rowMeans(immune_logical)
  tumor_frac  <- rowMeans(tumor_logical)
  out <- m_s %>%
    select(cell_id, Condition, treatment, timepoint_h, cell_type_major,
           necrosis_zone, x_slide_mm, y_slide_mm) %>%
    mutate(immune_neighbor_frac = immune_frac,
           tumor_neighbor_frac  = tumor_frac,
           Slide = s)
  nn_results[[as.character(s)]] <- out
  cat(sprintf("  Slide %s: %d cells, %.1f min elapsed\n",
              s, nrow(m_s),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
nn_df <- bind_rows(nn_results)
cat(sprintf("Total cells with NN: %d\n", nrow(nn_df)))

# Apply necrosis-zone exclusion at day2/day6 for aggregation
nn_df_for_agg <- nn_df %>%
  filter(!(timepoint_h %in% c(48, 144) & necrosis_zone %in% TRUE))
cat(sprintf("After necrosis-zone exclusion at day2/day6: %d cells\n", nrow(nn_df_for_agg)))

# Aggregate per (cell_type_major x timepoint x condition x treatment)
agg <- nn_df_for_agg %>%
  group_by(cell_type_major, timepoint_h, treatment, Condition) %>%
  summarise(mean_immune_frac = mean(immune_neighbor_frac, na.rm = TRUE),
            sd_immune_frac   = sd(immune_neighbor_frac, na.rm = TRUE),
            mean_tumor_frac  = mean(tumor_neighbor_frac, na.rm = TRUE),
            sd_tumor_frac    = sd(tumor_neighbor_frac, na.rm = TRUE),
            n_cells = n(),
            .groups = "drop")
write_tsv(agg, file.path(DATA_DIR, "spatial_nn_kinetics.tsv"))
cat(sprintf("Wrote: %s (%d rows)\n",
            file.path(DATA_DIR, "spatial_nn_kinetics.tsv"), nrow(agg)))

# Per-cell long output (for downstream re-use)
write_tsv(nn_df %>% select(cell_id, Condition, cell_type_major,
                           immune_neighbor_frac, tumor_neighbor_frac, necrosis_zone),
          file.path(DATA_DIR, "spatial_nn_per_cell.tsv.gz"))

# === Plots ===
TREATMENT_COLORS <- c("MBRT" = "#E74C3C", "SBRT" = "#3498DB", "NT" = "#7F8C8D")

flank_conditions <- c("Control",
                      "MBRT_1h", "MBRT_4h", "MBRT_day2", "MBRT_day6",
                      "SBRT_4h", "SBRT_day2", "SBRT_day6")

# 1. Immune fraction kinetics per cell type (line plots)
p_imm <- agg %>%
  filter(Condition %in% flank_conditions, treatment %in% c("MBRT","SBRT","NT")) %>%
  ggplot(aes(x = factor(timepoint_h), y = mean_immune_frac,
             color = treatment, group = treatment)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_immune_frac - sd_immune_frac,
                    ymax = mean_immune_frac + sd_immune_frac),
                width = 0.2, alpha = 0.4) +
  facet_wrap(~cell_type_major, scales = "free_y", ncol = 3) +
  scale_color_manual(values = TREATMENT_COLORS) +
  labs(title = "Spatial NN: immune-neighbor fraction (k=20)",
       x = "Hours post-RT", y = "Mean fraction immune neighbors",
       caption = "n=1 per condition; necrosis-zone cells excluded at day2/day6.") +
  theme_bw()
ggsave(file.path(PLOT_DIR, "nn_immune_fraction_kinetics.png"),
       plot = p_imm, width = 12, height = 9, dpi = 150)

# 2. Density per timepoint × cell type
p_dens <- nn_df_for_agg %>%
  filter(Condition %in% flank_conditions, treatment %in% c("MBRT","SBRT","NT")) %>%
  mutate(timepoint_label = case_when(
    timepoint_h == 0 ~ "Control",
    timepoint_h == 1 ~ "1h",
    timepoint_h == 4 ~ "4h",
    timepoint_h == 48 ~ "day2",
    timepoint_h == 144 ~ "day6")) %>%
  ggplot(aes(x = immune_neighbor_frac, fill = treatment)) +
  geom_density(alpha = 0.4) +
  facet_grid(cell_type_major ~ timepoint_label, scales = "free_y") +
  scale_fill_manual(values = TREATMENT_COLORS) +
  labs(title = "Distribution of immune-neighbor fraction (k=20)",
       x = "Fraction immune neighbors", y = "Density") +
  theme_bw() +
  theme(strip.text.y = element_text(size = 7))
ggsave(file.path(PLOT_DIR, "nn_immune_density.png"),
       plot = p_dens, width = 14, height = 12, dpi = 130)

# 3. Tumor neighbor fraction (mutual exclusivity check)
p_tum <- agg %>%
  filter(Condition %in% flank_conditions, treatment %in% c("MBRT","SBRT","NT")) %>%
  ggplot(aes(x = factor(timepoint_h), y = mean_tumor_frac,
             color = treatment, group = treatment)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  facet_wrap(~cell_type_major, scales = "free_y", ncol = 3) +
  scale_color_manual(values = TREATMENT_COLORS) +
  labs(title = "Spatial NN: tumor-neighbor fraction (k=20)",
       x = "Hours post-RT", y = "Mean fraction tumor neighbors") +
  theme_bw()
ggsave(file.path(PLOT_DIR, "nn_tumor_fraction_kinetics.png"),
       plot = p_tum, width = 12, height = 9, dpi = 150)

cat("\n=== Sanity check: tumor immune-neighbor fraction by condition ===\n")
print(agg %>%
        filter(cell_type_major == "tumor_epithelial",
               Condition %in% flank_conditions) %>%
        select(Condition, treatment, timepoint_h, mean_immune_frac, n_cells) %>%
        arrange(timepoint_h, treatment))

cat("\n=== 05_spatial_nn.R complete ===\n")
