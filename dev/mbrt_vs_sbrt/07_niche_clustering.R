## 07_niche_clustering.R
## Data-driven niche clustering: each cell is described by the cell-type
## composition of its k=20 nearest neighbors. k-means on those vectors
## identifies spatial niches. Compare niche frequency MBRT vs SBRT over time.

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(RANN)
  library(patchwork)
})

DATA_DIR <- "data"
PLOT_DIR <- "plots"
OBJECTS_DIR <- "/mnt/data/projects/spatial-rads/analysis/objects/mbrt_vs_sbrt"
tmp_data <- "/mnt/data/projects/spatial-rads/analysis/tmp"
if (!dir.exists(tmp_data)) dir.create(tmp_data, recursive = TRUE)
Sys.setenv(TMPDIR = tmp_data, TMP = tmp_data, TEMP = tmp_data)

obj <- readRDS(file.path(OBJECTS_DIR, "seurat_filtered.rds"))
cat(sprintf("Loaded: %d cells\n", ncol(obj)))

slide_col <- if ("Slide" %in% colnames(obj@meta.data)) "Slide" else "slide_ID_numeric"
slides <- unique(obj@meta.data[[slide_col]])
MAJOR_TYPES <- sort(unique(obj$cell_type_major))
cat(sprintf("Major types: %s\n", paste(MAJOR_TYPES, collapse=", ")))

meta <- as_tibble(obj@meta.data) %>% mutate(cell_id = colnames(obj))

# Build neighborhood composition vector per cell
cat("\nBuilding NN compositions per slide...\n")
comp_list <- list()
t0 <- Sys.time()
for (s in slides) {
  m_s <- meta %>% filter(!!sym(slide_col) == s)
  if (nrow(m_s) < 50) next
  coords <- as.matrix(m_s[, c("x_slide_mm", "y_slide_mm")])
  k_use <- min(21, nrow(m_s))
  nn <- RANN::nn2(coords, k = k_use)
  nn_idx <- nn$nn.idx[, -1, drop = FALSE]
  ct_neigh_mat <- matrix(m_s$cell_type_major[nn_idx],
                         nrow = nrow(nn_idx), ncol = ncol(nn_idx))
  # Composition vector per cell (length = #MAJOR_TYPES)
  comp_mat <- sapply(MAJOR_TYPES, function(ct) {
    rowMeans(matrix(ct_neigh_mat == ct, nrow = nrow(ct_neigh_mat), ncol = ncol(ct_neigh_mat)))
  })
  comp_list[[as.character(s)]] <- list(cells = m_s$cell_id, comp = comp_mat)
  cat(sprintf("  Slide %s: %d cells, %.1f min\n", s, nrow(m_s),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

# Combine into single matrix
all_cells <- do.call(c, lapply(comp_list, `[[`, "cells"))
all_comp  <- do.call(rbind, lapply(comp_list, `[[`, "comp"))
rownames(all_comp) <- all_cells
cat(sprintf("Total cells with composition vectors: %d\n", nrow(all_comp)))

# k-means cluster
set.seed(42)
K <- 6
km <- kmeans(all_comp, centers = K, nstart = 10, iter.max = 50)
cat(sprintf("k-means done. Cluster sizes:\n"))
print(table(km$cluster))

# Niche label per cell
niche_df <- tibble(cell_id = all_cells, niche = paste0("N", km$cluster))

# Niche centroids (mean composition vector per niche)
centroid_df <- as_tibble(km$centers) %>%
  mutate(niche = paste0("N", row_number())) %>%
  pivot_longer(-niche, names_to = "cell_type", values_to = "mean_frac")
write_tsv(centroid_df, file.path(DATA_DIR, "niche_centroids.tsv"))

# Niche frequency per (Condition × niche)
niche_meta <- meta %>%
  inner_join(niche_df, by = "cell_id")
niche_freq <- niche_meta %>%
  count(Condition, treatment, timepoint_h, niche) %>%
  group_by(Condition) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup()
write_tsv(niche_freq, file.path(DATA_DIR, "niche_frequency.tsv"))

# Per-cell niche assignments
write_tsv(niche_meta %>% select(cell_id, Condition, niche, cell_type_major,
                                 x_slide_mm, y_slide_mm),
          file.path(DATA_DIR, "niche_per_cell.tsv.gz"))

# Plot 1: niche centroid heatmap (cell-type composition profile of each niche)
p1 <- centroid_df %>%
  ggplot(aes(x = niche, y = cell_type, fill = mean_frac)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", mean_frac)), size = 3) +
  scale_fill_viridis_c(option = "magma", name = "Mean\nfraction") +
  labs(title = "Niche centroid composition",
       subtitle = "Each niche characterized by its k=20 neighbor cell-type fractions",
       x = NULL, y = NULL) +
  theme_bw() + theme(panel.grid = element_blank())
ggsave(file.path(PLOT_DIR, "niche_centroids.png"),
       plot = p1, width = 9, height = 6, dpi = 150)

# Plot 2: niche frequency over time, MBRT vs SBRT vs Control
TREATMENT_COLORS <- c("MBRT" = "#E74C3C", "SBRT" = "#3498DB", "NT" = "#7F8C8D")
flank_conds <- c("Control", "MBRT_1h", "MBRT_4h", "MBRT_day2", "MBRT_day6",
                 "SBRT_4h", "SBRT_day2", "SBRT_day6")
p2 <- niche_freq %>%
  filter(Condition %in% flank_conds, treatment %in% c("MBRT","SBRT","NT")) %>%
  ggplot(aes(x = factor(timepoint_h), y = frac,
             color = treatment, group = treatment)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  facet_wrap(~niche, scales = "free_y", ncol = 3) +
  scale_color_manual(values = TREATMENT_COLORS) +
  labs(title = "Spatial niche frequency over time",
       x = "Hours post-RT", y = "Fraction of cells in niche") +
  theme_bw()
ggsave(file.path(PLOT_DIR, "niche_frequency_kinetics.png"),
       plot = p2, width = 11, height = 7, dpi = 150)

# Plot 3: niche delta MBRT - SBRT at each timepoint
delta_df <- niche_freq %>%
  filter(treatment %in% c("MBRT","SBRT"), Condition %in% flank_conds) %>%
  select(niche, timepoint_h, treatment, frac) %>%
  pivot_wider(names_from = treatment, values_from = frac) %>%
  mutate(delta = MBRT - SBRT) %>%
  filter(!is.na(delta))

p3 <- delta_df %>%
  ggplot(aes(x = factor(timepoint_h), y = niche, fill = delta)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%+.3f", delta)), size = 3) +
  scale_fill_gradient2(low = "#3498DB", mid = "white", high = "#E74C3C",
                       midpoint = 0, name = "MBRT - SBRT") +
  labs(title = "Niche-frequency delta MBRT - SBRT (4T1 flank)",
       x = "Hours post-RT", y = NULL) +
  theme_bw() + theme(panel.grid = element_blank())
ggsave(file.path(PLOT_DIR, "niche_delta_mbrt_vs_sbrt.png"),
       plot = p3, width = 9, height = 5, dpi = 150)

# Plot 4: spatial niche maps (4 representative conditions)
p4_list <- list()
for (cd in c("Control", "MBRT_4h", "MBRT_day2", "SBRT_day2")) {
  d <- niche_meta %>% filter(Condition == cd)
  if (nrow(d) == 0) next
  pp <- d %>%
    ggplot(aes(x = x_slide_mm, y = y_slide_mm, color = niche)) +
    geom_point(size = 0.05, alpha = 0.4) +
    scale_color_brewer(palette = "Set2") +
    labs(title = cd) +
    theme_void() + theme(legend.position = "none", aspect.ratio = 1)
  p4_list[[cd]] <- pp
}
if (length(p4_list) > 0) {
  combined <- wrap_plots(p4_list, ncol = 2)
  ggsave(file.path(PLOT_DIR, "niche_spatial_maps.png"),
         plot = combined, width = 12, height = 9, dpi = 130)
}

cat("\n=== 07_niche_clustering.R complete ===\n")
