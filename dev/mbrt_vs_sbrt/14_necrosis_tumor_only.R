## 14_necrosis_tumor_only.R
## Recompute the necrosis story restricted to tumor_epithelial cells:
## - % necrosis = fraction of tumor cells in necrosis_zone (per Condition)
## - Connected components computed on tumor-only spatial grid

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

# Recompute Control-calibrated necrosis_calibrated flag (same threshold as before, applied to all cells)
ctrl_d20 <- obj@meta.data$local_density_d20[obj$Condition == "Control"]
ctrl_thresh <- quantile(ctrl_d20, 0.95, na.rm = TRUE)
cat(sprintf("Control 95th-pct local_density_d20 threshold: %.4f mm\n", ctrl_thresh))
obj$necrosis_calibrated <- obj$local_density_d20 > ctrl_thresh

slide_col <- if ("Slide" %in% colnames(obj@meta.data)) "Slide" else "slide_ID_numeric"

# === Tumor-only metrics ===
tumor_meta <- as_tibble(obj@meta.data) %>%
  mutate(cell_id = colnames(obj)) %>%
  filter(cell_type_major == "tumor_epithelial",
         !is.na(necrosis_calibrated))

cat(sprintf("Tumor cells with calibrated flag: %d\n", nrow(tumor_meta)))

# % necrosis among tumor cells per Condition
tumor_necr <- tumor_meta %>%
  group_by(Condition, treatment, timepoint_h) %>%
  summarise(n_tumor = n(),
            n_necr  = sum(necrosis_calibrated),
            pct_tumor_necr = 100 * n_necr / n_tumor,
            .groups = "drop")
write_tsv(tumor_necr, file.path(DATA_DIR, "necrosis_tumor_only.tsv"))
cat("\n=== Tumor-only necrosis fraction ===\n"); print(tumor_necr)

# Connected components on tumor-only necrosis bins (per slide × condition)
bin_size_mm <- 0.1
flood_fill_max <- function(bins_xy) {
  if (nrow(bins_xy) == 0) return(0L)
  gx <- round(bins_xy$x_bin, 4); gy <- round(bins_xy$y_bin, 4)
  keys <- paste(gx, gy, sep = "_")
  visited <- rep(FALSE, length(keys))
  max_size <- 0L
  for (i in seq_along(keys)) {
    if (visited[i]) next
    queue <- i; size <- 0L
    while (length(queue) > 0) {
      cur <- queue[1]; queue <- queue[-1]
      if (visited[cur]) next
      visited[cur] <- TRUE; size <- size + 1L
      x <- gx[cur]; y <- gy[cur]
      for (dx in c(-bin_size_mm, 0, bin_size_mm)) {
        for (dy in c(-bin_size_mm, 0, bin_size_mm)) {
          if (dx == 0 && dy == 0) next
          k <- paste(round(x+dx,4), round(y+dy,4), sep = "_")
          m <- match(k, keys)
          if (!is.na(m) && !visited[m]) queue <- c(queue, m)
        }
      }
    }
    if (size > max_size) max_size <- size
  }
  max_size
}

slides <- unique(obj@meta.data[[slide_col]])
cc_results <- list()
for (s in slides) {
  m_s <- as_tibble(obj@meta.data) %>%
    mutate(cell_id = colnames(obj)) %>%
    filter(.data[[slide_col]] == s,
           cell_type_major == "tumor_epithelial",
           !is.na(necrosis_calibrated)) %>%
    mutate(x_bin = round(x_slide_mm / bin_size_mm) * bin_size_mm,
           y_bin = round(y_slide_mm / bin_size_mm) * bin_size_mm)
  if (nrow(m_s) < 50) next
  for (cd in unique(m_s$Condition)) {
    bins <- m_s %>% filter(Condition == cd) %>%
      group_by(x_bin, y_bin) %>%
      summarise(n_cells = n(), n_necr = sum(necrosis_calibrated), .groups = "drop") %>%
      filter(n_cells >= 3, n_necr / n_cells >= 0.5)
    largest <- if (nrow(bins) == 0) 0 else flood_fill_max(bins %>% select(x_bin, y_bin))
    cc_results[[length(cc_results)+1]] <- tibble(
      Slide = s, Condition = cd,
      necrotic_bins = nrow(bins),
      largest_region_bins = largest,
      largest_region_mm2 = largest * bin_size_mm^2,
      pct_in_largest = if (nrow(bins) == 0) NA else 100 * largest / nrow(bins))
  }
}
cc_df <- bind_rows(cc_results) %>%
  mutate(treatment = case_when(Condition == "Control" ~ "NT",
                               grepl("MBRT", Condition) ~ "MBRT",
                               grepl("SBRT", Condition) ~ "SBRT"),
         timepoint_h = case_when(Condition == "Control" ~ 0,
                                 grepl("1h", Condition) ~ 1,
                                 grepl("4h", Condition) ~ 4,
                                 grepl("day2", Condition) ~ 48,
                                 grepl("day6", Condition) ~ 144))
write_tsv(cc_df, file.path(DATA_DIR, "necrosis_tumor_only_connected_components.tsv"))
cat("\n=== Tumor-only connected components ===\n"); print(cc_df)

# === Build best figure (tumor-only version) ===
TREATMENT_COLORS <- c("MBRT" = "#E74C3C", "SBRT" = "#3498DB", "NT" = "#7F8C8D")
flank_conds <- c("Control","MBRT_4h","MBRT_day2","MBRT_day6",
                 "SBRT_4h","SBRT_day2","SBRT_day6")

make_tumor_panel <- function(cond) {
  m <- tumor_meta %>% filter(Condition == cond) %>%
    mutate(zone_label = ifelse(necrosis_calibrated, "necrosis", "viable"))
  if (nrow(m) == 0) return(ggplot() + theme_void())

  q <- tumor_necr %>% filter(Condition == cond)
  c <- cc_df %>% filter(Condition == cond)
  pct_n  <- if (nrow(q) > 0) round(q$pct_tumor_necr[1], 1) else NA
  largest_mm2 <- if (nrow(c) > 0) round(c$largest_region_mm2[1], 2) else NA
  pct_largest <- if (nrow(c) > 0) round(c$pct_in_largest[1], 0) else NA
  ann <- sprintf("%.1f%% of tumor in necrosis-zone · largest %.2f mm² · %d%% in single cluster",
                 pct_n, largest_mm2, pct_largest)
  ggplot(m, aes(x = x_slide_mm, y = y_slide_mm, color = zone_label)) +
    geom_point(size = 0.25, alpha = 0.6) +
    scale_color_manual(values = c("viable" = "grey85", "necrosis" = "#E74C3C"),
                       guide = "none") +
    labs(title = cond, subtitle = ann) +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(size = 10, color = "grey20"),
          aspect.ratio = 1,
          plot.background = element_rect(fill = "white", color = NA))
}

p_m2 <- make_tumor_panel("MBRT_day2")
p_s2 <- make_tumor_panel("SBRT_day2")
p_m6 <- make_tumor_panel("MBRT_day6")
p_s6 <- make_tumor_panel("SBRT_day6")
p_ctrl <- make_tumor_panel("Control")

# Kinetics
p_pct <- tumor_necr %>%
  filter(Condition %in% flank_conds, Condition != "MBRT_1h") %>%
  ggplot(aes(x = factor(timepoint_h), y = pct_tumor_necr,
             color = treatment, group = treatment)) +
  geom_line(linewidth = 1.2) + geom_point(size = 4) +
  scale_color_manual(values = TREATMENT_COLORS, name = NULL) +
  labs(title = "% of tumor cells in necrosis-zone",
       subtitle = "Tumor-only denominator (excludes stroma/immune)",
       x = "Hours post-RT", y = "% tumor cells in sparse-density regions") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"))

p_frag <- cc_df %>%
  filter(Condition %in% flank_conds, Condition != "MBRT_1h") %>%
  ggplot(aes(x = factor(timepoint_h), y = pct_in_largest,
             color = treatment, group = treatment)) +
  geom_line(linewidth = 1.2) + geom_point(size = 4) +
  scale_color_manual(values = TREATMENT_COLORS, name = NULL) +
  labs(title = "Tumor-necrosis consolidation",
       subtitle = "% of tumor-necrosis bins in single largest cluster",
       x = "Hours post-RT", y = "% in largest cluster") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"))

top <- p_m2 | p_s2
middle <- p_ctrl | p_m6 | p_s6
bottom <- p_pct | p_frag

composite <- top / middle / bottom +
  plot_layout(heights = c(2.0, 1.0, 1.0)) +
  plot_annotation(
    title = "Tumor-cell central necrosis: MBRT vs SBRT (4T1 flank, tumor-only denominator)",
    subtitle = "Calibrated to Control: cells with k=20 NN-distance > 95th percentile of Control. Restricted to cell_type_major == tumor_epithelial.",
    caption = "n=1 per condition; descriptive only.",
    theme = theme(plot.title = element_text(face = "bold", size = 16),
                  plot.subtitle = element_text(size = 11),
                  plot.caption = element_text(size = 9, hjust = 0)))

ggsave(file.path(PLOT_DIR, "necrosis_tumor_only_figure.png"),
       plot = composite, width = 17, height = 18, dpi = 150)

cat("\n=== 14_necrosis_tumor_only.R complete ===\n")
