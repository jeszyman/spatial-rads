if (!exists("PATHWAY_COLS")) source("00_load_data.R")
library(RANN)

pv <- read_tsv(file.path(DATA_DIR, "mbrt4h_peak_valley.tsv"), show_col_types = FALSE)

# ---- k=20 NN within each zone ----
cat("Computing spatial neighborhoods by zone...\n")
nn_by_zone <- map_dfr(c("peak", "valley"), function(z) {
  zone_cells <- pv %>% filter(zone == z)
  coords <- as.matrix(zone_cells[, c("x_slide_mm", "y_slide_mm")])
  ct <- zone_cells$cell_type_validated

  nn <- nn2(coords, k = min(21, nrow(coords)))
  nn_idx <- nn$nn.idx[, -1]

  immune_frac <- rowMeans(matrix(ct[nn_idx] != "tumor_epithelial", nrow = nrow(nn_idx)))

  tibble(zone = z, cell_type = ct, immune_neighbor_frac = immune_frac)
})

# Density plot
p_nn <- ggplot(nn_by_zone, aes(x = immune_neighbor_frac, fill = zone)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = ZONE_COLORS) +
  facet_wrap(~cell_type, scales = "free_y") +
  labs(x = "Fraction immune neighbors (k=20)",
       title = "Spatial immune neighborhoods: Peak vs Valley (MBRT_4h)") +
  theme_bw()
ggsave(file.path(PLOT_DIR, "pv_spatial_nn.png"), plot = p_nn, width = 14, height = 10, dpi = 150)

# Summary stats
nn_summary <- nn_by_zone %>%
  group_by(zone, cell_type) %>%
  summarise(mean_immune_frac = mean(immune_neighbor_frac), n = n(), .groups = "drop") %>%
  filter(n > 100) %>%
  pivot_wider(names_from = zone, values_from = c(mean_immune_frac, n))
cat("=== Mean immune neighbor fraction by zone ===\n")
print(nn_summary)
write_tsv(nn_summary, file.path(DATA_DIR, "pv_spatial_nn_summary.tsv"))

cat("Spatial neighborhood analysis complete.\n")
