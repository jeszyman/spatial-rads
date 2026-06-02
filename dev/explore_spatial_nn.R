library(Seurat)
library(RANN)
library(tidyverse)

obj <- readRDS("/mnt/data/projects/spatial-rads/analysis/objects/seurat_clustered.rds")

# --- Nearest-neighbor immune fraction at 4h ---
cat("Computing spatial nearest neighbors at 4h...\n")
nn_results <- map_dfr(c("MBRT_4h", "SBRT_4h"), function(cond) {
  cat(sprintf("  %s...\n", cond))
  meta_sub <- obj@meta.data %>% as.data.frame() %>% filter(Condition == cond)
  coords <- as.matrix(meta_sub[, c("x_slide_mm", "y_slide_mm")])
  ct <- meta_sub$cell_type_validated

  nn <- nn2(coords, k = 21)
  nn_idx <- nn$nn.idx[, -1]

  # Fraction of neighbors that are immune (not tumor_epithelial)
  immune_frac <- rowMeans(matrix(ct[nn_idx] != "tumor_epithelial", nrow = nrow(nn_idx)))

  tibble(condition = cond, cell_type = ct, immune_neighbor_frac = immune_frac)
})

nn_results %>%
  ggplot(aes(x = immune_neighbor_frac, fill = condition)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~cell_type, scales = "free_y") +
  labs(x = "Fraction immune neighbors (k=20)",
       title = "Spatial immune neighborhoods: MBRT vs SBRT at 4h") +
  theme_bw()
ggsave("/mnt/data/projects/spatial-rads/analysis/figures/spatial_nn_immune_frac.pdf",
       width = 14, height = 10)
cat("Spatial NN plot saved.\n")

# Summary stats
cat("\nMean immune neighbor fraction by cell type and condition:\n")
nn_results %>%
  group_by(condition, cell_type) %>%
  summarise(mean_immune_frac = mean(immune_neighbor_frac),
            n = n(), .groups = "drop") %>%
  filter(n > 100) %>%
  pivot_wider(names_from = condition, values_from = c(mean_immune_frac, n)) %>%
  print(n = 30)
