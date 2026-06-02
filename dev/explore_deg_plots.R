library(tidyverse)

all_degs <- read.csv("/mnt/data/projects/spatial-rads/analysis/tables/deg_all_cells.csv")
cat(sprintf("Loaded %d DEG rows\n", nrow(all_degs)))

# --- Volcano: MBRT vs SBRT ---
all_degs %>%
  filter(comparison_type == "mbrt_vs_sbrt") %>%
  ggplot(aes(x = avg_log2FC, y = -log10(p_val), color = abs(avg_log2FC) > 0.5)) +
  geom_point(alpha = 0.5, size = 0.5) +
  facet_wrap(~comparison) +
  scale_color_manual(values = c("grey60", "red3"), guide = "none") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey40") +
  labs(title = "MBRT vs SBRT DEGs (effect size ranked)",
       x = "log2FC", y = "-log10(p)",
       caption = "n=1 per condition; p-values for visualization only") +
  theme_bw()
ggsave("/mnt/data/projects/spatial-rads/analysis/figures/volcano_mbrt_vs_sbrt.pdf",
       width = 12, height = 5)
cat("Volcano plot saved.\n")

# --- Volcano: MBRT vs Control ---
all_degs %>%
  filter(comparison_type == "mbrt_vs_ctrl") %>%
  ggplot(aes(x = avg_log2FC, y = -log10(p_val), color = abs(avg_log2FC) > 0.5)) +
  geom_point(alpha = 0.5, size = 0.5) +
  facet_wrap(~comparison) +
  scale_color_manual(values = c("grey60", "red3"), guide = "none") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey40") +
  labs(title = "MBRT vs Control DEGs",
       x = "log2FC", y = "-log10(p)",
       caption = "n=1 per condition; p-values for visualization only") +
  theme_bw()
ggsave("/mnt/data/projects/spatial-rads/analysis/figures/volcano_mbrt_vs_ctrl.pdf",
       width = 14, height = 5)
cat("MBRT vs Control volcano saved.\n")

# --- Top DEGs heatmap data ---
top_genes <- all_degs %>%
  filter(comparison_type == "mbrt_vs_sbrt") %>%
  group_by(comparison) %>%
  slice_max(abs(avg_log2FC), n = 10) %>%
  pull(gene) %>% unique()
cat(sprintf("Top genes for heatmap: %s\n", paste(top_genes, collapse = ", ")))

# --- Effect size comparison across timepoints ---
all_degs %>%
  filter(comparison_type == "mbrt_vs_sbrt", abs(avg_log2FC) > 0.5) %>%
  count(comparison) %>%
  print()
