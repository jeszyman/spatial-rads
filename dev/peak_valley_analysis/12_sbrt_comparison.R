if (!exists("PATHWAY_COLS")) source("00_load_data.R")

obj <- readRDS(file.path(ANALYSIS_DIR, "objects", "seurat_clustered.rds"))
pv <- read_tsv(file.path(DATA_DIR, "mbrt4h_peak_valley.tsv"), show_col_types = FALSE)
bulk_degs <- read_tsv(file.path(DATA_DIR, "pv_degs_bulk.tsv"), show_col_types = FALSE)

# ---- Get SBRT_4h cells ----
sbrt4h <- obj@meta.data %>%
  as.data.frame() %>%
  filter(Condition == "SBRT_4h") %>%
  select(cell_id, all_of(PATHWAY_COLS), cell_type_validated)

# Get expression for top peak/valley genes in SBRT
top_peak_genes <- bulk_degs %>% filter(avg_log2FC > 0) %>% head(20) %>% pull(gene)
top_valley_genes <- bulk_degs %>% filter(avg_log2FC < 0) %>% head(20) %>% pull(gene)
sig_genes <- unique(c(top_peak_genes, top_valley_genes))
sig_genes <- sig_genes[sig_genes %in% rownames(obj)]

# Extract expression for MBRT peak, MBRT valley, and SBRT
norm_data <- GetAssayData(obj, layer = "data")
peak_cells <- pv %>% filter(zone == "peak") %>% pull(cell_id)
valley_cells <- pv %>% filter(zone == "valley") %>% pull(cell_id)
sbrt_cells <- sbrt4h$cell_id

if (length(sig_genes) > 0) {
  expr_matrix <- as.matrix(norm_data[sig_genes, c(peak_cells, valley_cells, sbrt_cells)])

  # Mean expression per group
  group_means <- tibble(
    gene = sig_genes,
    peak = rowMeans(expr_matrix[, peak_cells, drop = FALSE]),
    valley = rowMeans(expr_matrix[, valley_cells, drop = FALSE]),
    sbrt = rowMeans(expr_matrix[, sbrt_cells, drop = FALSE])
  )
  write_tsv(group_means, file.path(DATA_DIR, "pv_sbrt_gene_comparison.tsv"))

  # Heatmap-style comparison
  group_long <- group_means %>%
    pivot_longer(cols = c(peak, valley, sbrt), names_to = "group", values_to = "mean_expr") %>%
    mutate(group = factor(group, levels = c("peak", "valley", "sbrt")))

  p_heat <- ggplot(group_long, aes(x = group, y = gene, fill = mean_expr)) +
    geom_tile() +
    scale_fill_viridis_c(option = "inferno") +
    labs(title = "Peak/Valley signature genes in SBRT_4h",
         subtitle = "Top 20 peak-up and valley-up genes from MBRT_4h",
         fill = "Mean expr", x = NULL) +
    theme_bw(base_size = 12) +
    theme(axis.text.y = element_text(size = 8))
  ggsave(file.path(PLOT_DIR, "pv_sbrt_heatmap.png"), plot = p_heat, width = 8, height = 10, dpi = 150)
}

# ---- Pathway score comparison: MBRT peak vs MBRT valley vs SBRT ----
sbrt_pathways <- sbrt4h %>%
  pivot_longer(cols = all_of(PATHWAY_COLS), names_to = "pathway", values_to = "score") %>%
  mutate(group = "SBRT_4h")

mbrt_pathways <- pv %>%
  pivot_longer(cols = all_of(PATHWAY_COLS), names_to = "pathway", values_to = "score") %>%
  mutate(group = paste0("MBRT_", zone))

combined <- bind_rows(mbrt_pathways %>% select(group, pathway, score),
                      sbrt_pathways %>% select(group, pathway, score)) %>%
  mutate(pathway = PATHWAY_LABELS[pathway],
         group = factor(group, levels = c("MBRT_peak", "MBRT_valley", "SBRT_4h")))

p_compare <- ggplot(combined, aes(x = group, y = score, fill = group)) +
  geom_boxplot(outlier.size = 0.1, outlier.alpha = 0.1) +
  scale_fill_manual(values = c("MBRT_peak" = "#E74C3C", "MBRT_valley" = "#3498DB", "SBRT_4h" = "#9B59B6")) +
  facet_wrap(~pathway, scales = "free_y") +
  labs(title = "Pathway scores: MBRT peak vs MBRT valley vs SBRT_4h",
       y = "Score", x = NULL) +
  theme_bw() + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")
ggsave(file.path(PLOT_DIR, "pv_sbrt_pathway_comparison.png"), plot = p_compare, width = 12, height = 8, dpi = 150)

cat("SBRT comparison complete.\n")
