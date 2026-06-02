if (!exists("PATHWAY_COLS")) source("00_load_data.R")

pv <- read_tsv(file.path(DATA_DIR, "mbrt4h_peak_valley.tsv"), show_col_types = FALSE)

# ---- Overall pathway comparison ----
pathway_by_zone <- pv %>%
  pivot_longer(cols = all_of(PATHWAY_COLS), names_to = "pathway", values_to = "score") %>%
  mutate(pathway = PATHWAY_LABELS[pathway])

p_box <- ggplot(pathway_by_zone, aes(x = pathway, y = score, fill = zone)) +
  geom_boxplot(outlier.size = 0.1, outlier.alpha = 0.1) +
  scale_fill_manual(values = ZONE_COLORS) +
  labs(title = "Pathway scores: Peak vs Valley (MBRT_4h)",
       subtitle = "NanoString Mouse UCC standard modules (DDR/IFN-I/IFN-II) + custom STING",
       y = "Score", x = NULL) +
  theme_bw() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(PLOT_DIR, "pv_pathway_boxplots.png"), plot = p_box, width = 10, height = 6, dpi = 150)

# ---- Summary stats ----
pathway_stats <- pathway_by_zone %>%
  group_by(pathway, zone) %>%
  summarise(n = n(), mean = mean(score, na.rm = TRUE), sd = sd(score, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = zone, values_from = c(n, mean, sd))
cat("=== Pathway scores by zone ===\n")
print(pathway_stats)
write_tsv(pathway_stats, file.path(DATA_DIR, "pv_pathway_stats.tsv"))

# ---- Stratified by cell type family ----
obj <- readRDS(file.path(ANALYSIS_DIR, "objects", "seurat_clustered.rds"))
mbrt4h_meta <- obj@meta.data %>%
  as.data.frame() %>%
  filter(Condition == "MBRT_4h") %>%
  select(cell_id, cell_family = ImmuneAtlas_ImmGen_cellFamily_Main_cell_Types)

pv_with_family <- pv %>%
  left_join(mbrt4h_meta, by = "cell_id")

major_families <- pv_with_family %>%
  count(cell_family) %>%
  filter(n >= 500) %>%
  pull(cell_family)

p_strat <- pv_with_family %>%
  filter(cell_family %in% major_families) %>%
  pivot_longer(cols = all_of(PATHWAY_COLS), names_to = "pathway", values_to = "score") %>%
  mutate(pathway = PATHWAY_LABELS[pathway]) %>%
  ggplot(aes(x = zone, y = score, fill = zone)) +
  geom_boxplot(outlier.size = 0.1, outlier.alpha = 0.1) +
  scale_fill_manual(values = ZONE_COLORS) +
  facet_grid(pathway ~ cell_family, scales = "free_y") +
  labs(title = "Pathway scores by zone and cell family", y = "Score") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(PLOT_DIR, "pv_pathway_by_celltype.png"), plot = p_strat, width = 16, height = 12, dpi = 150)

cat("Pathway comparison complete.\n")
