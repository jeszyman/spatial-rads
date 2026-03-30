if (!exists("PATHWAY_COLS")) source("00_load_data.R")
library(RANN)

obj <- readRDS(file.path(ANALYSIS_DIR, "objects", "seurat_clustered.rds"))
cat(sprintf("Loaded: %d genes x %d cells\n", nrow(obj), ncol(obj)))

# ---- DEG function (effect-size ranked) ----
run_deg <- function(obj, ident.1, ident.2, group.by = "Condition") {
  Idents(obj) <- group.by
  markers <- FindMarkers(obj, ident.1 = ident.1, ident.2 = ident.2,
                         min.pct = 0.1, logfc.threshold = 0.1)
  markers$gene <- rownames(markers)
  markers$comparison <- paste0(ident.1, "_vs_", ident.2)
  markers %>% arrange(desc(abs(avg_log2FC)))
}

# ---- Bulk DEGs: MBRT/SBRT vs Control, MBRT vs SBRT ----
cat("Running bulk DEGs...\n")
mbrt_tp <- c("MBRT_1h", "MBRT_4h", "MBRT_day2", "MBRT_day6")
deg_mbrt_ctrl <- map_dfr(mbrt_tp, function(tp) {
  sub <- subset(obj, cells = colnames(obj)[obj$Condition %in% c(tp, "Control")])
  run_deg(sub, tp, "Control")
})

sbrt_tp <- c("SBRT_4h", "SBRT_day2", "SBRT_day6")
deg_sbrt_ctrl <- map_dfr(sbrt_tp, function(tp) {
  sub <- subset(obj, cells = colnames(obj)[obj$Condition %in% c(tp, "Control")])
  run_deg(sub, tp, "Control")
})

matched <- list(c("MBRT_4h", "SBRT_4h"), c("MBRT_day2", "SBRT_day2"), c("MBRT_day6", "SBRT_day6"))
deg_mbrt_sbrt <- map_dfr(matched, function(pair) {
  sub <- subset(obj, cells = colnames(obj)[obj$Condition %in% pair])
  run_deg(sub, pair[1], pair[2])
})

all_degs <- bind_rows(
  mbrt_vs_ctrl = deg_mbrt_ctrl,
  sbrt_vs_ctrl = deg_sbrt_ctrl,
  mbrt_vs_sbrt = deg_mbrt_sbrt,
  .id = "comparison_type"
)
write_tsv(all_degs, file.path(DATA_DIR, "deg_all_cells.tsv"))
cat(sprintf("Total bulk DEG rows: %d\n", nrow(all_degs)))

# ---- Volcano: MBRT vs SBRT ----
p_vol <- all_degs %>%
  filter(comparison_type == "mbrt_vs_sbrt") %>%
  ggplot(aes(x = avg_log2FC, y = -log10(p_val), color = abs(avg_log2FC) > 0.5)) +
  geom_point(alpha = 0.5, size = 0.5) +
  facet_wrap(~comparison) +
  scale_color_manual(values = c("grey60", "red3"), guide = "none") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey40") +
  labs(title = "MBRT vs SBRT DEGs (effect-size ranked)", x = "log2FC", y = "-log10(p)",
       caption = "n=1 per condition; p-values for visualization only") +
  theme_bw()
ggsave(file.path(PLOT_DIR, "volcano_mbrt_vs_sbrt.png"), plot = p_vol, width = 12, height = 5, dpi = 150)

# ---- Pathway kinetics ----
kinetics <- obj@meta.data %>%
  as.data.frame() %>%
  filter(model == "flank") %>%
  group_by(timepoint_h, treatment) %>%
  summarise(across(all_of(PATHWAY_COLS), mean, na.rm = TRUE), .groups = "drop") %>%
  pivot_longer(cols = all_of(PATHWAY_COLS), names_to = "pathway", values_to = "score")

p_kin <- kinetics %>%
  filter(treatment %in% c("MBRT", "SBRT")) %>%
  mutate(pathway = PATHWAY_LABELS[pathway]) %>%
  ggplot(aes(x = factor(timepoint_h), y = score, color = treatment, group = treatment)) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  facet_wrap(~pathway, scales = "free_y") +
  scale_color_manual(values = TREATMENT_COLORS) +
  labs(x = "Hours post-RT", y = "Mean pathway score", caption = "n=1 per condition") +
  theme_bw()
ggsave(file.path(PLOT_DIR, "pathway_kinetics.png"), plot = p_kin, width = 10, height = 8, dpi = 150)
write_tsv(kinetics, file.path(DATA_DIR, "pathway_kinetics.tsv"))

# ---- Composition ----
composition <- obj@meta.data %>%
  as.data.frame() %>%
  filter(model == "flank") %>%
  count(Condition, timepoint_h, treatment, cell_type = cell_type_validated) %>%
  group_by(Condition) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

p_comp <- composition %>%
  filter(Condition %in% CONDITION_LEVELS[1:8]) %>%
  mutate(Condition = factor(Condition, levels = CONDITION_LEVELS)) %>%
  ggplot(aes(x = Condition, y = prop, fill = cell_type)) +
  geom_col() +
  labs(x = NULL, y = "Proportion", fill = "Cell type", title = "Cell type composition (flank)") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(PLOT_DIR, "composition_stacked.png"), plot = p_comp, width = 10, height = 6, dpi = 150)
write_tsv(composition, file.path(DATA_DIR, "composition.tsv"))

# ---- Spatial NN immune fraction (MBRT vs SBRT at 4h) ----
cat("Computing spatial nearest neighbors at 4h...\n")
nn_results <- map_dfr(c("MBRT_4h", "SBRT_4h"), function(cond) {
  meta_sub <- obj@meta.data %>% as.data.frame() %>% filter(Condition == cond)
  coords <- as.matrix(meta_sub[, c("x_slide_mm", "y_slide_mm")])
  ct <- meta_sub$cell_type_validated
  nn <- nn2(coords, k = 21)
  nn_idx <- nn$nn.idx[, -1]
  immune_frac <- rowMeans(matrix(ct[nn_idx] != "tumor_epithelial", nrow = nrow(nn_idx)))
  tibble(condition = cond, cell_type = ct, immune_neighbor_frac = immune_frac)
})

p_nn <- nn_results %>%
  ggplot(aes(x = immune_neighbor_frac, fill = condition)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~cell_type, scales = "free_y") +
  scale_fill_manual(values = c("MBRT_4h" = TREATMENT_COLORS["MBRT"], "SBRT_4h" = TREATMENT_COLORS["SBRT"])) +
  labs(x = "Fraction immune neighbors (k=20)", title = "Spatial immune neighborhoods: MBRT vs SBRT at 4h") +
  theme_bw()
ggsave(file.path(PLOT_DIR, "spatial_nn_immune_frac.png"), plot = p_nn, width = 14, height = 10, dpi = 150)

cat("Layer 2 condensed analysis complete.\n")
