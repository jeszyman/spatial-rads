#!/usr/bin/env Rscript
# Presentation figures for MBRT peak/valley analysis
# Generates: cdkn1a_evidence_panel, pv_signature_dotplot, pv_celltype_heatmap

library(tidyverse)
library(patchwork)
library(ComplexHeatmap)
library(circlize)

DATA_DIR <- "data"
PLOT_DIR <- "plots"
ZONE_COLORS <- c("peak" = "#E74C3C", "valley" = "#3498DB")

# ============================================================
# Figure 2: Cdkn1a / p21 Spatial Evidence Panel
# ============================================================
cat("=== Figure 2: Cdkn1a evidence panel ===\n")

pv <- read_tsv(file.path(DATA_DIR, "mbrt4h_peak_valley.tsv"), show_col_types = FALSE)
stripe_model <- readRDS(file.path(DATA_DIR, "stripe_model.rds"))

# Panel A: spatial scatter colored by p21 expression
p_a <- ggplot(pv, aes(x = x_slide_mm, y = y_slide_mm, color = p21)) +
  geom_point(size = 0.01, alpha = 0.4) +
  scale_color_viridis_c(option = "inferno", name = "p21\n(Cdkn1a)") +
  coord_fixed() +
  labs(title = "A. p21 expression across tissue",
       x = "x (mm)", y = "y (mm)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", size = 12))

# Add beam center lines
rad <- stripe_model$tilt_deg * pi / 180
for (bc in stripe_model$beam_centers) {
  p_a <- p_a + geom_abline(intercept = bc, slope = -tan(rad),
                            color = "yellow", linewidth = 0.5, alpha = 0.7)
}

# Panel B: p21 profile along corrected y-axis with beam centers
library(zoo)
y_profile <- pv %>%
  mutate(y_bin = round(y_corr / 0.1) * 0.1) %>%
  group_by(y_bin) %>%
  summarise(mean_p21 = mean(p21, na.rm = TRUE), n_cells = n(), .groups = "drop") %>%
  filter(n_cells >= 50) %>%
  arrange(y_bin) %>%
  mutate(p21_smooth = rollmean(mean_p21, k = 5, fill = NA, align = "center"))

p_b <- ggplot(y_profile, aes(x = y_bin)) +
  geom_point(aes(y = mean_p21), alpha = 0.3, size = 1.5) +
  geom_line(aes(y = p21_smooth), color = "#E74C3C", linewidth = 1.2, na.rm = TRUE) +
  geom_vline(xintercept = stripe_model$beam_centers,
             color = "#D4A017", linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = mean(stripe_model$beam_centers), y = max(y_profile$mean_p21, na.rm = TRUE) * 1.02,
           label = sprintf("1.02mm spacing\n(matches known beam geometry)"),
           fontface = "bold", size = 3.5, color = "#8B4513") +
  labs(title = "B. p21 profile along beam axis",
       subtitle = "Tilt-corrected y-axis; dashed lines = beam centers",
       x = "Corrected y position (mm)", y = "Mean p21 (Cdkn1a)") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12))

# Panel C: peak vs valley p21 violin
p_c <- ggplot(pv, aes(x = zone, y = p21, fill = zone)) +
  geom_violin(alpha = 0.7, draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_fill_manual(values = ZONE_COLORS) +
  labs(title = "C. p21 by zone",
       subtitle = sprintf("Peak: %.3f vs Valley: %.3f",
                          mean(pv$p21[pv$zone == "peak"], na.rm = TRUE),
                          mean(pv$p21[pv$zone == "valley"], na.rm = TRUE)),
       x = NULL, y = "p21 (Cdkn1a)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 12))

# Assemble
p_fig2 <- (p_a | (p_b / p_c + plot_layout(heights = c(2, 1)))) +
  plot_layout(widths = c(1.2, 1)) +
  plot_annotation(
    title = "Cdkn1a (p21) recovers MBRT beam geometry from gene expression",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

ggsave(file.path(PLOT_DIR, "cdkn1a_evidence_panel.png"),
       plot = p_fig2, width = 16, height = 10, dpi = 200)
cat("Saved: cdkn1a_evidence_panel.png\n")

# ============================================================
# Figure 4: Peak/Valley Signature Dot Plot
# ============================================================
cat("\n=== Figure 4: Signature dot plot ===\n")

degs <- read_tsv(file.path(DATA_DIR, "pv_degs_bulk.tsv"), show_col_types = FALSE)

# Assign biological categories
gene_categories <- tribble(
  ~gene,     ~category,
  "Cdkn1a",  "DDR / Cell cycle",
  "Gadd45b", "DDR / Cell cycle",
  "Ccl8",    "Myeloid / Immune",
  "Cd163",   "Myeloid / Immune",
  "Mrc1",    "Myeloid / Immune",
  "C1qa",    "Myeloid / Immune",
  "Fcer1g",  "Myeloid / Immune",
  "Aif1",    "Myeloid / Immune",
  "Arg1",    "Myeloid / Immune",
  "Lag3",    "Immune checkpoint",
  "Clu",     "Stress response",
  "Mt1",     "Stress response",
  "Mt2",     "Stress response",
  "Vegfa",   "Angiogenesis",
  "Vegfb",   "Angiogenesis",
  "Ccnd1",   "DDR / Cell cycle",
  "Krt14",   "Structural",
  "Vim",     "Structural"
)

# Get top genes by |FC|
top_genes <- degs %>%
  arrange(desc(abs(avg_log2FC))) %>%
  head(30)

top_genes <- top_genes %>%
  left_join(gene_categories, by = "gene") %>%
  mutate(
    category = replace_na(category, "Other"),
    direction = ifelse(avg_log2FC > 0, "Peak-enriched", "Valley-enriched"),
    gene = fct_reorder(gene, avg_log2FC)
  )

category_colors <- c(
  "DDR / Cell cycle" = "#D4A017",
  "Myeloid / Immune" = "#E65100",
  "Immune checkpoint" = "#2E7D32",
  "Stress response" = "#7B1FA2",
  "Angiogenesis" = "#1565C0",
  "Structural" = "#888888",
  "Other" = "#CCCCCC"
)

p_fig4 <- ggplot(top_genes, aes(x = avg_log2FC, y = gene)) +
  geom_vline(xintercept = 0, color = "grey50", linetype = "dashed") +
  geom_segment(aes(x = 0, xend = avg_log2FC, y = gene, yend = gene, color = direction),
               linewidth = 0.5, alpha = 0.5) +
  geom_point(aes(size = pct.1, fill = category, color = direction),
             shape = 21, stroke = 1.2) +
  scale_color_manual(values = c("Peak-enriched" = "#E74C3C", "Valley-enriched" = "#3498DB"),
                     name = "Direction") +
  scale_fill_manual(values = category_colors, name = "Function") +
  scale_size_continuous(name = "% cells\nexpressing", range = c(2, 8)) +
  labs(title = "Peak vs Valley Signature: Top Differentially Expressed Genes",
       subtitle = "MBRT_4h, 149K cells, Wilcoxon rank-sum (effect size ranking, n=1)",
       x = "log2 Fold Change (peak / valley)",
       y = NULL) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14),
        legend.position = "right")

ggsave(file.path(PLOT_DIR, "pv_signature_dotplot.png"),
       plot = p_fig4, width = 12, height = 10, dpi = 200)
cat("Saved: pv_signature_dotplot.png\n")


# ============================================================
# Figure 5: Cell-Type Response Heatmap
# ============================================================
cat("\n=== Figure 5: Cell-type response heatmap ===\n")

ct_degs <- read_tsv(file.path(DATA_DIR, "pv_degs_by_celltype.tsv"), show_col_types = FALSE)

# Summary stats per cell type
ct_summary <- ct_degs %>%
  group_by(cell_type) %>%
  summarise(
    n_degs = n(),
    max_abs_fc = max(abs(avg_log2FC)),
    .groups = "drop"
  ) %>%
  arrange(desc(max_abs_fc))

# Get top 15 cell types and top 20 genes (by max |FC| across any cell type)
top_cts <- ct_summary %>% head(15) %>% pull(cell_type)

top_ct_genes <- ct_degs %>%
  filter(cell_type %in% top_cts) %>%
  group_by(gene) %>%
  summarise(max_fc = max(abs(avg_log2FC)), .groups = "drop") %>%
  arrange(desc(max_fc)) %>%
  head(20) %>%
  pull(gene)

# Build matrix
hm_data <- ct_degs %>%
  filter(cell_type %in% top_cts, gene %in% top_ct_genes) %>%
  select(cell_type, gene, avg_log2FC) %>%
  pivot_wider(names_from = gene, values_from = avg_log2FC, values_fill = 0)

hm_mat <- as.matrix(hm_data[, -1])
rownames(hm_mat) <- hm_data$cell_type

# Order rows by max |FC|
row_order <- ct_summary %>%
  filter(cell_type %in% top_cts) %>%
  arrange(desc(max_abs_fc)) %>%
  pull(cell_type)
hm_mat <- hm_mat[row_order, ]

# Tier annotations
tier_labels <- ifelse(
  ct_summary$max_abs_fc[match(row_order, ct_summary$cell_type)] >= 0.72,
  "Tier 1\n(strong)",
  ifelse(ct_summary$max_abs_fc[match(row_order, ct_summary$cell_type)] < 0.45,
         "Zone-\nblind", "Tier 2")
)

tier_colors <- c("Tier 1\n(strong)" = "#E74C3C", "Tier 2" = "#FFA726", "Zone-\nblind" = "#90CAF9")

# Row annotations
n_degs_vec <- ct_summary$n_degs[match(row_order, ct_summary$cell_type)]
max_fc_vec <- ct_summary$max_abs_fc[match(row_order, ct_summary$cell_type)]

row_ha <- rowAnnotation(
  `Response\ntier` = tier_labels,
  `DEGs` = anno_barplot(n_degs_vec, width = unit(2, "cm"),
                        gp = gpar(fill = "#FFF3CC", col = "#D4A017")),
  `Max |FC|` = anno_barplot(max_fc_vec, width = unit(2, "cm"),
                            gp = gpar(fill = "#FFCCCC", col = "#E74C3C")),
  col = list(`Response\ntier` = tier_colors),
  annotation_name_gp = gpar(fontsize = 9)
)

# Color scale
col_fun <- colorRamp2(c(-0.8, 0, 0.8), c("#3498DB", "white", "#E74C3C"))

png(file.path(PLOT_DIR, "pv_celltype_heatmap.png"), width = 14, height = 10, units = "in", res = 200)
ht <- Heatmap(
  hm_mat,
  name = "log2FC\n(peak/valley)",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  column_names_rot = 45,
  right_annotation = row_ha,
  column_title = "Peak vs Valley: Cell-Type Specific Response",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  row_title = "Cell type (ordered by max |log2FC|)",
  row_title_gp = gpar(fontsize = 11),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 10),
    labels_gp = gpar(fontsize = 9)
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    val <- hm_mat[i, j]
    if (abs(val) > 0.5) {
      grid.text(sprintf("%.2f", val), x, y, gp = gpar(fontsize = 7, col = "white"))
    }
  }
)
draw(ht)
dev.off()
cat("Saved: pv_celltype_heatmap.png\n")

cat("\nAll presentation figures complete.\n")
