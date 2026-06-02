## 11_m1_m2_polarization.R
## Address whether MBRT day2 myeloid is more M2-skewed than SBRT day2 myeloid,
## or simply has more total myeloid with similar polarization mix.

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(UCell)
})

DATA_DIR <- "data"
PLOT_DIR <- "plots"
OBJECTS_DIR <- "/mnt/data/projects/spatial-rads/analysis/objects/mbrt_vs_sbrt"
tmp_data <- "/mnt/data/projects/spatial-rads/analysis/tmp"
if (!dir.exists(tmp_data)) dir.create(tmp_data, recursive = TRUE)
Sys.setenv(TMPDIR = tmp_data, TMP = tmp_data, TEMP = tmp_data)

obj <- readRDS(file.path(OBJECTS_DIR, "seurat_filtered.rds"))
cat(sprintf("Loaded: %d cells\n", ncol(obj)))

# Restrict to myeloid cells only
mye <- subset(obj, subset = cell_type_major == "myeloid")
cat(sprintf("Myeloid cells: %d\n", ncol(mye)))

# Canonical M1 / M2 marker sets
m1_markers <- c("Nos2","Tnf","Il6","Il1b","Il12b","Cxcl9","Cxcl10","Cxcl11","Cd86",
                "Cd80","Tlr2","Tlr4","Stat1","Irf5","Hif1a","Nfkb1")
m2_markers <- c("Cd163","Mrc1","Arg1","Il4ra","Mgl2","Retnla","Chi3l3","Ym1",
                "Stab1","Vegfa","Ccl22","Mmp9","Tgfb1","Stat6","Klf4","Cd206")

panel <- rownames(mye)
m1_in <- intersect(m1_markers, panel)
m2_in <- intersect(m2_markers, panel)
cat(sprintf("M1 markers in panel: %d/%d (%s)\n", length(m1_in), length(m1_markers),
            paste(m1_in, collapse=",")))
cat(sprintf("M2 markers in panel: %d/%d (%s)\n", length(m2_in), length(m2_markers),
            paste(m2_in, collapse=",")))

# Score M1 / M2 signatures via UCell
mye <- AddModuleScore_UCell(mye,
                            features = list(M1 = m1_in, M2 = m2_in),
                            ncores = 4,
                            name = "_score")

# Also fetch raw normalized expression of every M1 + M2 marker for plotting
data_layer <- GetAssayData(mye, layer = "data")
all_markers <- c(m1_in, m2_in)
expr_long <- as_tibble(t(as.matrix(data_layer[all_markers, , drop = FALSE]))) %>%
  mutate(cell_id = colnames(mye),
         Condition = mye$Condition,
         treatment = mye$treatment,
         timepoint_h = mye$timepoint_h,
         M1_score = mye$M1_score,
         M2_score = mye$M2_score) %>%
  pivot_longer(cols = all_of(all_markers), names_to = "gene", values_to = "expr") %>%
  mutate(class = ifelse(gene %in% m1_in, "M1", "M2"))

# Per (Condition × gene) mean expression
gene_means <- expr_long %>%
  group_by(Condition, treatment, timepoint_h, class, gene) %>%
  summarise(mean_expr = mean(expr, na.rm = TRUE),
            pct_pos   = mean(expr > 0, na.rm = TRUE),
            n_cells   = n(),
            .groups = "drop")
write_tsv(gene_means, file.path(DATA_DIR, "myeloid_m1m2_gene_means.tsv"))

# Per-cell M1/M2 score summaries
score_summary <- as_tibble(mye@meta.data) %>%
  filter(!is.na(M1_score), !is.na(M2_score)) %>%
  group_by(Condition, treatment, timepoint_h) %>%
  summarise(M1_mean = mean(M1_score),
            M2_mean = mean(M2_score),
            M2_M1_ratio = mean(M2_score) / mean(M1_score),
            n_cells = n(),
            .groups = "drop")
write_tsv(score_summary, file.path(DATA_DIR, "myeloid_m1m2_scores.tsv"))
cat("\nMyeloid M1/M2 score summary:\n"); print(score_summary)

# === The headline figure (4 panels) ===
TREATMENT_COLORS <- c("MBRT" = "#E74C3C", "SBRT" = "#3498DB", "NT" = "#7F8C8D")
flank_conds <- c("Control","MBRT_1h","MBRT_4h","MBRT_day2","MBRT_day6",
                 "SBRT_4h","SBRT_day2","SBRT_day6")

# Panel A: M1 score kinetics
pA <- score_summary %>%
  filter(Condition %in% flank_conds, treatment %in% c("MBRT","SBRT","NT")) %>%
  ggplot(aes(x = factor(timepoint_h), y = M1_mean,
             color = treatment, group = treatment)) +
  geom_line(linewidth = 1.0) + geom_point(size = 3) +
  scale_color_manual(values = TREATMENT_COLORS) +
  labs(title = "A. Myeloid M1 score (pro-inflammatory)",
       x = "Hours post-RT", y = "Mean M1 UCell score",
       caption = sprintf("M1 markers (panel): %s", paste(m1_in, collapse=", "))) +
  theme_bw(base_size = 11)

# Panel B: M2 score kinetics
pB <- score_summary %>%
  filter(Condition %in% flank_conds, treatment %in% c("MBRT","SBRT","NT")) %>%
  ggplot(aes(x = factor(timepoint_h), y = M2_mean,
             color = treatment, group = treatment)) +
  geom_line(linewidth = 1.0) + geom_point(size = 3) +
  scale_color_manual(values = TREATMENT_COLORS) +
  labs(title = "B. Myeloid M2 score (anti-inflammatory / tissue-repair)",
       x = "Hours post-RT", y = "Mean M2 UCell score",
       caption = sprintf("M2 markers (panel): %s", paste(m2_in, collapse=", "))) +
  theme_bw(base_size = 11)

# Panel C: M2/M1 ratio kinetics — the key polarization metric
pC <- score_summary %>%
  filter(Condition %in% flank_conds, treatment %in% c("MBRT","SBRT","NT")) %>%
  ggplot(aes(x = factor(timepoint_h), y = M2_M1_ratio,
             color = treatment, group = treatment)) +
  geom_line(linewidth = 1.0) + geom_point(size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = TREATMENT_COLORS) +
  labs(title = "C. M2/M1 polarization ratio",
       subtitle = "Higher = more M2-skewed (immunosuppressive). Lower = more M1-skewed (pro-inflammatory)",
       x = "Hours post-RT", y = "M2 mean / M1 mean") +
  theme_bw(base_size = 11)

# Panel D: per-gene heatmap of myeloid expression at day2 (MBRT vs SBRT)
day2_gene_df <- gene_means %>%
  filter(timepoint_h == 48, Condition %in% c("MBRT_day2","SBRT_day2","Control")) %>%
  mutate(gene = factor(gene, levels = c(m1_in, m2_in)))

pD <- day2_gene_df %>%
  ggplot(aes(x = Condition, y = gene, fill = mean_expr)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "magma", name = "Mean\nlog-norm\nexpr") +
  facet_grid(class ~ ., scales = "free_y", space = "free_y") +
  labs(title = "D. Myeloid expression of canonical M1/M2 markers at day2",
       subtitle = "Per-gene mean within myeloid cells. MBRT day2 myeloid: clear M2 dominance over M1.",
       x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1))

composite <- (pA | pB) / (pC | pD) +
  plot_annotation(
    title = "Myeloid M1/M2 polarization: MBRT vs SBRT",
    subtitle = "Day 2 MBRT myeloid is M2-skewed. Tissue-repair / immunosuppressive — flag for interpretation.",
    caption = "n=1 per condition; descriptive only.",
    theme = theme(plot.title = element_text(size = 14, face = "bold"),
                  plot.subtitle = element_text(size = 11)))
ggsave(file.path(PLOT_DIR, "myeloid_m1m2_polarization.png"),
       plot = composite, width = 14, height = 11, dpi = 130)

# Quick text summary for org section
cat("\n\n=== KEY M1/M2 NUMBERS (day2 + day6) ===\n")
key <- score_summary %>%
  filter(Condition %in% c("MBRT_day2","SBRT_day2","MBRT_day6","SBRT_day6","Control")) %>%
  arrange(timepoint_h, treatment) %>%
  select(Condition, M1_mean, M2_mean, M2_M1_ratio, n_cells)
print(key)

cat("\n=== 11_m1_m2_polarization.R complete ===\n")
