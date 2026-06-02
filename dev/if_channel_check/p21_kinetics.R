suppressPackageStartupMessages({
  library(arrow); library(data.table); library(ggplot2); library(dplyr)
})

OUT <- "/home/jeszyman/repos/spatial-rads/dev/if_channel_check"

# ---- Load metadata and counts ----
meta <- as.data.table(read_parquet("/tmp/mutter01_meta.parquet"))
counts_df <- as.data.table(read_parquet("/tmp/mutter01_counts.parquet"))

# ---- Join slide/block to treatment/timepoint via xlsx metadata ----
design <- as.data.table(readxl::read_excel("/home/jeszyman/repos/spatial-rads/data/metadata.xlsx"))
design <- design[, .(slide, block_id, condition, model, treatment, timepoint_h)]

# CosMx metadata has Slide (with 20250529 prefix) and Block (like Block_21).
# Design uses slide = same format and block_id = "Block_21".
meta <- merge(meta, design, by.x = c("Slide", "Block"), by.y = c("slide", "block_id"))
cat("After design merge:", nrow(meta), "cells across",
    length(unique(meta$condition)), "conditions\n")

# ---- Subset to flank-model conditions only (exclude tongue) ----
meta <- meta[is.na(model) | model == "flank"]
cat("Flank-only:", nrow(meta), "cells\n")

# ---- Compute per-cell Cdkn1a log-normalized expression ----
counts_df <- counts_df[cell_id %in% meta$cell_id]
setkeyv(meta, "cell_id"); setkeyv(counts_df, "cell_id")
meta <- meta[cell_id %in% counts_df$cell_id]

# Extract Cdkn1a column + library size
lib <- rowSums(counts_df[, !c("Slide", "fov", "cell_id"), with = FALSE])
cdkn1a <- counts_df$Cdkn1a
stopifnot(nrow(counts_df) == length(lib))
norm_cdkn1a <- log2((cdkn1a / pmax(lib, 1)) * 1e4 + 1)

meta[, Cdkn1a := norm_cdkn1a[match(cell_id, counts_df$cell_id)]]
cat("Cells with Cdkn1a expression > 0:", sum(meta$Cdkn1a > 0),
    "/", nrow(meta), "\n")

# ---- Pick well-defined cell types by radiosensitivity ----
# Use Main_cell_Types column (NOT cellFamily) — tumor cells are "a" in Main
setnames(meta, "ImmuneAtlas_ImmGen_Main_cell_Types", "main_type")

classify <- function(x) {
  out <- rep(NA_character_, length(x))
  out[x == "a"] <- "Tumor cells (a)"
  out[x %in% c("B.cell", "Memory.B", "Plasma", "Plasmablast",
               "GC_centroblasts", "GC_centrocyes", "Spleen.CD19")] <- "B cells"
  out[x %in% c("CD8.T.cell", "Spleen.Naive.CD4", "Spleen.Naive.CD8",
               "Spleen.CD4Act.48hrs", "Spleen.LN.Naive.CD4",
               "Spleen.Treg", "Colon.Treg.Nrplo")] <- "T cells"
  out[x %in% c("Macrophage", "PerC.macrophage", "spleen.red.pulp.macs",
               "small.peritoneal.macs", "Microglia")] <- "Macrophages"
  out[x %in% c("Ly6Clo.blood.monocytes", "Ly6Chi.blood.monocytes",
               "Dendritic")] <- "Monocytes/DC"
  out[x %in% c("Fibroblastic.reticular", "Pericyte")] <- "Stroma"
  out
}
meta[, cell_class := classify(main_type)]
meta <- meta[!is.na(cell_class)]
meta[, cell_class := factor(cell_class,
                            levels = c("B cells", "T cells", "Monocytes/DC",
                                       "Macrophages", "Tumor cells (a)", "Stroma"))]

# Encode treatment neatly. Controls get treatment = "Control".
meta[, treat := fifelse(is.na(treatment), "Control",
                        fifelse(treatment == "NT", "Control", treatment))]
meta[, timepoint_h := as.numeric(timepoint_h)]
# Set control at x = 0 for visualization with a slight offset per treatment arm
meta[condition == "Control", treat := "Control"]

cat("\nCells per cell_class x condition:\n")
print(dcast(meta[, .N, by = .(cell_class, condition)],
            cell_class ~ condition, value.var = "N", fill = 0))

# ---- Summaries: mean and % expressed per (condition, cell_class) ----
summ <- meta[, .(
  n_cells = .N,
  mean_expr = mean(Cdkn1a, na.rm = TRUE),
  pct_expr = 100 * mean(Cdkn1a > 0, na.rm = TRUE)
), by = .(cell_class, treat, timepoint_h, condition)]
setorder(summ, cell_class, treat, timepoint_h)

# ---- Add MBRT 4h peak / valley overlay from manual FOV labels ----
peak_fovs   <- unique(c(224, 209, 201, 186, 172, 225, 229, 215, 175,
                        176, 206, 181, 168, 167))
valley_fovs <- c(199, 213, 227, 188, 204, 218, 178, 165, 170, 231, 189,
                 173, 154)
pv_meta <- meta[condition == "MBRT_4h"]
pv_meta[, pv_label := fifelse(fov %in% peak_fovs, "MBRT peak",
                              fifelse(fov %in% valley_fovs, "MBRT valley",
                                      NA_character_))]
pv_summ <- pv_meta[!is.na(pv_label), .(
  n_cells = .N,
  mean_expr = mean(Cdkn1a, na.rm = TRUE),
  pct_expr = 100 * mean(Cdkn1a > 0, na.rm = TRUE),
  treat = pv_label,
  timepoint_h = 4,
  condition = pv_label
), by = .(cell_class, pv_label)]
pv_summ[, pv_label := NULL]
summ <- rbind(summ, pv_summ, fill = TRUE)
tp_levels <- c("0 Gy", "1h", "4h", "2d", "6d")
tp_map <- c("0" = "0 Gy", "1" = "1h", "4" = "4h", "48" = "2d", "144" = "6d")
summ[, tp := factor(tp_map[as.character(timepoint_h)], levels = tp_levels)]
meta[, tp := factor(tp_map[as.character(timepoint_h)], levels = tp_levels)]
fwrite(summ, file.path(OUT, "cdkn1a_kinetics_summary.tsv"), sep = "\t")
cat("\nSummary:\n"); print(summ)

# ---- Plot 1: mean expression kinetics with per-sample points ----
treat_levels <- c("Control", "MBRT", "SBRT", "MBRT peak", "MBRT valley")
summ[, treat := factor(treat, levels = treat_levels)]
treat_colors <- c("Control" = "gray40", "MBRT" = "#e41a1c",
                  "SBRT" = "#377eb8", "MBRT peak" = "#b30000",
                  "MBRT valley" = "#fc9272")
treat_shapes <- c("Control" = 16, "MBRT" = 17, "SBRT" = 15,
                  "MBRT peak" = 8, "MBRT valley" = 4)
line_series <- summ[treat %in% c("Control", "MBRT", "SBRT")]
pv_points   <- summ[treat %in% c("MBRT peak", "MBRT valley")]

p1 <- ggplot() +
  geom_line(data = line_series,
            aes(x = tp, y = mean_expr, color = treat, group = treat),
            linewidth = 0.6) +
  geom_point(data = line_series,
             aes(x = tp, y = mean_expr, color = treat, shape = treat),
             size = 3) +
  geom_point(data = pv_points,
             aes(x = tp, y = mean_expr, color = treat, shape = treat),
             size = 4, stroke = 1.2) +
  facet_wrap(~ cell_class, scales = "free_y") +
  scale_color_manual(values = treat_colors) +
  scale_shape_manual(values = treat_shapes) +
  labs(x = "Timepoint",
       y = "Cdkn1a (p21) — log-normalized mean",
       title = "p21 / Cdkn1a kinetics per cell class",
       subtitle = "Mutter_01 flank tumors — n=1 per condition") +
  theme_minimal() +
  theme(legend.position = "bottom", strip.text = element_text(size = 10))
ggsave(file.path(OUT, "cdkn1a_kinetics_mean.png"), p1,
       width = 11, height = 7, dpi = 150)
cat("Saved cdkn1a_kinetics_mean.png\n")

# ---- Plot 2: % expressing cells (robust to sparsity) ----
p2 <- ggplot() +
  geom_line(data = line_series,
            aes(x = tp, y = pct_expr, color = treat, group = treat),
            linewidth = 0.6) +
  geom_point(data = line_series,
             aes(x = tp, y = pct_expr, color = treat, shape = treat),
             size = 3) +
  geom_point(data = pv_points,
             aes(x = tp, y = pct_expr, color = treat, shape = treat),
             size = 4, stroke = 1.2) +
  facet_wrap(~ cell_class, scales = "free_y") +
  scale_color_manual(values = treat_colors) +
  scale_shape_manual(values = treat_shapes) +
  labs(x = "Timepoint", y = "% cells with Cdkn1a > 0",
       title = "Fraction of cells expressing p21",
       subtitle = "Mutter_01 flank tumors") +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(file.path(OUT, "cdkn1a_kinetics_pctexpr.png"), p2,
       width = 11, height = 7, dpi = 150)
cat("Saved cdkn1a_kinetics_pctexpr.png\n")

# ---- Plot 3: boxplot of per-cell expression distribution ----
p3 <- ggplot(meta[Cdkn1a > 0], aes(x = tp, y = Cdkn1a, fill = treat)) +
  geom_boxplot(outlier.size = 0.3, position = position_dodge(0.7), width = 0.6) +
  facet_wrap(~ cell_class) +
  scale_fill_manual(values = c("Control" = "gray40",
                               "MBRT" = "#e41a1c",
                               "SBRT" = "#377eb8")) +
  labs(x = "Timepoint (h)", y = "Cdkn1a log-normalized expression (cells with >0)",
       title = "p21 distribution in expressing cells",
       subtitle = "Mutter_01 flank tumors") +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(file.path(OUT, "cdkn1a_kinetics_boxplot.png"), p3,
       width = 12, height = 7, dpi = 150)
cat("Saved cdkn1a_kinetics_boxplot.png\n")
