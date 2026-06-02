library(arrow)
library(tidyverse)
library(cowplot)

INPUT_DIR <- "/mnt/data/projects/spatial-rads/inputs"
PLOT_DIR  <- "plots"

meta <- read_parquet(
  file.path(INPUT_DIR, "Analysis_Mutter_01_CosMmR_Mutter_updated_metadata.parquet"),
  col_select = c("cell_id", "Condition", "x_slide_mm", "y_slide_mm",
                  "qcFlagsCell", "ImmuneAtlas_ImmGen_Main_cell_Types")
) %>% as.data.frame() %>% filter(qcFlagsCell == "Pass")

flank_conds <- c("Control", "MBRT_1h", "MBRT_4h", "SBRT_4h",
                  "MBRT_day2", "SBRT_day2", "MBRT_day6", "SBRT_day6")
meta <- meta %>% filter(Condition %in% flank_conds)

tumor_types <- c("a", "Stem.Prog")
immune_lymphoid <- c("B.cell", "Memory.B", "CD8.T.cell", "NKT", "gdT", "ILC",
                      "Thymic.preT.DN1", "Thymic.4.8.CD3lo.DP", "Thymic.CD4SP",
                      "Thymic.CD8SP", "Colon.Treg.Nrplo", "Plasmablast")
immune_myeloid <- c("Macrophage", "PerC.macrophage", "spleen.red.pulp.macs",
                     "Dendritic", "Microglia", "Ly6Clo.blood.monocytes",
                     "Ly6Chi.blood.monocytes", "BM.Neutrophil",
                     "Thio.induced.Peritoneal.Neutrophil")
stromal <- c("Pericyte", "Blood.endothelial", "Lymphatic.endothelial",
             "Fibroblastic.reticular")

expr <- read_parquet(
  file.path(INPUT_DIR, "projects_Mutter_01_CosMmR_app_data_raw_tx_counts_matrix.parquet"),
  col_select = c("cell_id", "Krt8", "Krt18", "Epcam", "Cdh1", "Pecam1", "Cdh5", "Vwf")
) %>% as.data.frame()

meta <- left_join(meta, expr, by = "cell_id")
meta$epithelial_score <- rowMeans(cbind(
  replace_na(meta$Krt8, 0), replace_na(meta$Krt18, 0),
  replace_na(meta$Epcam, 0), replace_na(meta$Cdh1, 0)), na.rm = TRUE)
meta$endothelial_score <- rowMeans(cbind(
  replace_na(meta$Pecam1, 0), replace_na(meta$Cdh5, 0),
  replace_na(meta$Vwf, 0)), na.rm = TRUE)
meta$is_4t1 <- meta$ImmuneAtlas_ImmGen_Main_cell_Types %in% tumor_types &
               meta$epithelial_score > 0.3 & meta$endothelial_score < 0.1

meta$cell_class <- case_when(
  meta$is_4t1 ~ "4T1 Tumor",
  meta$ImmuneAtlas_ImmGen_Main_cell_Types %in% tumor_types ~ "Non-tumor 'a'",
  meta$ImmuneAtlas_ImmGen_Main_cell_Types %in% immune_lymphoid ~ "Lymphoid",
  meta$ImmuneAtlas_ImmGen_Main_cell_Types %in% immune_myeloid ~ "Myeloid",
  meta$ImmuneAtlas_ImmGen_Main_cell_Types %in% stromal ~ "Stromal",
  TRUE ~ "Other"
)

meta$cell_class <- factor(meta$cell_class,
                          levels = c("4T1 Tumor", "Non-tumor 'a'", "Myeloid",
                                     "Lymphoid", "Stromal", "Other"))

meta$Condition <- factor(meta$Condition, levels = flank_conds)

class_colors <- c(
  "4T1 Tumor"      = "#E74C3C",
  "Non-tumor 'a'"  = "#F5B7B1",
  "Myeloid"        = "#3498DB",
  "Lymphoid"       = "#2ECC71",
  "Stromal"        = "#F39C12",
  "Other"          = "#BDC3C7"
)

cat(sprintf("4T1 tumor cells: %d / %d total 'a' bucket (%.0f%%)\n",
            sum(meta$is_4t1), sum(meta$ImmuneAtlas_ImmGen_Main_cell_Types %in% tumor_types),
            100 * mean(meta$is_4t1[meta$ImmuneAtlas_ImmGen_Main_cell_Types %in% tumor_types])))

plot_one <- function(df, cond_label) {
  ggplot(df, aes(x = x_slide_mm, y = y_slide_mm, color = cell_class)) +
    geom_point(size = 0.05, alpha = 0.4) +
    scale_color_manual(values = class_colors, name = "Cell class", drop = FALSE) +
    coord_fixed() +
    labs(title = cond_label) +
    theme_void(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      legend.position = "none"
    )
}

plots <- lapply(levels(meta$Condition), function(cond) {
  plot_one(meta %>% filter(Condition == cond) %>% arrange(desc(cell_class)), cond)
})

legend_plot <- cowplot::get_legend(
  ggplot(meta[1:10, ], aes(x = x_slide_mm, y = y_slide_mm, color = cell_class)) +
    geom_point(size = 4) +
    scale_color_manual(values = class_colors, name = "Cell class", drop = FALSE) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 13, face = "bold")) +
    guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)))
)

p <- cowplot::plot_grid(
  cowplot::plot_grid(plotlist = plots, ncol = 4),
  legend_plot, ncol = 1, rel_heights = c(1, 0.08)
)

ggsave(file.path(PLOT_DIR, "spatial_cellclass_flank.png"), plot = p,
       width = 16, height = 8, dpi = 200)
cat("Saved: plots/spatial_cellclass_flank.png\n")
