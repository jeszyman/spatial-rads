library(Seurat)
library(tidyverse)

obj <- readRDS("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/objects/seurat_clustered.rds")
pathway_cols <- c("TypeI_interferon", "TypeII_interferon", "STING", "DNA_Damage_Repair")

# --- By cell family ---
kinetics <- obj@meta.data %>%
  as.data.frame() %>%
  filter(model == "flank") %>%
  group_by(timepoint_h, treatment,
           cell_family = ImmuneAtlas_ImmGen_cellFamily_Main_cell_Types) %>%
  summarise(n_cells = n(),
            across(all_of(pathway_cols), mean, na.rm = TRUE),
            .groups = "drop") %>%
  pivot_longer(cols = all_of(pathway_cols), names_to = "pathway", values_to = "score")

kinetics %>%
  filter(treatment %in% c("MBRT", "SBRT", "NT")) %>%
  ggplot(aes(x = factor(timepoint_h), y = score, color = treatment, group = treatment)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  facet_grid(pathway ~ cell_family, scales = "free_y") +
  labs(x = "Hours post-RT", y = "Mean pathway score", caption = "n=1 per condition") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/figures/pathway_kinetics_by_celltype.pdf",
       width = 20, height = 12)
cat("Pathway kinetics by cell type saved.\n")

# --- All cells (simplified) ---
obj@meta.data %>%
  as.data.frame() %>%
  filter(model == "flank", treatment %in% c("MBRT", "SBRT")) %>%
  group_by(timepoint_h, treatment) %>%
  summarise(across(all_of(pathway_cols), mean, na.rm = TRUE), .groups = "drop") %>%
  pivot_longer(cols = all_of(pathway_cols), names_to = "pathway", values_to = "score") %>%
  ggplot(aes(x = factor(timepoint_h), y = score, color = treatment, group = treatment)) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  facet_wrap(~pathway, scales = "free_y") +
  scale_color_manual(values = c("MBRT" = "#E74C3C", "SBRT" = "#3498DB")) +
  labs(x = "Hours post-RT", y = "Mean pathway score", caption = "n=1 per condition") +
  theme_bw()
ggsave("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/figures/pathway_kinetics_all_cells.pdf",
       width = 10, height = 8)
cat("Pathway kinetics all cells saved.\n")

write.csv(kinetics, "/mnt/gcs/jeszyman/projects/spatial-rads/analysis/tables/pathway_kinetics.csv", row.names = FALSE)
cat("Kinetics table saved.\n")
