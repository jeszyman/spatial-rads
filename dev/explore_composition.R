library(Seurat)
library(tidyverse)

obj <- readRDS("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/objects/seurat_clustered.rds")

# --- Proportions ---
composition <- obj@meta.data %>%
  as.data.frame() %>%
  filter(model == "flank") %>%
  count(Condition, timepoint_h, treatment,
        cell_type = cell_type_validated) %>%
  group_by(Condition) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# --- Stacked bar ---
cond_levels <- c("Control", "MBRT_1h", "MBRT_4h", "SBRT_4h",
                 "MBRT_day2", "SBRT_day2", "MBRT_day6", "SBRT_day6")
composition %>%
  filter(Condition %in% cond_levels) %>%
  mutate(Condition = factor(Condition, levels = cond_levels)) %>%
  ggplot(aes(x = Condition, y = prop, fill = cell_type)) +
  geom_col() +
  labs(x = NULL, y = "Proportion", fill = "Cell type",
       title = "Cell type composition (flank)") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/figures/composition_stacked.pdf",
       width = 10, height = 6)
cat("Stacked bar saved.\n")

# --- Trajectories for key cell types ---
composition %>%
  filter(treatment %in% c("MBRT", "SBRT")) %>%
  ggplot(aes(x = factor(timepoint_h), y = prop, color = treatment, group = treatment)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  facet_wrap(~cell_type, scales = "free_y") +
  labs(x = "Hours post-RT", y = "Proportion", color = "Treatment",
       caption = "n=1 per condition") +
  theme_bw()
ggsave("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/figures/composition_trajectories.pdf",
       width = 14, height = 10)
cat("Composition trajectories saved.\n")

write.csv(composition, "/mnt/gcs/jeszyman/projects/spatial-rads/analysis/tables/composition.csv", row.names = FALSE)
cat("Composition table saved.\n")

# --- Key findings ---
cat("\nTop cell types by proportion across conditions:\n")
composition %>%
  group_by(cell_type) %>%
  summarise(mean_prop = mean(prop), .groups = "drop") %>%
  arrange(desc(mean_prop)) %>%
  head(10) %>%
  print()
