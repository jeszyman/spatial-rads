## 12_mixing_score.R
## Keren et al. 2018 (Cell) mixing score for tumor-immune spatial architecture.
## MS = tumor-immune adjacencies / immune-immune adjacencies. Higher = mixed; lower = compartmentalized.
## Per slide × condition × timepoint. Compare MBRT vs SBRT.

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(RANN)
  library(patchwork)
})

DATA_DIR <- "data"
PLOT_DIR <- "plots"
OBJECTS_DIR <- "/mnt/data/projects/spatial-rads/analysis/objects/mbrt_vs_sbrt"
tmp_data <- "/mnt/data/projects/spatial-rads/analysis/tmp"
if (!dir.exists(tmp_data)) dir.create(tmp_data, recursive = TRUE)
Sys.setenv(TMPDIR = tmp_data, TMP = tmp_data, TEMP = tmp_data)

obj <- readRDS(file.path(OBJECTS_DIR, "seurat_filtered.rds"))
cat(sprintf("Loaded: %d cells\n", ncol(obj)))

slide_col <- if ("Slide" %in% colnames(obj@meta.data)) "Slide" else "slide_ID_numeric"
slides <- unique(obj@meta.data[[slide_col]])
immune_types <- c("T.cell", "B.cell", "myeloid", "NK")
tumor_type   <- "tumor_epithelial"

# Build k=20 NN per slide and tabulate edges
all_per_cell <- list()
all_per_condition <- list()

for (s in slides) {
  m_s <- as_tibble(obj@meta.data) %>%
    mutate(cell_id = colnames(obj),
           x = obj$x_slide_mm,
           y = obj$y_slide_mm) %>%
    filter(.data[[slide_col]] == s)
  if (nrow(m_s) < 50) next

  coords <- as.matrix(m_s[, c("x","y")])
  k_use <- min(21, nrow(m_s))
  nn <- RANN::nn2(coords, k = k_use)
  nn_idx <- nn$nn.idx[, -1, drop = FALSE]

  # Per-cell: count immune neighbors and tumor neighbors
  ct <- m_s$cell_type_major
  ct_neigh <- matrix(ct[nn_idx], nrow = nrow(nn_idx), ncol = ncol(nn_idx))
  is_immune_neigh <- matrix(ct_neigh %in% immune_types,
                            nrow = nrow(ct_neigh), ncol = ncol(ct_neigh))
  is_tumor_neigh  <- matrix(ct_neigh == tumor_type,
                            nrow = nrow(ct_neigh), ncol = ncol(ct_neigh))
  immune_n <- rowSums(is_immune_neigh)
  tumor_n  <- rowSums(is_tumor_neigh)
  total_n  <- ncol(ct_neigh)

  # Per-cell record
  per_cell <- m_s %>%
    select(cell_id, Condition, treatment, timepoint_h, cell_type_major, x, y, necrosis_zone) %>%
    mutate(immune_neighbor_n = immune_n,
           tumor_neighbor_n  = tumor_n,
           total_neighbors   = total_n,
           immune_frac = immune_n / total_n,
           Slide = s)
  all_per_cell[[as.character(s)]] <- per_cell

  # Per-condition: aggregate edges (treat NN as directed; a cell's neighbors count as outgoing edges)
  per_cond <- per_cell %>%
    group_by(Condition, treatment, timepoint_h) %>%
    summarise(
      n_total_cells = n(),
      n_tumor       = sum(cell_type_major == tumor_type),
      n_immune      = sum(cell_type_major %in% immune_types),
      tumor_immune_edges = sum(immune_neighbor_n[cell_type_major == tumor_type]),
      immune_immune_edges = sum(immune_neighbor_n[cell_type_major %in% immune_types]),
      tumor_tumor_edges  = sum(tumor_neighbor_n[cell_type_major == tumor_type]),
      .groups = "drop") %>%
    mutate(Slide = s)
  all_per_condition[[as.character(s)]] <- per_cond
}

per_cell_df <- bind_rows(all_per_cell)
per_cond_df <- bind_rows(all_per_condition)

# Apply day2/day6 necrosis-zone exclusion to per-cell calc (default behavior)
per_cell_main <- per_cell_df %>%
  filter(!(timepoint_h %in% c(48, 144) & necrosis_zone %in% TRUE))

# Recompute per-condition aggregates from filtered cells
per_cond_main <- per_cell_main %>%
  group_by(Condition, treatment, timepoint_h, Slide) %>%
  summarise(
    n_total_cells = n(),
    n_tumor       = sum(cell_type_major == tumor_type),
    n_immune      = sum(cell_type_major %in% immune_types),
    tumor_immune_edges = sum(immune_neighbor_n[cell_type_major == tumor_type]),
    immune_immune_edges = sum(immune_neighbor_n[cell_type_major %in% immune_types]),
    tumor_tumor_edges  = sum(tumor_neighbor_n[cell_type_major == tumor_type]),
    .groups = "drop") %>%
  mutate(
    mixing_score = tumor_immune_edges / pmax(immune_immune_edges, 1),
    # Also: tumor-immune-edges normalized by tumor count (same as mean immune_neighbor_n per tumor cell)
    immune_per_tumor = tumor_immune_edges / pmax(n_tumor, 1),
    # global immune fraction (random expectation)
    global_immune_frac = n_immune / n_total_cells,
    # Per-tumor expected immune neighbors under random = 20 * global_immune_frac
    expected_immune_per_tumor = 20 * global_immune_frac,
    enrichment_over_random = immune_per_tumor / pmax(expected_immune_per_tumor, 0.001)
  )

write_tsv(per_cond_main, file.path(DATA_DIR, "mixing_score_per_condition.tsv"))
cat("\n=== Mixing Score (Keren 2018) by condition ===\n")
print(per_cond_main %>%
        select(Condition, treatment, timepoint_h, n_tumor, n_immune,
               tumor_immune_edges, immune_immune_edges,
               mixing_score, immune_per_tumor, expected_immune_per_tumor,
               enrichment_over_random) %>%
        arrange(timepoint_h, treatment))

# Per-cell enrichment for spatial visualization
per_cell_main <- per_cell_main %>%
  group_by(Condition) %>%
  mutate(global_immune_frac = mean(cell_type_major %in% immune_types),
         per_cell_enrichment = (immune_frac / pmax(global_immune_frac, 0.01))) %>%
  ungroup()

write_tsv(per_cell_main %>% select(cell_id, Condition, cell_type_major,
                                    immune_frac, per_cell_enrichment),
          file.path(DATA_DIR, "mixing_per_cell.tsv.gz"))

# === The headline mixing figure ===
TREATMENT_COLORS <- c("MBRT" = "#E74C3C", "SBRT" = "#3498DB", "NT" = "#7F8C8D")
flank_conds <- c("Control","MBRT_1h","MBRT_4h","MBRT_day2","MBRT_day6",
                 "SBRT_4h","SBRT_day2","SBRT_day6")

# Panel A: mixing score kinetics
pA <- per_cond_main %>%
  filter(Condition %in% flank_conds, treatment %in% c("MBRT","SBRT","NT")) %>%
  ggplot(aes(x = factor(timepoint_h), y = mixing_score,
             color = treatment, group = treatment)) +
  geom_line(linewidth = 1.0) + geom_point(size = 3) +
  scale_color_manual(values = TREATMENT_COLORS) +
  labs(title = "A. Mixing Score (Keren 2018, Cell)",
       subtitle = "Tumor↔immune edges / Immune↔immune edges. Higher = mixed; lower = compartmentalized.",
       x = "Hours post-RT", y = "Mixing score") +
  theme_bw(base_size = 11)

# Panel B: enrichment over random
pB <- per_cond_main %>%
  filter(Condition %in% flank_conds, treatment %in% c("MBRT","SBRT","NT")) %>%
  ggplot(aes(x = factor(timepoint_h), y = enrichment_over_random,
             color = treatment, group = treatment)) +
  geom_line(linewidth = 1.0) + geom_point(size = 3) +
  geom_hline(yintercept = 1, color = "grey50", linetype = "dashed") +
  scale_color_manual(values = TREATMENT_COLORS) +
  labs(title = "B. Tumor immune-neighbor enrichment over random",
       subtitle = "Observed immune neighbors per tumor / expected under random spatial mixing. >1 = preferentially mixed.",
       x = "Hours post-RT", y = "Enrichment ratio (obs/exp)") +
  theme_bw(base_size = 11)

# Panel C: spatial scatter of tumor cells colored by enrichment, faceted by condition
spatial_df <- per_cell_main %>%
  filter(cell_type_major == tumor_type, Condition %in% flank_conds)

pC <- spatial_df %>%
  ggplot(aes(x = x, y = y, color = pmin(per_cell_enrichment, 5))) +
  geom_point(size = 0.05, alpha = 0.5) +
  scale_color_gradient2(low = "#3498DB", mid = "white", high = "#E74C3C",
                       midpoint = 1, name = "Tumor cell\nimmune-enrich\n(obs/exp)",
                       limits = c(0, 5)) +
  facet_wrap(~Condition, ncol = 4, scales = "free") +
  labs(title = "C. Spatial map: per-tumor-cell immune-neighbor enrichment over random expectation",
       subtitle = "Each tumor cell colored by how many more immune neighbors it has vs random expectation",
       x = NULL, y = NULL) +
  theme_void(base_size = 9) +
  theme(strip.text = element_text(size = 9), aspect.ratio = 1)

# Panel D: per-condition density plot of per-tumor-cell enrichment
pD <- spatial_df %>%
  filter(treatment %in% c("MBRT","SBRT","NT")) %>%
  mutate(timepoint_label = case_when(
    timepoint_h == 0 ~ "Control",
    timepoint_h == 1 ~ "1h",
    timepoint_h == 4 ~ "4h",
    timepoint_h == 48 ~ "day2",
    timepoint_h == 144 ~ "day6")) %>%
  ggplot(aes(x = pmin(per_cell_enrichment, 5), fill = treatment)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 1, color = "grey50", linetype = "dashed") +
  facet_wrap(~factor(timepoint_label, levels=c("Control","1h","4h","day2","day6")), nrow = 1, scales = "free_y") +
  scale_fill_manual(values = TREATMENT_COLORS) +
  labs(title = "D. Distribution of per-tumor-cell mixing enrichment",
       x = "Enrichment over random (obs/exp)", y = "Density") +
  theme_bw(base_size = 10)

composite <- (pA / pB / pD) | pC +
  plot_annotation(
    title = "Tumor-immune spatial mixing: MBRT vs SBRT (4T1 flank)",
    subtitle = "Keren mixing score, enrichment over random, spatial heterogeneity",
    theme = theme(plot.title = element_text(size = 14, face = "bold")))
ggsave(file.path(PLOT_DIR, "mixing_score_composite.png"),
       plot = composite, width = 22, height = 12, dpi = 130)

# Also save individual panel A (the headline)
ggsave(file.path(PLOT_DIR, "mixing_score_kinetics.png"),
       plot = pA / pB, width = 10, height = 8, dpi = 150)

# Spatial maps standalone
ggsave(file.path(PLOT_DIR, "mixing_score_spatial.png"),
       plot = pC, width = 14, height = 9, dpi = 130)

cat("\n=== 12_mixing_score.R complete ===\n")
