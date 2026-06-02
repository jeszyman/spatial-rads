## 08b_set2_validation_fast.R
## Fast Mutter_02 day2 concordance: reuses existing typed RDS, subsamples
## per stratum to 3000 cells for DEG, runs Wilcoxon. Avoids the slow
## FindMarkers-on-millions-of-cells bottleneck of script 08.

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
})

DATA_DIR <- "data"
PLOT_DIR <- "plots"
OBJECTS_DIR <- "/mnt/data/projects/spatial-rads/analysis/objects/mbrt_vs_sbrt"

# Set TMPDIR to data disk (rule from feedback_no_local_tmp.md)
tmp_data <- "/mnt/data/projects/spatial-rads/analysis/tmp"
if (!dir.exists(tmp_data)) dir.create(tmp_data, recursive = TRUE)
Sys.setenv(TMPDIR = tmp_data, TMP = tmp_data, TEMP = tmp_data)

# Load already-typed Mutter_02
typed_rds <- file.path(OBJECTS_DIR, "seurat_mutter02_typed.rds")
stopifnot(file.exists(typed_rds))
obj <- readRDS(typed_rds)
cat(sprintf("Loaded typed Mutter_02: %d genes x %d cells\n", nrow(obj), ncol(obj)))
cat("Cell type counts:\n"); print(table(obj$cell_type_major))

# Harmonize treatment + Condition labels for DEG group identifiers
obj$treatment_clean <- case_when(
  toupper(obj$treatment) %in% c("MBRT") ~ "MBRT",
  toupper(obj$treatment) %in% c("SBRT") ~ "SBRT",
  toupper(obj$treatment) %in% c("CONTROL","CTRL","NT") ~ "NT",
  TRUE ~ as.character(obj$treatment))
obj$Condition_m02 <- ifelse(obj$treatment_clean == "NT",
                            "Control",
                            paste0(obj$treatment_clean, "_day2"))
cat("Condition_m02 distribution:\n"); print(table(obj$Condition_m02))

# Apply default necrosis exclusion at day2
keep_main <- !(obj$timepoint_h == 48 & obj$necrosis_zone %in% TRUE)
obj <- subset(obj, cells = colnames(obj)[keep_main])
cat(sprintf("After day2 necrosis-zone exclusion: %d cells\n", ncol(obj)))

# Subsample per (Condition_m02 x cell_type_major) to max 3000 cells
set.seed(42)
MAJOR_TYPES <- c("tumor_epithelial", "endothelial", "fibroblast", "smooth_muscle",
                 "T.cell", "B.cell", "myeloid", "NK", "plasma")

meta <- as_tibble(obj@meta.data) %>%
  mutate(cell_id = colnames(obj)) %>%
  filter(cell_type_major %in% MAJOR_TYPES, !is.na(Condition_m02))

cells_to_keep <- meta %>%
  group_by(Condition_m02, cell_type_major) %>%
  slice_sample(n = 3000) %>%
  ungroup() %>%
  pull(cell_id)
cat(sprintf("After subsampling to 3000 per stratum: %d cells\n", length(cells_to_keep)))

obj <- subset(obj, cells = cells_to_keep)

# DEGs at day2: 3 contrasts
day2_contrasts <- list(
  list(name = "MBRT_vs_Control", ident.1 = "MBRT_day2", ident.2 = "Control"),
  list(name = "SBRT_vs_Control", ident.1 = "SBRT_day2", ident.2 = "Control"),
  list(name = "MBRT_vs_SBRT",    ident.1 = "MBRT_day2", ident.2 = "SBRT_day2")
)

run_deg <- function(obj_sub, ident.1, ident.2, group.by) {
  Idents(obj_sub) <- group.by
  out <- tryCatch(
    FindMarkers(obj_sub, ident.1 = ident.1, ident.2 = ident.2,
                min.pct = 0.1, logfc.threshold = 0.1, verbose = FALSE),
    error = function(e) {
      cat(sprintf("  FindMarkers error: %s\n", conditionMessage(e)))
      NULL
    })
  if (is.null(out) || nrow(out) == 0) return(NULL)
  out$gene <- rownames(out)
  out
}

m02_strata <- as_tibble(obj@meta.data) %>% count(Condition_m02, cell_type_major)
deg_jobs <- list()
for (ce in day2_contrasts) {
  for (ct in MAJOR_TYPES) {
    n1 <- m02_strata %>% filter(Condition_m02 == ce$ident.1, cell_type_major == ct) %>% pull(n)
    n2 <- m02_strata %>% filter(Condition_m02 == ce$ident.2, cell_type_major == ct) %>% pull(n)
    if (length(n1) == 0) n1 <- 0; if (length(n2) == 0) n2 <- 0
    if (n1 >= 100 && n2 >= 100) {
      deg_jobs[[length(deg_jobs) + 1]] <- c(ce, list(cell_type = ct, n1 = n1, n2 = n2))
    }
  }
  deg_jobs[[length(deg_jobs) + 1]] <- c(ce, list(cell_type = "BULK", n1 = NA, n2 = NA))
}
cat(sprintf("Total DEG jobs: %d\n", length(deg_jobs)))

m02_degs <- list()
t0 <- Sys.time()
for (i in seq_along(deg_jobs)) {
  job <- deg_jobs[[i]]
  ct <- job$cell_type
  cells_keep <- if (ct == "BULK") {
    colnames(obj)[obj$Condition_m02 %in% c(job$ident.1, job$ident.2)]
  } else {
    colnames(obj)[obj$Condition_m02 %in% c(job$ident.1, job$ident.2) &
                  obj$cell_type_major == ct]
  }
  if (length(cells_keep) < 50) next
  obj_sub <- subset(obj, cells = cells_keep)
  d <- run_deg(obj_sub, job$ident.1, job$ident.2, "Condition_m02")
  if (!is.null(d)) {
    d$contrast <- job$name
    d$cell_type_major <- ct
    d$ident.1 <- job$ident.1; d$ident.2 <- job$ident.2
    m02_degs[[length(m02_degs) + 1]] <- d
  }
  cat(sprintf("  [%d/%d] %s %s: %.0fs elapsed\n", i, length(deg_jobs),
              job$name, ct,
              as.numeric(difftime(Sys.time(), t0, units="secs"))))
}
m02_degs_df <- bind_rows(m02_degs) %>%
  rename(log2FC = avg_log2FC) %>%
  arrange(contrast, cell_type_major, desc(abs(log2FC)))
write_tsv(m02_degs_df, file.path(DATA_DIR, "mutter02_degs_day2.tsv"))
cat(sprintf("\nWrote %d DEG rows\n", nrow(m02_degs_df)))

# Concordance with Mutter_01 day2
m01_degs <- read_tsv(file.path(DATA_DIR, "degs_kinetics.tsv"),
                     show_col_types = FALSE) %>%
  filter(timepoint_h == 48) %>%
  rename(log2FC_m01 = log2FC, p_val_adj_m01 = p_val_adj)

concord <- m01_degs %>%
  inner_join(m02_degs_df %>%
               rename(log2FC_m02 = log2FC, p_val_adj_m02 = p_val_adj),
             by = c("gene", "contrast", "cell_type_major")) %>%
  select(gene, contrast, cell_type_major,
         log2FC_m01, log2FC_m02, p_val_adj_m01, p_val_adj_m02)
write_tsv(concord, file.path(DATA_DIR, "set2_concordance.tsv"))

rho_table <- concord %>%
  group_by(contrast, cell_type_major) %>%
  summarise(rho = if (sum(!is.na(log2FC_m01) & !is.na(log2FC_m02)) >= 10)
              cor(log2FC_m01, log2FC_m02, method = "spearman", use = "pairwise") else NA_real_,
            n_genes = n(),
            .groups = "drop") %>%
  arrange(contrast, desc(rho))
write_tsv(rho_table, file.path(DATA_DIR, "set2_concordance_rho.tsv"))
cat("\n=== Concordance rho table ===\n"); print(rho_table)

# Plot: concordance scatter (MBRT vs SBRT only, per cell type)
TREATMENT_COLORS <- c("MBRT" = "#E74C3C", "SBRT" = "#3498DB", "NT" = "#7F8C8D")

p_conc <- concord %>%
  filter(contrast == "MBRT_vs_SBRT", cell_type_major != "BULK") %>%
  ggplot(aes(x = log2FC_m01, y = log2FC_m02)) +
  geom_point(alpha = 0.4, size = 0.6) +
  geom_abline(slope = 1, intercept = 0, color = "red3", linetype = "dashed") +
  geom_hline(yintercept = 0, color = "grey60") +
  geom_vline(xintercept = 0, color = "grey60") +
  facet_wrap(~cell_type_major, ncol = 3) +
  labs(title = "Mutter_01 day2 vs Mutter_02 day2 — MBRT vs SBRT log2FC",
       subtitle = "Sub-sampled to 3000 cells per stratum",
       x = "Mutter_01 log2FC", y = "Mutter_02 log2FC") +
  theme_bw()
ggsave(file.path(PLOT_DIR, "set2_concordance_scatter.png"),
       plot = p_conc, width = 12, height = 9, dpi = 130)

# Spearman heatmap
p_rho <- rho_table %>%
  ggplot(aes(x = contrast, y = cell_type_major, fill = rho)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("rho=%.2f\nn=%d", rho, n_genes)), size = 3) +
  scale_fill_gradient2(low = "#3498DB", mid = "white", high = "#E74C3C",
                       midpoint = 0, limits = c(-1, 1), na.value = "grey90") +
  labs(title = "Mutter_01 vs Mutter_02 day2 effect-size concordance (Spearman)",
       x = "Contrast", y = NULL) +
  theme_bw() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(PLOT_DIR, "set2_concordance_rho.png"),
       plot = p_rho, width = 8, height = 6, dpi = 150)

# Composition concordance (uses full obj, not subsampled — composition needs full counts)
obj_full <- readRDS(typed_rds)
m02_comp <- as_tibble(obj_full@meta.data) %>%
  filter(timepoint_h == 48,
         toupper(treatment) %in% c("MBRT","SBRT"),
         cell_type_major %in% MAJOR_TYPES) %>%
  mutate(treatment_clean = toupper(treatment)) %>%
  count(treatment_clean, cell_type_major) %>%
  group_by(treatment_clean) %>% mutate(prop = n/sum(n)) %>% ungroup() %>%
  pivot_wider(names_from = treatment_clean, values_from = c(n, prop)) %>%
  mutate(delta_m02 = prop_MBRT - prop_SBRT)

m01_comp <- read_tsv(file.path(DATA_DIR, "composition_delta_mbrt_vs_sbrt.tsv"),
                     show_col_types = FALSE) %>%
  filter(timepoint_h == 48) %>%
  rename(delta_m01 = delta_MBRT_SBRT) %>%
  select(cell_type_major, delta_m01)

comp_concord <- m02_comp %>% left_join(m01_comp, by = "cell_type_major")
write_tsv(comp_concord, file.path(DATA_DIR, "set2_composition_concordance.tsv"))
cat("\n=== Composition concordance day2 ===\n"); print(comp_concord)

# Composition concordance plot
p_comp <- comp_concord %>%
  pivot_longer(c(delta_m01, delta_m02), names_to = "dataset", values_to = "delta") %>%
  mutate(dataset = ifelse(dataset == "delta_m01", "Mutter_01 (n=1)", "Mutter_02 (n=4)")) %>%
  ggplot(aes(x = cell_type_major, y = delta, fill = dataset)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 0, color = "grey50") +
  scale_fill_manual(values = c("Mutter_01 (n=1)" = "#7F8C8D", "Mutter_02 (n=4)" = "#16A085")) +
  labs(title = "Composition delta (MBRT − SBRT) at day2: Mutter_01 vs Mutter_02",
       y = "ΔProp (MBRT − SBRT)", x = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(PLOT_DIR, "set2_composition_concordance.png"),
       plot = p_comp, width = 10, height = 5, dpi = 150)

cat("\n=== 08b_set2_validation_fast.R complete ===\n")
