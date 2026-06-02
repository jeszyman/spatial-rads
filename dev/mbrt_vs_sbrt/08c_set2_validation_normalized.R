## 08c_set2_validation_normalized.R
## Bugfix: Mutter_02 RDS "data" layer was identical to counts (never normalized).
## Re-normalize, then redo DEGs + concordance with Mutter_01.

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
})

DATA_DIR <- "data"
PLOT_DIR <- "plots"
OBJECTS_DIR <- "/mnt/data/projects/spatial-rads/analysis/objects/mbrt_vs_sbrt"

tmp_data <- "/mnt/data/projects/spatial-rads/analysis/tmp"
if (!dir.exists(tmp_data)) dir.create(tmp_data, recursive = TRUE)
Sys.setenv(TMPDIR = tmp_data, TMP = tmp_data, TEMP = tmp_data)

typed_rds <- file.path(OBJECTS_DIR, "seurat_mutter02_typed.rds")
stopifnot(file.exists(typed_rds))
obj <- readRDS(typed_rds)
cat(sprintf("Loaded typed Mutter_02: %d genes x %d cells\n", nrow(obj), ncol(obj)))

# Verify data layer issue, then re-normalize
counts_first <- GetAssayData(obj, layer = "counts")[1:3, 1:3]
data_first   <- GetAssayData(obj, layer = "data")[1:3, 1:3]
cat("Before normalize: counts == data?\n"); print(all.equal(counts_first, data_first))

cat("\nRunning NormalizeData (LogNormalize, scale.factor=10000)...\n")
t0 <- Sys.time()
obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                     scale.factor = 10000, verbose = FALSE)
cat(sprintf("NormalizeData done: %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# Verify normalization succeeded
data_first_post <- GetAssayData(obj, layer = "data")[1:3, 1:3]
cat("After normalize: counts == data?\n"); print(all.equal(counts_first, data_first_post))
cat("Sample data range post-normalize:", range(data_first_post), "\n")

# Save normalized object
saveRDS(obj, file.path(OBJECTS_DIR, "seurat_mutter02_typed_normalized.rds"))

# Harmonize labels and apply default day2 necrosis exclusion
obj$treatment_clean <- case_when(
  toupper(obj$treatment) %in% c("MBRT") ~ "MBRT",
  toupper(obj$treatment) %in% c("SBRT") ~ "SBRT",
  toupper(obj$treatment) %in% c("CONTROL","CTRL","NT") ~ "NT",
  TRUE ~ as.character(obj$treatment))
obj$Condition_m02 <- ifelse(obj$treatment_clean == "NT", "Control",
                            paste0(obj$treatment_clean, "_day2"))
keep_main <- !(obj$timepoint_h == 48 & obj$necrosis_zone %in% TRUE)
obj <- subset(obj, cells = colnames(obj)[keep_main])
cat(sprintf("After day2 necrosis-zone exclusion: %d cells\n", ncol(obj)))

# Subsample to 3000 per stratum for fast DEG
set.seed(42)
MAJOR_TYPES <- c("tumor_epithelial", "endothelial", "fibroblast", "smooth_muscle",
                 "T.cell", "B.cell", "myeloid", "NK", "plasma")
meta <- as_tibble(obj@meta.data) %>%
  mutate(cell_id = colnames(obj)) %>%
  filter(cell_type_major %in% MAJOR_TYPES, !is.na(Condition_m02))
cells_to_keep <- meta %>%
  group_by(Condition_m02, cell_type_major) %>%
  slice_sample(n = 3000) %>% ungroup() %>% pull(cell_id)
cat(sprintf("After subsampling to 3000 per stratum: %d cells\n", length(cells_to_keep)))
obj <- subset(obj, cells = cells_to_keep)

# DEGs
day2_contrasts <- list(
  list(name = "MBRT_vs_Control", ident.1 = "MBRT_day2", ident.2 = "Control"),
  list(name = "SBRT_vs_Control", ident.1 = "SBRT_day2", ident.2 = "Control"),
  list(name = "MBRT_vs_SBRT",    ident.1 = "MBRT_day2", ident.2 = "SBRT_day2"))

run_deg <- function(obj_sub, ident.1, ident.2, group.by) {
  Idents(obj_sub) <- group.by
  out <- tryCatch(
    FindMarkers(obj_sub, ident.1 = ident.1, ident.2 = ident.2,
                min.pct = 0.1, logfc.threshold = 0.1, verbose = FALSE),
    error = function(e) NULL)
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
    if (n1 >= 100 && n2 >= 100) deg_jobs[[length(deg_jobs)+1]] <- c(ce, list(cell_type=ct, n1=n1, n2=n2))
  }
  deg_jobs[[length(deg_jobs)+1]] <- c(ce, list(cell_type="BULK", n1=NA, n2=NA))
}
cat(sprintf("Total DEG jobs: %d\n", length(deg_jobs)))

m02_degs <- list()
t0 <- Sys.time()
for (i in seq_along(deg_jobs)) {
  job <- deg_jobs[[i]]; ct <- job$cell_type
  cells_keep <- if (ct == "BULK") {
    colnames(obj)[obj$Condition_m02 %in% c(job$ident.1, job$ident.2)]
  } else {
    colnames(obj)[obj$Condition_m02 %in% c(job$ident.1, job$ident.2) & obj$cell_type_major == ct]
  }
  if (length(cells_keep) < 50) next
  obj_sub <- subset(obj, cells = cells_keep)
  d <- run_deg(obj_sub, job$ident.1, job$ident.2, "Condition_m02")
  if (!is.null(d)) {
    d$contrast <- job$name; d$cell_type_major <- ct
    d$ident.1 <- job$ident.1; d$ident.2 <- job$ident.2
    m02_degs[[length(m02_degs)+1]] <- d
  }
}
m02_degs_df <- bind_rows(m02_degs) %>%
  rename(log2FC = avg_log2FC) %>%
  arrange(contrast, cell_type_major, desc(abs(log2FC)))
write_tsv(m02_degs_df, file.path(DATA_DIR, "mutter02_degs_day2.tsv"))

cat(sprintf("\nWrote %d DEG rows. log2FC range: [%.2f, %.2f]\n",
            nrow(m02_degs_df), min(m02_degs_df$log2FC), max(m02_degs_df$log2FC)))

# Concordance
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
            n_genes = n(), .groups = "drop") %>%
  arrange(contrast, desc(rho))
write_tsv(rho_table, file.path(DATA_DIR, "set2_concordance_rho.tsv"))
cat("\n=== Concordance rho (post-normalization) ===\n"); print(rho_table)

# Plots
TREATMENT_COLORS <- c("MBRT" = "#E74C3C", "SBRT" = "#3498DB", "NT" = "#7F8C8D")

p_conc <- concord %>%
  filter(contrast == "MBRT_vs_SBRT", cell_type_major != "BULK") %>%
  ggplot(aes(x = log2FC_m01, y = log2FC_m02)) +
  geom_point(alpha = 0.4, size = 0.6) +
  geom_abline(slope = 1, intercept = 0, color = "red3", linetype = "dashed") +
  geom_hline(yintercept = 0, color = "grey60") +
  geom_vline(xintercept = 0, color = "grey60") +
  facet_wrap(~cell_type_major, ncol = 3) +
  labs(title = "Mutter_01 day2 vs Mutter_02 day2 (post-normalization) — MBRT vs SBRT log2FC",
       subtitle = "Sub-sampled to 3000 cells per stratum; LogNormalize re-applied",
       x = "Mutter_01 log2FC", y = "Mutter_02 log2FC") +
  theme_bw()
ggsave(file.path(PLOT_DIR, "set2_concordance_scatter.png"),
       plot = p_conc, width = 12, height = 9, dpi = 130)

p_rho <- rho_table %>%
  ggplot(aes(x = contrast, y = cell_type_major, fill = rho)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("rho=%.2f\nn=%d", rho, n_genes)), size = 3) +
  scale_fill_gradient2(low = "#3498DB", mid = "white", high = "#E74C3C",
                       midpoint = 0, limits = c(-1, 1), na.value = "grey90") +
  labs(title = "Mutter_01 vs Mutter_02 day2 effect-size concordance (post-normalization)",
       subtitle = "Spearman rank correlation",
       x = "Contrast", y = NULL) +
  theme_bw() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(PLOT_DIR, "set2_concordance_rho.png"),
       plot = p_rho, width = 8, height = 6, dpi = 150)

cat("\n=== 08c_set2_validation_normalized.R complete ===\n")
