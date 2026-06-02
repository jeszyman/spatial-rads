## 01_degs_kinetics.R
## DEGs by cell type x timepoint x contrast (3 contrasts)
##   - MBRT vs Control (1h, 4h, day2, day6)
##   - SBRT vs Control (4h, day2, day6)
##   - MBRT vs SBRT  (4h, day2, day6)
## Per (timepoint x contrast x cell_type_major) where both groups have >= 100 cells.
## Runs FindMarkers (Wilcox), effect-size ranked. Necrosis_zone cells excluded
## at day2 and day6 (default downstream behavior).

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(future)
  library(future.apply)
})

DATA_DIR <- "data"
PLOT_DIR <- "plots"
OBJECTS_DIR <- "/mnt/data/projects/spatial-rads/analysis/objects/mbrt_vs_sbrt"

obj <- readRDS(file.path(OBJECTS_DIR, "seurat_filtered.rds"))
cat(sprintf("Loaded filtered object: %d genes x %d cells\n", nrow(obj), ncol(obj)))

# Apply default necrosis-zone exclusion at day2 (48h) and day6 (144h)
exclude_necr_at <- c(48, 144)
keep_mask <- !(obj$timepoint_h %in% exclude_necr_at & isTRUE(obj$necrosis_zone))
keep_mask <- !(obj$timepoint_h %in% exclude_necr_at & obj$necrosis_zone %in% TRUE)
n_drop <- sum(!keep_mask, na.rm = TRUE)
obj <- subset(obj, cells = colnames(obj)[keep_mask])
cat(sprintf("Dropped %d necrosis-zone cells at day2/day6; %d cells remaining\n", n_drop, ncol(obj)))

# DEG helper
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

# Define contrasts
contrasts_def <- list(
  list(name = "MBRT_vs_Control",  ident.1 = "MBRT_1h",   ident.2 = "Control",   timepoint_h = 1),
  list(name = "MBRT_vs_Control",  ident.1 = "MBRT_4h",   ident.2 = "Control",   timepoint_h = 4),
  list(name = "MBRT_vs_Control",  ident.1 = "MBRT_day2", ident.2 = "Control",   timepoint_h = 48),
  list(name = "MBRT_vs_Control",  ident.1 = "MBRT_day6", ident.2 = "Control",   timepoint_h = 144),
  list(name = "SBRT_vs_Control",  ident.1 = "SBRT_4h",   ident.2 = "Control",   timepoint_h = 4),
  list(name = "SBRT_vs_Control",  ident.1 = "SBRT_day2", ident.2 = "Control",   timepoint_h = 48),
  list(name = "SBRT_vs_Control",  ident.1 = "SBRT_day6", ident.2 = "Control",   timepoint_h = 144),
  list(name = "MBRT_vs_SBRT",     ident.1 = "MBRT_4h",   ident.2 = "SBRT_4h",   timepoint_h = 4),
  list(name = "MBRT_vs_SBRT",     ident.1 = "MBRT_day2", ident.2 = "SBRT_day2", timepoint_h = 48),
  list(name = "MBRT_vs_SBRT",     ident.1 = "MBRT_day6", ident.2 = "SBRT_day6", timepoint_h = 144)
)

MAJOR_TYPES <- c("tumor_epithelial", "endothelial", "fibroblast", "smooth_muscle",
                 "T.cell", "B.cell", "myeloid", "NK", "plasma")

# Pre-compute strata sizes for filtering
strata_n <- as_tibble(obj@meta.data) %>% count(Condition, cell_type_major)

# Build job list: (contrast, cell_type) tuples passing n>=100 in both groups
jobs <- list()
for (ce in contrasts_def) {
  for (ct in MAJOR_TYPES) {
    n1 <- strata_n %>% filter(Condition == ce$ident.1, cell_type_major == ct) %>% pull(n)
    n2 <- strata_n %>% filter(Condition == ce$ident.2, cell_type_major == ct) %>% pull(n)
    if (length(n1) == 0) n1 <- 0
    if (length(n2) == 0) n2 <- 0
    if (n1 >= 100 && n2 >= 100) {
      jobs[[length(jobs) + 1]] <- c(ce, list(cell_type = ct, n1 = n1, n2 = n2))
    }
  }
}
# Also: bulk DEGs (all cells per condition, no celltype split)
bulk_jobs <- lapply(contrasts_def, function(ce) c(ce, list(cell_type = "BULK", n1 = NA, n2 = NA)))
jobs <- c(jobs, bulk_jobs)
cat(sprintf("Total DEG jobs: %d (%d per-celltype + %d bulk)\n",
            length(jobs), length(jobs) - length(bulk_jobs), length(bulk_jobs)))

# Run jobs sequentially (FindMarkers internal loops; future complicates with large Seurat objects)
all_degs <- list()
t0 <- Sys.time()
for (i in seq_along(jobs)) {
  job <- jobs[[i]]
  ct <- job$cell_type
  cells_keep <- if (ct == "BULK") {
    colnames(obj)[obj$Condition %in% c(job$ident.1, job$ident.2)]
  } else {
    colnames(obj)[obj$Condition %in% c(job$ident.1, job$ident.2) & obj$cell_type_major == ct]
  }
  if (length(cells_keep) < 50) next
  obj_sub <- subset(obj, cells = cells_keep)
  deg <- run_deg(obj_sub, job$ident.1, job$ident.2, "Condition")
  if (!is.null(deg)) {
    deg$contrast <- job$name
    deg$timepoint_h <- job$timepoint_h
    deg$cell_type_major <- ct
    deg$ident.1 <- job$ident.1
    deg$ident.2 <- job$ident.2
    all_degs[[length(all_degs) + 1]] <- deg
  }
  if (i %% 5 == 0) {
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    cat(sprintf("  [%d/%d] elapsed %.0fs\n", i, length(jobs), elapsed))
  }
}
all_degs_df <- bind_rows(all_degs)
cat(sprintf("\nTotal DEG rows: %d\n", nrow(all_degs_df)))

# Add convenience cols
all_degs_df <- all_degs_df %>%
  rename(log2FC = avg_log2FC) %>%
  mutate(timepoint_label = case_when(
    timepoint_h == 1 ~ "1h",
    timepoint_h == 4 ~ "4h",
    timepoint_h == 48 ~ "day2",
    timepoint_h == 144 ~ "day6"
  )) %>%
  arrange(contrast, timepoint_h, cell_type_major, desc(abs(log2FC)))

write_tsv(all_degs_df, file.path(DATA_DIR, "degs_kinetics.tsv"))
cat(sprintf("Wrote: %s\n", file.path(DATA_DIR, "degs_kinetics.tsv")))

# Volcano plots retired: -log10(p_val) here is a cell-level Wilcoxon p (FindMarkers),
# which pseudoreplicates the n=1-per-condition design -- cells are not biological
# replicates. The effect-size tables (degs_kinetics.tsv, degs_top10_per_stratum.tsv)
# stay valid descriptively. Correct significance testing is deferred until Mutter_02
# (2d) replicates: pseudobulk per sample -> DESeq2/edgeR/limma -> MA-plot/volcano with
# real p-values.

# Top-table summary: top 10 |log2FC| per (contrast x timepoint x cell_type_major)
top_table <- all_degs_df %>%
  group_by(contrast, timepoint_label, cell_type_major) %>%
  slice_max(abs(log2FC), n = 10) %>%
  ungroup() %>%
  select(contrast, timepoint_label, cell_type_major, gene, log2FC, pct.1, pct.2, p_val_adj)
write_tsv(top_table, file.path(DATA_DIR, "degs_top10_per_stratum.tsv"))

cat("\n=== Sanity check: 4h MBRT-vs-Control tumor DDR genes ===\n")
ddr_check <- all_degs_df %>%
  filter(contrast == "MBRT_vs_Control",
         timepoint_h == 4,
         cell_type_major == "tumor_epithelial",
         gene %in% c("Cdkn1a", "Mdm2", "Bax", "Gadd45a", "Gadd45b", "Trp53", "Ccng1")) %>%
  select(gene, log2FC, pct.1, pct.2, p_val_adj) %>%
  arrange(desc(log2FC))
print(ddr_check)

cat("\n=== 01_degs_kinetics.R complete ===\n")
