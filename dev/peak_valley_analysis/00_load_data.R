library(data.table)
library(Seurat)
library(arrow)
library(tidyverse)

# --- Paths ---
INPUT_DIR <- "/mnt/gcs/jeszyman/projects/spatial-rads/inputs"
ANALYSIS_DIR <- "/mnt/gcs/jeszyman/projects/spatial-rads/analysis"
DATA_DIR <- "data"
PLOT_DIR <- "plots"
for (d in c(DATA_DIR, PLOT_DIR)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# --- Load pre-computed object if available, otherwise build from parquet ---
cached_rds <- file.path(ANALYSIS_DIR, "objects", "seurat_clustered.rds")
if (file.exists(cached_rds)) {
  cat("Loading pre-computed clustered object...\n")
  obj <- readRDS(cached_rds)
  cat(sprintf("Loaded: %d genes x %d cells\n", nrow(obj), ncol(obj)))
} else {
  cat("No cached object — building from parquet...\n")
  counts_df <- read_parquet(file.path(INPUT_DIR, "projects_Mutter_01_CosMmR_app_data_raw_tx_counts_matrix.parquet"))
  counts_dt <- as.data.table(counts_df)
  non_gene_cols <- c("Slide", "fov", "cell_id")
  gene_cols <- setdiff(colnames(counts_dt), non_gene_cols)
  cat(sprintf("Panel size: %d genes\n", length(gene_cols)))
  count_mat <- t(as.matrix(counts_dt[, ..gene_cols]))
  colnames(count_mat) <- counts_dt$cell_id
  rownames(count_mat) <- gene_cols
  meta_df <- read_parquet(file.path(INPUT_DIR, "Analysis_Mutter_01_CosMmR_Mutter_updated_metadata.parquet"))
  meta_for_seurat <- as.data.frame(meta_df)
  rownames(meta_for_seurat) <- meta_for_seurat$cell_id
  obj <- CreateSeuratObject(counts = count_mat, meta.data = meta_for_seurat, project = "CosMx_Mutter")
}

# --- Ensure parsed metadata columns exist ---
if (!"model" %in% colnames(obj@meta.data)) {
  obj$model <- ifelse(grepl("^Tongue", obj$Condition), "tongue", "flank")
}
if (!"treatment" %in% colnames(obj@meta.data)) {
  obj$treatment <- case_when(
    obj$Condition == "Control" ~ "NT",
    grepl("MBRT", obj$Condition) ~ "MBRT",
    grepl("SBRT", obj$Condition) ~ "SBRT"
  )
}
if (!"timepoint_h" %in% colnames(obj@meta.data)) {
  obj$timepoint_h <- case_when(
    obj$Condition == "Control" ~ 0,
    grepl("1h", obj$Condition) ~ 1,
    grepl("4h", obj$Condition) ~ 4,
    grepl("day2", obj$Condition) ~ 48,
    grepl("day6", obj$Condition) ~ 144,
    grepl("day8", obj$Condition) ~ 192,
    grepl("day10", obj$Condition) ~ 240
  )
}

# --- Verify pathway columns ---
PATHWAY_COLS <- c("TypeI_interferon", "TypeII_interferon", "STING", "DNA_Damage_Repair")
stopifnot(all(PATHWAY_COLS %in% colnames(obj@meta.data)))
cat("Pathway score columns confirmed present.\n")

# --- Canonical color palettes ---
TREATMENT_COLORS <- c("MBRT" = "#E74C3C", "SBRT" = "#3498DB", "NT" = "#7F8C8D")

CONDITION_LEVELS <- c("Control", "MBRT_1h", "MBRT_4h", "SBRT_4h",
                      "MBRT_day2", "SBRT_day2", "MBRT_day6", "SBRT_day6",
                      "Tongue_MBRT_day8", "Tongue_MBRT_day10", "Tongue_SBRT_day10")

ZONE_COLORS <- c("peak" = "#E74C3C", "valley" = "#3498DB")

PATHWAY_LABELS <- c(
  "TypeI_interferon" = "Type I IFN",
  "TypeII_interferon" = "Type II IFN",
  "STING" = "STING",
  "DNA_Damage_Repair" = "DDR"
)

cat(sprintf("Raw object ready: %d genes x %d cells\n", nrow(obj), ncol(obj)))
