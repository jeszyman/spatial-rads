if (!exists("PATHWAY_COLS")) source("00_load_data.R")
library(UCell)

# Prefer local copy over GCS mount (FUSE reads are slow for 1.7G object)
local_rds <- file.path(Sys.getenv("HOME"), "data/spatial-rads/seurat_clustered.rds")
gcs_rds <- file.path(ANALYSIS_DIR, "objects", "seurat_clustered.rds")
rds_path <- if (file.exists(local_rds)) local_rds else gcs_rds
cat(sprintf("Loading Seurat object from: %s\n", rds_path))
obj <- readRDS(rds_path)
cat(sprintf("Loaded: %d genes x %d cells\n", nrow(obj), ncol(obj)))

# --- Load signature gene lists ---
peak_bulk   <- read_tsv(file.path(DATA_DIR, "peak_signature_bulk.tsv"), show_col_types = FALSE)
valley_bulk <- read_tsv(file.path(DATA_DIR, "valley_signature_bulk.tsv"), show_col_types = FALSE)

peak_genes   <- peak_bulk$gene
valley_genes <- valley_bulk$gene
cat(sprintf("Signature genes: %d peak-up, %d valley-up\n", length(peak_genes), length(valley_genes)))

# --- UCell scoring (primary) ---
cat("Scoring with UCell...\n")
signatures <- list(peak_up = peak_genes, valley_up = valley_genes)
obj <- AddModuleScore_UCell(obj, features = signatures, name = NULL)
# UCell adds columns: peak_up, valley_up

# --- AddModuleScore scoring (secondary) ---
# ctrl must be smaller than panel size (~1000 genes); default 100 is fine but nbin can cause
# sampling issues with small panels. Use ctrl=50, nbin=12 to be safe.
cat("Scoring with AddModuleScore...\n")
n_genes <- nrow(obj)
ctrl_size <- min(50, floor(n_genes / 10))
cat(sprintf("  Panel: %d genes, ctrl=%d, nbin=12\n", n_genes, ctrl_size))
obj <- AddModuleScore(obj, features = list(peak_genes), name = "peak_up_AMS",
                      seed = 42, ctrl = ctrl_size, nbin = 12)
obj <- AddModuleScore(obj, features = list(valley_genes), name = "valley_up_AMS",
                      seed = 42, ctrl = ctrl_size, nbin = 12)
# Rename the auto-numbered columns
obj$peak_up_AMS <- obj$peak_up_AMS1
obj$peak_up_AMS1 <- NULL
obj$valley_up_AMS <- obj$valley_up_AMS1
obj$valley_up_AMS1 <- NULL

# --- Method correlation ---
cat("Computing method correlation...\n")
meta <- obj@meta.data %>% as.data.frame()
cor_peak <- cor(meta$peak_up, meta$peak_up_AMS, use = "complete.obs", method = "spearman")
cor_valley <- cor(meta$valley_up, meta$valley_up_AMS, use = "complete.obs", method = "spearman")
cor_df <- tibble(
  signature = c("peak_up", "valley_up"),
  spearman_r = c(cor_peak, cor_valley),
  n_cells = c(sum(!is.na(meta$peak_up)), sum(!is.na(meta$valley_up)))
)
write_tsv(cor_df, file.path(DATA_DIR, "scoring_method_correlation.tsv"))
cat(sprintf("UCell vs AMS correlation: peak r=%.3f, valley r=%.3f\n", cor_peak, cor_valley))

# --- Per-condition × cell-type summary (UCell scores) ---
cat("Summarizing scores by condition and cell type...\n")
score_summary <- meta %>%
  group_by(Condition, cell_type_validated) %>%
  summarise(
    mean_peak_score = mean(peak_up, na.rm = TRUE),
    sd_peak_score = sd(peak_up, na.rm = TRUE),
    mean_valley_score = mean(valley_up, na.rm = TRUE),
    sd_valley_score = sd(valley_up, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )
write_tsv(score_summary, file.path(DATA_DIR, "signature_scores_all_cells.tsv"))

# --- Save scored object for downstream scripts ---
saveRDS(obj, file.path(DATA_DIR, "seurat_scored.rds"))
cat(sprintf("Scored object saved: %s (%.1f GB)\n",
            file.path(DATA_DIR, "seurat_scored.rds"),
            file.size(file.path(DATA_DIR, "seurat_scored.rds")) / 1e9))

cat("Signature scoring complete.\n")
