library(data.table)
library(Seurat)
library(arrow)
library(dplyr)

# --- Load counts ---
counts_df <- read_parquet("/mnt/gcs/jeszyman/projects/spatial-rads/inputs/projects_Mutter_01_CosMmR_app_data_raw_tx_counts_matrix.parquet")
counts_dt <- as.data.table(counts_df)

non_gene_cols <- c("Slide", "fov", "cell_id")
gene_cols <- setdiff(colnames(counts_dt), non_gene_cols)
cat(sprintf("Panel size: %d genes\n", length(gene_cols)))
cat(sprintf("Total cells in counts: %d\n", nrow(counts_dt)))
stopifnot(!anyDuplicated(counts_dt$cell_id))

# --- Load metadata ---
meta_df <- read_parquet("/mnt/gcs/jeszyman/projects/spatial-rads/inputs/Analysis_Mutter_01_CosMmR_Mutter_updated_metadata.parquet")
meta_dt <- as.data.table(meta_df)
cat(sprintf("Total cells in metadata: %d\n", nrow(meta_dt)))

# --- Build count matrix (genes x cells) ---
count_mat <- t(as.matrix(counts_dt[, ..gene_cols]))
colnames(count_mat) <- counts_dt$cell_id
rownames(count_mat) <- gene_cols

meta_for_seurat <- as.data.frame(meta_dt)
rownames(meta_for_seurat) <- meta_for_seurat$cell_id
stopifnot(all(colnames(count_mat) %in% rownames(meta_for_seurat)))

# --- Create Seurat object ---
obj <- CreateSeuratObject(
  counts    = count_mat,
  meta.data = meta_for_seurat,
  project   = "CosMx_Mutter"
)
cat(sprintf("Raw Seurat: %d genes x %d cells\n", nrow(obj), ncol(obj)))

# --- QC filtering ---
pre_qc <- table(obj$Condition)
obj <- subset(obj, subset = qcFlagsCell == "Pass")
obj <- subset(obj, subset = nCount_RNA > 20 & nFeature_RNA > 10 & propNegative < 0.5)
post_qc <- table(obj$Condition)

qc_summary <- data.frame(
  condition = names(pre_qc),
  pre_filter = as.integer(pre_qc),
  post_filter = as.integer(post_qc[names(pre_qc)]),
  pct_retained = round(as.integer(post_qc[names(pre_qc)]) / as.integer(pre_qc) * 100, 1)
)
print(qc_summary)
write.csv(qc_summary, "/mnt/gcs/jeszyman/projects/spatial-rads/analysis/tables/qc_summary.csv", row.names = FALSE)
cat(sprintf("Post-QC: %d cells retained\n", ncol(obj)))

# --- Add parsed metadata ---
obj$model <- ifelse(grepl("^Tongue", obj$Condition), "tongue", "flank")
obj$treatment <- case_when(
  obj$Condition == "Control" ~ "NT",
  grepl("MBRT", obj$Condition) ~ "MBRT",
  grepl("SBRT", obj$Condition) ~ "SBRT"
)
obj$timepoint_h <- case_when(
  obj$Condition == "Control" ~ 0,
  grepl("1h", obj$Condition) ~ 1,
  grepl("4h", obj$Condition) ~ 4,
  grepl("day2", obj$Condition) ~ 48,
  grepl("day6", obj$Condition) ~ 144,
  grepl("day8", obj$Condition) ~ 192,
  grepl("day10", obj$Condition) ~ 240
)

# --- Verify pathway columns exist ---
pathway_cols <- c("TypeI_interferon", "TypeII_interferon", "STING", "DNA_Damage_Repair")
stopifnot(all(pathway_cols %in% colnames(obj@meta.data)))
cat("Pathway score columns confirmed present.\n")

# --- Save ---
saveRDS(obj, "/mnt/gcs/jeszyman/projects/spatial-rads/analysis/objects/seurat_qc.rds")
cat("QC'd Seurat object saved.\n")
