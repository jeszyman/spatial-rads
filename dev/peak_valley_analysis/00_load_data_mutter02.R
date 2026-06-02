library(Seurat)
library(tidyverse)
library(data.table)

# --- Paths ---
INPUT_DIR <- "/mnt/data/projects/spatial-rads/inputs/mutter02"
OUTPUT_DIR <- "data"
PLOT_DIR <- "plots"
for (d in c(OUTPUT_DIR, PLOT_DIR)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# --- RDS file paths ---
rds_files <- list(
  slide1 = file.path(INPUT_DIR, "seuratObject_01_Mutter_02_CosMmR.RDS"),
  slide2 = file.path(INPUT_DIR, "seuratObject_02_Mutter_02_CosMmR.RDS"),
  slide3 = file.path(INPUT_DIR, "seuratObject_03_Mutter_02_CosMmR.RDS"),
  slide4 = file.path(INPUT_DIR, "seuratObject_04_Mutter_02_CosMmR.RDS")
)
stopifnot(all(file.exists(unlist(rds_files))))

# --- Pathway gene lists (NanoString Mouse UCC panel modules) ---
# Type I IFN (19 genes), Type II IFN (17 genes), DDR (18 genes) from Annotations sheet
# STING is custom: Irf3, Ifna, Isg15, Cxcl10
PATHWAY_GENES <- list(
  TypeI_interferon = c("Ifit1", "Ifit2", "Ifit3", "Isg15", "Mx1", "Mx2", "Oas1a", "Oas2",
                       "Oasl1", "Rsad2", "Ifitm1", "Ifitm3", "Irf7", "Stat1", "Stat2",
                       "Isg20", "Usp18", "Ifi44", "Zbp1"),
  TypeII_interferon = c("Gbp2", "Gbp3", "Gbp4", "Gbp5", "Gbp6", "Gbp7", "Irf1", "Irf8",
                        "Stat1", "Cxcl9", "Cxcl10", "Cxcl11", "Ciita", "H2-Aa", "H2-Ab1",
                        "H2-Eb1", "Cd74"),
  DNA_Damage_Repair = c("Cdkn1a", "Gadd45a", "Gadd45b", "Mdm2", "Bax", "Pmaip1", "Bbc3",
                        "Ddit3", "Atm", "Atr", "Chek1", "Chek2", "Rad51", "Brca1", "Brca2",
                        "Xrcc5", "Xrcc6", "Parp1"),
  STING = c("Irf3", "Isg15", "Cxcl10")  # Ifna not typically single-gene, skip
)

# --- Sample layout per Jenn Fazzari (2026-04-09) + spatial assignment ---
# Each slide has 3 samples arranged top-to-bottom (highest to lowest y_slide_mm):
#   top = Control (0 Gy), middle = MBRT 2d, bottom = SBRT 20 Gy 2d
# Verified against CosMx Control Center FOV maps in Mutter02_Slide Powerpoint.pptx
SAMPLE_LAYOUT <- tribble(
  ~slide, ~rank, ~sample_id, ~treatment,
  1, 1, 1,  "Control",  1, 2, 34, "MBRT",  1, 3, 28, "SBRT",
  2, 1, 2,  "Control",  2, 2, 35, "MBRT",  2, 3, 30, "SBRT",
  3, 1, 3,  "Control",  3, 2, 19, "MBRT",  3, 3, 16, "SBRT",
  4, 1, 1,  "Control",  4, 2, 20, "MBRT",  4, 3, 18, "SBRT"
)

SAMPLE_MAP <- tribble(
  ~slide, ~sample_id, ~treatment, ~timepoint_h,
  1, 1, "Control", 48,
  1, 34, "MBRT", 48,
  1, 28, "SBRT", 48,
  2, 2, "Control", 48,
  2, 35, "MBRT", 48,
  2, 30, "SBRT", 48,
  3, 3, "Control", 48,
  3, 19, "MBRT", 48,
  3, 16, "SBRT", 48,
  4, 1, "Control", 48,
  4, 20, "MBRT", 48,
  4, 18, "SBRT", 48
)

# --- Helper: compute pathway scores using AddModuleScore ---
compute_pathway_scores <- function(obj, pathway_genes) {
  available_genes <- rownames(obj)
  n_genes <- length(available_genes)
  # Adjust ctrl to avoid sampling error with small gene panels
  ctrl_size <- min(50, floor(n_genes / 10))

  for (pathway_name in names(pathway_genes)) {
    genes <- pathway_genes[[pathway_name]]
    genes_present <- genes[genes %in% available_genes]
    if (length(genes_present) < 2) {
      warning(sprintf("Pathway %s: only %d/%d genes found, skipping",
                      pathway_name, length(genes_present), length(genes)))
      obj[[pathway_name]] <- NA
      next
    }
    cat(sprintf("Computing %s score (%d/%d genes, ctrl=%d)...\n",
                pathway_name, length(genes_present), length(genes), ctrl_size))
    obj <- AddModuleScore(obj, features = list(genes_present), name = pathway_name,
                          ctrl = ctrl_size, nbin = 12)
    obj[[pathway_name]] <- obj[[paste0(pathway_name, "1")]]
    obj[[paste0(pathway_name, "1")]] <- NULL
  }
  obj
}

# --- Load and process each slide ---
cat("Loading Mutter_02 slides...\n")
objects <- list()

for (slide_name in names(rds_files)) {
  slide_num <- as.integer(gsub("slide", "", slide_name))
  cat(sprintf("\n=== Processing %s (slide %d) ===\n", slide_name, slide_num))

  # Load and update to current Seurat version
  obj <- readRDS(rds_files[[slide_name]])
  obj <- UpdateSeuratObject(obj)
  cat(sprintf("Loaded: %d genes x %d cells\n", nrow(obj), ncol(obj)))

  # Add slide identifier
  obj$slide <- slide_num
  obj$dataset <- "Mutter_02"

  # Inspect metadata columns
  cat("Metadata columns:", paste(head(colnames(obj@meta.data), 20), collapse = ", "), "...\n")

  # --- QC filtering (apply Mutter_01 thresholds) ---
  pre_qc <- ncol(obj)

  # Check which QC columns exist
  qc_cols <- c("qcFlagsCell", "nCount_RNA", "nFeature_RNA", "propNegative")
  existing_qc <- qc_cols[qc_cols %in% colnames(obj@meta.data)]
  cat("QC columns available:", paste(existing_qc, collapse = ", "), "\n")

  # Build filter expression based on available columns
  keep <- rep(TRUE, ncol(obj))
  if ("qcFlagsCell" %in% colnames(obj@meta.data)) {
    keep <- keep & (obj$qcFlagsCell == "Pass")
  }
  if ("nCount_RNA" %in% colnames(obj@meta.data)) {
    keep <- keep & (obj$nCount_RNA > 20)
  }
  if ("nFeature_RNA" %in% colnames(obj@meta.data)) {
    keep <- keep & (obj$nFeature_RNA > 10)
  }
  if ("propNegative" %in% colnames(obj@meta.data)) {
    keep <- keep & (obj$propNegative < 0.5)
  }

  obj <- obj[, keep]
  post_qc <- ncol(obj)
  cat(sprintf("QC: %d -> %d cells (%.1f%% retained)\n", pre_qc, post_qc, 100*post_qc/pre_qc))

  # --- Assign samples by y-coordinate clustering ---
  km <- kmeans(obj$y_slide_mm, centers = 3, nstart = 10)
  rank_map <- order(km$centers[, 1], decreasing = TRUE)
  cluster_rank <- match(km$cluster, rank_map)
  layout <- SAMPLE_LAYOUT %>% filter(slide == slide_num)
  obj$spatial_rank <- cluster_rank
  obj$sample_id    <- layout$sample_id[cluster_rank]
  obj$treatment    <- layout$treatment[cluster_rank]
  obj$timepoint_h  <- 48
  cat(sprintf("Sample assignment: %s\n",
    paste(names(table(obj$treatment)), table(obj$treatment), sep = "=", collapse = ", ")))

  # --- Compute pathway scores ---
  obj <- compute_pathway_scores(obj, PATHWAY_GENES)

  # --- Store ---
  objects[[slide_name]] <- obj
}

# --- Summary stats ---
cat("\n=== Summary ===\n")
summary_df <- data.frame(
  slide = names(objects),
  cells = sapply(objects, ncol),
  genes = sapply(objects, nrow)
)
print(summary_df)
cat(sprintf("Total cells: %d\n", sum(summary_df$cells)))

# --- Save per-slide objects ---
for (slide_name in names(objects)) {
  out_path <- file.path(OUTPUT_DIR, paste0("seurat_mutter02_", slide_name, "_qc.rds"))
  saveRDS(objects[[slide_name]], out_path)
  cat(sprintf("Saved: %s\n", out_path))
}

# --- Write QC summary ---
write_tsv(summary_df, file.path(OUTPUT_DIR, "qc_summary_mutter02.tsv"))
cat("Wrote: data/qc_summary_mutter02.tsv\n")

cat("\nMutter_02 preprocessing complete. Per-slide objects saved with sample assignments.\n")
