## 08d_insitutype_mutter02.R (v2)
## Properly cell-type Mutter_02 by:
## 1. Building a reference profile from Mutter_01 labeled cells: gene × cell-type
##    matrix of mean expression (skip InSituType::Estep, use direct mean per type).
## 2. Running InSituType supervised on Mutter_02.

suppressPackageStartupMessages({
  library(Seurat)
  library(InSituType)
  library(tidyverse)
  library(Matrix)
})

DATA_DIR <- "data"
OBJECTS_DIR <- "/mnt/data/projects/spatial-rads/analysis/objects/mbrt_vs_sbrt"
tmp_data <- "/mnt/data/projects/spatial-rads/analysis/tmp"
if (!dir.exists(tmp_data)) dir.create(tmp_data, recursive = TRUE)
Sys.setenv(TMPDIR = tmp_data, TMP = tmp_data, TEMP = tmp_data)

# ============================================================
# Step 1: Build reference profile from Mutter_01 labeled cells
# ============================================================
cat("Loading Mutter_01...\n")
m01 <- readRDS("/mnt/data/projects/spatial-rads/analysis/objects/seurat_clustered.rds")

label_col <- "ImmuneAtlas_ImmGen_Main_cell_Types"
labels <- m01@meta.data[[label_col]]
cat(sprintf("Mutter_01 unique cell types: %d\n", length(unique(labels))))
cat("Per-type counts (top 15):\n")
print(head(sort(table(labels), decreasing=TRUE), 15))

# Build reference: gene x celltype = mean(counts of cells with that label)
counts_m01 <- GetAssayData(m01, layer = "counts")
genes <- rownames(counts_m01)

# Filter to types with >= 100 cells (drop ultra-rare types that won't have stable profiles)
type_counts <- table(labels)
keep_types <- names(type_counts[type_counts >= 100])
cat(sprintf("Types with >= 100 cells: %d (dropping %d rare types)\n",
            length(keep_types), length(type_counts) - length(keep_types)))

ref_profile <- sapply(keep_types, function(ct) {
  idx <- which(labels == ct & !is.na(labels))
  if (length(idx) == 0) return(rep(0, length(genes)))
  rowMeans(counts_m01[, idx, drop = FALSE])
})
rownames(ref_profile) <- genes
cat(sprintf("Reference profile built: %d genes x %d cell types\n",
            nrow(ref_profile), ncol(ref_profile)))

# Save for re-use
write.csv(ref_profile,
          file.path(OBJECTS_DIR, "ImmuneAtlas_ImmGen_derived_from_M01.csv"))

# Free Mutter_01 from memory
rm(m01, counts_m01, labels); invisible(gc(verbose = FALSE))

# ============================================================
# Step 2: Run InSituType on Mutter_02
# ============================================================
cat("\nLoading Mutter_02 typed RDS (already has marker-based labels — will replace)...\n")
m02 <- readRDS(file.path(OBJECTS_DIR, "seurat_mutter02_typed.rds"))
cat(sprintf("Mutter_02: %d cells\n", ncol(m02)))

counts_m02 <- t(as.matrix(GetAssayData(m02, layer = "counts")))
common_genes <- intersect(colnames(counts_m02), rownames(ref_profile))
cat(sprintf("Common genes M02 panel ∩ M01 reference: %d\n", length(common_genes)))
counts_m02 <- counts_m02[, common_genes]
ref_profile <- ref_profile[common_genes, ]

# Negative probes — fall back to small constant if not found
neg_count_col <- intersect(c("nCount_negprobes","nCount_NegativeProbes"), colnames(m02@meta.data))[1]
if (!is.na(neg_count_col)) {
  # Approx mean neg probe count per cell: nCount_negprobes / number_of_neg_probes_in_panel
  # Mutter_02 typically has ~50 neg probes. Use 50 as denominator if features col missing.
  neg_n <- 50
  neg_m02 <- m02@meta.data[[neg_count_col]] / neg_n
  neg_m02[is.na(neg_m02) | !is.finite(neg_m02)] <- 0.01
} else {
  cat("No neg probe count column — using small constant\n")
  neg_m02 <- rep(0.01, ncol(m02))
}
cat("M02 neg_mean range:", range(neg_m02), "\n")

cat("\nRunning insitutype supervised on Mutter_02...\n")
t0 <- Sys.time()
res <- insitutype(x = counts_m02,
                  neg = neg_m02,
                  reference_profiles = ref_profile,
                  bg = NULL,
                  cohort = NULL,
                  n_clusts = 0,
                  update_reference_profiles = FALSE,
                  n_phase1 = 1000,
                  n_phase2 = 5000,
                  n_phase3 = 20000,
                  n_starts = 1,
                  max_iters = 5)
cat(sprintf("InSituType done: %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

m02$insitutype_label <- res$clust
m02$insitutype_prob  <- res$prob

cat("\nMutter_02 cell type counts (insitutype):\n")
print(head(sort(table(m02$insitutype_label), decreasing=TRUE), 20))

saveRDS(m02, file.path(OBJECTS_DIR, "seurat_mutter02_insitutype.rds"))

out_df <- tibble(cell_id = colnames(m02),
                 sample_id = m02$sample_id,
                 treatment = m02$treatment,
                 timepoint_h = m02$timepoint_h,
                 insitutype_label = m02$insitutype_label,
                 insitutype_prob  = m02$insitutype_prob)
write_tsv(out_df, file.path(DATA_DIR, "mutter02_insitutype_labels.tsv.gz"))

# Map insitutype labels to major-8 categories (same mapping as 00_load_and_filter.R)
ct_map <- c(
  "tumor_epithelial" = "tumor_epithelial",
  "Blood.endothelial" = "endothelial",
  "Lymphatic.endothelial" = "endothelial",
  "Fibroblastic.reticular" = "fibroblast",
  "Pericyte" = "smooth_muscle",
  "B.cell" = "B.cell",
  "Memory.B" = "B.cell",
  "GC_centroblasts" = "B.cell",
  "GC_centrocyes" = "B.cell",
  "Plasmablast" = "B.cell",
  "Spleen.CD19" = "B.cell",
  "Plasma" = "B.cell",
  "NK" = "NK",
  "ILC" = "NK",
  "CD8.T.cell" = "T.cell", "gdT" = "T.cell", "NKT" = "T.cell",
  "Spleen.CD4Act.48hrs" = "T.cell", "Spleen.LN.Naive.CD4" = "T.cell",
  "Spleen.Naive.CD4" = "T.cell", "Spleen.Naive.CD8" = "T.cell",
  "Spleen.Treg" = "T.cell", "Colon.Treg.Nrplo" = "T.cell",
  "Thymic.4.8.CD3lo.DP" = "T.cell", "Thymic.CD4SP" = "T.cell",
  "Thymic.CD8SP" = "T.cell", "Thymic.preT.DN1" = "T.cell",
  "ISP" = "T.cell", "DN2a" = "T.cell", "DN2b" = "T.cell",
  "DN3" = "T.cell", "DN4" = "T.cell",
  "Macrophage" = "myeloid", "Microglia" = "myeloid",
  "PerC.macrophage" = "myeloid", "small.peritoneal.macs" = "myeloid",
  "spleen.red.pulp.macs" = "myeloid", "BM.Neutrophil" = "myeloid",
  "Thio.induced.Peritoneal.Neutrophil" = "myeloid",
  "Ly6Chi.blood.monocytes" = "myeloid", "Ly6Clo.blood.monocytes" = "myeloid",
  "Dendritic" = "myeloid"
)
m02$cell_type_major_insitutype <- unname(ct_map[m02$insitutype_label])
m02$cell_type_major_insitutype[is.na(m02$cell_type_major_insitutype)] <- "other"
cat("\nMutter_02 major-8 mapping (insitutype-derived):\n")
print(table(m02$cell_type_major_insitutype))

# Compare insitutype-derived vs marker-gate-derived
cat("\nConcordance: marker-gate vs insitutype labels (Mutter_02):\n")
tab <- table(marker = m02$cell_type_major,
             insitutype = m02$cell_type_major_insitutype)
print(tab)
agree <- sum(diag(tab[intersect(rownames(tab),colnames(tab)),
                      intersect(rownames(tab),colnames(tab))]))
total <- sum(tab)
cat(sprintf("Marker vs insitutype agreement: %.1f%% (%d/%d)\n",
            100*agree/total, agree, total))

saveRDS(m02, file.path(OBJECTS_DIR, "seurat_mutter02_insitutype.rds"))

cat("\n=== 08d_insitutype_mutter02.R complete ===\n")
