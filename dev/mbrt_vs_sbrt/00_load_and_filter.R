## 00_load_and_filter.R
## Layer 2 Extended — comprehensive MBRT vs SBRT (4T1 flank)
## Loads Mutter_01 cached object (local disk), applies 3 filters, caches result.
##   Filter 1: model == "flank"            (drops tongue)
##   Filter 2: 4T1 tumor purity            (relabels low-purity tumor cells)
##   Filter 3: necrosis-zone density flag  (per slide; flagged for downstream exclusion)
## Run from dev/mbrt_vs_sbrt/.

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(data.table)
  library(RANN)
})

# --- Paths (LARGE objects on data disk; small TSVs/plots in repo) ---
ANALYSIS_DIR <- "/mnt/data/projects/spatial-rads/analysis"
OBJECTS_DIR  <- file.path(ANALYSIS_DIR, "objects", "mbrt_vs_sbrt")  # cached RDS go here (data disk)
DATA_DIR <- "data"   # small TSVs in repo (root partition is fine)
PLOT_DIR <- "plots"  # PNGs in repo
for (d in c(OBJECTS_DIR, DATA_DIR, PLOT_DIR)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# --- Load Mutter_01 cached object ---
cached_rds <- file.path(ANALYSIS_DIR, "objects", "seurat_clustered.rds")
stopifnot(file.exists(cached_rds))
obj <- readRDS(cached_rds)
cat(sprintf("Loaded Mutter_01: %d genes x %d cells\n", nrow(obj), ncol(obj)))

# --- Ensure parsed metadata (idempotent; matches dev/peak_valley_analysis/00_load_data.R) ---
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

n_pre <- ncol(obj)
cat(sprintf("Pre-filter: %d cells\n", n_pre))
cat("Pre-filter Condition counts:\n"); print(table(obj$Condition))

# ============================================================
# FILTER 1: 4T1 flank only
# ============================================================
cat("\n=== Filter 1: model == 'flank' ===\n")
obj <- subset(obj, subset = model == "flank")
cat(sprintf("After flank filter: %d cells (dropped %d tongue)\n", ncol(obj), n_pre - ncol(obj)))

# ============================================================
# FILTER 2: 4T1 tumor purity (matches pv_simple_dge.R:31-38)
#   epith_score = log2-normalized mean of Krt8/Krt18/Epcam/Cdh1
#   endo_score  = log2-normalized mean of Pecam1/Cdh5/Vwf
#   keep tumor_epithelial label only if epith_score > 0.3 & endo_score < 0.1
# ============================================================
cat("\n=== Filter 2: 4T1 tumor purity ===\n")
counts_mat <- GetAssayData(obj, layer = "counts")
lib <- Matrix::colSums(counts_mat); lib[lib == 0] <- 1
epith_genes <- intersect(c("Krt8", "Krt18", "Epcam", "Cdh1"), rownames(counts_mat))
endo_genes  <- intersect(c("Pecam1", "Cdh5", "Vwf"), rownames(counts_mat))
cat(sprintf("Epithelial markers in panel (%d/4): %s\n", length(epith_genes), paste(epith_genes, collapse = ",")))
cat(sprintf("Endothelial markers in panel (%d/3): %s\n", length(endo_genes), paste(endo_genes, collapse = ",")))

normed_epi <- log2(t(t(as.matrix(counts_mat[epith_genes, , drop = FALSE])) / lib) * 1e4 + 1)
normed_end <- log2(t(t(as.matrix(counts_mat[endo_genes,  , drop = FALSE])) / lib) * 1e4 + 1)
obj$epith_score <- colMeans(normed_epi)
obj$endo_score  <- colMeans(normed_end)
obj$is_4t1_pure <- obj$epith_score > 0.3 & obj$endo_score < 0.1

obj$cell_type_strict <- as.character(obj$cell_type_validated)
n_tumor_pre <- sum(obj$cell_type_strict == "tumor_epithelial", na.rm = TRUE)
mask <- obj$cell_type_strict == "tumor_epithelial" & !obj$is_4t1_pure
obj$cell_type_strict[mask] <- "tumor_uncertain"
n_tumor_post <- sum(obj$cell_type_strict == "tumor_epithelial", na.rm = TRUE)
cat(sprintf("Tumor_epithelial purity: %d / %d = %.1f%% retained\n",
            n_tumor_post, n_tumor_pre, 100 * n_tumor_post / n_tumor_pre))

# Free memory
rm(counts_mat, normed_epi, normed_end); invisible(gc(verbose = FALSE))

# ============================================================
# FILTER 3: necrosis-zone density flag (per slide)
# ============================================================
cat("\n=== Filter 3: necrosis-zone density flag ===\n")
slide_col <- if ("Slide" %in% colnames(obj@meta.data)) "Slide" else "slide_ID_numeric"
cat(sprintf("Slide column: %s\n", slide_col))
slides <- unique(obj@meta.data[[slide_col]])
cat(sprintf("Number of slides: %d\n", length(slides)))

obj$local_density_d20 <- NA_real_
obj$necrosis_zone <- NA
xy <- as.data.frame(obj@meta.data[, c("x_slide_mm", "y_slide_mm")])
slide_vec <- obj@meta.data[[slide_col]]

for (s in slides) {
  idx <- which(slide_vec == s)
  if (length(idx) < 50) next
  coords <- as.matrix(xy[idx, , drop = FALSE])
  k_use <- min(21, length(idx))
  nn <- RANN::nn2(coords, k = k_use)
  d20 <- rowMeans(nn$nn.dists[, -1, drop = FALSE])
  obj$local_density_d20[idx] <- d20
  thresh <- quantile(d20, 0.90, na.rm = TRUE)
  obj$necrosis_zone[idx] <- d20 > thresh
}
n_necr <- sum(obj$necrosis_zone, na.rm = TRUE)
cat(sprintf("Necrosis-zone cells: %d / %d = %.1f%% (top 10%% NN-distance per slide)\n",
            n_necr, ncol(obj), 100 * n_necr / ncol(obj)))

# Optional secondary criterion: low DAPI
if ("Mean.DAPI" %in% colnames(obj@meta.data)) {
  obj$necrosis_zone_strict <- NA
  for (s in slides) {
    idx <- which(slide_vec == s)
    if (length(idx) < 50) next
    dapi <- obj@meta.data$Mean.DAPI[idx]
    dapi_thresh <- quantile(dapi, 0.10, na.rm = TRUE)
    obj$necrosis_zone_strict[idx] <- obj$necrosis_zone[idx] & (dapi < dapi_thresh)
  }
  cat(sprintf("Necrosis-zone-strict cells (also low DAPI): %d / %d = %.1f%%\n",
              sum(obj$necrosis_zone_strict, na.rm = TRUE), ncol(obj),
              100 * sum(obj$necrosis_zone_strict, na.rm = TRUE) / ncol(obj)))
} else {
  cat("Mean.DAPI not in metadata — secondary criterion skipped.\n")
  obj$necrosis_zone_strict <- obj$necrosis_zone
}

# ============================================================
# Map granular cell_type_strict labels to 9 major categories
# ============================================================
cat("\n=== Cell-type major-class mapping ===\n")
ct_map <- c(
  "tumor_epithelial" = "tumor_epithelial",
  "tumor_uncertain"  = "tumor_uncertain",
  "Blood.endothelial" = "endothelial",
  "Lymphatic.endothelial" = "endothelial",
  "Fibroblastic.reticular" = "fibroblast",
  "Pericyte" = "smooth_muscle",
  "B.cell" = "B.cell",
  "Memory.B" = "B.cell",
  "GC_centroblasts" = "B.cell",
  "GC_centrocyes" = "B.cell",
  "Plasmablast" = "B.cell",   # plasmablast also folded under B-lineage
  "Spleen.CD19" = "B.cell",
  "Plasma" = "B.cell",   # plasma cells folded into B-lineage (fixed tissue, low cell counts)
  "NK" = "NK",
  "ILC" = "NK",
  "CD8.T.cell" = "T.cell",
  "gdT" = "T.cell",
  "NKT" = "T.cell",
  "Spleen.CD4Act.48hrs" = "T.cell",
  "Spleen.LN.Naive.CD4" = "T.cell",
  "Spleen.Naive.CD4" = "T.cell",
  "Spleen.Naive.CD8" = "T.cell",
  "Spleen.Treg" = "T.cell",
  "Colon.Treg.Nrplo" = "T.cell",
  "Thymic.4.8.CD3lo.DP" = "T.cell",
  "Thymic.CD4SP" = "T.cell",
  "Thymic.CD8SP" = "T.cell",
  "Thymic.preT.DN1" = "T.cell",
  "ISP" = "T.cell",
  "DN2a" = "T.cell",
  "DN2b" = "T.cell",
  "DN3" = "T.cell",
  "DN4" = "T.cell",
  "Macrophage" = "myeloid",
  "Microglia" = "myeloid",
  "PerC.macrophage" = "myeloid",
  "small.peritoneal.macs" = "myeloid",
  "spleen.red.pulp.macs" = "myeloid",
  "BM.Neutrophil" = "myeloid",
  "Thio.induced.Peritoneal.Neutrophil" = "myeloid",
  "Ly6Chi.blood.monocytes" = "myeloid",
  "Ly6Clo.blood.monocytes" = "myeloid",
  "Dendritic" = "myeloid"
)
cell_type_major_vec <- unname(ct_map[obj$cell_type_strict])
cell_type_major_vec[is.na(cell_type_major_vec)] <- "other"
obj$cell_type_major <- cell_type_major_vec
cat("cell_type_major counts:\n"); print(table(obj$cell_type_major))

# Persist universe (pre-major-class subset) for composition denominators
counts_universe <- as_tibble(obj@meta.data) %>%
  count(Condition, treatment, timepoint_h, cell_type_major)
write_tsv(counts_universe, file.path(DATA_DIR, "cell_counts_universe.tsv"))

# ============================================================
# Restrict to 9 major types for analytical scripts
# ============================================================
MAJOR_TYPES <- c("tumor_epithelial", "endothelial", "fibroblast", "smooth_muscle",
                 "T.cell", "B.cell", "myeloid", "NK")  # 8 types (plasma folded into B.cell)
obj_major <- subset(obj, subset = cell_type_major %in% MAJOR_TYPES)
cat(sprintf("\nAfter major-type filter: %d cells\n", ncol(obj_major)))

strata <- as_tibble(obj_major@meta.data) %>%
  count(Condition, treatment, timepoint_h, cell_type_major) %>%
  mutate(keep_for_deg = n >= 100)
strata_necr <- as_tibble(obj_major@meta.data) %>%
  filter(!is.na(necrosis_zone)) %>%
  group_by(Condition, treatment, timepoint_h, cell_type_major) %>%
  summarise(n_necrosis_zone = sum(necrosis_zone), .groups = "drop")
strata <- strata %>% left_join(strata_necr,
                               by = c("Condition", "treatment", "timepoint_h", "cell_type_major"))
write_tsv(strata, file.path(DATA_DIR, "cell_counts_by_strata.tsv"))
cat("Strata table written: data/cell_counts_by_strata.tsv\n")

# ============================================================
# Diagnostic: necrosis zone overlay (flank conditions only, post-major-filter)
# ============================================================
cat("\nGenerating necrosis-zone diagnostic overlay...\n")
diag_meta <- as_tibble(obj_major@meta.data) %>%
  filter(!is.na(necrosis_zone)) %>%
  mutate(necrosis_zone_chr = ifelse(necrosis_zone, "TRUE", "FALSE"),
         Condition = factor(Condition, levels = sort(unique(Condition))))
p_necr <- ggplot(diag_meta, aes(x = x_slide_mm, y = y_slide_mm, color = necrosis_zone_chr)) +
  geom_point(size = 0.05, alpha = 0.4) +
  scale_color_manual(values = c("FALSE" = "grey80", "TRUE" = "red3"),
                     name = "necrosis zone") +
  theme_void() +
  facet_wrap(~Condition, ncol = 4, scales = "free") +
  labs(title = "Necrosis-zone diagnostic (top 10% NN-distance per slide)",
       caption = "Sparse-region cells likely correspond to true necrosis or scar tissue") +
  theme(legend.position = "top",
        strip.text = element_text(size = 8),
        aspect.ratio = 1)
ggsave(file.path(PLOT_DIR, "necrosis_zone_overlay.png"),
       plot = p_necr, width = 16, height = 10, dpi = 120)

# ============================================================
# Cache filtered object on DATA DISK (not root partition)
# ============================================================
out_rds <- file.path(OBJECTS_DIR, "seurat_filtered.rds")
saveRDS(obj_major, out_rds)
cat(sprintf("\nCached filtered object: %s (%.1f MB)\n",
            out_rds, file.info(out_rds)$size / 1e6))

cat("\n=== 00_load_and_filter.R complete ===\n")
