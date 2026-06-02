## 08_set2_validation.R
## Validate Mutter_01 day2 findings against Mutter_02 (4 slides, 12 samples).
## Mutter_02 lacks cell type labels; assign marker-based then run same analyses.
## All operations on local disk — no GCS reads.

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(UCell)
  library(RANN)
  library(Matrix)
})

DATA_DIR <- "data"
PLOT_DIR <- "plots"
OBJECTS_DIR <- "/mnt/data/projects/spatial-rads/analysis/objects/mbrt_vs_sbrt"
if (!dir.exists(OBJECTS_DIR)) dir.create(OBJECTS_DIR, recursive = TRUE)

# Set TMPDIR to data disk so R scratch (merge intermediates) doesn't fill root partition
tmp_data <- "/mnt/data/projects/spatial-rads/analysis/tmp"
if (!dir.exists(tmp_data)) dir.create(tmp_data, recursive = TRUE)
Sys.setenv(TMPDIR = tmp_data, TMP = tmp_data, TEMP = tmp_data)

# ============================================================
# Load 4 Mutter_02 RDS files (LOCAL DATA DISK)
# ============================================================
m02_dir <- "/home/jeszyman/repos/spatial-rads/dev/peak_valley_analysis/data"
slide_files <- file.path(m02_dir, sprintf("seurat_mutter02_slide%d_qc.rds", 1:4))
stopifnot(all(file.exists(slide_files)))

cat("Loading Mutter_02 slides (4 RDS files, ~1 GB total)...\n")
slide_objs <- lapply(seq_along(slide_files), function(i) {
  cat(sprintf("  Loading slide %d...\n", i))
  o <- readRDS(slide_files[i])
  if (!"slide_idx" %in% colnames(o@meta.data)) o$slide_idx <- i
  o
})

cat(sprintf("Slide cell counts: %s\n",
            paste(sapply(slide_objs, ncol), collapse=", ")))

# Merge into single object
obj <- merge(slide_objs[[1]], slide_objs[-1])
rm(slide_objs); invisible(gc(verbose = FALSE))
cat(sprintf("Merged Mutter_02: %d genes x %d cells\n", nrow(obj), ncol(obj)))
cat(sprintf("Layers after merge: %s\n", paste(Layers(obj), collapse=", ")))

# For Seurat v5: JoinLayers consolidates split counts/data. For v4 Assay (this dataset), skip.
assay_class <- class(obj[["RNA"]])[1]
cat(sprintf("Assay class: %s\n", assay_class))
if (assay_class == "Assay5" && length(grep("^counts\\.", Layers(obj))) > 0) {
  cat("JoinLayers (Seurat v5 split-layer consolidation)...\n")
  t_join <- Sys.time()
  obj <- JoinLayers(obj)
  cat(sprintf("JoinLayers done in %.1f min. Layers now: %s\n",
              as.numeric(difftime(Sys.time(), t_join, units = "mins")),
              paste(Layers(obj), collapse=", ")))
} else {
  cat("Assay is v4 or already consolidated. Skipping JoinLayers.\n")
}

# Mutter_02 may already have data layer; if not, normalize
if (!"data" %in% Layers(obj)) {
  cat("Normalizing Mutter_02 (LogNormalize)...\n")
  obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                       scale.factor = 10000, verbose = FALSE)
}

# Verify treatment + timepoint columns
stopifnot("treatment" %in% colnames(obj@meta.data))
stopifnot("timepoint_h" %in% colnames(obj@meta.data))
cat("\nTreatment x timepoint_h:\n"); print(table(obj$treatment, obj$timepoint_h))

# Harmonize treatment label: spatial-rads project uses MBRT, SBRT, NT (Control)
# Mutter_02 may use lowercase or different — verify
cat("\nUnique treatment values:\n"); print(unique(obj$treatment))

# ============================================================
# Marker-based cell typing
# ============================================================
cat("\n=== Marker-based cell typing for Mutter_02 ===\n")
counts_mat <- GetAssayData(obj, layer = "counts")
lib <- Matrix::colSums(counts_mat); lib[lib == 0] <- 1

panel <- rownames(counts_mat)
need <- c("Krt8","Krt18","Epcam","Cdh1","Pecam1","Cdh5","Vwf",
          "Cd3e","Cd3d","Cd19","Ms4a1","Cd68","Adgre1","Itgam",
          "Ncr1","Klrb1c","Col1a1","Pdgfra","Acta2","Myh11",
          "Jchain","Mzb1")
present <- intersect(need, panel)
cat(sprintf("Markers in Mutter_02 panel (%d/%d): %s\n",
            length(present), length(need), paste(present, collapse=",")))

# Compute log2-normalized scores for marker groups
score_group <- function(mat, lib, genes) {
  g <- intersect(genes, rownames(mat))
  if (length(g) == 0) return(rep(0, ncol(mat)))
  norm <- log2(t(t(as.matrix(mat[g, , drop = FALSE])) / lib) * 1e4 + 1)
  if (length(g) == 1) as.numeric(norm) else colMeans(norm)
}

epi_s   <- score_group(counts_mat, lib, c("Krt8","Krt18","Epcam","Cdh1"))
endo_s  <- score_group(counts_mat, lib, c("Pecam1","Cdh5","Vwf"))
t_s     <- score_group(counts_mat, lib, c("Cd3e","Cd3d"))
b_s     <- score_group(counts_mat, lib, c("Cd19","Ms4a1"))
mye_s   <- score_group(counts_mat, lib, c("Cd68","Adgre1","Itgam"))
nk_s    <- score_group(counts_mat, lib, c("Ncr1","Klrb1c"))
fib_s   <- score_group(counts_mat, lib, c("Col1a1","Pdgfra"))
sm_s    <- score_group(counts_mat, lib, c("Acta2","Myh11"))
plasma_s<- score_group(counts_mat, lib, c("Jchain","Mzb1"))

# Cell type assignment in priority order
cell_type_strict <- rep("unassigned", ncol(obj))
# 4T1 first (epithelial-positive AND endothelial-negative)
is_4t1 <- epi_s > 0.3 & endo_s < 0.1
cell_type_strict[is_4t1] <- "tumor_epithelial"
# Then immune (in priority: T, B, myeloid, NK, plasma)
unassign <- cell_type_strict == "unassigned"
cell_type_strict[unassign & t_s > 0.3]      <- "T.cell"
cell_type_strict[cell_type_strict == "unassigned" & b_s > 0.3]      <- "B.cell"
cell_type_strict[cell_type_strict == "unassigned" & mye_s > 0.3]    <- "myeloid"
cell_type_strict[cell_type_strict == "unassigned" & nk_s > 0.3]     <- "NK"
cell_type_strict[cell_type_strict == "unassigned" & plasma_s > 0.3] <- "plasma"
# Stromal
cell_type_strict[cell_type_strict == "unassigned" & endo_s > 0.3 & !is_4t1] <- "endothelial"
cell_type_strict[cell_type_strict == "unassigned" & fib_s > 0.3]    <- "fibroblast"
cell_type_strict[cell_type_strict == "unassigned" & sm_s > 0.3]     <- "smooth_muscle"

obj$cell_type_strict <- cell_type_strict
obj$cell_type_major <- cell_type_strict
obj$epith_score <- epi_s; obj$endo_score <- endo_s
obj$is_4t1_pure <- is_4t1
cat("Cell type counts:\n"); print(table(obj$cell_type_major))
cat(sprintf("Coverage (assigned / total): %.1f%%\n",
            100 * mean(obj$cell_type_major != "unassigned")))

rm(counts_mat); invisible(gc(verbose = FALSE))

# ============================================================
# Necrosis-zone flag per slide
# ============================================================
cat("\n=== Necrosis-zone flag per slide ===\n")
slide_col <- "slide"
if (!slide_col %in% colnames(obj@meta.data)) {
  slide_col <- if ("slide_idx" %in% colnames(obj@meta.data)) "slide_idx" else "slide_ID_numeric"
}
cat(sprintf("Slide column: %s\n", slide_col))
slides <- unique(obj@meta.data[[slide_col]])
cat(sprintf("Slides: %d\n", length(slides)))

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
  obj$necrosis_zone[idx] <- d20 > quantile(d20, 0.90, na.rm = TRUE)
}
cat(sprintf("Necrosis-zone cells: %.1f%%\n",
            100 * mean(obj$necrosis_zone, na.rm = TRUE)))

# ============================================================
# Restrict to 9 major types (drop unassigned for analyses)
# ============================================================
MAJOR_TYPES <- c("tumor_epithelial", "endothelial", "fibroblast", "smooth_muscle",
                 "T.cell", "B.cell", "myeloid", "NK", "plasma")
obj_major <- subset(obj, subset = cell_type_major %in% MAJOR_TYPES)
cat(sprintf("\nAfter major-type filter: %d cells\n", ncol(obj_major)))

# Save typed object for re-use (data disk)
saveRDS(obj_major, file.path(OBJECTS_DIR, "seurat_mutter02_typed.rds"))

# ============================================================
# DEGs at day2: MBRT vs SBRT, MBRT vs Control, SBRT vs Control
# ============================================================
cat("\n=== Mutter_02 DEGs (day2 only) ===\n")

# Establish a Condition label compatible with Mutter_01: e.g., "MBRT_day2", "SBRT_day2", "Control"
# Mutter_02 has timepoint_h column. We need to map treatment values too.
treatments <- unique(obj_major$treatment)
cat(sprintf("Treatments: %s\n", paste(treatments, collapse=", ")))
# Use string matching to harmonize
obj_major$treatment_clean <- case_when(
  toupper(obj_major$treatment) %in% c("MBRT") ~ "MBRT",
  toupper(obj_major$treatment) %in% c("SBRT") ~ "SBRT",
  toupper(obj_major$treatment) %in% c("CONTROL","CTRL","NT") ~ "NT",
  TRUE ~ as.character(obj_major$treatment)
)
cat("treatment_clean distribution:\n"); print(table(obj_major$treatment_clean))
obj_major$Condition_m02 <- ifelse(obj_major$treatment_clean == "NT",
                                  "Control",
                                  paste0(obj_major$treatment_clean, "_day2"))

# Apply default necrosis exclusion at day2 for DEGs
keep_main <- !(obj_major$timepoint_h == 48 & obj_major$necrosis_zone %in% TRUE)
obj_main <- subset(obj_major, cells = colnames(obj_major)[keep_main])
cat(sprintf("After day2 necrosis-zone exclusion: %d cells\n", ncol(obj_main)))

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

day2_contrasts <- list(
  list(name = "MBRT_vs_Control", ident.1 = "MBRT_day2", ident.2 = "Control"),
  list(name = "SBRT_vs_Control", ident.1 = "SBRT_day2", ident.2 = "Control"),
  list(name = "MBRT_vs_SBRT",    ident.1 = "MBRT_day2", ident.2 = "SBRT_day2")
)

m02_strata <- as_tibble(obj_main@meta.data) %>%
  count(Condition_m02, cell_type_major)

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
  # bulk too
  deg_jobs[[length(deg_jobs) + 1]] <- c(ce, list(cell_type = "BULK", n1 = NA, n2 = NA))
}
cat(sprintf("Total DEG jobs: %d\n", length(deg_jobs)))

m02_degs <- list()
t0 <- Sys.time()
for (i in seq_along(deg_jobs)) {
  job <- deg_jobs[[i]]
  ct <- job$cell_type
  cells_keep <- if (ct == "BULK") {
    colnames(obj_main)[obj_main$Condition_m02 %in% c(job$ident.1, job$ident.2)]
  } else {
    colnames(obj_main)[obj_main$Condition_m02 %in% c(job$ident.1, job$ident.2) &
                       obj_main$cell_type_major == ct]
  }
  if (length(cells_keep) < 50) next
  obj_sub <- subset(obj_main, cells = cells_keep)
  d <- run_deg(obj_sub, job$ident.1, job$ident.2, "Condition_m02")
  if (!is.null(d)) {
    d$contrast <- job$name
    d$cell_type_major <- ct
    d$ident.1 <- job$ident.1; d$ident.2 <- job$ident.2
    m02_degs[[length(m02_degs) + 1]] <- d
  }
  if (i %% 5 == 0) cat(sprintf("  [%d/%d] elapsed %.0fs\n", i, length(deg_jobs),
                               as.numeric(difftime(Sys.time(), t0, units="secs"))))
}
m02_degs_df <- bind_rows(m02_degs) %>%
  rename(log2FC = avg_log2FC) %>%
  arrange(contrast, cell_type_major, desc(abs(log2FC)))
write_tsv(m02_degs_df, file.path(DATA_DIR, "mutter02_degs_day2.tsv"))
cat(sprintf("Wrote %d DEG rows\n", nrow(m02_degs_df)))

# ============================================================
# Concordance with Mutter_01 day2
# ============================================================
cat("\n=== Concordance with Mutter_01 day2 ===\n")
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

# Spearman by (contrast x cell_type_major)
rho_table <- concord %>%
  group_by(contrast, cell_type_major) %>%
  summarise(rho = if (sum(!is.na(log2FC_m01) & !is.na(log2FC_m02)) >= 10)
              cor(log2FC_m01, log2FC_m02, method = "spearman", use = "pairwise") else NA_real_,
            n_genes = n(),
            .groups = "drop") %>%
  arrange(contrast, desc(rho))
write_tsv(rho_table, file.path(DATA_DIR, "set2_concordance_rho.tsv"))
cat("Concordance rho table:\n"); print(rho_table)

# Plots
TREATMENT_COLORS <- c("MBRT" = "#E74C3C", "SBRT" = "#3498DB", "NT" = "#7F8C8D")

# Concordance scatter (MBRT vs SBRT only)
p_conc <- concord %>%
  filter(contrast == "MBRT_vs_SBRT", cell_type_major != "BULK") %>%
  ggplot(aes(x = log2FC_m01, y = log2FC_m02)) +
  geom_point(alpha = 0.4, size = 0.6) +
  geom_abline(slope = 1, intercept = 0, color = "red3", linetype = "dashed") +
  geom_hline(yintercept = 0, color = "grey60") +
  geom_vline(xintercept = 0, color = "grey60") +
  facet_wrap(~cell_type_major, ncol = 3) +
  labs(title = "Mutter_01 day2 vs Mutter_02 day2 — MBRT vs SBRT log2FC",
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

# Composition concordance: MBRT vs SBRT proportion deltas at day2
m02_comp <- as_tibble(obj_major@meta.data) %>%
  filter(timepoint_h == 48, treatment_clean %in% c("MBRT","SBRT")) %>%
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

cat("\nComposition concordance day2:\n"); print(comp_concord)

cat("\n=== 08_set2_validation.R complete ===\n")
