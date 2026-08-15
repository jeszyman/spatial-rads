#!/usr/bin/env Rscript
# Genome-wide per-cell smiDE differential expression, following the smiDE
# authors' published protocol (Vasconcelos et al., Genome Biology 2026, and the
# package vignettes at the installed commit).
#
# Protocol points implemented here:
#   1. Genes are PREFILTERED per cell type on the overlap-ratio metric (keep
#      ratio < 1). Covariate adjustment alone is the weaker half of the
#      authors' filter-and-adjust recommendation.
#   2. assay_matrix is the FULL cell population of the cohort; metadata is the
#      cell-type subset under analysis. Passing only the subset trips a package
#      warning and truncates the neighbor lookup.
#   3. pre_de(adjacencies_only = TRUE) once over the cohort; neighbor expression
#      is measured on the fly inside smi_de(). adjacencies_only = FALSE is
#      deprecated and materializes a per-cell-type neighbor-expression matrix
#      for every gene.
#   4. Inverse-distance neighbor weighting (the `weight` column that
#      fast_make_all_neighbors() puts on the adjacency table).
#   5. Neighbor expression is totalcount-normalized against a cohort-wide scale
#      factor, so the covariate is not a per-cell depth proxy.
#   6/7. Two model modes, selected by the `mode` argument:
#      screen  -- NB GLMM (nebula) with a sample random intercept:
#                 count ~ RankNorm(otherct_expr) + offset(log(totalcounts))
#                         + condition + (1|sample_id)
#                 Fast, genome-wide, and anti-conservative for a contrast that
#                 varies between samples. Every row is tagged mode="screen" so
#                 the results tier can hold it out of confirmatory inference.
#      spatial -- one fit per spatial unit with a Matern Gaussian-process random
#                 effect (spaMM) and no sample random intercept, combined across
#                 units by inverse-variance fixed-effect meta-analysis. The
#                 spatial random effect is the authors' fix for the inflated
#                 type I error of a per-cell model on spatially autocorrelated
#                 data.
#
# Spatial fitting unit: the authors fit per sample and meta-analyze. In this
# study `condition` is constant within sample_id (one tumor, one arm), so an arm
# contrast is not estimable inside a sample; it is estimable inside a slide,
# which carries a Control/MBRT/SBRT block set. The unit is therefore chosen as
# the finest metadata column in which the contrast levels co-occur, recorded in
# the `fit_unit` output column.
#
# Output is the engine schema (comparison/contrast/unit/feature_type/feature_id/
# estimate/se/df/stat/p) plus mode/sample_id/fit_unit/n_fits. `unit` is the cell
# type (engine convention); `sample_id` is the spatial fitting unit's value and
# is NA for screen and meta rows.
#
# Cohort membership is the cohort_samples.tsv whitelist intersected with the
# registry's contrast condition levels, the same cross-dataset-leakage guard
# used by scripts/engines/lm_engine.R.
#
# Args: <merged.rds> <full_labels.parquet> <coords_necrosis.parquet>
#       <obs.parquet> <comparisons.tsv> <samples.tsv> <overlap_ratio_qc.tsv>
#       <cohort> <mode> <out_tsv>
suppressPackageStartupMessages({
  library(data.table); library(arrow); library(Seurat); library(Matrix); library(smiDE)
})
a <- commandArgs(trailingOnly = TRUE)
if (length(a) != 10) {
  stop("usage: smide_de.R <merged.rds> <labels.parquet> <coords.parquet> <obs.parquet> ",
       "<comparisons.tsv> <samples.tsv> <overlap_ratio_qc.tsv> <cohort> <mode> <out.tsv>")
}
rds_path <- a[1]; labels_path <- a[2]; coords_path <- a[3]; obs_path <- a[4]
comp_path <- a[5]; samples_path <- a[6]; overlap_path <- a[7]
cohort_name <- a[8]; de_mode <- a[9]; out_tsv <- a[10]
stopifnot(`mode must be "screen" or "spatial"` = de_mode %in% c("screen", "spatial"))

RADIUS_MM <- 0.05
MIN_CELLS <- 50
# Fewer surviving genes than this leaves nothing worth a per-cell-type BH family.
MIN_TARGETS <- 25
# Overlap ratio >= 1 means neighbor-other-type expression matches or exceeds the
# gene's own expression in that cell type -- the authors' prefilter cut.
ORM_MAX <- 1
N_CORES <- 16
# Matern spatial clusters as a fraction of the cells in a fit; the spaMM fit
# cost grows steeply in the resulting number of levels.
SPATIAL_K_PROP <- 0.25
# Candidate fitting units for spatial mode, finest first.
FIT_UNIT_CANDIDATES <- c("sample_id", "slide_id")

dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)

OUT_COLS <- c("comparison", "contrast", "unit", "feature_type", "feature_id",
              "estimate", "se", "df", "stat", "p", "mode", "sample_id",
              "fit_unit", "n_fits")

EMPTY_OUT <- function() {
  data.table(comparison=character(), contrast=character(), unit=character(),
             feature_type=character(), feature_id=character(),
             estimate=numeric(), se=numeric(), df=numeric(),
             stat=numeric(), p=numeric(), mode=character(),
             sample_id=character(), fit_unit=character(), n_fits=integer())
}

write_out <- function(dt) {
  fwrite(dt[, ..OUT_COLS], out_tsv, sep = "\t")
}

comp <- fread(comp_path)
cohort_comp <- comp[cohort == cohort_name & !is.na(contrast_num_level)]
if (nrow(cohort_comp) == 0) {
  cat(sprintf("smide_de: no %s comparisons found\n", cohort_name))
  write_out(EMPTY_OUT())
  quit(save = "no")
}
contrast_defs <- unique(cohort_comp[, .(contrast = name, num_level = contrast_num_level,
                                        den_level = contrast_den_level)])
conditions <- unique(c(cohort_comp$contrast_num_level, cohort_comp$contrast_den_level))

## cohort sample whitelist (prevents cross-dataset condition-name leakage --
## e.g. Mutter_01's Control/MBRT_day2/SBRT_day2 samples sharing condition names
## with mutter02_day2); same lookup as scripts/engines/lm_engine.R lines 36-40.
cohort_samples_path <- file.path(dirname(comp_path), "cohort_samples.tsv")
COHORT_SAMPLES <- if (file.exists(cohort_samples_path)) {
  fread(cohort_samples_path)[cohort == cohort_name, sample_id]
} else character()

# =============================================================================
# SECTION: LOAD + ALIGN
# =============================================================================

seu <- readRDS(rds_path)
counts <- GetAssayData(seu, assay = "RNA", layer = "counts")
lab <- as.data.table(read_parquet(labels_path))[, .(cell, cell_subtype)]
co  <- as.data.table(read_parquet(coords_path))[, .(cell, sample_id, x_slide_mm, y_slide_mm)]
ob  <- as.data.table(read_parquet(obs_path))[, .(cell, slide_id)]
meta <- lab[co, on = "cell"]
meta <- ob[meta, on = "cell"]
meta <- meta[!is.na(cell_subtype)]

## samples.tsv is the single source of truth for sample-level metadata; join
## condition from it rather than trust obs.parquet's baked-in copy (same
## pattern as pathway_scores.R lines 51-63).
ss <- fread(samples_path)
smeta <- unique(ss[, .(sample_id, condition)])
meta[, condition := smeta$condition[match(sample_id, smeta$sample_id)]]
stopifnot(!anyNA(meta$condition))

## cohort membership = sample whitelist AND registry condition levels. The
## surviving cells are the "full population" for smi_de(): all cell types of the
## cohort's samples. Cells outside the cohort are dropped rather than kept as
## potential neighbors because the neighbor graph never crosses a sample
## boundary (split_neighbors_by_colname = "sample_id").
meta <- meta[sample_id %in% COHORT_SAMPLES & condition %in% conditions]

cells_keep <- intersect(colnames(counts), meta$cell)
counts <- counts[, cells_keep]
meta <- meta[cell %in% cells_keep]
setkey(meta, cell)
meta <- meta[colnames(counts)]
meta[, totalcounts := Matrix::colSums(counts)]

genes <- rownames(counts)
subtypes <- sort(unique(meta$cell_subtype))
sample_ids <- unique(meta$sample_id)

cat(sprintf("smide_de: cohort=%s mode=%s | %d cells, %d genes, %d subtypes, %d samples\n",
            cohort_name, de_mode, ncol(counts), length(genes), length(subtypes),
            length(sample_ids)))

if (nrow(meta) == 0) {
  cat("smide_de: no cells after cohort filtering\n")
  write_out(EMPTY_OUT())
  quit(save = "no")
}

# =============================================================================
# SECTION: OVERLAP-RATIO GENE PREFILTER
# =============================================================================

# The authors filter AND adjust: genes whose neighbor-other-type expression
# matches their own expression in a cell type carry a segmentation-bleed signal
# the covariate cannot fully separate, so they are dropped from that cell type's
# target list rather than tested and caveated.
orm <- fread(overlap_path)
stopifnot(all(c("gene", "cell_subtype", "ratio") %in% names(orm)))
orm_pass <- orm[is.finite(ratio) & ratio < ORM_MAX]
targets_by_ct <- split(orm_pass$gene, orm_pass$cell_subtype)
targets_by_ct <- lapply(targets_by_ct, function(g) intersect(genes, unique(g)))

# =============================================================================
# SECTION: NEIGHBOR ADJACENCIES (once over the cohort)
# =============================================================================

# adjacencies_only = TRUE is the supported path: pre_de() returns only the
# cell-cell adjacency table (with per-edge inverse-distance `weight`), and
# smi_de() measures each target gene's neighbor expression on the fly. The
# deprecated adjacencies_only = FALSE path precomputes a neighbor-expression
# matrix per cell type over all genes, which is what forced the previous
# per-sample chunking here.
pde <- pre_de(
  adjacencies_only = TRUE,
  metadata = meta,
  ref_celltype = subtypes,
  cell_type_metadata_colname = "cell_subtype",
  cellid_colname = "cell",
  sdimx_colname = "x_slide_mm",
  sdimy_colname = "y_slide_mm",
  split_neighbors_by_colname = "sample_id",
  mm_radius = RADIUS_MM,
  weight_colname = "weight"
)

# weight = 1/distance, so two cells sharing a centroid coordinate produce a
# non-finite edge weight that would dominate the neighbor covariate.
n_bad_weight <- pde$cell_adjacency_dt[from != to & !is.finite(weight), .N]
cat(sprintf("smide_de: %d adjacency edges, %d non-finite inverse-distance weights\n",
            nrow(pde$cell_adjacency_dt), n_bad_weight))

# Cohort-wide scale factor for the on-the-fly neighbor expression, matching
# smiDE's own convention (mean depth / cell depth, empty cells held at 1). Must
# be computed over the full population passed as assay_matrix, not per subset,
# or the covariate shifts by a per-subset constant.
tc <- Matrix::colSums(counts)
mean_tc <- mean(tc)
tc[tc == 0] <- 1
nb_scalefactor <- setNames(as.numeric(mean_tc / tc), colnames(counts))

# =============================================================================
# SECTION: FIT HELPERS
# =============================================================================

# emmeans returns response-scale ratios for a log-link fit; convert to log2 fold
# change and propagate the SE by the delta method.
map_pairwise <- function(de, ct, mode_tag, sample_tag, fit_unit_tag) {
  pw <- as.data.table(results(de, comparisons = "pairwise", variable = "condition")$pairwise)
  if (nrow(pw) == 0 || !all(c("target", "contrast", "ratio", "SE") %in% names(pw))) return(NULL)
  pw <- pw[!is.na(ratio) & !is.na(SE) & ratio > 0]
  if (nrow(pw) == 0) return(NULL)
  pw[, target := as.character(target)]
  pw[, c("levelA", "levelB") := tstrsplit(as.character(contrast), " / ", fixed = TRUE)]
  pw[, levelA := trimws(levelA)]
  pw[, levelB := trimws(levelB)]
  pw_df <- if ("df" %in% names(pw)) pw$df else rep(NA_real_, nrow(pw))

  mapped <- vector("list", nrow(contrast_defs))
  for (i in seq_len(nrow(contrast_defs))) {
    cdef <- contrast_defs[i]
    sel <- (pw$levelA == cdef$num_level & pw$levelB == cdef$den_level) |
           (pw$levelA == cdef$den_level & pw$levelB == cdef$num_level)
    if (!any(sel)) next
    hit <- pw[sel]
    hit_df <- pw_df[sel]
    dir_sign <- ifelse(hit$levelA == cdef$num_level, 1, -1)
    mapped[[i]] <- data.table(
      comparison = cohort_name,
      contrast = cdef$contrast,
      unit = ct,
      feature_type = "gene",
      feature_id = hit$target,
      estimate = dir_sign * log2(hit$ratio),
      se = hit$SE / (hit$ratio * log(2)),
      df = as.numeric(hit_df),
      stat = NA_real_,
      p = hit$p.value,
      mode = mode_tag,
      sample_id = sample_tag,
      fit_unit = fit_unit_tag,
      n_fits = NA_integer_
    )
  }
  mapped <- rbindlist(mapped, fill = TRUE)
  if (nrow(mapped) == 0) NULL else mapped
}

# Fixed-effect inverse-variance meta-analysis across the per-unit spatial fits.
# metafor is not installed in the spatial-rads env, and the closed form below is
# the same estimator metafor's rma(method="FE") returns.
meta_combine <- function(dt) {
  ok <- dt[is.finite(estimate) & is.finite(se) & se > 0]
  if (nrow(ok) == 0) return(NULL)
  m <- ok[, {
    w <- 1 / se^2
    est <- sum(w * estimate) / sum(w)
    sse <- sqrt(1 / sum(w))
    .(estimate = est, se = sse, n_fits = .N)
  }, by = .(comparison, contrast, unit, feature_type, feature_id)]
  m[, stat := estimate / se]
  m[, p := 2 * pnorm(-abs(stat))]
  m[, df := NA_real_]
  m[, mode := "spatial_meta"]
  m[, sample_id := NA_character_]
  m[, fit_unit := NA_character_]
  m[]
}

# =============================================================================
# SECTION: SCREEN MODE -- nebula NB GLMM with a sample random intercept
# =============================================================================

results_list <- list()

if (de_mode == "screen") {
  for (ct in subtypes) {
    ct_cells <- meta[cell_subtype == ct, cell]
    tg <- targets_by_ct[[ct]]
    if (is.null(tg)) tg <- character()
    cat(sprintf("smide_de: %s -- %d cells, %d/%d panel genes pass ORM ratio < %g\n",
                ct, length(ct_cells), length(tg), length(genes), ORM_MAX))
    if (length(ct_cells) < MIN_CELLS) {
      cat(sprintf("smide_de: SKIP %s (%d cells < %d)\n", ct, length(ct_cells), MIN_CELLS))
      next
    }
    if (length(tg) < MIN_TARGETS) {
      cat(sprintf("smide_de: SKIP %s (%d genes pass ORM filter < %d)\n",
                  ct, length(tg), MIN_TARGETS))
      next
    }
    if (uniqueN(meta[cell %in% ct_cells, sample_id]) < 2) {
      cat(sprintf("smide_de: SKIP %s (< 2 samples)\n", ct))
      next
    }

    ct_meta <- meta[cell %in% ct_cells]
    setkey(ct_meta, cell)
    ct_meta <- ct_meta[ct_cells]

    tryCatch({
      de <- smi_de(
        assay_matrix = counts,
        metadata = ct_meta,
        formula = ~ RankNorm(otherct_expr) + offset(log(totalcounts)) + condition + (1|sample_id),
        pre_de_obj = pde,
        groupVar = "condition",
        family = "nbinom2",
        targets = tg,
        cellid_colname = "cell",
        nCores = N_CORES,
        neighbor_expr_overlap_agg = "sum",
        neighbor_expr_overlap_weight_colname = "weight",
        neighbor_expr_totalcount_normalize = TRUE,
        neighbor_expr_totalcount_scalefactor = nb_scalefactor,
        neighbor_expr_cell_type_metadata_colname = "cell_subtype"
      )
      mapped <- map_pairwise(de, ct, "screen", NA_character_, NA_character_)
      if (!is.null(mapped)) results_list[[ct]] <- mapped
      cat(sprintf("smide_de: %s -- %d genes tested, %d contrast rows\n",
                  ct, length(tg), if (is.null(mapped)) 0L else nrow(mapped)))
    }, error = function(e) {
      cat(sprintf("smide_de: ERROR %s -- %s\n", ct, conditionMessage(e)))
    })
  }
}

# =============================================================================
# SECTION: SPATIAL MODE -- per-unit GP_Matern fits + inverse-variance meta
# =============================================================================

if (de_mode == "spatial") {
  # Finest metadata column in which the arm contrast is estimable.
  fit_unit <- NA_character_
  for (cand in FIT_UNIT_CANDIDATES) {
    if (!cand %in% names(meta)) next
    if (meta[, uniqueN(condition), by = cand][, any(V1 >= 2)]) { fit_unit <- cand; break }
  }
  if (is.na(fit_unit)) {
    cat("smide_de: no metadata column carries >= 2 conditions; spatial mode cannot fit\n")
    write_out(EMPTY_OUT())
    quit(save = "no")
  }
  cat(sprintf("smide_de: spatial fitting unit = %s\n", fit_unit))

  # xy_kmeans_clusters() names its centroid columns <coord>_cluster and
  # smi_de() never rewrites the default GP_Matern random-effect formula, which
  # references CosMx's sdimx/sdimy/Run_Tissue_name names. The formula is passed
  # explicitly so it matches the columns actually created.
  spatial_spec <- list(
    name = "GP_Matern",
    k_prop_n = SPATIAL_K_PROP,
    x_coord_col = "x_slide_mm",
    y_coord_col = "y_slide_mm",
    split_neighbors_by_colname = "sample_id",
    spatial_random_effect = ~ Matern(1 | x_slide_mm_cluster + y_slide_mm_cluster %in% sample_id)
  )

  unit_values <- sort(unique(meta[[fit_unit]]))
  for (ct in subtypes) {
    tg <- targets_by_ct[[ct]]
    if (is.null(tg)) tg <- character()
    ct_meta_all <- meta[cell_subtype == ct]
    cat(sprintf("smide_de: %s -- %d cells, %d/%d panel genes pass ORM ratio < %g\n",
                ct, nrow(ct_meta_all), length(tg), length(genes), ORM_MAX))
    if (length(tg) < MIN_TARGETS) {
      cat(sprintf("smide_de: SKIP %s (%d genes pass ORM filter < %d)\n",
                  ct, length(tg), MIN_TARGETS))
      next
    }

    per_unit <- list()
    for (u in unit_values) {
      u_meta <- ct_meta_all[get(fit_unit) == u]
      if (nrow(u_meta) < MIN_CELLS) {
        cat(sprintf("smide_de: SKIP %s / %s=%s (%d cells < %d)\n",
                    ct, fit_unit, u, nrow(u_meta), MIN_CELLS))
        next
      }
      if (uniqueN(u_meta$condition) < 2) {
        cat(sprintf("smide_de: SKIP %s / %s=%s (single condition)\n", ct, fit_unit, u))
        next
      }
      setkey(u_meta, cell)

      tryCatch({
        de <- smi_de(
          assay_matrix = counts,
          metadata = u_meta,
          formula = ~ RankNorm(otherct_expr) + offset(log(totalcounts)) + condition,
          pre_de_obj = pde,
          groupVar = "condition",
          family = "nbinom2",
          targets = tg,
          cellid_colname = "cell",
          nCores = N_CORES,
          spatial_model = spatial_spec,
          neighbor_expr_overlap_agg = "sum",
          neighbor_expr_overlap_weight_colname = "weight",
          neighbor_expr_totalcount_normalize = TRUE,
          neighbor_expr_totalcount_scalefactor = nb_scalefactor,
          neighbor_expr_cell_type_metadata_colname = "cell_subtype"
        )
        mapped <- map_pairwise(de, ct, "spatial", as.character(u), fit_unit)
        if (!is.null(mapped)) per_unit[[as.character(u)]] <- mapped
        cat(sprintf("smide_de: %s / %s=%s -- %d genes tested, %d contrast rows\n",
                    ct, fit_unit, u, length(tg),
                    if (is.null(mapped)) 0L else nrow(mapped)))
      }, error = function(e) {
        cat(sprintf("smide_de: ERROR %s / %s=%s -- %s\n", ct, fit_unit, u,
                    conditionMessage(e)))
      })
    }

    if (length(per_unit) == 0) next
    ct_rows <- rbindlist(per_unit, fill = TRUE)
    ct_meta_rows <- meta_combine(ct_rows)
    results_list[[ct]] <- rbindlist(list(ct_rows, ct_meta_rows), use.names = TRUE, fill = TRUE)
    cat(sprintf("smide_de: %s -- %d per-unit rows, %d meta rows over %d unit(s)\n",
                ct, nrow(ct_rows),
                if (is.null(ct_meta_rows)) 0L else nrow(ct_meta_rows), length(per_unit)))
  }
}

# =============================================================================
# SECTION: WRITE
# =============================================================================

if (length(results_list) == 0) {
  cat("smide_de: all cell types failed or skipped\n")
  write_out(EMPTY_OUT())
} else {
  out <- rbindlist(results_list, use.names = TRUE, fill = TRUE)
  write_out(out)
  cat(sprintf("smide_de: %d rows | %d genes x %d cell types x %d contrasts | modes: %s\n",
              nrow(out), uniqueN(out$feature_id), uniqueN(out$unit), uniqueN(out$contrast),
              paste(sort(unique(out$mode)), collapse = ", ")))
}
