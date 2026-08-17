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
#      spatial -- one fit per spatial unit with a Gaussian-process spatial
#                 random effect and no sample random intercept, combined across
#                 units by inverse-variance fixed-effect meta-analysis. The
#                 spatial random effect is the authors' fix for the inflated
#                 type I error of a per-cell model on spatially autocorrelated
#                 data.
#
# Spatial backend: GP_INLA (SPDE Gaussian process fitted by INLA) is the
# default. On this cohort it is ~12x faster than GP_Matern on the spatial term
# at matched fit size and completes mid-sized subtype fits (>10k cells) that
# GP_Matern cannot complete at all. GP_Matern stays reachable through
# SPATIAL_BACKEND for cross-checking. INLA returns posterior summaries, not
# frequentist tests: the p reported for a spatial row is a normal approximation
# to the posterior credible interval (z = estimate/se, two-sided), and the
# credible bounds themselves are carried in ci_lo/ci_hi.
#
# Spatial fitting unit: the authors fit per sample and meta-analyze. In this
# study `condition` is constant within sample_id (one tumor, one arm), so an arm
# contrast is not estimable inside a sample; it is estimable inside a slide,
# which carries a Control/MBRT/SBRT block set. The unit is therefore chosen as
# the finest metadata column in which the contrast levels co-occur, recorded in
# the `fit_unit` output column.
#
# Output is the engine schema (comparison/contrast/unit/feature_type/feature_id/
# estimate/se/df/stat/p) plus ci_lo/ci_hi/mode/sample_id/fit_unit/n_fits/
# n_cells_fit and the per-arm evidence behind the fit (n_expr_num/n_expr_den =
# expressing cells, n_cells_num/n_cells_den = cells, oriented to the registry
# contrast's numerator and denominator). `unit` is the cell type (engine
# convention); `sample_id` is the spatial fitting unit's value and is NA for
# screen and meta rows.
#
# Cohort membership is the cohort_samples.tsv whitelist intersected with the
# registry's contrast condition levels, the same cross-dataset-leakage guard
# used by scripts/engines/lm_engine.R.
#
# Args: <merged.rds|cohort counts cache> <full_labels.parquet>
#       <coords_necrosis.parquet> <obs.parquet> <comparisons.tsv> <samples.tsv>
#       <overlap_ratio_qc.tsv> <cohort> <mode> <out_tsv> <out_skipped_tsv>
suppressPackageStartupMessages({
  library(data.table); library(arrow); library(Seurat); library(Matrix); library(smiDE)
})
a <- commandArgs(trailingOnly = TRUE)
if (length(a) != 11) {
  stop("usage: smide_de.R <merged.rds|cohort.counts.rds> <labels.parquet> <coords.parquet> <obs.parquet> ",
       "<comparisons.tsv> <samples.tsv> <overlap_ratio_qc.tsv> <cohort> <mode> <out.tsv> <out_skipped.tsv>")
}
rds_path <- a[1]; labels_path <- a[2]; coords_path <- a[3]; obs_path <- a[4]
comp_path <- a[5]; samples_path <- a[6]; overlap_path <- a[7]
cohort_name <- a[8]; de_mode <- a[9]; out_tsv <- a[10]; out_skipped <- a[11]
stopifnot(`mode must be "screen" or "spatial"` = de_mode %in% c("screen", "spatial"))

RADIUS_MM <- 0.05
MIN_CELLS <- 50
# Minimum expressing (non-zero) cells an arm must contribute to a spatial fit for
# that fit's contrasts to be trusted. MIN_CELLS is a floor on cells present, which
# says nothing about whether a given gene was seen in each arm; an arm with no
# expressing cell makes the NB log-link coefficient unidentified and INLA answers
# with a Gaussian approximation centred wherever the optimiser stopped, at a
# finite and often small posterior sd, so the runaway estimate then wins the
# inverse-variance meta-analysis. Measured on mutter02_day2: with no floor, 21%
# of per-slide fits had an arm with zero expressing cells and 100% of the 226 rows
# with |log2FC| > 100 were of that kind; at a floor of 3 the implausible rows fall
# to 14 of 34,962 while 62% of rows are kept. The floor applies to EVERY arm in
# the fit, not only the two levels of the contrast being read off it: the arms
# share one linear predictor, so a single separated arm carries the whole fit off.
MIN_EXPR_PER_ARM <- 3
# Real effects on a 950-gene panel run about +/- 1 to 3 log2 units. Nothing above
# this is a biological result, so a surviving meta estimate past it means the gate
# above has failed to catch a degenerate fit and the run must not be believed.
MAX_PLAUSIBLE_LOG2FC <- 10
# Fewer surviving genes than this leaves nothing worth a per-cell-type BH family.
MIN_TARGETS <- 25
# Overlap ratio >= 1 means neighbor-other-type expression matches or exceeds the
# gene's own expression in that cell type -- the authors' prefilter cut.
ORM_MAX <- 1
# nebula parallelizes across genes and each worker is light.
N_CORES <- 16
# The GP backends show ~1.0x speedup from extra workers while each forked worker
# holds a 14-22 GB copy of the cohort adjacency machinery, so the spatial path
# runs serial inside a call and parallelizes across calls instead.
N_CORES_SPATIAL <- 1
# Spatial fits are capped at this many cells per fit unit, subsampled
# proportionally within sample_id so both contrast arms survive. The cap keeps
# every fit inside the size range that was timed directly; the neighbor
# covariate is unaffected because assay_matrix and the adjacency table remain
# the full cohort.
SPATIAL_MAX_CELLS <- 2500
SPATIAL_SUBSAMPLE_SEED <- 1
# Matern spatial clusters as a fraction of the cells in a fit, and their hard
# ceiling; the spaMM fit cost grows steeply in the resulting number of levels.
SPATIAL_K_PROP <- 0.25
K_MAX <- 700
SPATIAL_BACKEND <- "GP_INLA"
SPATIAL_BACKEND_MATERN <- "GP_Matern"
# Candidate fitting units for spatial mode, finest first.
FIT_UNIT_CANDIDATES <- c("sample_id", "slide_id")

dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)

OUT_COLS <- c("comparison", "contrast", "unit", "feature_type", "feature_id",
              "estimate", "se", "df", "stat", "p", "ci_lo", "ci_hi", "mode",
              "sample_id", "fit_unit", "n_fits", "n_cells_fit",
              "n_expr_num", "n_expr_den", "n_cells_num", "n_cells_den")

# Gated fits, written alongside the results so the gate is auditable without
# re-deriving expressing counts from the counts matrix.
SKIP_COLS <- c("comparison", "contrast", "unit", "feature_id", "sample_id",
               "fit_unit", "estimate", "se", "p", "n_expr_num", "n_expr_den",
               "n_cells_num", "n_cells_den", "min_expr_fit", "reason")

EMPTY_OUT <- function() {
  data.table(comparison=character(), contrast=character(), unit=character(),
             feature_type=character(), feature_id=character(),
             estimate=numeric(), se=numeric(), df=numeric(),
             stat=numeric(), p=numeric(), ci_lo=numeric(), ci_hi=numeric(),
             mode=character(), sample_id=character(), fit_unit=character(),
             n_fits=integer(), n_cells_fit=integer(),
             n_expr_num=integer(), n_expr_den=integer(),
             n_cells_num=integer(), n_cells_den=integer())
}

EMPTY_SKIPPED <- function() {
  data.table(comparison=character(), contrast=character(), unit=character(),
             feature_id=character(), sample_id=character(), fit_unit=character(),
             estimate=numeric(), se=numeric(), p=numeric(),
             n_expr_num=integer(), n_expr_den=integer(),
             n_cells_num=integer(), n_cells_den=integer(),
             min_expr_fit=integer(), reason=character())
}

skipped_list <- list()

write_skipped <- function() {
  dt <- if (length(skipped_list)) rbindlist(skipped_list, fill = TRUE) else EMPTY_SKIPPED()
  fwrite(dt[, ..SKIP_COLS], out_skipped, sep = "\t")
}

write_out <- function(dt) {
  fwrite(dt[, ..OUT_COLS], out_tsv, sep = "\t")
  write_skipped()
}

# An empty spatial table is never a valid result: a cohort is on the spatial run
# manifest because it is expected to produce fits, and a silent empty file is
# indistinguishable downstream from a genuine null.
abort_or_empty <- function(reason) {
  if (de_mode == "spatial")
    stop(sprintf("smide_de: cohort %s spatial mode produced no results -- %s",
                 cohort_name, reason))
  cat(sprintf("smide_de: %s\n", reason))
  write_out(EMPTY_OUT())
  quit(save = "no")
}

comp <- fread(comp_path)
cohort_comp <- comp[cohort == cohort_name & !is.na(contrast_num_level)]
if (nrow(cohort_comp) == 0) {
  abort_or_empty(sprintf("no %s comparisons found in the registry", cohort_name))
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

## The first argument is either the merged object (all 20 samples) or the cohort
## counts cache built by build_cohort_counts.R, detected by class. The cache holds
## only this cohort's cells, so the peak below is the cohort rather than the whole
## study; the merged path is kept working because it is what in-flight runs pass.
seu <- readRDS(rds_path)
counts <- if (inherits(seu, "cohort_counts_cache")) seu$counts else
  GetAssayData(seu, assay = "RNA", layer = "counts")
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
# The subset above is an independent matrix, so the loaded object can go. It is
# several GB per process and the spatial stage runs many processes at once against a
# fixed memory budget.
rm(seu); invisible(gc())
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
  abort_or_empty("no cells survived cohort filtering")
}

# =============================================================================
# SECTION: SPATIAL FITTING UNIT + BORROWED-CONTROL FEASIBILITY
# =============================================================================

# Resolved before any fitting so an infeasible cohort fails in seconds rather
# than after hours of fits that map to nothing.
fit_unit <- NA_character_
unit_values <- character()
if (de_mode == "spatial") {
  for (cand in FIT_UNIT_CANDIDATES) {
    if (!cand %in% names(meta)) next
    if (meta[, uniqueN(condition), by = cand][, any(V1 >= 2)]) { fit_unit <- cand; break }
  }
  if (is.na(fit_unit)) {
    abort_or_empty("no metadata column carries >= 2 conditions, so no arm contrast is estimable within a fitting unit")
  }
  cat(sprintf("smide_de: spatial fitting unit = %s\n", fit_unit))

  # A unit contributes only if some registry contrast has BOTH its levels
  # present in that unit; a unit holding two conditions that are never
  # contrasted against each other (e.g. a treated-only slide in a cohort whose
  # contrasts are all treated-vs-Control) fits fine and maps to zero rows.
  unit_conds <- meta[, .(conds = list(unique(condition))), by = c(fit_unit)]
  unit_conds[, estimable := vapply(conds, function(cc)
    any(contrast_defs[, num_level %in% cc & den_level %in% cc]), logical(1))]
  unit_values <- sort(unit_conds[estimable == TRUE][[fit_unit]])
  dead_units  <- sort(unit_conds[estimable == FALSE][[fit_unit]])

  # A dropped treated unit only costs power inside the units that remain. A
  # dropped CONTROL unit means the cohort's borrowed controls contributed
  # nothing at all, while the output still carries the pooled cohort's name --
  # a result labelled pooled that rests entirely on the co-resident units.
  if (length(dead_units) > 0) {
    dead_samples <- meta[get(fit_unit) %in% dead_units,
                         .(conds = paste(sort(unique(condition)), collapse = "/")),
                         by = sample_id]
    cat(sprintf("smide_de: %s %s carry no estimable contrast and are dropped (%s)\n",
                fit_unit, paste(dead_units, collapse = ", "),
                paste(sprintf("%s=%s", dead_samples$sample_id, dead_samples$conds),
                      collapse = ", ")))
    dropped_den <- meta[get(fit_unit) %in% dead_units &
                        condition %in% contrast_defs$den_level, unique(sample_id)]
    if (length(dropped_den) > 0)
      stop(sprintf(paste0(
        "smide_de: cohort %s cannot run in spatial mode. Reference-arm sample(s) %s sit on ",
        "%s %s, which carry no estimable contrast, so the spatial fits would use none of ",
        "the borrowed controls while the output was still labelled %s. Run this cohort in ",
        "screen mode only."),
        cohort_name, paste(sort(dropped_den), collapse = ", "), fit_unit,
        paste(dead_units, collapse = ", "), cohort_name))
  }
  if (length(unit_values) == 0)
    abort_or_empty(sprintf("no %s carries both levels of any registry contrast", fit_unit))
  cat(sprintf("smide_de: %d estimable %s value(s): %s\n", length(unit_values), fit_unit,
              paste(unit_values, collapse = ", ")))
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
# non-finite edge weight. Left in place it becomes NaN in the weighted neighbor
# sum and poisons the covariate for every gene. Capping at the largest finite
# edge weight keeps the (genuinely adjacent) pair in the graph while bounding
# its leverage; edges are dropped only if no finite weight exists to cap to.
adj <- pde$cell_adjacency_dt
off_diag <- adj[["from"]] != adj[["to"]]
bad_w <- off_diag & !is.finite(adj[["weight"]])
n_bad_weight <- sum(bad_w)
if (n_bad_weight > 0) {
  finite_w <- adj[["weight"]][off_diag & is.finite(adj[["weight"]])]
  if (length(finite_w) > 0) {
    w_cap <- max(finite_w)
    data.table::set(adj, which(bad_w), "weight", w_cap)
    cat(sprintf("smide_de: capped %d non-finite inverse-distance weight(s) at %.4g\n",
                n_bad_weight, w_cap))
  } else {
    pde$cell_adjacency_dt <- adj[!bad_w]
    cat(sprintf("smide_de: dropped %d non-finite-weight edge(s) (no finite weight to cap to)\n",
                n_bad_weight))
  }
}
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

# Two backends, two result schemas, one engine schema.
#
#   emmeans-backed (nebula, spaMM): response-scale `ratio` + `SE` for a log-link
#     fit. estimate = log2(ratio); SE propagated to the log2 scale by the delta
#     method, SE / (ratio * log(2)).
#   INLA: posterior summaries of the same log-link contrast. smiDE exponentiates
#     `mean`, `mode` and every `*quant` column for an nbinom2 fit but leaves `sd`
#     on the natural-log scale (postprocess_inla_contrasts), so estimate =
#     log2(mean) and se = sd / log(2) land on the identical log2 scale as the
#     emmeans branch. There is no p-value: the reported p is the two-sided
#     normal approximation to the credible interval, z = estimate / se, and the
#     credible bounds are carried through in ci_lo/ci_hi.
#
# `dir_sign` flips both the estimate and (with a swap) the interval when the
# package emits the contrast in the opposite order to the registry's.
map_pairwise <- function(de, ct, mode_tag, sample_tag, fit_unit_tag, n_cells) {
  pw <- as.data.table(results(de, comparisons = "pairwise", variable = "condition")$pairwise)
  if (nrow(pw) == 0)
    stop(sprintf("smide_de: results() returned zero pairwise rows for %s", ct))
  if (!all(c("target", "contrast") %in% names(pw)))
    stop(sprintf("smide_de: results() output for %s has no target/contrast columns (got: %s)",
                 ct, paste(names(pw), collapse = ", ")))

  quant_cols <- grep("quant$", names(pw), value = TRUE)
  is_emmeans <- all(c("ratio", "SE") %in% names(pw))
  is_inla    <- all(c("mean", "sd") %in% names(pw)) && length(quant_cols) >= 2
  if (!is_emmeans && !is_inla)
    stop(sprintf(paste0("smide_de: unrecognized results() schema for %s -- expected an emmeans ",
                        "(ratio/SE) or INLA (mean/sd/*quant) summary, got: %s"),
                 ct, paste(names(pw), collapse = ", ")))

  if (is_emmeans) {
    pw <- pw[is.finite(ratio) & is.finite(SE) & ratio > 0]
    if (nrow(pw) == 0) return(NULL)
    est_v <- log2(pw$ratio)
    se_v  <- pw$SE / (pw$ratio * log(2))
    lo_v  <- rep(NA_real_, nrow(pw))
    hi_v  <- rep(NA_real_, nrow(pw))
    df_v  <- if ("df" %in% names(pw)) as.numeric(pw$df) else rep(NA_real_, nrow(pw))
    p_v   <- pw$p.value
    stat_v <- rep(NA_real_, nrow(pw))
  } else {
    q_lev <- suppressWarnings(as.numeric(sub("quant$", "", quant_cols)))
    stopifnot(`INLA quantile column names are not numeric levels` = !anyNA(q_lev))
    lo_col <- quant_cols[which.min(q_lev)]
    hi_col <- quant_cols[which.max(q_lev)]
    pw <- pw[is.finite(mean) & is.finite(sd) & mean > 0 & sd > 0]
    if (nrow(pw) == 0) return(NULL)
    est_v <- log2(pw$mean)
    se_v  <- pw$sd / log(2)
    lo_v  <- log2(pw[[lo_col]])
    hi_v  <- log2(pw[[hi_col]])
    df_v  <- rep(NA_real_, nrow(pw))
    stat_v <- est_v / se_v
    p_v   <- 2 * pnorm(-abs(stat_v))
  }

  pw[, target := as.character(target)]
  pw[, c("levelA", "levelB") := tstrsplit(as.character(contrast), " / ", fixed = TRUE)]
  pw[, levelA := trimws(levelA)]
  pw[, levelB := trimws(levelB)]

  # smiDE attaches the evidence behind each contrast row: propnz is the fraction
  # of that arm's cells carrying a non-zero count (mean(y > 0) in summarize_model),
  # so propnz * ncells is the number of cells that actually informed the fit for
  # this gene. The _1 columns describe the left level of the "A / B" contrast
  # string and the _2 columns the right. Every arm of the fit shows up as levelA
  # or levelB of some pairwise row, so the long form gives the whole fit's
  # per-arm evidence, not just the two levels of one contrast.
  arm_cols <- c("propnz_1", "ncells_1", "propnz_2", "ncells_2")
  if (all(arm_cols %in% names(pw))) {
    pw[, n_expr_A  := as.integer(round(propnz_1 * ncells_1))]
    pw[, n_expr_B  := as.integer(round(propnz_2 * ncells_2))]
    pw[, n_cells_A := as.integer(ncells_1)]
    pw[, n_cells_B := as.integer(ncells_2)]
    arm_long <- rbind(pw[, .(target, n_expr = n_expr_A)],
                      pw[, .(target, n_expr = n_expr_B)])
    fit_min <- arm_long[, .(m = min(n_expr)), by = target]
    pw[, min_expr_fit := fit_min$m[match(target, fit_min$target)]]
  } else {
    pw[, c("n_expr_A", "n_expr_B", "n_cells_A", "n_cells_B",
           "min_expr_fit") := NA_integer_]
  }

  mapped <- vector("list", nrow(contrast_defs))
  for (i in seq_len(nrow(contrast_defs))) {
    cdef <- contrast_defs[i]
    sel <- (pw$levelA == cdef$num_level & pw$levelB == cdef$den_level) |
           (pw$levelA == cdef$den_level & pw$levelB == cdef$num_level)
    if (!any(sel)) next
    dir_sign <- ifelse(pw$levelA[sel] == cdef$num_level, 1, -1)
    fwd   <- dir_sign > 0
    est_i <- dir_sign * est_v[sel]
    lo_i  <- ifelse(dir_sign > 0, lo_v[sel], -hi_v[sel])
    hi_i  <- ifelse(dir_sign > 0, hi_v[sel], -lo_v[sel])
    mapped[[i]] <- data.table(
      comparison = cohort_name,
      contrast = cdef$contrast,
      unit = ct,
      feature_type = "gene",
      feature_id = pw$target[sel],
      estimate = est_i,
      se = se_v[sel],
      df = df_v[sel],
      stat = if (is_emmeans) stat_v[sel] else dir_sign * stat_v[sel],
      p = p_v[sel],
      ci_lo = lo_i,
      ci_hi = hi_i,
      mode = mode_tag,
      sample_id = sample_tag,
      fit_unit = fit_unit_tag,
      n_fits = NA_integer_,
      n_cells_fit = as.integer(n_cells),
      n_expr_num  = ifelse(fwd, pw$n_expr_A[sel],  pw$n_expr_B[sel]),
      n_expr_den  = ifelse(fwd, pw$n_expr_B[sel],  pw$n_expr_A[sel]),
      n_cells_num = ifelse(fwd, pw$n_cells_A[sel], pw$n_cells_B[sel]),
      n_cells_den = ifelse(fwd, pw$n_cells_B[sel], pw$n_cells_A[sel]),
      min_expr_fit = pw$min_expr_fit[sel]
    )
  }
  mapped <- rbindlist(mapped, fill = TRUE)
  if (nrow(mapped) == 0) NULL else mapped
}

# Fixed-effect inverse-variance meta-analysis across the per-unit spatial fits.
# metafor is not installed in the spatial-rads env, and the closed form below is
# the same estimator metafor's rma(method="FE") returns. n_cells_fit on a meta
# row is the total cells fitted across the contributing units. Its per-arm columns
# summarize the contributing fits: cells add up across units, while the expressing
# counts carry the weakest contributing fit, which is the binding evidence.
meta_combine <- function(dt, fit_unit_tag) {
  ok <- dt[is.finite(estimate) & is.finite(se) & se > 0]
  if (nrow(ok) == 0) return(NULL)
  m <- ok[, {
    w <- 1 / se^2
    est <- sum(w * estimate) / sum(w)
    sse <- sqrt(1 / sum(w))
    .(estimate = est, se = sse, n_fits = .N,
      n_cells_fit = sum(unique(data.table(sample_id, n_cells_fit))$n_cells_fit),
      n_expr_num = min(n_expr_num), n_expr_den = min(n_expr_den),
      n_cells_num = sum(n_cells_num), n_cells_den = sum(n_cells_den))
  }, by = .(comparison, contrast, unit, feature_type, feature_id)]
  m[, stat := estimate / se]
  m[, p := 2 * pnorm(-abs(stat))]
  m[, df := NA_real_]
  m[, ci_lo := NA_real_]
  m[, ci_hi := NA_real_]
  m[, mode := "spatial_meta"]
  m[, sample_id := NA_character_]
  m[, fit_unit := fit_unit_tag]
  m[]
}

# =============================================================================
# SECTION: SCREEN MODE -- nebula NB GLMM with a sample random intercept
# =============================================================================

results_list <- list()
n_fit_ok <- 0L
n_fit_fail <- 0L

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

    fit <- tryCatch({
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
      list(ok = TRUE,
           mapped = map_pairwise(de, ct, "screen", NA_character_, NA_character_,
                                 nrow(ct_meta)))
    }, error = function(e) list(ok = FALSE, msg = conditionMessage(e)))

    if (isTRUE(fit$ok)) {
      n_fit_ok <- n_fit_ok + 1L
      if (!is.null(fit$mapped)) results_list[[ct]] <- fit$mapped
      cat(sprintf("smide_de: %s -- %d genes tested, %d contrast rows\n",
                  ct, length(tg), if (is.null(fit$mapped)) 0L else nrow(fit$mapped)))
    } else {
      n_fit_fail <- n_fit_fail + 1L
      cat(sprintf("smide_de: ERROR %s -- %s\n", ct, fit$msg))
    }
  }
}

# =============================================================================
# SECTION: SPATIAL MODE -- per-unit GP fits + inverse-variance meta
# =============================================================================

if (de_mode == "spatial") {
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

    # INLA reports credible intervals rather than tests, so the vignette's
    # multiplicity control is a Bonferroni-widened interval over the gene list
    # actually fitted in this call.
    q_tail <- 0.025 / length(tg)

    per_unit <- list()
    for (u in unit_values) {
      u_meta <- ct_meta_all[get(fit_unit) == u]
      n_cells_unit <- nrow(u_meta)
      if (n_cells_unit < MIN_CELLS) {
        cat(sprintf("smide_de: SKIP %s / %s=%s (%d cells < %d)\n",
                    ct, fit_unit, u, n_cells_unit, MIN_CELLS))
        next
      }
      if (uniqueN(u_meta$condition) < 2) {
        cat(sprintf("smide_de: SKIP %s / %s=%s (single condition)\n", ct, fit_unit, u))
        next
      }
      if (n_cells_unit > SPATIAL_MAX_CELLS) {
        set.seed(SPATIAL_SUBSAMPLE_SEED)
        frac <- SPATIAL_MAX_CELLS / n_cells_unit
        u_meta <- u_meta[, .SD[sample(.N, max(1L, round(.N * frac)))], by = sample_id]
      }
      setkey(u_meta, cell)
      n_cells_fit <- nrow(u_meta)

      # xy_kmeans_clusters() takes k as a total and splits it proportionally
      # across sample_id; it warns and ignores k_prop_n whenever k is also
      # given, so the absolute count is what actually binds. Consumed by
      # GP_Matern only -- GP_INLA builds an SPDE mesh from the coordinates and
      # has no k-means step, and passing k there would leak into INLA::inla().
      k_abs <- min(floor(n_cells_fit * SPATIAL_K_PROP), K_MAX)
      spatial_spec <- if (SPATIAL_BACKEND == "GP_INLA") {
        list(
          name = "GP_INLA",
          x_coord_col = "x_slide_mm",
          y_coord_col = "y_slide_mm",
          quantiles = c(q_tail, 0.5, 1 - q_tail),
          num.threads = 1
        )
      } else {
        # xy_kmeans_clusters() names its centroid columns <coord>_cluster and
        # smi_de() never rewrites the default GP_Matern random-effect formula,
        # which references CosMx's sdimx/sdimy/Run_Tissue_name names. The
        # formula is passed explicitly so it matches the columns actually created.
        list(
          name = "GP_Matern",
          k = k_abs,
          k_prop_n = NULL,
          x_coord_col = "x_slide_mm",
          y_coord_col = "y_slide_mm",
          split_neighbors_by_colname = "sample_id",
          spatial_random_effect = ~ Matern(1 | x_slide_mm_cluster + y_slide_mm_cluster %in% sample_id)
        )
      }
      cat(sprintf("smide_de: FIT %s / %s=%s | %d of %d cells (seed %d) | %d genes | %s | k=%d\n",
                  ct, fit_unit, u, n_cells_fit, n_cells_unit, SPATIAL_SUBSAMPLE_SEED,
                  length(tg), SPATIAL_BACKEND, k_abs))

      fit <- tryCatch({
        de <- smi_de(
          assay_matrix = counts,
          metadata = u_meta,
          formula = ~ RankNorm(otherct_expr) + offset(log(totalcounts)) + condition,
          pre_de_obj = pde,
          groupVar = "condition",
          family = "nbinom2",
          targets = tg,
          cellid_colname = "cell",
          nCores = N_CORES_SPATIAL,
          spatial_model = spatial_spec,
          neighbor_expr_overlap_agg = "sum",
          neighbor_expr_overlap_weight_colname = "weight",
          neighbor_expr_totalcount_normalize = TRUE,
          neighbor_expr_totalcount_scalefactor = nb_scalefactor,
          neighbor_expr_cell_type_metadata_colname = "cell_subtype"
        )
        list(ok = TRUE,
             mapped = map_pairwise(de, ct, "spatial", as.character(u), fit_unit,
                                   n_cells_fit))
      }, error = function(e) list(ok = FALSE, msg = conditionMessage(e)))

      if (isTRUE(fit$ok)) {
        n_fit_ok <- n_fit_ok + 1L
        n_gated <- 0L
        if (!is.null(fit$mapped)) {
          mp <- fit$mapped
          if (anyNA(mp$min_expr_fit))
            stop(sprintf(paste0("smide_de: %s / %s=%s -- results() carried no per-arm ",
                                "counts (propnz_*/ncells_*), so the expressing-cell gate ",
                                "cannot be applied and a separated fit would go unchecked"),
                         ct, fit_unit, u))
          gate_pass <- mp$min_expr_fit >= MIN_EXPR_PER_ARM
          if (any(!gate_pass)) {
            g <- mp[!gate_pass]
            g[, reason := sprintf("min_expr_per_arm: weakest arm has %d expressing cell(s) in this fit (need %d)",
                                  min_expr_fit, MIN_EXPR_PER_ARM)]
            skipped_list[[length(skipped_list) + 1L]] <- g
            n_gated <- nrow(g)
          }
          if (any(gate_pass)) per_unit[[as.character(u)]] <- mp[gate_pass]
        }
        cat(sprintf("smide_de: %s / %s=%s -- %d genes tested, %d contrast rows, %d gated below %d expressing cells/arm\n",
                    ct, fit_unit, u, length(tg),
                    if (is.null(fit$mapped)) 0L else nrow(fit$mapped),
                    n_gated, MIN_EXPR_PER_ARM))
      } else {
        n_fit_fail <- n_fit_fail + 1L
        cat(sprintf("smide_de: ERROR %s / %s=%s -- %s\n", ct, fit_unit, u, fit$msg))
      }
    }

    if (length(per_unit) == 0) next
    ct_rows <- rbindlist(per_unit, fill = TRUE)
    ct_meta_rows <- meta_combine(ct_rows, fit_unit)
    results_list[[ct]] <- rbindlist(list(ct_rows, ct_meta_rows), use.names = TRUE, fill = TRUE)
    cat(sprintf("smide_de: %s -- %d per-unit rows, %d meta rows over %d unit(s)\n",
                ct, nrow(ct_rows),
                if (is.null(ct_meta_rows)) 0L else nrow(ct_meta_rows), length(per_unit)))
  }

  if (n_fit_ok == 0L)
    stop(sprintf("smide_de: cohort %s spatial mode -- all %d attempted fit(s) failed; no results produced",
                 cohort_name, n_fit_fail))
}

cat(sprintf("smide_de: %d fit(s) succeeded, %d failed\n", n_fit_ok, n_fit_fail))

# =============================================================================
# SECTION: WRITE
# =============================================================================

if (length(results_list) == 0) {
  abort_or_empty("every cell type failed or was skipped")
} else {
  out <- rbindlist(results_list, use.names = TRUE, fill = TRUE)
  if (de_mode == "spatial" && out[mode == "spatial_meta", .N] == 0)
    stop(sprintf("smide_de: cohort %s spatial mode produced per-unit rows but no meta rows",
                 cohort_name))
  write_out(out)
  cat(sprintf("smide_de: %d rows | %d genes x %d cell types x %d contrasts | modes: %s\n",
              nrow(out), uniqueN(out$feature_id), uniqueN(out$unit), uniqueN(out$contrast),
              paste(sort(unique(out$mode)), collapse = ", ")))
  if (de_mode == "spatial") {
    n_gated_total <- if (length(skipped_list)) sum(vapply(skipped_list, nrow, integer(1))) else 0L
    cat(sprintf("smide_de: %d per-unit row(s) gated below %d expressing cells per arm -> %s\n",
                n_gated_total, MIN_EXPR_PER_ARM, out_skipped))
    # Checked after the write so a hand-run keeps its tables for inspection; under
    # snakemake the job still fails and the outputs are discarded, which is the
    # intent -- an estimate this large means a degenerate fit reached the results.
    bad <- out[mode == "spatial_meta" & abs(estimate) > MAX_PLAUSIBLE_LOG2FC]
    if (nrow(bad) > 0) {
      worst <- bad[order(-abs(estimate))][seq_len(min(5L, .N))]
      stop(sprintf(paste0("smide_de: cohort %s spatial mode -- %d meta estimate(s) exceed %g log2 ",
                          "units after the >=%d expressing-cells-per-arm gate, which no real ",
                          "effect on this panel can reach. Worst: %s"),
                   cohort_name, nrow(bad), MAX_PLAUSIBLE_LOG2FC, MIN_EXPR_PER_ARM,
                   paste(sprintf("%s/%s/%s %.1f", worst$unit, worst$contrast,
                                 worst$feature_id, worst$estimate), collapse = "; ")))
    }
  }
}
