#!/usr/bin/env Rscript
# Negative-binomial pseudobulk count engine (modularization Task 3). Per-cell-type DE on
# the (sample x cell_type) pseudobulk SummarizedExperiment via DESeq2, reading the arm
# contrasts + their reference levels from the resolved comparison registry -- no embedded
# arm literals. Ports deg_pseudobulk.R's statistics verbatim, with three changes: contrasts
# from the registry, thresholds from engine_params.yaml, and NO padj emitted (correction is
# the results tier's job).
#
# apeglm two-fit pattern, preserved: lfcShrink(type="apeglm") accepts only coef=, not an
# arbitrary contrast vector, so a shrunk MBRT-vs-SBRT estimate needs SBRT as the DESeq2
# reference. The engine therefore builds ONE DESeq fit per distinct reference level named
# in the registry (here Control and SBRT_day2) and pulls each contrast from the fit whose
# reference is that contrast's denominator.
#
# The model TERMS come from the registry `formula` column, like lm_engine.R -- a blocked
# cohort (~0 + condition + slide_id) is fit with its block, an unblocked one
# (~0 + condition) without. Only the parameterization is method-internal: the registry
# writes the intercept-free form limma needs for makeContrasts, while DESeq2's
# results(name=)/lfcShrink(coef=) path needs treatment contrasts, so the same term set is
# re-expressed with an intercept (~ covariates + condition). Term set from the registry,
# parameterization from the method.
# Args: <pseudobulk_se.rds> <cohort> <comparisons.tsv> <engine_params.yaml>
#       <out_degs.tsv> <out_skipped.tsv>
suppressPackageStartupMessages({
  library(SummarizedExperiment); library(DESeq2); library(apeglm)
  library(data.table); library(yaml)
})
a <- commandArgs(trailingOnly = TRUE)
se_path <- a[1]; coh <- a[2]; comp_path <- a[3]; params_path <- a[4]
out_degs <- a[5]; out_skipped <- a[6]

## ---- thresholds (study-wide, from engine config) ----
P <- read_yaml(params_path)$count_engine
MIN_CELLS <- P$min_cells; MIN_SAMPLES <- P$min_samples
GENE_MIN_CT <- P$gene_min_count; GENE_MIN_SMP <- P$gene_min_samples

## ---- cohort sample whitelist (prevents cross-dataset condition-name leakage) ----
cohort_samples_path <- file.path(dirname(comp_path), "cohort_samples.tsv")
COHORT_SAMPLES <- if (file.exists(cohort_samples_path)) {
  fread(cohort_samples_path)[cohort == coh, sample_id]
} else character()

## ---- contrasts + model terms from the registry: name -> (num, den=ref) ----
comp <- fread(comp_path)
fam  <- comp[cohort == coh & kind == "sample" & !is.na(contrast_num_level)]
fam  <- unique(fam, by = c("name", "contrast_num_level", "contrast_den_level"))
stopifnot(nrow(fam) >= 1)
formula_str <- fam$formula[1]
stopifnot(!is.na(formula_str), nzchar(formula_str))
rhs_terms <- attr(terms(as.formula(formula_str)), "term.labels")
if (!"condition" %in% rhs_terms)
  stop(sprintf("registry formula for cohort %s has no condition term: %s", coh, formula_str))
COVARS <- setdiff(rhs_terms, "condition")
DESIGN <- as.formula(paste("~", paste(c(COVARS, "condition"), collapse = " + ")))
DESIGN_STR <- paste(deparse(DESIGN), collapse = "")
CONTRASTS <- lapply(seq_len(nrow(fam)), function(i)
  list(name = fam$name[i], num = fam$contrast_num_level[i], den = fam$contrast_den_level[i]))
# condition levels: primary reference (den shared by the vs-baseline contrasts) first, so the
# default fit's reference matches deg_pseudobulk.R (Control), then the remaining arms sorted.
den_counts  <- sort(table(vapply(CONTRASTS, `[[`, "", "den")), decreasing = TRUE)
primary_ref <- names(den_counts)[1]
all_levels  <- unique(c(vapply(CONTRASTS, `[[`, "", "num"), vapply(CONTRASTS, `[[`, "", "den")))
CONDS       <- c(primary_ref, sort(setdiff(all_levels, primary_ref)))
ref_levels  <- unique(vapply(CONTRASTS, `[[`, "", "den"))     # one DESeq fit per distinct ref

dir.create(dirname(out_degs), recursive = TRUE, showWarnings = FALSE)
se <- readRDS(se_path)
cd <- as.data.table(as.data.frame(colData(se)), keep.rownames = "col")
missing_covar <- setdiff(COVARS, names(cd))
if (length(missing_covar))
  stop(sprintf("registry formula term(s) absent from the pseudobulk colData: %s",
               paste(missing_covar, collapse = ", ")))

if (length(COHORT_SAMPLES) > 0) {
  keep_cols <- cd$col[cd$sample_id %in% COHORT_SAMPLES]
  se <- se[, keep_cols]
  cd <- cd[sample_id %in% COHORT_SAMPLES]
}

deg_rows <- list(); skip_rows <- list()
usable   <- cd[condition %in% CONDS & n_cells >= MIN_CELLS, .(n_pass = .N), by = .(cell_type, condition)]
all_types <- sort(unique(cd[condition %in% CONDS, cell_type]))

for (ct in all_types) {
  u  <- usable[cell_type == ct]
  np <- setNames(rep(0L, length(CONDS)), CONDS); np[u$condition] <- u$n_pass
  failing <- names(np)[np < MIN_SAMPLES]
  if (length(failing) > 0) {                                   # abundance floor -> skip type
    for (cc in CONTRASTS) {
      bad <- intersect(c(cc$num, cc$den), failing); if (!length(bad)) bad <- failing[1]
      ncells_bad <- cd[cell_type == ct & condition == bad[1], sum(n_cells)]
      n_avail    <- cd[cell_type == ct & condition == bad[1], .N]
      skip_rows[[length(skip_rows) + 1]] <- data.table(
        cell_type = ct, contrast = cc$name,
        reason = sprintf("abundance_floor: condition %s has %d/%d samples >=%d cells (need %d)",
                         bad[1], np[[bad[1]]], n_avail, MIN_CELLS, MIN_SAMPLES),
        n_cells_in_failed_group = ncells_bad)
    }
    next
  }

  cols <- cd[cell_type == ct & condition %in% CONDS & n_cells >= MIN_CELLS, col]
  sub  <- se[, cols]
  colData(sub)$condition <- factor(colData(sub)$condition, levels = CONDS)
  for (v in COVARS) {                       # blocking terms enter as factors, measured ones as-is
    x <- colData(sub)[[v]]
    if (!is.numeric(x)) colData(sub)[[v]] <- factor(x)
  }

  mm <- model.matrix(DESIGN, data = as.data.frame(colData(sub)))
  if (qr(mm)$rank < ncol(mm)) {                                # design not full rank
    for (cc in CONTRASTS) skip_rows[[length(skip_rows) + 1]] <- data.table(
      cell_type = ct, contrast = cc$name,
      reason = sprintf("design not full rank after abundance filter (%s confounded)", DESIGN_STR),
      n_cells_in_failed_group = NA_integer_)
    next
  }

  cm     <- assay(sub)
  keep_g <- rowSums(cm >= GENE_MIN_CT) >= GENE_MIN_SMP
  sub    <- sub[keep_g, ]
  n_used <- ncol(sub)
  df_res <- n_used - qr(mm)$rank                               # residual df for the fit

  # one DESeq fit per distinct reference level (generalizes the Control/SBRT two-fit)
  fits <- setNames(lapply(ref_levels, function(rl) {
    d <- DESeqDataSetFromMatrix(assay(sub), as.data.frame(colData(sub)),
                                design = DESIGN)
    colData(d)$condition <- relevel(factor(colData(d)$condition, levels = CONDS), ref = rl)
    DESeq(d, quiet = TRUE)
  }), ref_levels)

  for (cc in CONTRASTS) {
    coef <- sprintf("condition_%s_vs_%s", cc$num, cc$den)
    obj  <- fits[[cc$den]]
    res  <- results(obj, name = coef, independentFiltering = TRUE)  # baseMean, stat, pvalue
    shr  <- lfcShrink(obj, coef = coef, type = "apeglm", quiet = TRUE)  # log2FC, lfcSE
    deg_rows[[length(deg_rows) + 1]] <- data.table(
      comparison = coh, contrast = cc$name, unit = ct, feature_type = "gene",
      feature_id = rownames(res),
      estimate = shr$log2FoldChange, se = shr$lfcSE, df = df_res,
      stat = res$stat, p = res$pvalue,
      baseMean = res$baseMean, n_samples_used = n_used)
  }
}

degs <- rbindlist(deg_rows)
setorder(degs, unit, contrast, p, na.last = TRUE)
fwrite(degs, out_degs, sep = "\t")
skipped <- if (length(skip_rows)) rbindlist(skip_rows) else data.table(
  cell_type = character(), contrast = character(),
  reason = character(), n_cells_in_failed_group = integer())
fwrite(skipped, out_skipped, sep = "\t")

cat(sprintf("count_engine[%s]: registry %s -> DESeq2 %s | %d cell types tested, %d skipped | %d (gene x ct x contrast) rows\n",
            coh, formula_str, DESIGN_STR,
            uniqueN(degs$unit), uniqueN(skipped$cell_type), nrow(degs)))
