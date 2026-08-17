#!/usr/bin/env Rscript
# limma-voom moderated-t pseudobulk engine. Same input, same cohort whitelist, same
# abundance/gene filters and the same registry design as count_engine.R, fit by
# limma-voom instead of DESeq2 so the reported test is a moderated t.
#
# Why a second engine rather than a switch inside count_engine.R: the two are not
# alternatives to choose between, they are two views of one fit that the results tier
# carries side by side. A blocked 4h cohort has 1-2 residual df, where a per-gene
# variance estimate is essentially unusable; empirical-Bayes shrinkage borrows df
# across genes and restores a usable denominator (Smyth 2004; Law 2014 for the voom
# mean-variance weights that make the borrowing valid on counts). DESeq2's Wald p
# stays available as a secondary column, so nothing is discarded.
#
# The filters are deliberately identical to count_engine.R (min_cells, min_samples,
# gene_min_count, gene_min_samples read from the same engine_params block): both
# engines must test exactly the same (gene x cell_type x contrast) set or the two
# columns in the results table would not be comparable row by row.
#
# Emits sufficient statistics only, in the fixed engine schema (comparison, contrast,
# unit, feature_type, feature_id, estimate, se, df, stat, p) plus the two df columns a
# moderated fit needs to be readable: df_residual (what the design realizes) and df
# (df.total = residual + prior, what the moderated t is actually tested against).
# No adjusted p -- correction is the results tier's job.
# Args: <pseudobulk_se.rds> <cohort> <comparisons.tsv> <engine_params.yaml>
#       <out_degs.tsv> <out_skipped.tsv>
suppressPackageStartupMessages({
  library(SummarizedExperiment); library(edgeR); library(limma)
  library(data.table); library(yaml)
})
a <- commandArgs(trailingOnly = TRUE)
se_path <- a[1]; coh <- a[2]; comp_path <- a[3]; params_path <- a[4]
out_degs <- a[5]; out_skipped <- a[6]

## ---- thresholds: shared with count_engine so both engines test the same gene set ----
PAR <- read_yaml(params_path)
P <- PAR$count_engine
MIN_CELLS <- P$min_cells; MIN_SAMPLES <- P$min_samples
GENE_MIN_CT <- P$gene_min_count; GENE_MIN_SMP <- P$gene_min_samples
ROBUST <- isTRUE(PAR$voom_engine$robust)

## ---- cohort sample whitelist (prevents cross-dataset condition-name leakage) ----
cohort_samples_path <- file.path(dirname(comp_path), "cohort_samples.tsv")
COHORT_SAMPLES <- if (file.exists(cohort_samples_path)) {
  fread(cohort_samples_path)[cohort == coh, sample_id]
} else character()

## ---- contrasts + model terms from the registry: name -> (num, den) ----
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
# The registry writes the cell-means form (~0 + condition + ...), which is exactly what
# makeContrasts needs -- limma uses the registry formula verbatim, no re-parameterization.
DESIGN <- as.formula(formula_str)
CONTRASTS <- lapply(seq_len(nrow(fam)), function(i)
  list(name = fam$name[i], num = fam$contrast_num_level[i], den = fam$contrast_den_level[i]))
all_levels <- unique(c(vapply(CONTRASTS, `[[`, "", "num"), vapply(CONTRASTS, `[[`, "", "den")))
den_counts  <- sort(table(vapply(CONTRASTS, `[[`, "", "den")), decreasing = TRUE)
CONDS <- c(names(den_counts)[1], sort(setdiff(all_levels, names(den_counts)[1])))

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
    if (!is.numeric(x)) colData(sub)[[v]] <- droplevels(factor(x))
  }

  mm <- model.matrix(DESIGN, data = as.data.frame(colData(sub)))
  colnames(mm) <- make.names(colnames(mm))
  if (qr(mm)$rank < ncol(mm)) {                                # design not full rank
    for (cc in CONTRASTS) skip_rows[[length(skip_rows) + 1]] <- data.table(
      cell_type = ct, contrast = cc$name,
      reason = sprintf("design not full rank after abundance filter (%s confounded)", formula_str),
      n_cells_in_failed_group = NA_integer_)
    next
  }

  cm     <- assay(sub)
  keep_g <- rowSums(cm >= GENE_MIN_CT) >= GENE_MIN_SMP
  sub    <- sub[keep_g, ]
  n_used <- ncol(sub)
  df_res <- n_used - qr(mm)$rank                               # residual df before moderation
  if (df_res < 1) {                        # nothing left for the prior to shrink toward
    for (cc in CONTRASTS) skip_rows[[length(skip_rows) + 1]] <- data.table(
      cell_type = ct, contrast = cc$name,
      reason = sprintf("zero residual df under %s (%d samples, rank %d)",
                       formula_str, n_used, qr(mm)$rank),
      n_cells_in_failed_group = NA_integer_)
    next
  }

  y <- DGEList(assay(sub)); y <- calcNormFactors(y)
  v <- voom(y, mm)
  fit <- lmFit(v, mm)

  for (cc in CONTRASTS) {
    num <- make.names(paste0("condition", cc$num)); den <- make.names(paste0("condition", cc$den))
    if (!all(c(num, den) %in% colnames(mm))) {
      skip_rows[[length(skip_rows) + 1]] <- data.table(
        cell_type = ct, contrast = cc$name,
        reason = sprintf("contrast level(s) absent from the design: %s",
                         paste(setdiff(c(num, den), colnames(mm)), collapse = ", ")),
        n_cells_in_failed_group = NA_integer_)
      next
    }
    cmat <- makeContrasts(contrasts = paste(num, "-", den), levels = mm)
    f2   <- eBayes(contrasts.fit(fit, cmat), robust = ROBUST)
    # Moderated SE = sqrt(posterior variance) x unscaled SE of the contrast; the
    # moderated t divides the same log2FC by this instead of the per-gene-only SE.
    se_mod <- sqrt(f2$s2.post) * f2$stdev.unscaled[, 1]
    deg_rows[[length(deg_rows) + 1]] <- data.table(
      comparison = coh, contrast = cc$name, unit = ct, feature_type = "gene",
      feature_id = rownames(f2),
      estimate = as.numeric(f2$coefficients[, 1]), se = as.numeric(se_mod),
      df = as.numeric(f2$df.total), stat = as.numeric(f2$t[, 1]), p = as.numeric(f2$p.value[, 1]),
      df_residual = df_res, df_prior = as.numeric(f2$df.prior),
      n_samples_used = n_used)
  }
}

degs <- rbindlist(deg_rows)
setorder(degs, unit, contrast, p, na.last = TRUE)
fwrite(degs, out_degs, sep = "\t")
skipped <- if (length(skip_rows)) rbindlist(skip_rows) else data.table(
  cell_type = character(), contrast = character(),
  reason = character(), n_cells_in_failed_group = integer())
fwrite(skipped, out_skipped, sep = "\t")

cat(sprintf("voom_engine[%s]: registry %s -> limma-voom moderated t (robust=%s) | %d cell types tested, %d skipped | %d (gene x ct x contrast) rows\n",
            coh, formula_str, ROBUST, uniqueN(degs$unit), uniqueN(skipped$cell_type), nrow(degs)))
if (nrow(degs) > 0)
  cat(sprintf("voom_engine[%s]: residual df %s -> moderated df.total %s\n", coh,
              paste(sort(unique(degs$df_residual)), collapse = "/"),
              paste(sprintf("%.1f", range(degs$df)), collapse = "-")))
