#!/usr/bin/env Rscript
# Arm test on InSituCor module scores. This is OUR extension: InSituCor's documented
# workflow stops at prioritization and never tests modules against a treatment variable.
# The test is admissible because discovery pooled all 4h cells across arms with no
# treatment variable in the model, so modules were defined blind to treatment.
# Machinery mirrors pathway_arm_test.R exactly: per-cell scores aggregated to one mean per
# (sample x cell_type x module), then a per-cell_type stratified limma moderated-t fit with
# the design and contrast levels resolved from the comparison registry. The moderated t is
# the reported inference for the 4h cohorts (the same low-residual-df policy declared in
# scripts/aggregate/inference_source.R for the pseudobulk layer); there is no Wald test on
# this path for any cohort. Exploratory tier only: BH within cohort, never part of the
# confirmatory family. Runs only for 4h cohorts, because per-cell scores exist only for
# the cells the (4h-pooled) discovery scored.
# Args: <scores_sc.parquet> <full_labels.parquet> <coords.parquet> <modules.tsv>
#       <comparisons.tsv> <cohort> <out.tsv>
suppressPackageStartupMessages({
  library(data.table); library(arrow); library(limma)
})
a <- commandArgs(trailingOnly = TRUE)
scores_path <- a[1]; labels_path <- a[2]; coords_path <- a[3]; modules_path <- a[4]
comp_path <- a[5]; coh <- a[6]; out_path <- a[7]

MIN_CELLS <- 10L; MIN_SAMPLES <- 2L; MIN_GENES <- 4L
# diffuse catch-all rule shared with insitucor_prioritize.R
CATCHALL_WEIGHT_FRAC <- 0.1
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

## ---- resolve the contrast family from the registry (lm_engine.R convention) ----
comp <- fread(comp_path)
fam  <- comp[cohort == coh & kind == "sample" & !is.na(contrast_num_level)]
fam  <- unique(fam, by = c("name", "contrast_num_level", "contrast_den_level"))
stopifnot(nrow(fam) >= 1)
formula_str <- fam$formula[1]
stopifnot(!is.na(formula_str))
CONDS <- unique(c(fam$contrast_num_level, fam$contrast_den_level))

## ---- cohort sample whitelist (prevents cross-dataset condition-name leakage) ----
cohort_samples_path <- file.path(dirname(comp_path), "cohort_samples.tsv")
COHORT_SAMPLES <- if (file.exists(cohort_samples_path)) {
  fread(cohort_samples_path)[cohort == coh, sample_id]
} else character()

## ---- tested module set: >= MIN_GENES genes, excluding diffuse background catch-alls
## (per-gene weights an order of magnitude below every real module) ----
mod <- fread(modules_path)
msize <- mod[, .(n_genes = .N, median_weight = median(abs(gene_weight))), by = module]
msize <- msize[n_genes >= MIN_GENES]
msize[, diffuse_catchall := median_weight < CATCHALL_WEIGHT_FRAC * median(median_weight)]
TESTED <- msize[diffuse_catchall == FALSE, module]
stopifnot(length(TESTED) >= 1)

## ---- per-cell scores -> per (sample x cell_type x module) means ----
sc  <- as.data.table(read_parquet(scores_path))
lab <- as.data.table(read_parquet(labels_path))[, .(cell, cell_subtype)]
co  <- as.data.table(read_parquet(coords_path))[, .(cell, sample_id)]
sc  <- lab[sc, on = "cell"][co, on = "cell", nomatch = NULL]
sc  <- sc[!is.na(cell_subtype)]

long <- melt(sc, id.vars = c("cell", "cell_subtype", "sample_id"),
             variable.name = "module", value.name = "score", variable.factor = FALSE)
long <- long[module %in% TESTED]
agg  <- long[, .(mean_score = mean(score), n_cells = .N),
             by = .(sample_id, cell_type = cell_subtype, module)]
agg  <- agg[n_cells >= MIN_CELLS]

## ---- sample-level metadata is authoritative from samples.tsv ----
ss    <- fread("results/data_model/samples.tsv")
scol  <- names(ss)[grepl("sample", names(ss), ignore.case = TRUE)][1]
smeta <- unique(ss[, .(sample_id = get(scol), condition, slide_id, dataset = name)])
agg   <- merge(agg, smeta, by = "sample_id", sort = FALSE)

d <- agg
if (length(COHORT_SAMPLES) > 0) d <- d[sample_id %in% COHORT_SAMPLES]
d <- d[condition %in% CONDS]
cat(sprintf("insitucor_arm_test[%s]: %d samples selected: %s\n", coh,
            length(unique(d$sample_id)), paste(sort(unique(d$sample_id)), collapse = ",")))
stopifnot(nrow(d) > 0)

## ---- one contrast expression per registry row, shared across every stratum's fit ----
con_exprs <- setNames(
  sprintf("condition%s - condition%s", fam$contrast_num_level, fam$contrast_den_level),
  fam$name)

## ---- per-cell_type stratified fit ----
test_rows <- list()
for (ct in sort(unique(d$cell_type))) {
  sub <- d[cell_type == ct]
  if (nrow(sub) == 0) next
  scnd <- unique(sub[, .(sample_id, condition, slide_id)])
  cnt  <- setNames(rep(0L, length(CONDS)), CONDS)
  cc   <- scnd[, .N, by = condition]; cnt[as.character(cc$condition)] <- cc$N
  nmin <- min(cnt); if (nmin < MIN_SAMPLES) next
  w   <- dcast(sub, module ~ sample_id, value.var = "mean_score")
  mat <- as.matrix(w[, -1]); rownames(mat) <- w$module
  mat <- mat[complete.cases(mat), , drop = FALSE]; if (nrow(mat) == 0) next
  samp <- scnd[match(colnames(mat), sample_id)]
  samp[, condition := factor(condition, levels = CONDS)]; samp[, slide_id := factor(slide_id)]
  design <- model.matrix(as.formula(formula_str), data = samp)
  colnames(design) <- make.names(colnames(design))
  if (qr(design)$rank < ncol(design)) next  # a CONDS level absent from this stratum -> skip
  cont <- makeContrasts(contrasts = con_exprs, levels = design)
  colnames(cont) <- names(con_exprs)
  fit2 <- eBayes(contrasts.fit(lmFit(mat, design), cont), robust = TRUE)
  df_res   <- setNames(fit2$df.residual, rownames(mat))
  df_total <- setNames(fit2$df.total, rownames(mat))
  for (cn in colnames(cont)) {
    est <- fit2$coefficients[, cn]
    se  <- fit2$stdev.unscaled[, cn] * sqrt(fit2$s2.post)
    test_rows[[length(test_rows) + 1]] <- data.table(
      cell_type = ct, module = rownames(mat), contrast = cn,
      estimate = est, se = se, t_stat = fit2$t[, cn], pvalue = fit2$p.value[, cn],
      df_residual = df_res[rownames(mat)], df_total = df_total[rownames(mat)],
      n_samples_per_group = nmin)
  }
}
stopifnot(length(test_rows) > 0)
test <- rbindlist(test_rows)
test[, padj_bh := p.adjust(pvalue, method = "BH")]   # BH within cohort, exploratory tier
test <- merge(test, msize, by = "module", sort = FALSE)
test[, `:=`(cohort = coh, tier = "exploratory", inference = "limma_moderated_t")]
test_out <- test[, .(cohort, cell_type, module, n_genes, contrast, estimate, se, t_stat,
                     df_residual, df_total, pvalue, padj_bh, n_samples_per_group,
                     tier, inference)]
setorder(test_out, padj_bh, pvalue, na.last = TRUE)
fwrite(test_out, out_path, sep = "\t")

cat(sprintf("insitucor_arm_test[%s]: %d modules x %d strata rows (%d padj<0.05)\n",
            coh, length(TESTED), nrow(test_out),
            test_out[!is.na(padj_bh) & padj_bh < 0.05, .N]))
