#!/usr/bin/env Rscript
# aggregate.smk pathway track -- ARM-TEST half (cross-timepoint parity, task 3). Reads the
# cached per-(sample x cell_type x pathway x score_type) summary from pathway_scores.R (the
# ~5h per-cell UCell/AddModuleScore compute, cohort-agnostic) and runs a per-(cell_type x
# score_type) stratified limma moderated-t test for ONE cohort's contrast family, with the
# design/contrast levels resolved from the comparison registry (lm_engine.R convention) --
# no hardcoded day-2 literals. Split from pathway_scores.R so a design/contrast change never
# forces the ~5h recompute, and parameterized by cohort so the test runs for mutter02_day2,
# combined_4h_treated, and combined_4h alike.
# Args: <pathway_scores_summary.tsv> <comparisons.tsv> <cohort> <out_test.tsv>
suppressPackageStartupMessages({
  library(limma); library(data.table)
})
args         <- commandArgs(trailingOnly = TRUE)
summary_path <- args[1]; comp_path <- args[2]; coh <- args[3]; out_test <- args[4]

# n=2/arm is the smallest per-arm n across the 3 registered cohorts (mutter02_day2 after the
# relabel correction; combined_4h's Control arm) -- 3 (the old day2-only threshold) would
# empty every mutter02_day2 stratum.
MIN_CELLS <- 10L; MIN_SAMPLES <- 2L
dir.create(dirname(out_test), recursive = TRUE, showWarnings = FALSE)

## ---- resolve the contrast family from the registry (lm_engine.R convention) ----
comp <- fread(comp_path)
fam  <- comp[cohort == coh & kind == "sample" & !is.na(contrast_num_level)]
fam  <- unique(fam, by = c("name", "contrast_num_level", "contrast_den_level"))
stopifnot(nrow(fam) >= 1)
formula_str <- fam$formula[1]
stopifnot(!is.na(formula_str))
CONDS <- unique(c(fam$contrast_num_level, fam$contrast_den_level))

## ---- cohort sample whitelist (prevents cross-dataset condition-name leakage, e.g. the
## shared "Control" condition name between mutter02_day2 and the M02 4h samples) ----
cohort_samples_path <- file.path(dirname(comp_path), "cohort_samples.tsv")
COHORT_SAMPLES <- if (file.exists(cohort_samples_path)) {
  fread(cohort_samples_path)[cohort == coh, sample_id]
} else character()

summ <- fread(summary_path)
## ---- sample-level metadata is authoritative from samples.tsv, not the summary cache.
## The summary's own condition/slide_id/dataset columns are whatever was baked into
## merged.rds at pathway_scores.R's last run and can drift out of sync with samples.tsv
## (the declared single source of truth for sample metadata); drop them and rejoin. ----
ss    <- fread("results/data_model/samples.tsv")
scol  <- names(ss)[grepl("sample", names(ss), ignore.case = TRUE)][1]
smeta <- unique(ss[, .(sample_id = get(scol), condition, slide_id, dataset = name)])
summ[, c("condition", "slide_id", "dataset") := NULL]
summ  <- merge(summ, smeta, by = "sample_id", sort = FALSE)

d <- summ
if (length(COHORT_SAMPLES) > 0) d <- d[sample_id %in% COHORT_SAMPLES]
d <- d[condition %in% CONDS]
cat(sprintf("pathway_arm_test[%s]: %d samples selected: %s\n", coh,
            length(unique(d$sample_id)), paste(sort(unique(d$sample_id)), collapse = ",")))
stopifnot(nrow(d) > 0)

pw_meta <- unique(summ[, .(pathway_name, pathway_source, tier, n_panel_genes, panel_coverage_frac)])

## ---- one contrast expression per registry row, shared across every stratum's fit ----
con_exprs <- setNames(
  sprintf("condition%s - condition%s", fam$contrast_num_level, fam$contrast_den_level),
  fam$name)

## ---- per-(cell_type x score_type) stratified fit ----
test_rows <- list()
for (ct in sort(unique(d$cell_type))) {
  for (st in c("UCell", "AMS")) {
    sub <- d[cell_type == ct & score_type == st & n_cells >= MIN_CELLS]
    if (nrow(sub) == 0) next
    sc  <- unique(sub[, .(sample_id, condition, slide_id)])
    cnt <- setNames(rep(0L, length(CONDS)), CONDS)
    cc  <- sc[, .N, by = condition]; cnt[as.character(cc$condition)] <- cc$N
    nmin <- min(cnt); if (nmin < MIN_SAMPLES) next
    w   <- dcast(sub, pathway_name ~ sample_id, value.var = "mean")
    mat <- as.matrix(w[, -1]); rownames(mat) <- w$pathway_name
    mat <- mat[complete.cases(mat), , drop = FALSE]; if (nrow(mat) == 0) next
    samp <- sc[match(colnames(mat), sample_id)]
    samp[, condition := factor(condition, levels = CONDS)]; samp[, slide_id := factor(slide_id)]
    design <- model.matrix(as.formula(formula_str), data = samp)
    colnames(design) <- make.names(colnames(design))
    if (qr(design)$rank < ncol(design)) next  # a CONDS level absent from this stratum -> skip
    cont <- makeContrasts(contrasts = con_exprs, levels = design)
    colnames(cont) <- names(con_exprs)
    fit2 <- eBayes(contrasts.fit(lmFit(mat, design), cont), robust = TRUE)
    for (cn in colnames(cont)) {
      est <- fit2$coefficients[, cn]
      se  <- fit2$stdev.unscaled[, cn] * sqrt(fit2$s2.post)
      test_rows[[length(test_rows) + 1]] <- data.table(
        cell_type = ct, pathway_name = rownames(mat), score_type = st, contrast = cn,
        estimate = est, se = se, t_stat = fit2$t[, cn], pvalue = fit2$p.value[, cn],
        n_samples_per_group = nmin)
    }
  }
}
stopifnot(length(test_rows) > 0)
test <- rbindlist(test_rows)
test[, padj_bh := p.adjust(pvalue, method = "BH")]
test <- merge(test, pw_meta, by = "pathway_name", sort = FALSE)
test[, dataset := paste(sort(unique(d$dataset)), collapse = ";")]
test_out <- test[, .(cell_type, pathway_name, pathway_source, tier, score_type, contrast,
                     estimate, se, t_stat, pvalue, padj_bh, n_samples_per_group,
                     n_panel_genes, panel_coverage_frac, dataset)]
setorder(test_out, cell_type, score_type, contrast, padj_bh, na.last = TRUE)
fwrite(test_out, out_test, sep = "\t")

cat(sprintf("pathway_arm_test[%s]: %d rows (%d padj<0.05)\n",
            coh, nrow(test_out), test_out[!is.na(padj_bh) & padj_bh < 0.05, .N]))
