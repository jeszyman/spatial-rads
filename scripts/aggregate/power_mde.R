#!/usr/bin/env Rscript
# Minimum-detectable-effect (MDE) yardstick for the sample-level comparisons in the
# resolved registry. Every null in the cross-dataset analysis is reported against the MDE
# so a flat result is read as "no effect larger than X was detectable", not "no effect".
# One row per (comparison, contrast) x readout x unit: per-arm n, residual df and the
# model formula whose residual SD is the noise estimate all come from the registry, so a
# cohort's floor reflects the design that cohort actually has (day-2 n=2/arm blocked on
# slide; the 4h cohorts n=2-3/arm unblocked). The detectable standardized effect at that
# n/arm, 80% power, two-sided alpha=0.05 is d = power.t.test()$delta (base stats, sd=1;
# pwr not needed). MDE on each readout's native scale = d * SD, where SD is the within-arm
# biological SD estimated as the residual SD under the cohort's registry formula.
#   (a) composition  -- SD of empirical-logit cell-type proportion across samples.
#   (b) pseudobulk DE -- per cell type, BCV = sqrt(median DESeq2 gene dispersion);
#                        SD on log2 scale = BCV/log(2); MDE is a log2FC. Abundance / gene
#                        floors come from config/engine_params.yaml -- the same constants
#                        count_engine.R applies -- so DE-MDE rows align with the cell
#                        types actually tested.
#   (c) program score -- SD of per-sample mean UCell score (project/primary sets).
# Args: <composition_by_sample.tsv> <pseudobulk_se.rds> <pathway_scores_summary.tsv>
#       <out_power_mde.tsv> [<comparisons.tsv>] [<engine_params.yaml>]
suppressPackageStartupMessages({
  library(data.table); library(SummarizedExperiment); library(DESeq2); library(yaml)
})

a         <- commandArgs(trailingOnly = TRUE)
comp_in   <- a[1]
se_in     <- a[2]
path_in   <- a[3]
out_tsv   <- a[4]
reg_in    <- if (length(a) >= 5) a[5] else "results/data_model/comparisons.tsv"
params_in <- if (length(a) >= 6) a[6] else "config/engine_params.yaml"
samples_in <- file.path(dirname(reg_in), "cohort_samples.tsv")

POWER <- 0.80
ALPHA <- 0.05

P <- read_yaml(params_in)$count_engine
MIN_CELLS    <- P$min_cells
MIN_SAMPLES  <- P$min_samples
GENE_MIN_CT  <- P$gene_min_count
GENE_MIN_SMP <- P$gene_min_samples

## ---- design parameters per (comparison, contrast), from the registry ----
REG <- fread(reg_in)
reg <- unique(
  REG[kind == "sample" & !is.na(contrast_num_level) & !is.na(formula) & formula != "",
      .(comparison = cohort, contrast = name, design_formula = formula,
        num = contrast_num_level, den = contrast_den_level,
        n_per_arm = pmin(n_group1, n_group2), residual_df = resid_df)],
  by = c("comparison", "contrast"))
stopifnot(nrow(reg) >= 1)

# detectable standardized effect (Cohen's d) at the contrast's own n/arm, 80% power, 0.05
reg[, cohens_d := vapply(n_per_arm, function(n) {
  if (!is.finite(n) || n < 2L) return(NA_real_)
  power.t.test(n = n, sd = 1, power = POWER, sig.level = ALPHA,
               type = "two.sample", alternative = "two.sided")$delta
}, numeric(1))]

COHORT_SAMPLES <- fread(samples_in)
covars_of <- function(f) setdiff(attr(terms(as.formula(f)), "term.labels"), "condition")

# within-arm residual SD under the cohort's registry model (condition + its covariates)
resid_sd <- function(d, covars) {
  d <- d[is.finite(value)]
  if (nrow(d) < 3L || uniqueN(d$condition) < 2L) return(NA_real_)
  f   <- as.formula(paste("value ~", paste(c("condition", covars), collapse = " + ")))
  fit <- tryCatch(lm(f, data = d), error = function(e) NULL)
  if (is.null(fit) || fit$df.residual < 1L) return(NA_real_)
  s <- summary(fit)$sigma
  if (!is.finite(s)) NA_real_ else s
}

comp_all <- fread(comp_in)
se       <- readRDS(se_in)
cd_all   <- as.data.table(as.data.frame(colData(se)), keep.rownames = "col")
ps_all   <- fread(path_in)

## ---- (b) helper: pseudobulk DE BCV for one cell type under one design ---------------
de_one <- function(cd, ct, covars) {
  sub_cd <- cd[cell_type == ct]
  usable <- sub_cd[n_cells >= MIN_CELLS]
  np <- usable[, .N, by = condition]
  if (nrow(np) < uniqueN(sub_cd$condition) || any(np$N < MIN_SAMPLES)) return(NULL)
  cols <- usable$col
  cm   <- assay(se)[, cols, drop = FALSE]
  keep <- rowSums(cm >= GENE_MIN_CT) >= GENE_MIN_SMP
  if (sum(keep) < 10L) return(NULL)
  sub  <- se[keep, cols]
  colData(sub)$condition <- factor(colData(sub)$condition)
  for (v in covars) {
    x <- colData(sub)[[v]]
    if (!is.numeric(x)) colData(sub)[[v]] <- factor(x)
  }
  design <- as.formula(paste("~", paste(c(covars, "condition"), collapse = " + ")))
  mm <- model.matrix(design, data = as.data.frame(colData(sub)))
  if (qr(mm)$rank < ncol(mm)) return(NULL)          # rank-deficient: the engine skips it too
  dds <- tryCatch({
    d <- DESeqDataSetFromMatrix(assay(sub), as.data.frame(colData(sub)), design = design)
    d <- estimateSizeFactors(d); estimateDispersions(d, quiet = TRUE)
  }, error = function(e) NULL)
  if (is.null(dds)) return(NULL)
  disp <- dispersions(dds)
  md   <- median(disp[is.finite(disp)], na.rm = TRUE)
  bcv  <- sqrt(md)
  data.table(cell_type = ct, observed_sd = bcv / log(2), n_samples = length(cols),
             detail = sprintf("median_dispersion=%.4f; BCV=%.3f", md, bcv))
}

## ---- per-cohort SD tables, crossed with that cohort's contrasts ---------------------
cohort_rows <- function(coh) {
  r      <- reg[comparison == coh]
  covars <- covars_of(r$design_formula[1])
  smp    <- COHORT_SAMPLES[cohort == coh, sample_id]
  conds  <- unique(c(r$num, r$den))

  # A readout table that carries fewer of the cohort's samples than the registry declares
  # (stale condition labels, a cohort-restricted upstream step) yields an SD estimated on
  # the wrong replicate count, so say so instead of letting it pass as a floor.
  note_coverage <- function(readout, matched) {
    if (length(matched) < length(smp))
      cat(sprintf("power_mde[%s/%s]: %d of %d cohort samples present in the input (%s)\n",
                  coh, readout, length(matched), length(smp),
                  paste(setdiff(smp, matched), collapse = ",")))
  }

  # (a) composition: empirical-logit proportion
  cmp <- comp_all[sample_id %in% smp & condition %in% conds]
  note_coverage("composition", unique(cmp$sample_id))
  comp_rows <- NULL
  if (nrow(cmp)) {
    samp <- unique(cmp[, .(sample_id, condition, slide_id)])
    tot  <- cmp[, .(sample_total = sum(n_cells)), by = sample_id]
    grid <- CJ(sample_id = samp$sample_id, cell_type = unique(cmp$cell_type), unique = TRUE)
    grid <- cmp[, .(sample_id, cell_type, n_cells)][grid, on = .(sample_id, cell_type)]
    grid[is.na(n_cells), n_cells := 0L]
    grid <- samp[grid, on = "sample_id"]
    grid <- tot[grid, on = "sample_id"]
    grid[, value := qlogis((n_cells + 0.5) / (sample_total + 1))]
    comp_rows <- grid[, .(observed_sd = resid_sd(.SD, covars), n_samples = .N),
                      by = cell_type, .SDcols = unique(c("value", "condition", covars))]
    comp_rows[, `:=`(readout_class = "composition", detail = NA_character_,
                     sd_scale = "empirical_logit_proportion",
                     mde_scale = "logit_proportion_difference")]
  }

  # (b) pseudobulk DE: BCV from DESeq2 dispersions
  cd <- cd_all[sample_id %in% smp & condition %in% conds]
  note_coverage("pseudobulk_DE", unique(cd$sample_id))
  de_rows <- rbindlist(lapply(sort(unique(cd$cell_type)), function(ct) de_one(cd, ct, covars)),
                       fill = TRUE)
  if (nrow(de_rows)) de_rows[, `:=`(readout_class = "pseudobulk_DE",
                                    sd_scale = "log2_expression",
                                    mde_scale = "log2_fold_change")]

  # (c) program score: per-sample mean UCell score, project/primary sets
  psx <- ps_all[sample_id %in% smp & condition %in% conds &
                score_type == "UCell" & tier == "primary"]
  note_coverage("program_score", unique(psx$sample_id))
  prog_rows <- NULL
  if (nrow(psx)) {
    psx <- copy(psx); setnames(psx, "mean", "value")
    prog_rows <- psx[, .(observed_sd = resid_sd(.SD, covars), n_samples = .N),
                     by = .(cell_type, pathway_name),
                     .SDcols = unique(c("value", "condition", covars))]
    setnames(prog_rows, "pathway_name", "detail")
    prog_rows[, `:=`(readout_class = "program_score",
                     sd_scale = "ucell_score", mde_scale = "ucell_score_delta")]
  }

  sd_rows <- rbindlist(list(comp_rows, de_rows, prog_rows), fill = TRUE)
  if (!nrow(sd_rows)) return(NULL)
  sd_rows[, k := 1L]
  rr <- r[, .(comparison, contrast, n_per_arm, residual_df, design_formula, cohens_d, k = 1L)]
  merge(rr, sd_rows, by = "k", allow.cartesian = TRUE)[, k := NULL][]
}

## ---- assemble -------------------------------------------------------------------
cols <- c("comparison", "contrast", "readout_class", "cell_type", "detail",
          "n_per_arm", "residual_df", "design_formula", "cohens_d",
          "observed_sd", "sd_scale", "mde", "mde_scale", "n_samples")
out <- rbindlist(lapply(unique(reg$comparison), cohort_rows), fill = TRUE)
stopifnot(nrow(out) > 0)
out[, mde := cohens_d * observed_sd]
setcolorder(out, cols)
setorder(out, comparison, contrast, readout_class, cell_type, detail)
fwrite(out, out_tsv, sep = "\t")

summ <- out[, .(readouts = uniqueN(readout_class), rows = .N,
                d = round(cohens_d[1], 3), n = n_per_arm[1], df = residual_df[1],
                na_sd = sum(!is.finite(mde))), by = comparison]
cat("power/MDE per comparison (d at that cohort's n/arm, 80% power, alpha=0.05):\n")
print(summ)
cat(sprintf("power/MDE: %d rows over %d comparisons x %d contrasts | %d NA-SD\n",
            nrow(out), uniqueN(out$comparison), uniqueN(out$contrast),
            sum(!is.finite(out$mde))))
