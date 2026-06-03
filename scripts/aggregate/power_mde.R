#!/usr/bin/env Rscript
# Minimum-detectable-effect (MDE) yardstick for the M02 day-2 blocked design
# (n=4/arm, ~slide_id+condition, 6 residual df). Every null in the cross-dataset
# analysis is reported against the MDE so a flat result is read as "no effect
# larger than X was detectable", not "no effect". The detectable standardized
# effect at n=4/arm, 80% power, two-sided alpha=0.05 is d = power.t.test()$delta
# (base stats, sd=1; pwr not needed). MDE on each readout's native scale = d * SD,
# where SD is the within-arm biological SD estimated as the residual SD after
# ~slide_id+condition (the same block model the formal tests use).
#   (a) composition  -- SD of empirical-logit cell-type proportion across samples.
#   (b) pseudobulk DE -- per cell type, BCV = sqrt(median DESeq2 gene dispersion);
#                        SD on log2 scale = BCV/log(2); MDE is a log2FC. The
#                        abundance/gene floors match deg_pseudobulk.R so DE-MDE
#                        rows align with the cell types actually tested.
#   (c) program score -- SD of per-sample mean UCell score (project/primary sets).
# Args: <composition_by_sample.tsv> <pseudobulk_se.rds> <pathway_scores_summary.tsv> <out_power_mde.tsv>
suppressPackageStartupMessages({
  library(data.table); library(SummarizedExperiment); library(DESeq2)
})

a        <- commandArgs(trailingOnly = TRUE)
comp_in  <- a[1]
se_in    <- a[2]
path_in  <- a[3]
out_tsv  <- a[4]

N_ARM    <- 4L
POWER    <- 0.80
ALPHA    <- 0.05
MIN_CELLS    <- 10L
MIN_SAMPLES  <- 3L
GENE_MIN_CT  <- 10L
GENE_MIN_SMP <- 4L

# detectable standardized effect (Cohen's d) at n=4/arm, 80% power, two-sided 0.05
d_n4 <- power.t.test(n = N_ARM, sd = 1, power = POWER, sig.level = ALPHA,
                     type = "two.sample", alternative = "two.sided")$delta

# within-arm residual SD under the ~slide_id+condition block model (6 resid df)
resid_sd <- function(value, slide, cond) {
  ok <- is.finite(value)
  if (sum(ok) < 4L || length(unique(cond[ok])) < 2L) return(NA_real_)
  fit <- tryCatch(lm(value[ok] ~ slide[ok] + cond[ok]), error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  s <- summary(fit)$sigma
  if (!is.finite(s)) NA_real_ else s
}

# ---- (a) composition: empirical-logit proportion --------------------------------
comp <- fread(comp_in)
m2c  <- comp[dataset == "Mutter_02" & timepoint_h == 48L]
samp <- unique(m2c[, .(sample_id, condition, slide_id)])
tot  <- m2c[, .(sample_total = sum(n_cells)), by = sample_id]
grid <- CJ(sample_id = samp$sample_id, cell_type = unique(m2c$cell_type), unique = TRUE)
grid <- m2c[, .(sample_id, cell_type, n_cells)][grid, on = .(sample_id, cell_type)]
grid[is.na(n_cells), n_cells := 0L]
grid <- samp[grid, on = "sample_id"]
grid <- tot[grid, on = "sample_id"]
grid[, emp_logit := qlogis((n_cells + 0.5) / (sample_total + 1))]

comp_rows <- grid[, .(
  observed_sd = resid_sd(emp_logit, slide_id, condition),
  n_samples   = .N), by = cell_type]
comp_rows[, `:=`(readout_class = "composition", detail = NA_character_,
                 sd_scale = "empirical_logit_proportion",
                 mde_scale = "logit_proportion_difference")]

# ---- (b) pseudobulk DE: BCV from DESeq2 dispersions ------------------------------
se <- readRDS(se_in)
cd <- as.data.table(as.data.frame(colData(se)), keep.rownames = "col")
de_one <- function(ct) {
  sub_cd <- cd[cell_type == ct]
  usable <- sub_cd[n_cells >= MIN_CELLS]
  np <- usable[, .N, by = condition]
  if (nrow(np) < length(unique(sub_cd$condition)) ||
      any(np$N < MIN_SAMPLES)) return(NULL)          # abundance floor (matches DE)
  cols <- usable$col
  cm   <- assay(se)[, cols, drop = FALSE]
  keep <- rowSums(cm >= GENE_MIN_CT) >= GENE_MIN_SMP
  if (sum(keep) < 10L) return(NULL)
  sub  <- se[keep, cols]
  dds  <- tryCatch({
    d <- DESeqDataSetFromMatrix(assay(sub), as.data.frame(colData(sub)),
                                design = ~ slide_id + condition)
    d <- estimateSizeFactors(d); estimateDispersions(d, quiet = TRUE)
  }, error = function(e) NULL)
  if (is.null(dds)) return(NULL)
  disp <- dispersions(dds)
  md   <- median(disp[is.finite(disp)], na.rm = TRUE)
  bcv  <- sqrt(md)
  data.table(cell_type = ct, observed_sd = bcv / log(2), n_samples = length(cols),
             detail = sprintf("median_dispersion=%.4f; BCV=%.3f", md, bcv))
}
de_rows <- rbindlist(lapply(sort(unique(cd$cell_type)), de_one), fill = TRUE)
if (nrow(de_rows)) de_rows[, `:=`(readout_class = "pseudobulk_DE",
                                  sd_scale = "log2_expression",
                                  mde_scale = "log2_fold_change")]

# ---- (c) program score: per-sample mean UCell score, project/primary sets --------
ps   <- fread(path_in)
m2p  <- ps[dataset == "Mutter_02" & timepoint_h == 48L &
           score_type == "UCell" & tier == "primary"]
prog_rows <- m2p[, .(
  observed_sd = resid_sd(mean, slide_id, condition),
  n_samples   = .N), by = .(cell_type, pathway_name)]
setnames(prog_rows, "pathway_name", "detail")
prog_rows[, `:=`(readout_class = "program_score",
                 sd_scale = "ucell_score", mde_scale = "ucell_score_delta")]

# ---- assemble -------------------------------------------------------------------
cols <- c("readout_class", "cell_type", "detail", "n_per_arm", "residual_df",
          "cohens_d", "observed_sd", "sd_scale", "mde", "mde_scale", "n_samples")
out  <- rbindlist(list(comp_rows, de_rows, prog_rows), fill = TRUE)
out[, `:=`(n_per_arm = N_ARM, residual_df = 6L, cohens_d = d_n4)]
out[, mde := cohens_d * observed_sd]
setcolorder(out, cols)
setorder(out, readout_class, cell_type, detail)
fwrite(out, out_tsv, sep = "\t")

cat(sprintf("power/MDE: d(n=%d,80%%,a=0.05)=%.3f | composition %d, DE %d, program %d rows | %d NA-SD\n",
            N_ARM, d_n4, nrow(comp_rows), nrow(de_rows), nrow(prog_rows),
            sum(!is.finite(out$mde))))
