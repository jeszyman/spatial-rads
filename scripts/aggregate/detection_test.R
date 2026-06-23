#!/usr/bin/env Rscript
# Detection vs level decomposition for the pseudobulk DE (plan-differential-robustness Fix 2).
# Per cell_type: limma empirical-Bayes moderated test on arcsin-sqrt per-sample DETECTION
# fractions, and (in parallel) on log2 per-expresser LEVEL. Sparsity-aware classifier
# separates regulation from a cell-state-fraction shift. n=4/arm; small-n is handled by
# eBayes borrowing strength across genes (the propeller/Phipson rationale), NOT a per-gene
# beta-binomial (too few residual df at n=4). Args: <pseudobulk_se.rds> <out.tsv>
suppressPackageStartupMessages({
  library(SummarizedExperiment); library(limma); library(data.table)
})
a  <- commandArgs(trailingOnly = TRUE)
se <- readRDS(a[1]); out <- a[2]
cd <- as.data.table(as.data.frame(colData(se)))
cd[, condition := factor(condition, levels = c("Control", "MBRT_day2", "SBRT_day2"))]
ln2         <- log(2)
LEVEL_FLOOR <- 2      # raw counts among expressers below this => per-cell level near-binary/uninterpretable => ambiguous (tunable)

# moderated fit on a (genes x samples) transformed matrix; cols index cd rows (global)
fit_one <- function(mat, cols) {
  samp <- cd[cols]
  d <- model.matrix(~ 0 + condition + slide_id, data = samp)
  colnames(d) <- make.names(colnames(d))
  cm <- makeContrasts(
    MBRT_vs_Ctrl = conditionMBRT_day2 - conditionControl,
    SBRT_vs_Ctrl = conditionSBRT_day2 - conditionControl,
    MBRT_vs_SBRT = conditionMBRT_day2 - conditionSBRT_day2, levels = d)
  eBayes(contrasts.fit(lmFit(mat, d), cm), robust = TRUE)
}
CONTR <- c("MBRT_vs_Ctrl", "SBRT_vs_Ctrl", "MBRT_vs_SBRT")

res <- list()
for (ct in unique(cd$cell_type)) {
  cols <- which(cd$cell_type == ct)
  if (length(unique(cd$condition[cols])) < 3 || length(cols) < 6) next  # need all arms + replicates
  pe  <- assay(se, "pct_expr")[, cols, drop = FALSE]              # genes x n(this type's samples)
  mae <- assay(se, "mean_among_expr")[, cols, drop = FALSE]
  keep <- rowSums(pe > 0) >= 3                                    # gene detected in >=3 samples
  if (sum(keep) < 2) next
  fit <- tryCatch({
    fd <- fit_one(asin(sqrt(pe[keep, , drop = FALSE])), cols)     # detection model
    fl <- fit_one(log2(mae[keep, , drop = FALSE] + 1), cols)      # level model
    list(fd = fd, fl = fl)
  }, error = function(e) NULL)                                    # skip rank-deficient sparse types
  if (is.null(fit)) next
  arms    <- split(seq_along(cols), cd$condition[cols])          # local indices into pe/mae
  mae_arm <- sapply(arms, function(j) rowMeans(mae[keep, j, drop = FALSE]))
  mae_max <- apply(as.matrix(mae_arm), 1, max)                   # max per-arm mean-among-expr
  for (cn in CONTR) {
    td <- topTable(fit$fd, coef = cn, number = Inf, sort.by = "none")
    tl <- topTable(fit$fl, coef = cn, number = Inf, sort.by = "none")
    res[[paste(ct, cn)]] <- data.table(
      gene = rownames(td), cell_type = ct, contrast = cn,
      detection_log2fc = td$logFC / ln2, detection_p = td$P.Value,
      level_log2fc = tl$logFC / ln2,     level_p = tl$P.Value,
      mean_among_expr_max = round(mae_max, 2))
  }
}
dt <- rbindlist(res)
dt[, detection_padj := p.adjust(detection_p, "BH"), by = .(cell_type, contrast)]  # match DESeq2: BH within cell_type x contrast, not globally
dt[, level_padj     := p.adjust(level_p, "BH"),     by = .(cell_type, contrast)]
# Symmetric floor: below it the per-cell level component is uninterpretable, so neither a
# level-based "regulation" nor a flat-level "fraction_shift" call is trustworthy -> ambiguous.
dt[, call_class := fifelse(mean_among_expr_max < LEVEL_FLOOR, "ambiguous",
                    fifelse(level_padj < 0.05, "regulation",
                    fifelse(detection_padj < 0.05, "fraction_shift", "ambiguous")))]
setcolorder(dt, c("gene", "cell_type", "contrast", "detection_log2fc", "detection_padj",
                  "level_log2fc", "level_padj", "mean_among_expr_max", "call_class"))
fwrite(dt, out, sep = "\t")
cat(sprintf("detection_test: %d rows | %d regulation / %d fraction_shift / %d ambiguous\n",
            nrow(dt), dt[call_class == "regulation", .N],
            dt[call_class == "fraction_shift", .N], dt[call_class == "ambiguous", .N]))
