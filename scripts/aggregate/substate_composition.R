#!/usr/bin/env Rscript
# Sub-state composition adjudication (plan-differential-robustness Fix 3b). Propeller on
# the M02 day-2 labels with Fibroblast split into resting/activated. If the activated-
# fibroblast proportion shifts with SBRT and the collagen/Acta2 fibroblast DE tracks it,
# the fibrosis program is (partly) composition -- more activated cells -- not within-state
# regulation. Mirrors composition.R's propeller block on the relabelled clusters.
# Args: <obs.parquet> <full_labels.parquet> <fibroblast_substate.parquet> <out_test.tsv>
suppressPackageStartupMessages({ library(data.table); library(arrow); library(speckle); library(limma) })
a   <- commandArgs(trailingOnly = TRUE)
m   <- as.data.table(read_parquet(a[1]))
lab <- as.data.table(read_parquet(a[2]))[, .(cell, cell_type = cell_subtype)]
sub <- as.data.table(read_parquet(a[3]))                 # cell, substate
m[, cell_type := NULL]
m <- merge(m, lab, by = "cell", all.x = TRUE)[!is.na(cell_type)]
m <- merge(m, sub, by = "cell", all.x = TRUE)
m[cell_type == "Fibroblast", cell_type := paste0("Fibroblast_", substate)]
m[, condition := as.character(condition)]

m2 <- m[dataset == "Mutter_02"]
m2[, condition := factor(condition, levels = c("Control", "MBRT_day2", "SBRT_day2"))]
props <- getTransformedProps(clusters = m2$cell_type, sample = m2$sample_id, transform = "logit")
tp <- props$TransformedProps
tp <- tp[apply(tp, 1, function(r) all(is.finite(r))), , drop = FALSE]
samp <- unique(m2[, .(sample_id, condition, slide_id)]); setkey(samp, sample_id); samp <- samp[colnames(tp)]
design <- model.matrix(~ 0 + condition + slide_id, data = samp); colnames(design) <- make.names(colnames(design))
cm <- makeContrasts(
  MBRT_vs_Ctrl = conditionMBRT_day2 - conditionControl,
  SBRT_vs_Ctrl = conditionSBRT_day2 - conditionControl,
  MBRT_vs_SBRT = conditionMBRT_day2 - conditionSBRT_day2, levels = design)
fit2 <- eBayes(contrasts.fit(lmFit(tp, design), cm), robust = TRUE)
ln2  <- log(2)
mean_n_cells <- rowMeans(props$Counts)[rownames(tp)]
test <- rbindlist(lapply(colnames(cm), function(cn) {
  tt <- topTable(fit2, coef = cn, number = Inf, sort.by = "none", confint = TRUE)
  data.table(cell_type = rownames(tt), contrast = cn, log2FC_logit = tt$logFC / ln2,
             ci_low_log2 = tt$CI.L / ln2, ci_high_log2 = tt$CI.R / ln2,
             t_stat = tt$t, pvalue = tt$P.Value)
}))
test[, padj := p.adjust(pvalue, "BH")]
test[, mean_n_cells := round(mean_n_cells[cell_type], 1)]
setorder(test, contrast, -log2FC_logit)
fwrite(test, a[4], sep = "\t")
af <- function(ct, cn) test[cell_type == ct & contrast == cn,
                            sprintf("log2FC=%.2f padj=%.3g", log2FC_logit, padj)]
cat(sprintf("substate composition (SBRT_vs_Ctrl): activated %s | resting %s\n",
            af("Fibroblast_activated", "SBRT_vs_Ctrl"), af("Fibroblast_resting", "SBRT_vs_Ctrl")))
