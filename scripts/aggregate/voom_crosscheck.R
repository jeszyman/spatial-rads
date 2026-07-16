#!/usr/bin/env Rscript
# limma-voom DESeq2 concordance cross-check. Runs limma-voom on the SAME sample x
# cell_type pseudobulk SummarizedExperiment the count engine (DESeq2) consumes, with
# the same design (~0+condition+slide_id) and contrasts, and reports the per-gene
# log2FC Spearman concordance to DESeq2. DESeq2 stays primary; voom is a robustness
# anchor. Precedent for the voom idiom: dev/if_channel_check/pv_step2_fov_pseudobulk.R.
suppressPackageStartupMessages({
  library(SummarizedExperiment); library(edgeR); library(limma); library(data.table)
})
a <- commandArgs(trailingOnly = TRUE)
se_path <- if (length(a) >= 1) a[1] else "/mnt/data/projects/spatial-rads/aggregate/pseudobulk_se.rds"
de_path <- if (length(a) >= 2) a[2] else "results/aggregate/engine/de_engine.tsv"
out_p   <- if (length(a) >= 3) a[3] else "results/aggregate/deseq2_voom_concordance.tsv"

se <- readRDS(se_path)
de <- fread(de_path)                                     # DESeq2 log2FC per gene (estimate)
CONTRASTS <- list(
  MBRT_vs_Ctrl = c("conditionMBRT_day2", "conditionControl"),
  SBRT_vs_Ctrl = c("conditionSBRT_day2", "conditionControl"),
  MBRT_vs_SBRT = c("conditionMBRT_day2", "conditionSBRT_day2"))
MIN_CELLS <- 10L; MIN_SAMPLES <- 3L

rows <- list()
for (ct in unique(colData(se)$cell_type)) {
  sub <- se[, colData(se)$cell_type == ct & colData(se)$n_cells >= MIN_CELLS]
  cd  <- as.data.frame(colData(sub))
  if (length(unique(cd$condition)) < 2 || ncol(sub) < MIN_SAMPLES) next
  cd$condition <- factor(cd$condition, levels = c("Control","MBRT_day2","SBRT_day2"))
  cd$slide_id  <- factor(cd$slide_id)
  design <- model.matrix(~ 0 + condition + slide_id, data = cd)
  colnames(design) <- make.names(colnames(design))
  if (qr(design)$rank < ncol(design)) next
  y <- DGEList(assay(sub, "counts")); y <- calcNormFactors(y)
  v <- voom(y, design)
  fit <- lmFit(v, design)
  for (cn in names(CONTRASTS)) {
    num <- make.names(CONTRASTS[[cn]][1]); den <- make.names(CONTRASTS[[cn]][2])
    if (!all(c(num, den) %in% colnames(design))) next
    cm  <- makeContrasts(contrasts = paste(num, "-", den), levels = design)
    f2  <- eBayes(contrasts.fit(fit, cm))
    tt  <- topTable(f2, number = Inf, sort.by = "none")
    vv  <- data.table(cell_type = ct, contrast = cn, feature_id = rownames(tt),
                      voom_log2fc = tt$logFC)
    dd  <- de[unit == ct & contrast == cn, .(feature_id, deseq2_log2fc = estimate)]
    j   <- merge(vv, dd, by = "feature_id")
    if (nrow(j) >= 10)
      rows[[length(rows)+1]] <- data.table(cell_type = ct, contrast = cn,
        spearman_rho = cor(j$voom_log2fc, j$deseq2_log2fc, method = "spearman"),
        n_genes = nrow(j))
  }
}
conc <- rbindlist(rows)
fwrite(conc, out_p, sep = "\t")
cat(sprintf("voom concordance: %d strata, median rho %.3f (range %.3f-%.3f)\n",
            nrow(conc), median(conc$spearman_rho, na.rm=TRUE),
            min(conc$spearman_rho, na.rm=TRUE), max(conc$spearman_rho, na.rm=TRUE)))
