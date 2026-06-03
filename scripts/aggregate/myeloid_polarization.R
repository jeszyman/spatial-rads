#!/usr/bin/env Rscript
# aggregate.smk myeloid M1/M2 polarization. Re-tests the dev "MBRT day-2 M2 skew"
# finding (dev/mbrt_vs_sbrt/11_m1_m2_polarization.R) on the unified cross-dataset
# labels. Restricts to the Macrophages cell_subtype (the canonical M1/M2 substrate),
# UCell-scores the 16+16 canonical M1/M2 marker panels (panel-filtered, kept genes
# logged), and aggregates per-cell scores to per-sample means + M2/M1 ratio. M02
# day-2 gets a limma test on the per-sample metric matrix (M1, M2, M2/M1 ratio;
# ~0+condition+slide_id, 3 contrasts, 95% CIs, BH) - few outcome rows so eBayes
# moderation is negligible. M01 is descriptive (n=1).
# Args: <merged.rds> <full_labels.parquet> <out_scores.tsv> <out_test.tsv> <plot_m02>
suppressPackageStartupMessages({
  library(Seurat); library(arrow); library(UCell)
  library(limma); library(data.table); library(ggplot2)
})
a <- commandArgs(trailingOnly = TRUE)
merged_path <- a[1]; labels_path <- a[2]
out_scores <- a[3]; out_test <- a[4]; plot_m02 <- a[5]
SEED <- 42L; UCELL_CORES <- 8L
CONDS <- c("Control", "MBRT_day2", "SBRT_day2")

m1_markers <- c("Nos2","Tnf","Il6","Il1b","Il12b","Cxcl9","Cxcl10","Cxcl11","Cd86",
                "Cd80","Tlr2","Tlr4","Stat1","Irf5","Hif1a","Nfkb1")
m2_markers <- c("Cd163","Mrc1","Arg1","Il4ra","Mgl2","Retnla","Chi3l3","Ym1",
                "Stab1","Vegfa","Ccl22","Mmp9","Tgfb1","Stat6","Klf4","Cd206")

dir.create(dirname(out_scores), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(plot_m02), recursive = TRUE, showWarnings = FALSE)

o   <- readRDS(merged_path)
lab <- as.data.table(read_parquet(labels_path))[, .(cell, cell_subtype)]
lv  <- setNames(lab$cell_subtype, lab$cell)
o$cell_subtype <- unname(lv[colnames(o)])
mye <- subset(o, cells = colnames(o)[!is.na(o$cell_subtype) &
                                     o$cell_subtype == "Macrophages"])
rm(o); invisible(gc())

panel <- rownames(mye)
m1_in <- intersect(m1_markers, panel); m2_in <- intersect(m2_markers, panel)
cat(sprintf("Macrophages: %d cells | M1 %d/%d on panel (%s) | M2 %d/%d on panel (%s)\n",
            ncol(mye), length(m1_in), length(m1_markers), paste(m1_in, collapse=","),
            length(m2_in), length(m2_markers), paste(m2_in, collapse=",")))

set.seed(SEED)
mye <- AddModuleScore_UCell(mye, features = list(M1 = m1_in, M2 = m2_in),
                            assay = "RNA", slot = "data",
                            ncores = UCELL_CORES, name = "_score", force.gc = TRUE)

md <- as.data.table(mye@meta.data, keep.rownames = "cell")
score <- md[!is.na(M1_score) & !is.na(M2_score),
            .(M1_mean = mean(M1_score), M2_mean = mean(M2_score),
              M2_M1_ratio = mean(M2_score) / mean(M1_score), n_cells = .N),
            by = .(sample_id, condition, slide_id, timepoint_h, dataset)]
setorder(score, dataset, timepoint_h, sample_id)
fwrite(score, out_scores, sep = "\t")

# --- M02 day2 limma test on the per-sample metric matrix ------------------------
m2 <- score[dataset == "Mutter_02" & timepoint_h == 48L]
m2[, condition := factor(condition, levels = CONDS)]
setkey(m2, sample_id)
metrics <- c("M1_mean", "M2_mean", "M2_M1_ratio")
mat <- t(as.matrix(m2[, ..metrics])); colnames(mat) <- m2$sample_id
design <- model.matrix(~ 0 + condition + slide_id, data = m2)
colnames(design) <- make.names(colnames(design))
fit <- lmFit(mat, design)
cm  <- makeContrasts(
  MBRT_vs_Ctrl = conditionMBRT_day2 - conditionControl,
  SBRT_vs_Ctrl = conditionSBRT_day2 - conditionControl,
  MBRT_vs_SBRT = conditionMBRT_day2 - conditionSBRT_day2, levels = design)
fit2 <- eBayes(contrasts.fit(fit, cm))
test <- rbindlist(lapply(colnames(cm), function(cn) {
  tt <- topTable(fit2, coef = cn, number = Inf, sort.by = "none", confint = TRUE)
  data.table(metric = rownames(tt), contrast = cn, estimate = tt$logFC,
             ci_low = tt$CI.L, ci_high = tt$CI.R, t_stat = tt$t, pvalue = tt$P.Value)
}))
test[, padj := p.adjust(pvalue, method = "BH")]
test[, `:=`(method = "limma", n_samples_per_group = 4L, dataset = "Mutter_02")]
fwrite(test, out_test, sep = "\t")

# --- plot: M1 / M2 / ratio by arm (M02 day2) -----------------------------------
pd <- melt(m2[, c("sample_id", "condition", metrics), with = FALSE],
           id.vars = c("sample_id", "condition"),
           variable.name = "metric", value.name = "value")
p <- ggplot(pd, aes(condition, value, fill = condition)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_point(size = 1, position = position_jitter(width = 0.12)) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_fill_brewer(palette = "Set1") +
  labs(x = NULL, y = NULL,
       title = "M02 day2 macrophage M1/M2 polarization by arm (n=4/arm)") +
  theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(plot_m02, p, width = 8, height = 4.5, dpi = 150)

cat(sprintf("myeloid M1/M2: %d samples | M02 day2 test %d rows, %d padj<0.05 | M2/M1 range %.3f-%.3f\n",
            nrow(score), nrow(test), test[padj < 0.05, .N],
            min(m2$M2_M1_ratio), max(m2$M2_M1_ratio)))
