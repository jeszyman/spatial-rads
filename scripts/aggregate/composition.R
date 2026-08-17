#!/usr/bin/env Rscript
# aggregate.smk Track 1 -- cell-type composition.
# M02 day2 (n=4/condition, 4 slides crossed) gets a formal propeller test:
# getTransformedProps(logit) -> lmFit(~0+condition+slide_id) -> 3 contrasts ->
# eBayes(robust) -> topTable, with global BH across (cell_type x contrast).
# M01 (n=1 timecourse) is descriptive proportions only. Per plan-aggregate.md
# Track 1: no silent drops -- a cell type is logged to the dropped table only if
# the logit transform yields a non-finite row (genuine method failure).
# Args: <obs.parquet> <full_labels.parquet> <out_by_sample> <out_test>
#       <out_dropped> <plot_bars> <plot_forest> <plot_timecourse>
suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(speckle)
  library(limma)
  library(ggplot2)
})

args             <- commandArgs(trailingOnly = TRUE)
obs_path         <- args[1]
labels_path      <- args[2]
out_by_sample    <- args[3]
out_test         <- args[4]
out_dropped      <- args[5]
plot_bars        <- args[6]
plot_forest      <- args[7]
plot_timecourse  <- args[8]
out_sens         <- args[9]
out_lm_input     <- args[10]   # per-cell engine input; the shared lm_engine produces the
                               # canonical composition arm test for results_master. The inline
                               # propeller below is retained only for this producer's own forest
                               # plot + unassigned-sensitivity diagnostics (not the master family).

dir.create(dirname(out_test), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(plot_bars), recursive = TRUE, showWarnings = FALSE)

m <- as.data.table(read_parquet(obs_path))   # cell, sample_id, dataset, slide_id, condition, treatment, timepoint_h, cell_type(stale), neg
# Use the unified cross-dataset cell_subtype (tier-1/2 joined) in place of the
# stale per-sample cell_type, so composition reflects identical typing (plan-aggregate.md).
lab <- as.data.table(read_parquet(labels_path))[, .(cell, cell_type = cell_subtype)]
# `Epithelial cells` is the tier-2 immune tumour-contamination bucket (1,180 cells, 0.04%
# of the cohort) -- epithelial cells that clustered with immune cells, not a lineage of
# its own. Fold it into Tumor at the point of label read so the census and every fit below
# see one tumour stratum. Collapsing here is what makes the collapse exact: renaming an
# already-fitted stratum in the results tier instead yields two rows both claiming to be
# Tumor. verify_arm_tables.R gates the absence of the alias downstream.
lab[cell_type == "Epithelial cells", cell_type := "Tumor"]
m[, cell_type := NULL]
m <- merge(m, lab, by = "cell", all.x = TRUE)
m <- m[!is.na(cell_type)]
m[, condition := as.character(condition)]

# --- composition_by_sample: every (sample x cell_type) proportion ---------
by_sample <- m[, .(n_cells = .N), by = .(sample_id, cell_type, condition,
                                          timepoint_h, dataset, slide_id)]
by_sample[, fraction := n_cells / sum(n_cells), by = sample_id]
setcolorder(by_sample, c("sample_id", "cell_type", "n_cells", "fraction",
                         "condition", "timepoint_h", "dataset", "slide_id"))
setorder(by_sample, dataset, sample_id, -n_cells)
fwrite(by_sample, out_by_sample, sep = "\t")

# --- engine input: per-cell labels (cell, sample_id, label, condition, slide_id), M02 day2.
# The canonical composition arm test that enters results_master runs in the shared lm_engine
# (proportion/logit path, robust eBayes) on this; the inline propeller below stays only for
# this producer's forest plot + unassigned-sensitivity diagnostics. ---
lm_input <- m[!grepl("^Tongue", condition), .(cell, sample_id, label = cell_type, condition, slide_id)]
fwrite(lm_input, out_lm_input, sep = "\t")

# --- M02 day2 propeller test (inline, retained for forest plot only) -------
m2 <- m[dataset == "Mutter_02" & as.integer(timepoint_h) == 48L]
m2[, condition := factor(condition, levels = c("Control", "MBRT_day2", "SBRT_day2"))]

props <- getTransformedProps(clusters = m2$cell_type, sample = m2$sample_id,
                             transform = "logit")
tp <- props$TransformedProps                                   # celltype x sample

# Genuine method failures only: a cell type whose logit row is non-finite.
nonfinite <- apply(tp, 1, function(r) any(!is.finite(r)))
mean_cells <- rowMeans(props$Counts)                           # per-sample mean count
dropped <- data.table(
  cell_type   = rownames(tp)[nonfinite],
  reason      = "non-finite logit-transformed proportion (zero cells in >=1 sample)",
  mean_n_cells = round(mean_cells[nonfinite], 1)
)
fwrite(dropped, out_dropped, sep = "\t")

tp <- tp[!nonfinite, , drop = FALSE]

samp <- unique(m2[, .(sample_id, condition, slide_id)])
setkey(samp, sample_id)
samp <- samp[colnames(tp)]                                     # align to tp columns
stopifnot(identical(samp$sample_id, colnames(tp)))
design <- model.matrix(~ 0 + condition + slide_id, data = samp)
colnames(design) <- make.names(colnames(design))

fit <- lmFit(tp, design)
cm  <- makeContrasts(
  MBRT_vs_Ctrl = conditionMBRT_day2 - conditionControl,
  SBRT_vs_Ctrl = conditionSBRT_day2 - conditionControl,
  MBRT_vs_SBRT = conditionMBRT_day2 - conditionSBRT_day2,
  levels = design)
fit2 <- eBayes(contrasts.fit(fit, cm), robust = TRUE)

mean_n_cells <- rowMeans(props$Counts)[rownames(tp)]
ln2 <- log(2)
test <- rbindlist(lapply(colnames(cm), function(cn) {
  tt <- topTable(fit2, coef = cn, number = Inf, sort.by = "none", confint = TRUE)
  data.table(
    cell_type    = rownames(tt),
    contrast     = cn,
    log2FC_logit = tt$logFC / ln2,
    ci_low_log2  = tt$CI.L  / ln2,
    ci_high_log2 = tt$CI.R  / ln2,
    t_stat       = tt$t,
    pvalue       = tt$P.Value)
}))
test[, padj := p.adjust(pvalue, method = "BH")]                # global across all rows
nmin <- min(table(unique(m2[, .(sample_id, condition)])$condition))
test[, `:=`(method = "propeller", n_samples_per_group = as.integer(nmin),
            mean_n_cells = round(mean_n_cells[cell_type], 1), dataset = "Mutter_02")]
setcolorder(test, c("cell_type", "contrast", "log2FC_logit", "ci_low_log2",
                    "ci_high_log2", "t_stat", "pvalue", "padj", "method",
                    "n_samples_per_group", "mean_n_cells", "dataset"))
fwrite(test, out_test, sep = "\t")

# --- unassigned sensitivity (Fix 1): rerun the propeller test with unassigned dropped ---
# Composition is a closed system, so a treatment-shifted unassigned fraction distorts every
# labelled type. Report each labelled type's effect with vs without unassigned; flag flips.
m2x   <- m2[cell_type != "unassigned"]
propx <- getTransformedProps(clusters = m2x$cell_type, sample = m2x$sample_id, transform = "logit")
tpx   <- propx$TransformedProps
tpx   <- tpx[apply(tpx, 1, function(r) all(is.finite(r))), , drop = FALSE]
sx    <- unique(m2x[, .(sample_id, condition, slide_id)]); setkey(sx, sample_id); sx <- sx[colnames(tpx)]
dx    <- model.matrix(~ 0 + condition + slide_id, data = sx); colnames(dx) <- make.names(colnames(dx))
cmx   <- makeContrasts(
  MBRT_vs_Ctrl = conditionMBRT_day2 - conditionControl,
  SBRT_vs_Ctrl = conditionSBRT_day2 - conditionControl,
  MBRT_vs_SBRT = conditionMBRT_day2 - conditionSBRT_day2, levels = dx)
fx    <- eBayes(contrasts.fit(lmFit(tpx, dx), cmx), robust = TRUE)
excl  <- rbindlist(lapply(colnames(cmx), function(cn) {
  tt <- topTable(fx, coef = cn, number = Inf, sort.by = "none")
  data.table(cell_type = rownames(tt), contrast = cn,
             log2FC_excl = tt$logFC / log(2), p_excl = tt$P.Value)
}))
excl[, padj_excl := p.adjust(p_excl, "BH")]
sens <- merge(test[, .(cell_type, contrast, log2FC_incl = log2FC_logit, padj_incl = padj)],
              excl[, .(cell_type, contrast, log2FC_excl, padj_excl)],
              by = c("cell_type", "contrast"))
sens[, flip := (padj_incl < 0.05) != (padj_excl < 0.05) | sign(log2FC_incl) != sign(log2FC_excl)]
setorder(sens, contrast, cell_type)
fwrite(sens, out_sens, sep = "\t")
cat(sprintf("unassigned sensitivity: %d labelled-type rows, %d flip inclusion regimes\n",
            nrow(sens), sens[flip == TRUE, .N]))

# --- plot 1: M02 stacked composition bars (top-15 types + Other) ----------
m2_bys <- by_sample[dataset == "Mutter_02" & as.integer(timepoint_h) == 48L]
top15  <- m2_bys[, .(tot = sum(n_cells)), by = cell_type][order(-tot)][1:15, cell_type]
m2_bys[, ct_lab := ifelse(cell_type %in% top15, cell_type, "Other")]
lvls   <- c(top15, "Other")
m2_bys[, ct_lab := factor(ct_lab, levels = lvls)]
pal    <- c(scales::hue_pal()(15), Other = "grey70"); names(pal) <- lvls

p_bars <- ggplot(m2_bys, aes(sample_id, fraction, fill = ct_lab)) +
  geom_col(width = 0.9) +
  facet_wrap(~ condition, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = pal, name = "cell type") +
  labs(x = NULL, y = "fraction of cells",
       title = "M02 day2 cell-type composition (top 15 + Other)") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0.35, "cm"))
ggsave(plot_bars, p_bars, width = 9, height = 5.5, dpi = 150)

# --- plot 2: forest of the three contrasts (all tested cell types) --------
ord <- test[contrast == "MBRT_vs_Ctrl"][order(mean_n_cells), cell_type]
test[, cell_type := factor(cell_type, levels = ord)]
test[, contrast  := factor(contrast,
       levels = c("MBRT_vs_Ctrl", "SBRT_vs_Ctrl", "MBRT_vs_SBRT"))]

p_forest <- ggplot(test, aes(log2FC_logit, cell_type, colour = contrast)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = ci_low_log2, xmax = ci_high_log2),
                orientation = "y", width = 0,
                position = position_dodge(width = 0.7)) +
  geom_point(size = 1.4, position = position_dodge(width = 0.7)) +
  scale_colour_brewer(palette = "Set1") +
  labs(x = "log2 fold-change (logit proportion), 95% CI", y = NULL,
       title = "M02 day2 composition contrasts (propeller, global BH)") +
  theme_bw(base_size = 9)
ggsave(plot_forest, p_forest, width = 7.5, height = 10, dpi = 150)

# --- plot 3: M01 descriptive timecourse (fraction vs timepoint) -----------
# Collapse to the top-15 most abundant M01 cell types + Other (mirrors the M02
# bars plot) so the panel is legible instead of ~40 thumbnail facets. Per
# sample_id, fractions sum within a label, so Other is the combined residual.
treat_lookup <- unique(m[, .(sample_id, treatment)])           # treatment not in by_sample schema
m1_bys <- merge(by_sample[dataset == "Mutter_01"], treat_lookup, by = "sample_id")
m1_top <- m1_bys[, .(tot = sum(n_cells)), by = cell_type][order(-tot)][1:15, cell_type]
m1_bys[, ct_lab := ifelse(cell_type %in% m1_top, cell_type, "Other")]
m1_tc  <- m1_bys[, .(fraction = sum(fraction)),
                 by = .(sample_id, timepoint_h, treatment, ct_lab)]
m1_tc[, ct_lab    := factor(ct_lab, levels = c(m1_top, "Other"))]
m1_tc[, treatment := factor(treatment)]
p_tc <- ggplot(m1_tc, aes(timepoint_h, fraction, colour = treatment)) +
  geom_line(aes(group = treatment), linewidth = 0.5) +
  geom_point(size = 1.2) +
  facet_wrap(~ ct_lab, scales = "free_y", ncol = 4) +
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "timepoint (h)", y = "fraction of cells",
       title = "M01 cell-type composition over time (top 15 + Other; n=1, descriptive)") +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(size = 9))
ggsave(plot_timecourse, p_tc, width = 11, height = 9, dpi = 150)

cat(sprintf("composition: %d samples, %d cell types | M02 test %d rows, %d dropped | %d padj<0.05\n",
            uniqueN(m$sample_id), uniqueN(m$cell_type), nrow(test),
            nrow(dropped), test[padj < 0.05, .N]))
