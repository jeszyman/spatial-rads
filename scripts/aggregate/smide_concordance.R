#!/usr/bin/env Rscript
# Pseudobulk vs smiDE concordance analysis. Reads results_master.tsv (which
# contains both readout_class="DE" and readout_class="smiDE_*" rows) and compares
# effect sizes, hit overlap, and contamination-ratio enrichment, scoped to one
# cohort (mutter02_day2 / combined_4h / combined_4h_treated).
# Args: <results_master.tsv> <overlap_ratio_qc.tsv> <cohort> <out_concordance.tsv> <out_plot.png>
suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
})
a <- commandArgs(trailingOnly = TRUE)
master_path <- a[1]; overlap_path <- a[2]; cohort <- a[3]; out_tsv <- a[4]; out_plot <- a[5]

dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_plot), recursive = TRUE, showWarnings = FALSE)

master_raw <- fread(master_path)
# Engine outputs (composition/DE/niche/mixing/myeloid/substate/smiDE/pathway/gsea)
# all carry a `comparison` column set to the cohort name; once results_master.tsv
# is assembled from all cohorts it carries the same column. A master table built
# before cross-cohort assembly has no `comparison` column and is implicitly
# mutter02_day2-only.
has_comparison <- "comparison" %in% names(master_raw)
if (has_comparison) {
  master <- master_raw[comparison == cohort]
} else if (cohort == "mutter02_day2") {
  master <- master_raw
} else {
  master <- master_raw[0]
}

if (nrow(master) == 0) {
  reason <- if (has_comparison) {
    sprintf("no rows with comparison == '%s'", cohort)
  } else {
    sprintf("results_master.tsv has no comparison column (single-cohort, mutter02_day2-only) but cohort '%s' was requested", cohort)
  }
  cat(sprintf("smide_concordance: no results_master.tsv data for cohort '%s' -- %s\n", cohort, reason))
  quit(save = "no", status = 1)
}

pb <- master[readout_class == "DE", .(unit, feature, contrast,
  pb_effect = effect, pb_pvalue = pvalue, pb_padj = padj_own)]

# smiDE enters the master under two readout classes: "smiDE_spatial" (Gaussian-process
# fits meta-analyzed across spatial units) and "smiDE_screen" (nebula NB GLMM with a
# sample random intercept, discovery only). Compare against the spatial fits when
# they exist and fall back to the screen otherwise, naming which was used.
SMIDE_CLASS <- if (master[readout_class == "smiDE_spatial", .N] > 0) "smiDE_spatial" else "smiDE_screen"
cat(sprintf("smide_concordance: comparing pseudobulk against readout_class '%s'\n", SMIDE_CLASS))
sm <- master[readout_class == SMIDE_CLASS, .(unit, feature, contrast,
  sm_effect = effect, sm_pvalue = pvalue, sm_padj = padj_own)]

merged <- pb[sm, on = .(unit, feature, contrast), nomatch = NULL]

if (nrow(merged) == 0) {
  cat("smide_concordance: 0 shared rows between pseudobulk and smiDE\n")
  fwrite(data.table(metric = character(), value = numeric()), out_tsv, sep = "\t")
  # Snakemake declares out_plot a rule output; a missing file is a
  # MissingOutputException even on an intentional early exit, so write a minimal
  # placeholder in place of the effect-size scatter.
  p <- ggplot() +
    annotate("text", x = 0, y = 0,
             label = sprintf("No shared pseudobulk/smiDE rows for cohort '%s'", cohort)) +
    theme_void()
  ggsave(out_plot, p, width = 10, height = 4, dpi = 150)
  quit(save = "no")
}

orm <- fread(overlap_path)
setnames(orm, old = c("gene", "cell_subtype", "ratio"),
         new = c("feature", "unit", "contamination_ratio"), skip_absent = TRUE)
merged <- orm[, .(feature, unit, contamination_ratio)][merged, on = .(feature, unit)]

# --- 1. effect-size concordance per contrast ---
concord <- merged[!is.na(pb_effect) & !is.na(sm_effect),
  .(n = .N,
    spearman = cor(pb_effect, sm_effect, method = "spearman", use = "complete.obs"),
    pearson = cor(pb_effect, sm_effect, method = "pearson", use = "complete.obs"),
    sign_agree_frac = mean(sign(pb_effect) == sign(sm_effect))),
  by = contrast]

# --- 2. hit classification at padj < 0.05 ---
ALPHA <- 0.05
merged[, `:=`(
  pb_sig = !is.na(pb_padj) & pb_padj < ALPHA,
  sm_sig = !is.na(sm_padj) & sm_padj < ALPHA
)]
merged[, hit_class := fcase(
  pb_sig & sm_sig, "both",
  pb_sig & !sm_sig, "pseudobulk_only",
  !pb_sig & sm_sig, "smide_only",
  default = "neither"
)]

hit_counts <- merged[, .N, by = .(contrast, hit_class)]

# --- 3. contamination-ratio enrichment in discordant hits ---
contam_by_class <- merged[!is.na(contamination_ratio),
  .(median_ratio = median(contamination_ratio),
    frac_ratio_ge1 = mean(contamination_ratio >= 1),
    n = .N),
  by = hit_class]

# wilcoxon: pseudobulk-only vs both (are pseudobulk-only hits more contaminated?)
pb_only_ratios <- merged[hit_class == "pseudobulk_only" & !is.na(contamination_ratio), contamination_ratio]
both_ratios <- merged[hit_class == "both" & !is.na(contamination_ratio), contamination_ratio]
if (length(pb_only_ratios) >= 3 & length(both_ratios) >= 3) {
  wt <- wilcox.test(pb_only_ratios, both_ratios, alternative = "greater")
  contam_enrichment_p <- wt$p.value
} else {
  contam_enrichment_p <- NA_real_
}

# --- 4. confirmatory hit survival ---
confirm <- master[tier == "confirmatory" & readout_class == "DE",
                  .(unit, feature, contrast, pb_effect = effect, pb_padj = padj_confirmatory)]
confirm_sm <- sm[confirm, on = .(unit, feature, contrast)]
confirm_sm[, survives := !is.na(sm_padj) & sm_padj < ALPHA &
             sign(sm_effect) == sign(pb_effect)]

# --- output ---
results <- rbindlist(list(
  concord[, .(metric = paste0("spearman_", contrast), value = spearman)],
  concord[, .(metric = paste0("pearson_", contrast), value = pearson)],
  concord[, .(metric = paste0("sign_agree_", contrast), value = sign_agree_frac)],
  concord[, .(metric = paste0("n_shared_", contrast), value = as.numeric(n))],
  data.table(metric = "n_total_shared", value = nrow(merged)),
  data.table(metric = "n_both_sig", value = sum(merged$hit_class == "both")),
  data.table(metric = "n_pseudobulk_only", value = sum(merged$hit_class == "pseudobulk_only")),
  data.table(metric = "n_smide_only", value = sum(merged$hit_class == "smide_only")),
  contam_by_class[, .(metric = paste0("contam_median_ratio_", hit_class), value = median_ratio)],
  contam_by_class[, .(metric = paste0("contam_frac_ge1_", hit_class), value = frac_ratio_ge1)],
  data.table(metric = "contam_enrichment_p_pbonly_vs_both", value = contam_enrichment_p),
  data.table(metric = "n_confirmatory_tested", value = nrow(confirm_sm)),
  data.table(metric = "n_confirmatory_survives_smide",
             value = as.numeric(sum(confirm_sm$survives, na.rm = TRUE)))
), fill = TRUE)
results[, smide_class := SMIDE_CLASS]
fwrite(results, out_tsv, sep = "\t")

# --- plot: effect-size scatter, colored by hit class ---
pd <- merged[!is.na(pb_effect) & !is.na(sm_effect)]
p <- ggplot(pd, aes(pb_effect, sm_effect, color = hit_class)) +
  geom_point(alpha = 0.3, size = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  facet_wrap(~ contrast, scales = "free") +
  scale_color_manual(values = c(both = "grey30", pseudobulk_only = "#E41A1C",
                                smide_only = "#377EB8", neither = "grey80")) +
  labs(x = "Pseudobulk DESeq2 log2FC", y = "smiDE per-cell log2FC",
       title = sprintf("Pseudobulk vs %s effect-size concordance (%s)", SMIDE_CLASS, cohort)) +
  theme_bw(base_size = 10)
ggsave(out_plot, p, width = 10, height = 4, dpi = 150)

cat(sprintf("smide_concordance: %d shared rows | effect Spearman %s | hits: both=%d, pb_only=%d, sm_only=%d\n",
            nrow(merged),
            paste(concord$contrast, round(concord$spearman, 3), sep = "=", collapse = ", "),
            sum(merged$hit_class == "both"),
            sum(merged$hit_class == "pseudobulk_only"),
            sum(merged$hit_class == "smide_only")))
cat(sprintf("confirmatory: %d / %d survive smiDE | contam enrichment p=%.4g\n",
            sum(confirm_sm$survives, na.rm = TRUE), nrow(confirm_sm), contam_enrichment_p))
