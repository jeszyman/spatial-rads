#!/usr/bin/env Rscript
# Cross-arm QC balance (report-only): the confound check on the day-2 composition (fraction-shift)
# result. Joins the per-sample technical metrics (SpatialQM sensitivity/SNR/specificityFDR) + the
# per-sample MECR contamination metric to the arm design, and asks whether each metric is BALANCED
# across Control/MBRT/SBRT for the M02 day-2 cohort (n=4/arm). If sensitivity/SNR/MECR overlap across
# arms, the "more cells expressing collagen/Acta2" signal is not a sensitivity/contamination artifact;
# if an arm's mean separates beyond within-arm spread, the composition claim is confounded on that axis.
# Effect-size-forward (n=4): per-arm mean +/- sd, between-arm gap vs within-arm sd -- no p-value theater.
# Tabular terminus: the per-metric verdict table + the tidy per-sample table the figure consumes;
# plotting lives in scripts/fig_qc_arm_balance.R.
# Args: <sample_tech_metrics.tsv> <contamination_qc.tsv> <samples.tsv> <out_verdict.tsv> <out_samples.tsv>
suppressPackageStartupMessages({library(tidyverse)})
a <- commandArgs(trailingOnly = TRUE)
TECH <- a[1]; CONTAM <- a[2]; SAMPLES <- a[3]; OUT_TSV <- a[4]; OUT_SAMPLES <- a[5]

METRICS <- c(tpc = "TPC (transcripts/cell)", med_nFeature = "genes/cell (median)",
             sparsity = "sparsity", snr = "SNR (log10)", specificity_fdr = "specificity FDR",
             mecr = "MECR (contamination)")

tech   <- read_tsv(TECH, show_col_types = FALSE)
contam <- read_tsv(CONTAM, show_col_types = FALSE) %>% select(sample_id, mecr)
ss     <- read_tsv(SAMPLES, show_col_types = FALSE) %>%
  transmute(sample_id, name, timepoint_h, model,
            arm = recode(treatment, NT = "Control"))

# M02 day-2 flank cohort: the only arm-replicated design (n=4/arm).
d <- tech %>%
  left_join(contam, by = "sample_id") %>%
  inner_join(ss, by = "sample_id") %>%
  filter(name == "Mutter_02", timepoint_h == 48, model == "flank") %>%
  mutate(arm = factor(arm, levels = c("Control", "MBRT", "SBRT")))
stopifnot(nrow(d) > 0, all(table(d$arm) > 1))

long <- d %>%
  select(sample_id, arm, all_of(names(METRICS))) %>%
  pivot_longer(all_of(names(METRICS)), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = names(METRICS)),
         metric_label = METRICS[as.character(metric)])

# Per-metric arm-balance readout: arm means/sds, largest pairwise arm gap vs the mean within-arm sd.
bal <- long %>%
  group_by(metric, arm) %>%
  summarise(mean = mean(value), sd = sd(value), .groups = "drop") %>%
  group_by(metric) %>%
  summarise(
    ctrl_mean = mean[arm == "Control"], mbrt_mean = mean[arm == "MBRT"], sbrt_mean = mean[arm == "SBRT"],
    within_sd    = mean(sd, na.rm = TRUE),
    max_arm_gap  = max(dist(c(ctrl_mean, mbrt_mean, sbrt_mean))),
    gap_sd       = max_arm_gap / within_sd,          # largest between-arm gap in within-arm sd units
    balanced     = gap_sd <= 1, .groups = "drop") %>% # heuristic: arm means within 1 sd = balanced
  mutate(metric_label = METRICS[as.character(metric)]) %>%
  arrange(desc(gap_sd))

write_tsv(bal, OUT_TSV)
write_tsv(long, OUT_SAMPLES)

cat("qc_arm_balance (M02 day-2):\n"); print(as.data.frame(bal), digits = 3)
cat(sprintf("balanced on %d/%d metrics; wrote %s + %s\n",
            sum(bal$balanced), nrow(bal), OUT_TSV, OUT_SAMPLES))
