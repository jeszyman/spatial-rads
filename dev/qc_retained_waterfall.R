#!/usr/bin/env Rscript
# Cells-retained QC funnel: per-sample cells before vs after the per-cell QC gates
# (counts / features / negative-probe fraction / area), one stacked bar per sample.
# Retained is colored by cohort, removed is grey; % retained labeled at the bar end.
# Complements qc_metrics_per_sample.R (which shows WHERE the gates sit) by showing the
# NET attrition per sample. Reads results/processing/qc_summary.tsv + samplesheet.tsv.
suppressMessages({library(data.table); library(ggplot2); library(scales); library(stringr)})
source("/home/jeszyman/repos/science/R/figure_schema.R")

QC  <- "results/processing/qc_summary.tsv"
SS  <- "results/processing/samplesheet.tsv"
OUT <- "results/aggregate/plots/qc_retained_waterfall"
stopifnot(file.exists(QC), file.exists(SS))

d  <- fread(QC)
ss <- unique(fread(SS)[, .(sample_id, dataset)])
d  <- merge(d, ss, by = "sample_id", all.x = TRUE)
d[, removed := pre_filter - post_filter]

n_samp <- nrow(d)
rng <- range(d$pct_retained); med <- median(d$pct_retained)
ds_lvls <- sort(unique(d$dataset))
ds_n <- d[, .N, by = dataset][order(dataset)]
ds_str <- paste(sprintf("%s n=%d", ds_n$dataset, ds_n$N), collapse = ", ")

# Order samples by retention (highest at top). Long form: retained + removed segments.
lvls <- d$sample_id[order(d$pct_retained)]
m <- melt(d, id.vars = c("sample_id", "dataset", "pct_retained", "pre_filter"),
          measure.vars = c("post_filter", "removed"),
          variable.name = "segment", value.name = "n_cells")
m[, sample_id := factor(sample_id, levels = lvls)]
# Fill key: retained segment carries the cohort hue; removed is always grey. Levels
# ordered so the dataset (retained) stacks at the base and Removed caps the bar.
m[, fillkey := fifelse(segment == "post_filter", dataset, "Removed")]
m[, fillkey := factor(fillkey, levels = c(ds_lvls, "Removed"))]

fill_cols <- c(setNames(c("#1b6ca8", "#e07b39")[seq_along(ds_lvls)], ds_lvls),
               Removed = "grey80")

lab <- d[, .(sample_id = factor(sample_id, levels = lvls), pre_filter,
             txt = paste0(round(pct_retained), "%"))]

key <- paste(
  sprintf("All %d samples retain %.0f-%.0f%% of cells after per-cell QC filtering (median %.0f%%); the lowest-retention samples sit at the bottom for review.",
          n_samp, rng[1], rng[2], med),
  sprintf("Mutter_01 + Mutter_02 (%s), CosMx 950-gene common panel.", ds_str),
  "Bars are cells per sample: retained (colored by cohort) stacked with removed (grey); % retained labeled at the bar end.",
  "Removed = cells failing the per-cell gates (counts, features, negative-probe fraction, cell area).")

p <- ggplot(m, aes(n_cells, sample_id, fill = fillkey)) +
  geom_col(width = 0.72, position = position_stack(reverse = TRUE)) +
  geom_text(data = lab, aes(pre_filter, sample_id, label = txt),
            inherit.aes = FALSE, hjust = -0.18, size = 2.8, colour = "gray25") +
  scale_fill_manual(values = fill_cols, breaks = c(ds_lvls, "Removed")) +
  scale_x_continuous(breaks = seq(0, 4e5, 1e5),
                     labels = scales::label_number(scale = 1e-3, suffix = "k"),
                     expand = expansion(mult = c(0, 0.13))) +
  labs(x = "Cells", y = NULL, fill = NULL, caption = str_wrap(key, 108)) +
  theme_scifig(base_size = 11) +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(margin = margin(t = 8)),
        plot.caption = element_text(size = 9, hjust = 0, colour = "gray30", margin = margin(t = 12)),
        plot.caption.position = "plot")

save_plot(p, OUT, w = 7, h = 6.4)
cat(sprintf("QC waterfall: %d samples | retention %.1f-%.1f%% (median %.1f) | %s\n",
            n_samp, rng[1], rng[2], med, ds_str))
