#!/usr/bin/env Rscript
# QC visualizations as pipeline outputs: (1) per-criterion removal attribution (cells failing
# counts / complexity / propNegative / area -- shows which criteria drive the cuts), and
# (2) metric distributions with cutoff lines. Cutoffs read from config (match scripts/qc_filter.R).
# Args: <config.yaml> <out_dir>
suppressMessages({library(Seurat); library(dplyr); library(tidyr); library(ggplot2); library(patchwork); library(readr); library(yaml)})
args <- commandArgs(trailingOnly = TRUE)
CFG <- args[1]; OUT <- args[2]
cfg <- yaml::read_yaml(CFG); q <- cfg$qc; datadir <- cfg$datadir
min_counts <- q$min_counts; min_features <- q$min_features
max_propneg <- q$max_prop_negative; area_nmads <- q$area_nmads
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

raw_files <- sort(list.files(file.path(datadir, "processing", "raw"), pattern = "\\.raw\\.rds$", full.names = TRUE))
stopifnot(length(raw_files) > 0)

set.seed(1)
attr_rows <- list(); metric_rows <- list()
for (f in raw_files) {
  m <- readRDS(f)@meta.data
  sid <- if ("sample_id" %in% names(m)) as.character(m$sample_id[1]) else sub("\\.raw\\.rds$", "", basename(f))
  ds  <- if ("dataset" %in% names(m)) as.character(m$dataset[1]) else NA_character_
  la  <- log10(m$Area); med <- median(la, na.rm = TRUE); md <- mad(la, na.rm = TRUE)
  attr_rows[[sid]] <- data.frame(
    sample_id = sid, dataset = ds, n_cells = nrow(m),
    counts       = sum(!(m$nCount_RNA   > min_counts)),
    complexity   = sum(!(m$nFeature_RNA > min_features)),
    propNegative = sum(!(m$propNegative < max_propneg)),
    area         = sum(is.na(la) | !(la > med - area_nmads * md & la < med + area_nmads * md)))
  idx <- sample(seq_len(nrow(m)), min(2000L, nrow(m)))
  metric_rows[[sid]] <- data.frame(sample_id = sid, dataset = ds,
    nCount_RNA = m$nCount_RNA[idx], nFeature_RNA = m$nFeature_RNA[idx],
    propNegative = m$propNegative[idx], Area = m$Area[idx])
  rm(m); gc()
}
attr_df <- bind_rows(attr_rows); met <- bind_rows(metric_rows)
write_tsv(attr_df, file.path(OUT, "qc_removal_attribution.tsv"))

## Plot A -- per-criterion removal attribution
aL <- attr_df %>%
  select(sample_id, dataset, counts, complexity, propNegative, area) %>%
  pivot_longer(c(counts, complexity, propNegative, area), names_to = "criterion", values_to = "n_fail")
pA <- ggplot(aL, aes(sample_id, n_fail, fill = criterion)) +
  geom_col() +
  facet_grid(~dataset, scales = "free_x", space = "free_x") +
  labs(title = "Cells failing each QC criterion (not mutually exclusive)",
       subtitle = "Count + complexity floors drive most removal", x = NULL, y = "cells flagged") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))

## Plot B -- metric distributions with cutoffs
mk <- function(col, cutoff, logy, lab) {
  p <- ggplot(met, aes(dataset, .data[[col]], fill = dataset)) + geom_violin(scale = "width")
  if (!is.na(cutoff)) p <- p + geom_hline(yintercept = cutoff, linetype = "dashed", color = "red")
  p <- p + labs(x = NULL, y = lab) + theme_bw() + theme(legend.position = "none")
  if (logy) p <- p + scale_y_log10()
  p
}
pB <- (mk("nCount_RNA", min_counts, TRUE, "nCount_RNA (log)") | mk("nFeature_RNA", min_features, TRUE, "nFeature_RNA (log)")) /
      (mk("propNegative", max_propneg, FALSE, "propNegative") | mk("Area", NA, TRUE, "Area (log, per-sample MAD gate)"))

ggsave(file.path(OUT, "qc_removal_attribution.png"), pA, width = 11, height = 5, dpi = 150)
ggsave(file.path(OUT, "qc_metric_distributions.png"), pB, width = 10, height = 7, dpi = 150)
cat(sprintf("wrote QC plots + attribution to %s\n", OUT)); print(attr_df)
