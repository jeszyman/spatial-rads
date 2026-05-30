#!/usr/bin/env Rscript
# Aggregate per-sample cell-type counts into one table. Args: <celltype_dir> <out.tsv>
suppressMessages({library(readr); library(dplyr); library(purrr); library(tidyr)})
args <- commandArgs(trailingOnly = TRUE)
DIR <- args[1]; OUT <- args[2]
files <- list.files(DIR, pattern = "\\.celltype\\.tsv$", full.names = TRUE)
stopifnot(length(files) > 0)
df <- map_dfr(files, read_tsv, show_col_types = FALSE)
wide <- df %>% group_by(dataset, sample_id, cell_type) %>% summarise(n = sum(n), .groups = "drop") %>%
  pivot_wider(names_from = cell_type, values_from = n, values_fill = 0) %>% arrange(dataset, sample_id)
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
write_tsv(wide, OUT)
cat(sprintf("aggregated cell types for %d samples -> %s\n", nrow(wide), OUT))
