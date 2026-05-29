#!/usr/bin/env Rscript
# Aggregate per-sample QC summaries into one table. Args: <qc_summary_dir> <out.tsv>
suppressMessages({library(readr); library(dplyr); library(purrr)})
args <- commandArgs(trailingOnly = TRUE)
DIR <- args[1]; OUT <- args[2]
files <- list.files(DIR, pattern = "\\.qcsummary\\.tsv$", full.names = TRUE)
stopifnot(length(files) > 0)
df <- map_dfr(files, read_tsv, show_col_types = FALSE) %>% arrange(sample_id)
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
write_tsv(df, OUT)
cat(sprintf("aggregated %d sample QC rows -> %s\n", nrow(df), OUT))
print(as.data.frame(df))
