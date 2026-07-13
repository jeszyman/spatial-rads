#!/usr/bin/env Rscript
# Task 1 TDD: the resolved registry must carry per-contrast levels + transform for the
# test engines. Asserts the three day-2 arm contrasts resolve to the correct
# numerator/denominator condition levels (denominator = DESeq2 reference).
suppressPackageStartupMessages(library(data.table))
tsv <- "results/data_model/comparisons.tsv"
stopifnot(file.exists(tsv))
d <- fread(tsv)
for (col in c("contrast_num_level", "contrast_den_level", "ref_level", "transform"))
  stopifnot(col %in% names(d))
row <- d[name == "MBRT_vs_SBRT" & cohort == "mutter02_day2" & resolution == "whole"]
stopifnot(nrow(row) == 1,
          row$contrast_num_level == "MBRT_day2",
          row$contrast_den_level == "SBRT_day2",
          row$ref_level == "SBRT_day2")
ctrl <- d[name == "SBRT_vs_Ctrl" & cohort == "mutter02_day2" & resolution == "whole"]
stopifnot(ctrl$contrast_den_level == "Control", ctrl$ref_level == "Control")
mbrt <- d[name == "MBRT_vs_Ctrl" & cohort == "mutter02_day2" & resolution == "whole"]
stopifnot(mbrt$contrast_num_level == "MBRT_day2", mbrt$contrast_den_level == "Control")
cat("PASS test_comparisons_resolve\n")
