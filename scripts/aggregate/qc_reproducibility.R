#!/usr/bin/env Rscript
# Replicate reproducibility QC (report-only): the n=4/arm M02 day-2 design is what the SpatialQM /
# Spatial Touchstone reproducibility metrics (Plummer et al., Nat Biotechnol 2025) were built for.
# Before trusting the composition inference, confirm no single outlier slide drives an arm's effect,
# two ways: (1) pseudobulk Spearman concordance of each sample to its arm-mates (SpatialQM getCorrelation
# -- per-sample gene pseudobulk, cell types collapsed); (2) a PCA of the scaled per-sample technical
# metrics (paper Fig 3a) -- a sample sitting apart from its arm is a candidate batch/outlier.
# Tabular terminus (per-sample concordance + PCA coords + PC variance-explained columns);
# plotting lives in scripts/fig_qc_reproducibility.R.
# Args: <pseudobulk_se.rds> <sample_tech_metrics.tsv> <out.tsv>
suppressPackageStartupMessages({library(SummarizedExperiment); library(tidyverse)})
a <- commandArgs(trailingOnly = TRUE)
SE <- a[1]; TECH <- a[2]; OUT_TSV <- a[3]

se  <- readRDS(SE)
cd  <- as.data.frame(colData(se))
arm_lab  <- recode(sub("_day2$", "", cd$condition), NT = "Control")   # MBRT_day2 -> MBRT, etc.
arm_of   <- setNames(arm_lab, cd$sample_id)
slide_of <- setNames(cd$slide_id, cd$sample_id)

# Collapse per-sample-per-celltype counts -> per-sample gene pseudobulk (950 x n_sample).
samp   <- cd$sample_id
pb     <- vapply(unique(samp), function(s) rowSums(assay(se, "counts")[, samp == s, drop = FALSE]),
                 numeric(nrow(se)))
arms   <- factor(arm_of[colnames(pb)], levels = c("Control", "MBRT", "SBRT"))

# Per-sample concordance to arm-mates: mean pairwise Spearman of gene pseudobulk within the arm.
rmat <- cor(pb, method = "spearman")
mate_r <- vapply(seq_len(ncol(pb)), function(i) {
  mates <- which(arms == arms[i]); mates <- setdiff(mates, i)
  mean(rmat[i, mates])
}, numeric(1))

# PCA of the scaled per-sample technical metrics (M02 day-2 only).
tech <- read_tsv(TECH, show_col_types = FALSE) %>% filter(sample_id %in% colnames(pb))
mets <- c("tpc", "med_nFeature", "sparsity", "snr", "specificity_fdr")
pc   <- prcomp(scale(as.matrix(tech[, mets])))
ve   <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
pcd  <- tibble(sample_id = tech$sample_id, PC1 = pc$x[, 1], PC2 = pc$x[, 2])

out <- tibble(sample_id = colnames(pb),
              arm = as.character(arms), slide_id = slide_of[colnames(pb)],
              mean_r_armmates = mate_r) %>%
  left_join(pcd, by = "sample_id") %>%
  mutate(r_outlier = mean_r_armmates < median(mean_r_armmates) - 2 * mad(mean_r_armmates),
         arm = factor(arm, levels = c("Control", "MBRT", "SBRT")),
         pc1_ve = ve[1], pc2_ve = ve[2]) %>%          # PC variance-explained for the figure axis labels
  arrange(mean_r_armmates)
write_tsv(out, OUT_TSV)

cat("per-arm mean within-replicate Spearman:\n")
print(tapply(out$mean_r_armmates, out$arm, mean), digits = 3)
cat(sprintf("r-outliers flagged: %d; wrote %s\n", sum(out$r_outlier), OUT_TSV))
