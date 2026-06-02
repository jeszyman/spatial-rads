#!/usr/bin/env Rscript
# Repair the broken Yi-T regex. Prior diag used grepl("^T\\.|^Tgd|NKT") which matched almost
# nothing because Yi codes T cells as CD8.T.cell / gdT / NKT / Thymic.* / Spleen.* / Treg / DN* / ISP.
# Step 1: dump the real label vocabulary in this sample. Step 2: classify each label as T / not-T with
# an explicit pattern and SHOW the assignment so it can be eyeballed. Step 3: recompute Yi-T fraction
# and gold-Cd3e+ concordance. Args: <norm.rds> <yi.tsv>
suppressMessages({library(Seurat); library(Matrix); library(readr)})
args <- commandArgs(trailingOnly = TRUE)
obj <- readRDS(args[1]); YI <- args[2]
yi   <- read_tsv(YI, show_col_types = FALSE)
ylab <- setNames(yi$ImmuneAtlas_ImmGen_Main_cell_Types, yi$cell_id)[colnames(obj)]

tab <- sort(table(ylab, useNA = "ifany"), decreasing = TRUE)
# Explicit T-lineage pattern over ImmGen Main nomenclature
Tpat <- "CD4|CD8|gdT|NKT|Treg|Tgd|\\bT\\.cell|T\\.cell|\\.DP|\\.DN[0-9]|preT|ISP|SP$|CD3|Tcon|Th1|Th2|Th17|Tfh|Tmem|Teff|Tnaive"
is_T <- grepl(Tpat, names(tab), ignore.case = TRUE)

cat("=== Yi ImmGen label vocabulary in this sample (label | n | T?) ===\n")
for (i in seq_along(tab)) {
  cat(sprintf("  %-28s %7d   %s\n", names(tab)[i], tab[i], ifelse(is_T[i], "T", "")))
}

# Recompute fractions with the corrected T set
Tlabels <- names(tab)[is_T]
yiT <- ylab %in% Tlabels
cnt <- GetAssayData(obj, layer = "counts")
cd3 <- (cnt["Cd3e", ] > 0) | (cnt["Cd3d", ] > 0)
cat(sprintf("\nCorrected: Yi-T = %.1f%% of cells (%d T labels)\n", 100*mean(yiT, na.rm=TRUE), length(Tlabels)))
cat(sprintf("Cd3e/Cd3d+ cells = %.1f%%\n", 100*mean(cd3)))
cat(sprintf("Gold Cd3+ cells: Yi-T concordance = %.1f%%\n", 100*mean(yiT[cd3], na.rm=TRUE)))
cat(sprintf("Yi-T cells: Cd3+ detection = %.1f%%  (panel sparsity: how often a true T even shows Cd3)\n",
            100*mean(cd3[yiT], na.rm=TRUE)))
