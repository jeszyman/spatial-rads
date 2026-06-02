#!/usr/bin/env Rscript
# Probe-vs-negative-control DIAGNOSTIC -- report-only, does NOT drop genes.
#
# NanoString good practice (Danaher, CosMx Analysis Scratch Space) advises AGAINST removing
# low-expression genes -- that loses rare-cell-type markers -- and instead models background
# per cell; InSituType accounts for background internally via its per-cell neg term. So we keep
# the full panel and only FLAG candidate failed probes (Ozirmak Lermi 2025, Nat Commun
# 10.1038/s41467-025-63414-1): genes sitting at/below the negative-control background in BOTH
# datasets (signal-to-background <= 1), i.e. behaving like a negative probe. Review, do not auto-drop.
# Args: <samplesheet.tsv> <common_genes.tsv> <out_report.tsv>
suppressMessages({library(readr); library(dplyr); library(arrow); library(data.table); library(Seurat); library(Matrix)})
args <- commandArgs(trailingOnly = TRUE)
SS <- args[1]; PANEL <- args[2]; OUT <- args[3]
ss <- read_tsv(SS, show_col_types = FALSE)
panel <- readLines(PANEL)

m02_paths <- unique(ss$raw_input_path[ss$name == "Mutter_02"])
o1 <- UpdateSeuratObject(readRDS(m02_paths[1]))
NEG_N <- nrow(GetAssayData(o1, assay = "negprobes", layer = "counts")); rm(o1); gc()

## Mutter_01: per-gene mean (counts parquet); per-negprobe background (metadata negprobe summary)
m01 <- ss %>% filter(name == "Mutter_01")
cdt <- as.data.table(read_parquet(unique(m01$counts_path)))
gcols <- intersect(setdiff(names(cdt), c("Slide", "fov", "cell_id")), panel)
m01_mean <- unlist(cdt[, lapply(.SD, mean), .SDcols = gcols]); rm(cdt); gc()
mdt <- as.data.table(read_parquet(unique(m01$metadata_path)))
m01_bg <- mean(mdt$nCount_NegativeProbes, na.rm = TRUE) / NEG_N; rm(mdt); gc()

## Mutter_02: per-gene mean (RNA assay); per-negprobe background (negprobes assay)
gm <- list(); bgs <- numeric(0)
for (p in m02_paths) {
  o <- UpdateSeuratObject(readRDS(p))
  rna <- GetAssayData(o, assay = "RNA", layer = "counts")
  gm[[p]] <- Matrix::rowMeans(rna[intersect(panel, rownames(rna)), , drop = FALSE])
  bgs <- c(bgs, mean(Matrix::rowMeans(GetAssayData(o, assay = "negprobes", layer = "counts")))); rm(o); gc()
}
g02 <- Reduce(intersect, lapply(gm, names))
m02_mean <- rowMeans(vapply(gm, function(x) x[g02], numeric(length(g02)))); m02_bg <- mean(bgs)

genes <- intersect(names(m01_mean), g02)
rep <- data.frame(gene = genes,
                  m01_mean = m01_mean[genes], m01_bg = m01_bg, m01_s2b = m01_mean[genes] / m01_bg,
                  m02_mean = m02_mean[genes], m02_bg = m02_bg, m02_s2b = m02_mean[genes] / m02_bg)
rep$flagged_failed <- rep$m01_s2b <= 1 & rep$m02_s2b <= 1   # at background in BOTH = candidate failed probe
rep <- rep[order(-rep$flagged_failed, rep$m02_s2b), ]
cat(sprintf("probe-QC diagnostic: %d genes | NEG_N=%d | M01 bg=%.4f M02 bg=%.4f | %d flagged at-background in BOTH (NOT dropped)\n",
            nrow(rep), NEG_N, m01_bg, m02_bg, sum(rep$flagged_failed)))
if (sum(rep$flagged_failed)) cat("flagged failed-probe candidates:", paste(rep$gene[rep$flagged_failed], collapse = ", "), "\n")
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
write_tsv(rep, OUT)
