#!/usr/bin/env Rscript
# Mutter_01 cross-check: our QC-pass cell set vs Yi's vendor QC (qcFlagsCell == Pass),
# per condition. Args: <yi_reference.tsv> <DATADIR> <out.tsv>
suppressMessages({library(data.table); library(Seurat); library(readr); library(dplyr)})
args <- commandArgs(trailingOnly = TRUE)
YI <- args[1]; DATADIR <- args[2]; OUT <- args[3]

yi <- fread(YI)
stopifnot(all(c("cell_id", "Condition", "qcFlagsCell") %in% names(yi)))

qc_files <- list.files(file.path(DATADIR, "processing", "qc"), pattern = "\\.qc\\.rds$", full.names = TRUE)
rows <- list()
for (f in qc_files) {
  o <- readRDS(f)
  if (!("dataset" %in% colnames(o@meta.data)) || o$dataset[1] != "Mutter_01") next
  cond <- as.character(o$condition[1]); our_cells <- colnames(o)
  yi_cells <- yi[Condition == cond & qcFlagsCell == "Pass", cell_id]
  rows[[length(rows) + 1]] <- data.frame(
    condition = cond, our_kept = length(our_cells), yi_pass = length(yi_cells),
    overlap = length(intersect(our_cells, yi_cells)),
    pct_of_yi = round(100 * length(our_cells) / max(length(yi_cells), 1), 1))
}
out <- bind_rows(rows) %>% arrange(condition)
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
write_tsv(out, OUT)
cat("Mutter_01 QC concordance vs Yi (qcFlagsCell == Pass):\n")
print(as.data.frame(out))
