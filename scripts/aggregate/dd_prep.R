#!/usr/bin/env Rscript
# Differential-detection prep (main spatial-rads env; the muscat container has no
# Seurat/arrow). Reads the merged cohort, restricts to the M02 day-2 arm samples,
# and writes a SingleCellExperiment (.rds) with a raw counts assay + colData
# (sample_id, cell_subtype, condition, slide_id) that dd_muscat.R consumes inside
# the muscat.sif container.
suppressPackageStartupMessages({
  library(Seurat); library(SingleCellExperiment); library(data.table); library(arrow)
})
a <- commandArgs(trailingOnly = TRUE)
merged_path <- if (length(a) >= 1) a[1] else "/mnt/data/projects/spatial-rads/aggregate/full/merged.rds"
samples_tsv <- if (length(a) >= 2) a[2] else "results/data_model/samples.tsv"
out_sce     <- if (length(a) >= 3) a[3] else "/mnt/data/projects/spatial-rads/aggregate/full/dd_sce.rds"

o <- readRDS(merged_path)

# M02 day-2 arm design: map sample_id -> condition/slide_id from the sample sheet.
ss <- fread(samples_tsv)
scol <- names(ss)[grepl("sample", names(ss), ignore.case = TRUE)][1]
day2 <- ss[grepl("Mutter_02|mutter02|dat0002", get(names(ss)[grepl("dataset", names(ss))][1])) &
           timepoint_h == 48]
keep_samples <- intersect(unique(o$sample_id), day2[[scol]])
o <- subset(o, cells = colnames(o)[o$sample_id %in% keep_samples & o$cell_subtype != "unassigned"])

# condition + slide_id per cell from the sample sheet
smeta <- day2[, .(sample_id = get(scol),
                  condition = fifelse(treatment == "NT", "Control", paste0(treatment, "_day2")),
                  slide_id)]
cd <- data.table(cell = colnames(o), sample_id = o$sample_id, cell_subtype = o$cell_subtype)
cd <- smeta[cd, on = "sample_id"]

sce <- SingleCellExperiment(
  assays = list(counts = LayerData(o, layer = "counts")),
  colData = DataFrame(sample_id = cd$sample_id, cell_subtype = cd$cell_subtype,
                      condition = cd$condition, slide_id = cd$slide_id))
saveRDS(sce, out_sce)
cat(sprintf("dd_sce: %d cells, %d samples, %d subtypes, conditions %s\n",
            ncol(sce), length(unique(sce$sample_id)), length(unique(sce$cell_subtype)),
            paste(unique(sce$condition), collapse = ",")))
