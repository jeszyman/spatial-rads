#!/usr/bin/env Rscript
# Stage B adapter -- Mutter_02: per-slide Seurat RDS -> per-sample common-format Seurat.
# Subset to common panel + shared meta; assign 3 y-bands by deterministic k-means; map
# bands to samples via extract_key (y-rank 1=top..3=bottom). Args: <samplesheet.tsv> <common_genes.tsv> <DATADIR>
suppressMessages({library(Seurat); library(readr); library(dplyr)})
args <- commandArgs(trailingOnly = TRUE)
SS <- args[1]; GENES <- args[2]; DATADIR <- args[3]
set.seed(1)

ss <- read_tsv(SS, show_col_types = FALSE) %>% filter(dataset == "Mutter_02")
common <- readLines(GENES)
shared <- c("nCount_RNA", "nFeature_RNA", "propNegative", "qcFlagsCell",
            "x_slide_mm", "y_slide_mm", "fov", "Area", "Mean.PanCK", "Mean.CD45",
            "Mean.CD298.B2M", "Mean.DAPI")
design <- c("sample_id", "dataset", "treatment", "condition", "timepoint_h", "model")
outdir <- file.path(DATADIR, "processing", "raw"); dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

slides <- ss %>% distinct(slide_id, raw_input_path)
for (j in seq_len(nrow(slides))) {
  sl <- slides[j, ]; this <- ss %>% filter(slide_id == sl$slide_id)
  obj <- UpdateSeuratObject(readRDS(sl$raw_input_path))
  obj$propNegative <- obj$nCount_negprobes / (obj$nCount_RNA + obj$nCount_negprobes)  # generate; vendor QC absent for Mutter_02
  DefaultAssay(obj) <- "RNA"
  for (a in setdiff(names(obj@assays), "RNA")) obj[[a]] <- NULL          # drop negprobes/falsecode
  obj <- subset(obj, features = intersect(common, rownames(obj)))

  set.seed(1)
  km <- kmeans(obj$y_slide_mm, centers = 3, nstart = 10)
  rank_of_cluster <- order(km$centers[, 1], decreasing = TRUE)           # 1 = top (highest y)
  cell_rank <- match(km$cluster, rank_of_cluster)
  stopifnot(length(unique(cell_rank)) == 3)                              # 3 separable bands
  obj$spatial_rank <- cell_rank

  for (r in 1:3) {
    srow <- this %>% filter(extract_key == as.character(r))
    stopifnot(nrow(srow) == 1)
    sub <- subset(obj, cells = colnames(obj)[cell_rank == r])
    sub$sample_id <- srow$sample_id; sub$dataset <- "Mutter_02"; sub$treatment <- srow$treatment
    sub$condition <- srow$condition; sub$timepoint_h <- srow$timepoint_h; sub$model <- srow$model
    keep <- intersect(c(shared, design), colnames(sub@meta.data))
    sub@meta.data <- sub@meta.data[, keep, drop = FALSE]
    saveRDS(sub, file.path(outdir, paste0(srow$sample_id, ".raw.rds")))
    cat(sprintf("%s (%s, slide %s rank %d): %d cells\n",
                srow$sample_id, srow$condition, sl$slide_id, r, ncol(sub)))
  }
}
cat("adapt_mutter02 done\n")
