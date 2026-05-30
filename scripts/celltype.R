#!/usr/bin/env Rscript
# Per-sample cell typing (label transfer against the M01-derived reference):
#   Mutter_01: carry Yi's ImmGen label as cell_type (yi_celltype = original); relabel
#              confident-epithelial cells in Yi's de-novo a/b buckets -> tumor_epithelial.
#   Mutter_02: Seurat anchor label transfer (FindTransferAnchors CCA + TransferData) against the
#              reference of pooled M01 cells -> cell_type (+ celltype_prob = prediction.score.max).
#              Batch-robust; replaces InSituType, whose count-likelihood collapsed M02's tumor
#              compartment into a Lymphatic.endothelial sink under the M01->M02 platform shift.
# Args: <in.norm.rds> <reference.rds> <yi_reference.tsv> <out.typed.rds> <out.summary.tsv>
suppressMessages({library(Seurat); library(Matrix); library(readr)})
args <- commandArgs(trailingOnly = TRUE)
IN <- args[1]; REF <- args[2]; YI <- args[3]; OUT_RDS <- args[4]; OUT_SUM <- args[5]
EPI <- c("Epcam", "Krt8", "Krt18", "Krt19")

obj <- readRDS(IN)
ds  <- as.character(obj$dataset[1]); sid <- as.character(obj$sample_id[1])

if (ds == "Mutter_01") {
  yi  <- read_tsv(YI, show_col_types = FALSE)
  lab <- setNames(yi$ImmuneAtlas_ImmGen_Main_cell_Types, yi$cell_id)
  ct  <- unname(lab[colnames(obj)])
  obj$yi_celltype <- ct
  epi <- intersect(EPI, rownames(obj))
  cnt <- GetAssayData(obj, layer = "counts")
  epi_detect <- Matrix::colSums(cnt[epi, , drop = FALSE] > 0)
  relabel <- (ct %in% c("a", "b") | is.na(ct)) & epi_detect >= 2
  ct[relabel] <- "tumor_epithelial"
  obj$cell_type <- ct
  cat(sprintf("%s (Mutter_01): %d cells; %d relabeled -> tumor_epithelial\n", sid, ncol(obj), sum(relabel)))
} else {
  refobj  <- readRDS(REF)
  feat    <- intersect(rownames(obj), rownames(refobj))
  anchors <- FindTransferAnchors(reference = refobj, query = obj, features = feat,
                                 dims = 1:30, reduction = "cca", verbose = FALSE)
  pred <- TransferData(anchorset = anchors, refdata = as.character(refobj$cell_type),
                       weight.reduction = "cca", dims = 1:30, verbose = FALSE)
  obj$cell_type     <- pred$predicted.id
  obj$celltype_prob <- pred$prediction.score.max
  cat(sprintf("%s (Mutter_02): Seurat label transfer typed %d cells\n", sid, ncol(obj)))
}

summ <- as.data.frame(table(obj$cell_type), responseName = "n")
names(summ)[1] <- "cell_type"; summ$sample_id <- sid; summ$dataset <- ds
dir.create(dirname(OUT_RDS), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(OUT_SUM), recursive = TRUE, showWarnings = FALSE)
saveRDS(obj, OUT_RDS); write_tsv(summ, OUT_SUM)
