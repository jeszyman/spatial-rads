#!/usr/bin/env Rscript
# Genome-wide per-cell smiDE differential expression with NB mixed model.
# Fits a per-cell negative-binomial GLMM (nebula backend) for every gene in each
# cell subtype, with a neighbor-expression covariate absorbing segmentation error
# and a sample-level random effect for pseudoreplication correction:
#   count ~ RankNorm(otherct_expr) + offset(log(totalcounts)) + condition + (1|sample_id)
# Output matches the engine schema (comparison/contrast/unit/feature_type/feature_id/
# estimate/se/df/stat/p) for direct incorporation into results_master.tsv as
# readout_class="smiDE". M02 day-2 cohort only (the inference cohort).
# Args: <merged.rds> <full_labels.parquet> <coords_necrosis.parquet>
#       <obs.parquet> <comparisons.tsv> <out_tsv>
suppressPackageStartupMessages({
  library(data.table); library(arrow); library(Seurat); library(Matrix); library(smiDE)
})
a <- commandArgs(trailingOnly = TRUE)
rds_path <- a[1]; labels_path <- a[2]; coords_path <- a[3]; obs_path <- a[4]
comp_path <- a[5]; out_tsv <- a[6]
RADIUS_MM <- 0.05
MIN_CELLS <- 50
N_CORES <- 4

dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)

EMPTY_OUT <- function() {
  data.table(comparison=character(), contrast=character(), unit=character(),
             feature_type=character(), feature_id=character(),
             estimate=numeric(), se=numeric(), df=numeric(),
             stat=numeric(), p=numeric())
}

comp <- fread(comp_path)
cohort_comp <- comp[name == "mutter02_day2" & !is.na(contrast_num_level)]
if (nrow(cohort_comp) == 0) {
  cat("smide_de: no mutter02_day2 comparisons found\n")
  fwrite(EMPTY_OUT(), out_tsv, sep = "\t")
  quit(save = "no")
}
contrast_defs <- unique(cohort_comp[, .(contrast = name, num_level = contrast_num_level,
                                        den_level = contrast_den_level)])
conditions <- unique(c(cohort_comp$contrast_num_level, cohort_comp$contrast_den_level))

seu <- readRDS(rds_path)
counts <- GetAssayData(seu, assay = "RNA", layer = "counts")
lab <- as.data.table(read_parquet(labels_path))[, .(cell, cell_subtype)]
co  <- as.data.table(read_parquet(coords_path))[, .(cell, sample_id, x_slide_mm, y_slide_mm)]
ob  <- as.data.table(read_parquet(obs_path))[, .(cell, condition, slide_id, dataset, timepoint_h)]
meta <- lab[co, on = "cell"]
meta <- ob[meta, on = "cell"]
meta <- meta[!is.na(cell_subtype)]
meta <- meta[dataset == "Mutter_02" & timepoint_h == 48]

cells_keep <- intersect(colnames(counts), meta$cell)
counts <- counts[, cells_keep]
meta <- meta[cell %in% cells_keep]
setkey(meta, cell)
meta <- meta[colnames(counts)]
meta[, totalcounts := Matrix::colSums(counts)]

genes <- rownames(counts)
subtypes <- sort(unique(meta$cell_subtype))

cat(sprintf("smide_de: %d cells, %d genes, %d subtypes, %d samples\n",
            ncol(counts), length(genes), length(subtypes), uniqueN(meta$sample_id)))

pde <- pre_de(
  adjacencies_only = FALSE,
  metadata = meta,
  ref_celltype = subtypes,
  cell_type_metadata_colname = "cell_subtype",
  cellid_colname = "cell",
  sdimx_colname = "x_slide_mm",
  sdimy_colname = "y_slide_mm",
  split_neighbors_by_colname = "sample_id",
  mm_radius = RADIUS_MM,
  counts = counts,
  normalized_data = NULL,
  aggregation = "sum"
)

results_list <- list()
for (ct in subtypes) {
  ct_cells <- meta[cell_subtype == ct, cell]
  if (length(ct_cells) < MIN_CELLS) {
    cat(sprintf("smide_de: SKIP %s (%d cells < %d)\n", ct, length(ct_cells), MIN_CELLS))
    next
  }
  if (uniqueN(meta[cell %in% ct_cells, sample_id]) < 2) {
    cat(sprintf("smide_de: SKIP %s (< 2 samples)\n", ct))
    next
  }

  ct_counts <- counts[, ct_cells, drop = FALSE]
  ct_meta <- meta[cell %in% ct_cells]
  setkey(ct_meta, cell)
  ct_meta <- ct_meta[ct_cells]

  tryCatch({
    de <- smi_de(
      assay_matrix = ct_counts,
      metadata = ct_meta,
      formula = ~ RankNorm(otherct_expr) + offset(log(totalcounts)) + condition + (1|sample_id),
      pre_de_obj = pde,
      groupVar = "condition",
      family = "nbinom2",
      targets = genes,
      cellid_colname = "cell",
      nCores = N_CORES,
      neighbor_expr_overlap_agg = "sum",
      neighbor_expr_totalcount_normalize = FALSE,
      neighbor_expr_cell_type_metadata_colname = "cell_subtype"
    )
    pw <- as.data.table(results(de, comparisons = "pairwise", variable = "condition")$pairwise)
    pw[, target := as.character(target)]
    pw[, c("levelA", "levelB") := tstrsplit(as.character(contrast), " / ", fixed = TRUE)]
    pw[, levelA := trimws(levelA)]
    pw[, levelB := trimws(levelB)]

    mapped <- vector("list", nrow(contrast_defs))
    for (i in seq_len(nrow(contrast_defs))) {
      cdef <- contrast_defs[i]
      hit <- pw[(levelA == cdef$num_level & levelB == cdef$den_level) |
                (levelA == cdef$den_level & levelB == cdef$num_level)]
      if (nrow(hit) == 0) next
      dir_sign <- ifelse(hit$levelA == cdef$num_level, 1, -1)
      mapped[[i]] <- data.table(
        comparison = "mutter02_day2",
        contrast = cdef$contrast,
        unit = ct,
        feature_type = "gene",
        feature_id = hit$target,
        estimate = dir_sign * log2(hit$ratio),
        se = hit$SE / (hit$ratio * log(2)),
        df = NA_real_,
        stat = NA_real_,
        p = hit$p.value
      )
    }
    mapped <- rbindlist(mapped, fill = TRUE)
    if (nrow(mapped) > 0) results_list[[ct]] <- mapped
    cat(sprintf("smide_de: %s — %d genes, %d contrast rows\n",
                ct, length(genes), nrow(mapped)))
  }, error = function(e) {
    cat(sprintf("smide_de: ERROR %s — %s\n", ct, conditionMessage(e)))
  })
}

if (length(results_list) == 0) {
  cat("smide_de: all cell types failed or skipped\n")
  fwrite(EMPTY_OUT(), out_tsv, sep = "\t")
} else {
  out <- rbindlist(results_list, fill = TRUE)
  fwrite(out, out_tsv, sep = "\t")
  cat(sprintf("smide_de: %d rows | %d genes x %d cell types x %d contrasts\n",
              nrow(out), uniqueN(out$feature_id), uniqueN(out$unit), uniqueN(out$contrast)))
}
