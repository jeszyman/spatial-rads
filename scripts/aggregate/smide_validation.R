#!/usr/bin/env Rscript
# smiDE per-cell negative-binomial DE (neighbor-expression covariate) as a
# targeted validation of the confirmatory-family pseudobulk DE hits, plus the
# Cdkn1a (p21) result across its cell types. Not a replacement for the
# pseudobulk pipeline -- pseudobulk stays the primary inference; this checks
# whether the same genes/cell types survive per-cell segmentation-error
# correction (RankNorm(otherct_expr) covariate, per overlap_ratio_qc.R's same
# spatial-neighbor logic). M02 day-2 only (the inference cohort).
#
# smiDE's results(comparisons="pairwise") returns ratio/SE/p.value on the
# emmeans contrast scale (log-link nbinom2 -> ratio, not a difference), and
# the contrast label ("levelA / levelB") ordering follows the groupVar
# factor's default level order, NOT the num/den direction of our named
# contrasts (e.g. "SBRT_vs_Ctrl"). comparisons.tsv supplies the actual
# num/den condition-level strings so each pairwise row can be matched to the
# correct named contrast and sign-corrected before converting ratio -> log2FC
# (same delta-method transform smiDE itself uses internally for its volcano
# plots, smiDE:::log2_ci).
# Args: <merged.rds> <full_labels.parquet> <coords_necrosis.parquet>
#       <obs.parquet> <results_master.tsv> <comparisons.tsv> <out_tsv>
suppressPackageStartupMessages({
  library(data.table); library(arrow); library(Seurat); library(Matrix); library(smiDE)
})
a <- commandArgs(trailingOnly = TRUE)
rds_path <- a[1]; labels_path <- a[2]; coords_path <- a[3]; obs_path <- a[4]
master_path <- a[5]; comp_path <- a[6]; out_tsv <- a[7]
RADIUS_MM <- 0.05  # neighbor search radius; matches overlap_ratio_qc.R / coords_necrosis.R
MIN_CELLS <- 20     # minimum cells of a subtype required to attempt a per-cell fit
N_CORES <- 4         # matches the smk rule's threads: 4 (smi_de forks via nCores)

dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)

EMPTY_OUT <- function() {
  data.table(gene = character(), cell_type = character(), contrast = character(),
             pseudobulk_estimate = numeric(), pseudobulk_pvalue = numeric(),
             smide_estimate = numeric(), smide_se = numeric(),
             smide_pvalue = numeric(), survives = logical())
}

# =============================================================================
# SECTION: TARGET LIST -- confirmatory DE hits + Cdkn1a across all its cell types
# =============================================================================

master <- fread(master_path)
confirm_de <- master[tier == "confirmatory" & readout_class == "DE",
                     .(gene = feature, cell_type = unit, contrast)]
cdkn1a <- master[readout_class == "DE" & feature == "Cdkn1a",
                 .(gene = feature, cell_type = unit, contrast)]
targets <- unique(rbind(confirm_de, cdkn1a))

if (nrow(targets) == 0) {
  cat("smide_validation: 0 targets to validate\n")
  fwrite(EMPTY_OUT(), out_tsv, sep = "\t")
  quit(save = "no")
}

# num/den condition-level strings for each named contrast (e.g. SBRT_vs_Ctrl ->
# num=SBRT_day2, den=Control) -- needed to sign-orient smiDE's pairwise output,
# which is oriented by groupVar factor level order, not by contrast name.
comp <- fread(comp_path)
contrast_defs <- unique(comp[name %in% unique(targets$contrast) & !is.na(contrast_num_level),
                             .(contrast = name, num_level = contrast_num_level, den_level = contrast_den_level)])

cat(sprintf("smide_validation: %d gene x cell_type x contrast targets (%d genes, %d cell types)\n",
            nrow(targets), uniqueN(targets$gene), uniqueN(targets$cell_type)))

# =============================================================================
# SECTION: LOAD + ALIGN (M02 day-2 inference cohort only)
# =============================================================================

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

# =============================================================================
# SECTION: PRE-COMPUTE NEIGHBOR EXPRESSION (once, reused across cell types)
# =============================================================================
# split_neighbors_by_colname="sample_id" keeps the neighbor search within a
# tissue section, matching overlap_ratio_qc.R. ref_celltype = the set of
# cell types we'll actually fit models for below (pre_de accepts a vector).
pde <- pre_de(
  adjacencies_only = FALSE,
  metadata = meta,
  ref_celltype = unique(targets$cell_type),
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

# =============================================================================
# SECTION: PER-CELL-TYPE smiDE DE + PAIRWISE CONTRAST EXTRACTION
# =============================================================================
# One smi_de() call per cell type fits all its target genes at once and
# returns pairwise contrasts between every pair of `condition` levels present
# (up to 3, for Control/MBRT_day2/SBRT_day2). Each pairwise row's `contrast`
# label ("levelA / levelB") is parsed and matched against contrast_defs to
# find which named contrast (if any) it corresponds to, sign-correcting so
# smide_estimate is always log2(num_level / den_level) -- the same convention
# count_engine.R uses for pseudobulk_estimate (DESeq2 log2FoldChange with
# den_level as the reference level).
results_list <- list()
for (ct in unique(targets$cell_type)) {
  ct_genes <- targets[cell_type == ct, unique(gene)]
  ct_genes <- ct_genes[ct_genes %in% rownames(counts)]
  if (length(ct_genes) == 0) next

  ct_cells <- meta[cell_subtype == ct, cell]
  if (length(ct_cells) < MIN_CELLS) next

  ct_counts <- counts[ct_genes, ct_cells, drop = FALSE]
  ct_meta <- meta[cell %in% ct_cells]
  setkey(ct_meta, cell)
  ct_meta <- ct_meta[ct_cells]

  tryCatch({
    de <- smi_de(
      assay_matrix = ct_counts,
      metadata = ct_meta,
      formula = ~ RankNorm(otherct_expr) + offset(log(totalcounts)) + condition,
      pre_de_obj = pde,
      groupVar = "condition",
      family = "nbinom2",
      targets = ct_genes,
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

    ct_targets <- targets[cell_type == ct]
    ct_defs <- contrast_defs[contrast %in% unique(ct_targets$contrast)]

    mapped <- vector("list", nrow(ct_defs))
    for (i in seq_len(nrow(ct_defs))) {
      cdef <- ct_defs[i]
      hit <- pw[(levelA == cdef$num_level & levelB == cdef$den_level) |
                (levelA == cdef$den_level & levelB == cdef$num_level)]
      if (nrow(hit) == 0) next
      dir_sign <- ifelse(hit$levelA == cdef$num_level, 1, -1)
      mapped[[i]] <- data.table(
        gene = hit$target,
        contrast = cdef$contrast,
        smide_estimate = dir_sign * log2(hit$ratio),
        smide_se = hit$SE / (hit$ratio * log(2)),  # delta-method log2 SE; sign-invariant
        smide_pvalue = hit$p.value
      )
    }
    mapped <- rbindlist(mapped, fill = TRUE)
    if (nrow(mapped) == 0) next
    mapped[, cell_type := ct]
    # restrict to the gene x contrast pairs actually needed for this cell type
    mapped <- mapped[ct_targets[, .(gene, contrast)], on = .(gene, contrast), nomatch = NULL]
    results_list[[ct]] <- mapped
  }, error = function(e) {
    cat(sprintf("smide_validation: SKIPPED %s (%s)\n", ct, conditionMessage(e)))
  })
}

if (length(results_list) == 0) {
  cat("smide_validation: all cell types skipped\n")
  fwrite(EMPTY_OUT(), out_tsv, sep = "\t")
  quit(save = "no")
}

smide_res <- rbindlist(results_list, fill = TRUE)

# =============================================================================
# SECTION: JOIN PSEUDOBULK EFFECT + SURVIVAL CALL
# =============================================================================

pb <- master[readout_class == "DE", .(gene = feature, cell_type = unit, contrast,
                                       pseudobulk_estimate = effect, pseudobulk_pvalue = pvalue)]
out <- pb[smide_res, on = .(gene, cell_type, contrast)]
out[, survives := !is.na(smide_pvalue) & smide_pvalue < 0.05 &
      sign(smide_estimate) == sign(pseudobulk_estimate)]

fwrite(out, out_tsv, sep = "\t")
n_surv <- sum(out$survives, na.rm = TRUE)
cat(sprintf("smide_validation: %d / %d gene x cell_type x contrast survive smiDE correction\n",
            n_surv, nrow(out)))
