#!/usr/bin/env Rscript
# Per-gene, per-cell_subtype segmentation-contamination QC via smiDE's
# overlap_ratio_metric: for each gene x cell_subtype pair, compares average
# expression in cells of that subtype ("self") against average expression
# attributable to neighboring cells of OTHER subtypes within a fixed spatial
# radius ("neighbor_othertype"). A ratio >= 1 means neighbor-contamination
# signal meets or exceeds the gene's own signal in that cell type -- i.e. the
# gene's expression in that subtype is not trustworthy at face value (CosMx
# segmentation bleeds transcripts across cell boundaries). Output is a flat QC
# annotation table; assemble_results.R left-joins `ratio` onto results_master
# DE rows (feature = gene, unit = cell_subtype) as `contamination_ratio`.
# Args: <merged.rds> <full_labels.parquet> <coords_necrosis.parquet> <out_overlap_ratio_qc.tsv>
suppressPackageStartupMessages({
  library(data.table); library(arrow); library(Seurat); library(Matrix); library(smiDE)
})
a <- commandArgs(trailingOnly = TRUE)
rds_path <- a[1]; labels_path <- a[2]; coords_path <- a[3]; out_tsv <- a[4]
RADIUS_MM <- 0.05  # neighbor search radius; matches the cell-cell contact scale used elsewhere (coords_necrosis.R kNN)

dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# SECTION: LOAD + ALIGN
# =============================================================================

seu <- readRDS(rds_path)
counts <- GetAssayData(seu, assay = "RNA", layer = "counts")
lab <- as.data.table(read_parquet(labels_path))[, .(cell, cell_subtype)]
co  <- as.data.table(read_parquet(coords_path))[, .(cell, sample_id, x_slide_mm, y_slide_mm)]
meta <- lab[co, on = "cell"]
meta <- meta[!is.na(cell_subtype)]

# Align counts columns and metadata rows to the same cell set/order -- smiDE
# indexes both by position, not by id.
cells_keep <- intersect(colnames(counts), meta$cell)
counts <- counts[, cells_keep]
meta <- meta[cell %in% cells_keep]
setkey(meta, cell)
meta <- meta[colnames(counts)]

# =============================================================================
# SECTION: OVERLAP RATIO METRIC
# =============================================================================

# split_neighbors_by_colname="sample_id" keeps the neighbor search within a
# tissue section (physical slides carry several non-overlapping regions; see
# coords_necrosis.R), so cross-sample coordinate coincidences never count as
# spatial neighbors.
orm <- overlap_ratio_metric(
  assay_matrix = counts,
  metadata     = meta,
  cluster_col  = "cell_subtype",
  cellid_col   = "cell",
  sdimx_col    = "x_slide_mm",
  sdimy_col    = "y_slide_mm",
  radius       = RADIUS_MM,
  split_neighbors_by_colname = "sample_id"
)
# smiDE's native names (verified against installed smiDE 0.0.2.5, see the Step 1b
# prototype in the task brief): target/avg_cluster/avg_neighbor_othercluster.
# `all_data` is a constant provenance flag ("all_cells"), dropped below -- one row
# per gene x cell_subtype already carries all the information it would add.
setnames(orm, old = c("target", "avg_cluster", "avg_neighbor_othercluster"),
         new = c("gene", "avg_self", "avg_neighbor_othertype"), skip_absent = TRUE)
orm <- orm[, .(gene, cell_subtype, avg_self, avg_neighbor_othertype, ratio)]
fwrite(orm, out_tsv, sep = "\t")

n_contaminated <- sum(orm$ratio >= 1, na.rm = TRUE)
n_total <- nrow(orm)
cat(sprintf("overlap_ratio_qc: %d / %d gene x cell_subtype pairs contamination-dominated (ratio >= 1)\n",
            n_contaminated, n_total))
