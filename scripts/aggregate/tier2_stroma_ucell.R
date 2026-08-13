#!/usr/bin/env Rscript
# Tier-2 stroma subtyping, step 2/2 (mirror of tier2_singler.R, UCell engine).
# Cluster-level lineage annotation of the stroma subclusters from step 1 using rank-based
# UCell scoring against the curated stromal marker sets in config/lineage_markers.yaml
# (Endothelial/Fibroblast/Pericyte/SmoothMuscle/Adipocyte). UCell ranks genes within each
# cell, so unlike a raw-count argmax it is invariant to the order-of-magnitude differences in
# marker counts across lineages (Fibroblast ECM genes otherwise dominate and collapse every
# subcluster to Fibroblast). One label per subcluster: mean per-cell UCell score per
# subcluster, argmax over the five lineages, with an unassigned rule (n>=N_MIN cells AND top
# mean score >= MIN_SCORE AND top-minus-second margin >= MARGIN). The margin also guards the
# documented panel-weak Pericyte vs SmoothMuscle separation. Per-cell scores are cached so a
# resolution re-run reuses them.
# Args: <stroma_mtx_dir> <stroma_subclusters.parquet> <lineage_markers.yaml> <outdir>
suppressPackageStartupMessages({
  library(Matrix); library(arrow); library(data.table); library(UCell); library(yaml)
})

N_MIN     <- 100   # min cells to annotate a subcluster (matches tier-1)
MIN_SCORE <- 0.10  # top-lineage mean UCell score floor
MARGIN    <- 0.02  # top-minus-second mean UCell score separation
NCORES    <- 8L
STROMA_LINEAGES <- c("Endothelial", "Fibroblast", "Pericyte", "SmoothMuscle", "Adipocyte")

a <- commandArgs(trailingOnly = TRUE)
mtxdir <- a[1]; sub_pq <- a[2]; markers_yaml <- a[3]; outdir <- a[4]

sub <- as.data.table(read_parquet(sub_pq))
sub[, stroma_subcluster := as.character(stroma_subcluster)]

PERCELL <- file.path(outdir, "stroma_ucell_percell.parquet")
if (file.exists(PERCELL)) {
  cat(sprintf("restoring cached per-cell UCell scores %s ...\n", PERCELL))
  sct <- as.data.table(read_parquet(PERCELL))
} else {
  cts <- as(Matrix::readMM(file.path(mtxdir, "counts.mtx")), "CsparseMatrix")
  feats <- readLines(file.path(mtxdir, "features.tsv"))
  bcs   <- readLines(file.path(mtxdir, "barcodes.tsv"))
  rownames(cts) <- feats; colnames(cts) <- bcs
  cat(sprintf("stroma counts: %d genes x %d cells\n", nrow(cts), ncol(cts)))

  mk <- yaml::read_yaml(markers_yaml)
  sigs <- lapply(STROMA_LINEAGES, function(l) intersect(mk[[l]], rownames(cts)))
  names(sigs) <- STROMA_LINEAGES
  for (l in STROMA_LINEAGES)
    cat(sprintf("  %s: %d/%d markers on panel\n", l, length(sigs[[l]]), length(mk[[l]])))

  # rank-based per-cell scoring; magnitude-invariant (UCell ranks genes within each cell)
  m <- ScoreSignatures_UCell(cts, features = sigs, name = "", ncores = NCORES)
  sct <- as.data.table(m, keep.rownames = "cell")
  write_parquet(sct, PERCELL)
  cat(sprintf("wrote %s\n", PERCELL))
}

# attach subcluster of the chosen resolution (merge on cell, order-independent)
sct <- sub[, .(cell, stroma_subcluster)][sct, on = "cell"]
stopifnot(!anyNA(sct$stroma_subcluster))

# cluster-level: mean UCell score per subcluster per lineage
clu <- sct[, lapply(.SD, mean), by = stroma_subcluster, .SDcols = STROMA_LINEAGES]
clu <- sub[, .(n_cells = .N), by = stroma_subcluster][clu, on = "stroma_subcluster"]

M <- as.matrix(clu[, ..STROMA_LINEAGES])
ord <- t(apply(M, 1, order, decreasing = TRUE))
top_score    <- M[cbind(seq_len(nrow(M)), ord[, 1])]
second_score <- M[cbind(seq_len(nrow(M)), ord[, 2])]
margin <- top_score - second_score
top_lin <- STROMA_LINEAGES[ord[, 1]]
keep <- clu$n_cells >= N_MIN & top_score >= MIN_SCORE & margin >= MARGIN
clu[, `:=`(top_lineage = top_lin,
           top_score = round(top_score, 4),
           second_score = round(second_score, 4),
           margin = round(margin, 4),
           stroma_subtype = fifelse(keep, top_lin, "stroma_unresolved"))]
setorder(clu, -n_cells)
fwrite(clu, file.path(outdir, "stroma_subcluster_ucell.tsv"), sep = "\t")
print(clu)

sub[clu, stroma_subtype := i.stroma_subtype, on = "stroma_subcluster"]
write_parquet(sub[, .(cell, stroma_subcluster, stroma_subtype, dataset, slide_id)],
              file.path(outdir, "stroma_subtypes.parquet"))

summ <- sub[, .(n_cells = .N), by = stroma_subtype][order(-n_cells)]
summ[, frac := round(n_cells / sum(n_cells), 4)]
fwrite(summ, file.path(outdir, "stroma_subtype_summary.tsv"), sep = "\t")
cat("\n--- stroma subtype composition ---\n"); print(summ)
cat("\n--- subtype x dataset ---\n")
print(dcast(sub[, .N, by = .(stroma_subtype, dataset)],
            stroma_subtype ~ dataset, value.var = "N", fill = 0))
peri <- summ[stroma_subtype == "Pericyte", sum(n_cells)]
smc  <- summ[stroma_subtype == "SmoothMuscle", sum(n_cells)]
cat(sprintf("\nQC note (panel-weak separation): Pericyte=%d SmoothMuscle=%d\n",
            ifelse(length(peri), peri, 0L), ifelse(length(smc), smc, 0L)))
cat(sprintf(paste0("\nwrote stroma_subtypes.parquet (%d cells) + ",
                   "stroma_subtype_summary.tsv + stroma_subcluster_ucell.tsv\n"),
            nrow(sub)))
