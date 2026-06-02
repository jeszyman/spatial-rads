#!/usr/bin/env Rscript
# aggregate.smk Stage 1a -- cell-type integration embedding (Harmony).
# ScaleData(all 950 panel genes) -> PCA(30) -> Harmony(slide_id; see below) ->
# FindNeighbors -> Leiden clustering (igraph) -> UMAP. Cell-type labels stay the
# per-sample assignments (Yi-ImmGen for M01, anchor TransferData for M02); the
# data-driven clusters are QC only, scored for marker concordance.
# Per plan-aggregate.md Stage 1a.
#
# Standard tools (spatial-rads conda env):
#   * Clustering: Leiden via igraph::cluster_leiden (Traag 2019), the field-
#     standard community-detection algorithm. (Earlier this script used Seurat's
#     Louvain because leidenalg was absent; the standard tool is now installed
#     and used. igraph's pure-C Leiden also vastly outpaces Seurat's single-
#     threaded Louvain on the 3.27M-cell graph -- minutes vs hours.)
#   * Batch mixing: LISI via lisi::compute_lisi (Korsunsky 2019), the published
#     perplexity-weighted inverse-Simpson diversity. (Earlier a hand-rolled
#     hard-kNN approximation stood in for the absent package.)
#
# One environment-gated deviation remains:
#   * Harmony groups on slide_id ALONE (one covariate), not the plan's
#     c("dataset","slide_id"). The env's r-harmony 2.0.2 binary throws
#     `inv(): use of LAPACK must be enabled` for >=2 covariates (reproducible;
#     single-covariate inv() works). slide_id nests dataset (each of the 7
#     slides belongs to exactly one dataset), so correcting slide_id also
#     corrects the coarser dataset axis -- the LISI(dataset) gate (now the real
#     package) validates dataset mixing. Restore c("dataset","slide_id") if the
#     gate underperforms and harmony is rebuilt against a LAPACK-enabled
#     Armadillo.
#
# Checkpoint: the expensive embedding (load -> scale -> PCA -> Harmony ->
# FindNeighbors, ~1.5 h on 3.27M cells) is written to `embed_checkpoint.rds`
# once complete, WITHOUT the dense 950 x 3.27M scale.data (only the ~50 marker
# rows needed downstream are stashed in o@misc$sd_markers). Re-invoking the
# script with the checkpoint present skips straight to clustering, so Leiden
# resolution / metrics / UMAP can be re-run in minutes. Delete the checkpoint to
# force a full re-embed.
#
# Integration gate (plan-aggregate.md): post-Harmony LISI(dataset) must reach
# `lisi_min` (production 1.8; interim floor 1.5) or the rule fails with a
# remediation hint. The gate runs BEFORE the (50-min, viz-only) UMAP so a badly
# mixed embedding fails fast. Scalar metrics are cat() to the log before the gate
# so the diagnostic survives even when snakemake removes the rule's outputs.
#
# Silhouette and LISI run on subsamples (full 3.27M-cell pairwise distance / kNN
# is intractable); the heavy embedding itself uses all cells.
#
# Args: <merged.rds> <out_celltype.rds> <out_qc.tsv> <out_metrics.tsv>
#       <out_umap.png> <lisi_min>
suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(data.table)
  library(igraph)
  library(lisi)
  library(cluster)
  library(ggplot2)
  library(patchwork)
  library(future)
})
set.seed(1)

# FindNeighbors routes the 3.27M-cell harmony embedding through future.apply,
# whose global-export size check (default 500 MiB) trips on the >2 GiB neighbor-
# search closure. Run sequentially (no workers -- consistent with threads:1 BLAS
# hygiene; the embedding is single-process anyway) and lift the cap.
plan("sequential")
options(future.globals.maxSize = 32 * 1024^3)

args        <- commandArgs(trailingOnly = TRUE)
merged      <- args[1]
out_rds     <- args[2]
out_qc      <- args[3]
out_metrics <- args[4]
out_umap    <- args[5]
lisi_min    <- as.numeric(args[6])

NPCS    <- 30
RES     <- 0.5
N_SIL   <- 5000      # silhouette distance subsample
N_LISI  <- 20000     # LISI subsample
PERPLEX <- 30        # lisi::compute_lisi perplexity (Korsunsky 2019 default)
N_PLOT  <- 150000    # UMAP scatter subsample
N_ITER  <- 10        # Leiden refinement iterations
for (d in unique(c(dirname(out_rds), dirname(out_qc), dirname(out_umap))))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
CKPT <- file.path(dirname(out_rds), "embed_checkpoint.rds")

# conda run fully buffers R stdout until exit, so the rule log stays empty for the
# whole multi-hour run. Mirror phase milestones to a sidecar that flushes each write
# (cat with file= opens/closes per call) so progress is observable live.
PROG <- file.path(dirname(out_metrics), "embed_progress.log")
if (!file.exists(CKPT)) cat("", file = PROG)
prog <- function(msg)
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg),
      file = PROG, append = TRUE)

# ---- marker sets for per-cluster QC concordance ------------------------------
marker_sets <- list(
  tumor_epithelial = c("Krt8","Krt18","Krt19","Epcam","Cdh1"),
  T_cell      = c("Cd3e","Cd3d","Cd3g","Cd8a","Cd4","Trbc2"),
  B_cell      = c("Cd19","Ms4a1","Cd79a","Cd79b","Ighm"),
  Plasma      = c("Jchain","Sdc1","Xbp1","Prdm1"),
  Macrophage  = c("Adgre1","Csf1r","Lyz2","Mrc1","C1qa","C1qb"),
  Mono_DC     = c("Itgax","Flt3","Cd209a","Ccr7","Itgam"),
  NK          = c("Ncr1","Klrb1c","Gzmb","Klrk1"),
  Neutrophil  = c("S100a8","S100a9","Retnlg"),
  Endothelial = c("Pecam1","Cdh5","Kdr","Vwf"),
  Fibroblast  = c("Col1a1","Col1a2","Col3a1","Pdgfra","Pdgfrb"),
  SmoothMuscle = c("Acta2","Myh11","Tagln"))

# ---- batch-mixing metrics ----------------------------------------------------
# silhouette: lower = better mixed. LISI (real package): higher = better mixed.
sil_batch <- function(emb, labels, n_sub = N_SIL) {
  i <- sample(nrow(emb), min(n_sub, nrow(emb)))
  L <- as.integer(factor(labels[i]))
  if (length(unique(L)) < 2) return(NA_real_)
  mean(silhouette(L, dist(emb[i, , drop = FALSE]))[, "sil_width"])
}
lisi_batch <- function(emb, labels, n_sub = N_LISI, perplexity = PERPLEX) {
  i <- sample(nrow(emb), min(n_sub, nrow(emb)))
  L <- labels[i]
  if (length(unique(L)) < 2) return(NA_real_)
  md_i <- data.frame(.batch = factor(L))
  mean(lisi::compute_lisi(emb[i, , drop = FALSE], md_i, ".batch",
                          perplexity = perplexity)[[".batch"]])
}

# ===== EMBED (full) or RESUME (from checkpoint) ===============================
if (file.exists(CKPT)) {
  prog("resuming from checkpoint (skip scale/PCA/Harmony/neighbors)")
  o          <- readRDS(CKPT)
  md         <- o@meta.data
  pca_emb    <- Embeddings(o, "pca")[, 1:NPCS]
  harm_emb   <- Embeddings(o, "harmony")[, 1:NPCS]
  sd_markers <- o@misc$sd_markers
  cat(sprintf("resumed %d cells x %d genes from %s\n", ncol(o), nrow(o), CKPT))
} else {
  # ---- load + ensure single normalized layer --------------------------------
  o  <- readRDS(merged)
  o  <- tryCatch(JoinLayers(o), error = function(e) o)   # no-op if already joined
  DefaultAssay(o) <- "RNA"
  dat <- tryCatch(LayerData(o, layer = "data"), error = function(e) NULL)
  if (is.null(dat) || length(dat@x) == 0) {
    cat("data layer empty -> NormalizeData(LogNormalize, 1e4)\n")
    o <- NormalizeData(o, normalization.method = "LogNormalize", scale.factor = 1e4,
                       verbose = FALSE)
  }
  md <- o@meta.data
  stopifnot(all(c("dataset", "slide_id", "cell_type") %in% colnames(md)))
  cat(sprintf("loaded %d cells x %d genes | datasets: %s\n",
              ncol(o), nrow(o), paste(table(md$dataset), collapse = "/")))
  prog(sprintf("loaded %d cells x %d genes", ncol(o), nrow(o)))

  # ---- scale -> PCA -> Harmony ----------------------------------------------
  o <- ScaleData(o, features = rownames(o), verbose = FALSE); prog("ScaleData done")
  o <- RunPCA(o, features = rownames(o), npcs = NPCS, verbose = FALSE); prog("PCA done")
  pca_emb <- Embeddings(o, "pca")[, 1:NPCS]
  o <- RunHarmony(o, group.by.vars = "slide_id", theta = 2,
                  reduction.use = "pca", dims.use = 1:NPCS, verbose = FALSE)
  harm_emb <- Embeddings(o, "harmony")[, 1:NPCS]; prog("Harmony done")

  # ---- neighbor graph (Leiden substrate) ------------------------------------
  o <- FindNeighbors(o, reduction = "harmony", dims = 1:NPCS, verbose = FALSE)
  prog("FindNeighbors done")

  # ---- stash marker scale.data rows, drop the dense 24 GB scale.data --------
  # Only the ~50 marker genes are needed downstream (per-cluster concordance);
  # dropping the full matrix frees ~24 GB before the igraph build + clustering.
  sd_full    <- LayerData(o, layer = "scale.data")
  mk_present <- intersect(unique(unlist(marker_sets)), rownames(sd_full))
  sd_markers <- sd_full[mk_present, , drop = FALSE]   # ~50 x 3.27M, ~1.3 GB
  rm(sd_full)
  o[["RNA"]]$scale.data <- NULL
  o@misc$sd_markers <- sd_markers
  invisible(gc())

  # ---- checkpoint (no full scale.data; uncompressed for fast I/O) ------------
  saveRDS(o, CKPT, compress = FALSE)
  prog("checkpoint written (embedding + graph + marker rows)")
}

# ===== CLUSTER -> METRICS -> GATE -> CONCORDANCE -> UMAP =======================
# ---- Leiden clustering via igraph (standard tool) ----------------------------
snn_name <- grep("_snn$", names(o@graphs), value = TRUE)[1]
ig <- graph_from_adjacency_matrix(as(o@graphs[[snn_name]], "dgCMatrix"),
                                  mode = "undirected", weighted = TRUE, diag = FALSE)
prog(sprintf("igraph built: %d vertices, %d edges", gorder(ig), gsize(ig)))
leiden_res <- cluster_leiden(ig, objective_function = "modularity",
                             weights = E(ig)$weight, resolution = RES,
                             n_iterations = N_ITER)
md$cluster <- as.integer(membership(leiden_res))   # vertex order == Cells(o)
nclust <- uniqueN(md$cluster)
rm(ig); invisible(gc())
cat(sprintf("clusters: %d (Leiden res %.1f, igraph)\n", nclust, RES))
prog(sprintf("Leiden done: %d clusters", nclust))

# ---- scalar integration metrics (long format) -------------------------------
metrics <- data.table(metric = character(), value = numeric())
add <- function(k, v) metrics <<- rbind(metrics, data.table(metric = k, value = v))
add("n_cells", ncol(o)); add("n_genes", nrow(o)); add("n_pcs", NPCS)
add("leiden_resolution", RES); add("n_clusters", nclust)
add("sil_dataset_pre",  sil_batch(pca_emb,  md$dataset))
add("sil_dataset_post", sil_batch(harm_emb, md$dataset))
add("sil_slide_pre",    sil_batch(pca_emb,  md$slide_id))
add("sil_slide_post",   sil_batch(harm_emb, md$slide_id))
add("lisi_dataset_pre",  lisi_batch(pca_emb,  md$dataset))
lisi_dataset_post <- lisi_batch(harm_emb, md$dataset)
add("lisi_dataset_post", lisi_dataset_post)
add("lisi_slide_pre",   lisi_batch(pca_emb,  md$slide_id))
add("lisi_slide_post",  lisi_batch(harm_emb, md$slide_id))
add("lisi_min_gate", lisi_min)
metrics[, value := round(value, 4)]
fwrite(metrics, out_metrics, sep = "\t")
cat("=== integration metrics ===\n"); print(metrics)
prog(sprintf("metrics written; LISI(dataset)_post=%.3f", lisi_dataset_post))

# ---- integration gate (BEFORE the 50-min UMAP; light outputs already written) ----
if (is.na(lisi_dataset_post) || lisi_dataset_post < lisi_min) {
  stop(sprintf(paste0("LISI(dataset) gate FAILED: post-Harmony %.3f < %.2f. ",
                      "Dataset axis under-mixed; cross-dataset concordance would be ",
                      "invalid. Remediation: re-run RunHarmony with theta=c(4,2) or ",
                      "revisit batch covariates (checkpoint preserved at %s)."),
               lisi_dataset_post, lisi_min, CKPT))
}

# ---- per-cluster marker concordance -----------------------------------------
present  <- lapply(marker_sets, intersect, rownames(sd_markers))
present  <- present[lengths(present) > 0]
# per-cell set score = mean scaled expression of the set's panel-present markers
set_score <- sapply(present, function(g) colMeans(sd_markers[g, , drop = FALSE]))  # cells x sets
clu <- sort(unique(md$cluster))
qc <- rbindlist(lapply(clu, function(cl) {
  ci   <- which(md$cluster == cl)
  tb   <- sort(table(md$cell_type[ci]), decreasing = TRUE)
  csc  <- colMeans(set_score[ci, , drop = FALSE])      # mean set score in cluster
  ord  <- order(csc, decreasing = TRUE)
  data.table(
    cluster                 = cl,
    n_cells                 = length(ci),
    modal_cell_type         = names(tb)[1],
    modal_frac              = round(as.numeric(tb[1]) / length(ci), 3),
    marker_lineage          = names(present)[ord[1]],
    marker_concordance_score = round(csc[ord[1]], 3),
    second_marker_lineage   = names(present)[ord[2]],
    second_marker_score     = round(csc[ord[2]], 3))
}))
setorder(qc, -n_cells)
fwrite(qc, out_qc, sep = "\t")

# ---- UMAP (viz only) + QC plot ----------------------------------------------
o <- RunUMAP(o, reduction = "harmony", dims = 1:NPCS,
             n.neighbors = 30, min.dist = 0.3, verbose = FALSE); prog("UMAP done")
um <- Embeddings(o, "umap")
pi <- sample(nrow(um), min(N_PLOT, nrow(um)))
top15 <- names(sort(table(md$cell_type), decreasing = TRUE))[1:15]
pdf <- data.table(
  UMAP1 = um[pi, 1], UMAP2 = um[pi, 2],
  dataset = md$dataset[pi],
  cell_type = ifelse(md$cell_type[pi] %in% top15, md$cell_type[pi], "Other"),
  cluster = factor(md$cluster[pi]))
base <- function() list(geom_point(size = 0.1, alpha = 0.4), coord_fixed(),
                        theme_void(base_size = 9),
                        theme(legend.key.size = unit(0.3, "cm")))
g1 <- ggplot(pdf, aes(UMAP1, UMAP2, color = dataset)) + base() +
  guides(color = guide_legend(override.aes = list(size = 2))) + labs(title = "dataset")
g2 <- ggplot(pdf, aes(UMAP1, UMAP2, color = cell_type)) + base() +
  guides(color = guide_legend(override.aes = list(size = 2), ncol = 1)) + labs(title = "cell type (top15)")
g3 <- ggplot(pdf, aes(UMAP1, UMAP2, color = cluster)) + base() +
  guides(color = "none") + labs(title = sprintf("Leiden clusters (n=%d)", nclust))
ggsave(out_umap, (g1 | g2 | g3) + plot_layout(widths = c(1, 1.3, 1)),
       width = 16, height = 5.5, dpi = 150)

# ---- persist (drop the stashed marker rows; downstream uses data/counts) -----
# Set both `cluster` and `seurat_clusters` (+ Idents) so the object is a drop-in
# for any consumer; active typing keys on cell_type, the clusters stay QC-only.
o$cluster         <- md$cluster
o$seurat_clusters <- factor(md$cluster)
Idents(o)         <- o$seurat_clusters
o@misc$sd_markers <- NULL; invisible(gc())
prog("gate passed; saving merged_celltype.rds")
saveRDS(o, out_rds)
prog("DONE")
cat(sprintf("embed_celltype: %d cells, %d clusters, LISI(dataset)_post=%.3f (gate %.2f) PASS\n",
            ncol(o), nclust, lisi_dataset_post, lisi_min))
