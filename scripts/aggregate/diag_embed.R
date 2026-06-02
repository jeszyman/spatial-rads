#!/usr/bin/env Rscript
# Diagnostic on the embed checkpoint -- no re-embedding.
#   (1) Leiden resolution calibration  : fix the 68,996-cluster over-partition.
#   (2) slide-only Harmony theta sweep  : can a harder diversity penalty clear the
#                                         LISI(dataset) gate without 2-covariate Harmony?
#   (3) 2-covariate Harmony reproduction: confirm the LAPACK error is real (the proper
#                                         tool for dataset mixing) so the fix is justified.
# Reuses the embed checkpoint at aggregate/embed_checkpoint.rds (PCA + Harmony + SNN
# graph + marker rows). Re-loads once; all three tests share the object.

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(data.table)
  library(igraph)
  library(lisi)
})
set.seed(1)

CKPT    <- "/mnt/data/projects/spatial-rads/aggregate/embed_checkpoint.rds"
LOG     <- "/home/jeszyman/repos/spatial-rads/results/aggregate/diag_embed.log"
OUT_TSV <- "/home/jeszyman/repos/spatial-rads/results/aggregate/diag_embed.tsv"
NPCS    <- 30
N_LISI  <- 20000
PERPLEX <- 30

cat("", file = LOG)
log <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg),
                         file = LOG, append = TRUE)
res <- data.table(test = character(), param = character(), metric = character(),
                  value = character())
add <- function(t, p, m, v) res <<- rbind(res, data.table(test = t, param = p,
                                                           metric = m, value = as.character(v)))

lisi_batch <- function(emb, labels, n_sub = N_LISI, perplexity = PERPLEX) {
  i <- sample(nrow(emb), min(n_sub, nrow(emb)))
  L <- labels[i]
  if (length(unique(L)) < 2) return(NA_real_)
  md_i <- data.frame(.batch = factor(L))
  mean(lisi::compute_lisi(emb[i, , drop = FALSE], md_i, ".batch",
                          perplexity = perplexity)[[".batch"]])
}

log("loading checkpoint")
o        <- readRDS(CKPT)
md       <- o@meta.data
pca_emb  <- Embeddings(o, "pca")[, 1:NPCS]
harm_emb <- Embeddings(o, "harmony")[, 1:NPCS]
log(sprintf("loaded %d cells | datasets %s | slides %d",
            ncol(o), paste(table(md$dataset), collapse = "/"),
            length(unique(md$slide_id))))

# ---- (1) Leiden resolution calibration --------------------------------------
snn_name <- grep("_snn$", names(o@graphs), value = TRUE)[1]
ig <- graph_from_adjacency_matrix(as(o@graphs[[snn_name]], "dgCMatrix"),
                                  mode = "undirected", weighted = TRUE, diag = FALSE)
log(sprintf("igraph: %d vertices, %d edges -- Leiden sweep", gorder(ig), gsize(ig)))
str_w <- strength(ig)                                  # weighted degree (modularity vtx wts)
for (obj in c("modularity", "CPM")) {
  for (r in c(1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 1e-4)) {
    cl <- tryCatch(
      cluster_leiden(ig, objective_function = obj, weights = E(ig)$weight,
                     vertex_weights = if (obj == "modularity") str_w else NULL,
                     resolution = r, n_iterations = 5),
      error = function(e) NULL)
    n <- if (is.null(cl)) NA_integer_ else length(unique(membership(cl)))
    log(sprintf("  leiden %s res=%g -> %s clusters", obj, r, n))
    add("leiden", sprintf("%s_res=%g", obj, r), "n_clusters", n)
  }
}
rm(ig, str_w); invisible(gc())

# ---- (2) slide-only Harmony theta sweep -------------------------------------
add("baseline", "theta=2_slide", "lisi_dataset", round(lisi_batch(harm_emb, md$dataset), 4))
add("baseline", "theta=2_slide", "lisi_slide",   round(lisi_batch(harm_emb, md$slide_id), 4))
add("baseline", "pca_uncorrected", "lisi_dataset", round(lisi_batch(pca_emb, md$dataset), 4))
log("slide-only Harmony theta sweep")
for (th in c(4, 8)) {
  o2 <- tryCatch(
    RunHarmony(o, group.by.vars = "slide_id", theta = th, reduction.use = "pca",
               dims.use = 1:NPCS, reduction.save = "harmony_t", verbose = FALSE),
    error = function(e) { log(sprintf("  harmony theta=%g ERROR: %s", th, conditionMessage(e))); NULL })
  if (!is.null(o2)) {
    he <- Embeddings(o2, "harmony_t")[, 1:NPCS]
    ld <- round(lisi_batch(he, md$dataset), 4)
    ls <- round(lisi_batch(he, md$slide_id), 4)
    log(sprintf("  theta=%g slide-only -> lisi_dataset=%.4f lisi_slide=%.4f", th, ld, ls))
    add("harmony_theta", sprintf("theta=%g_slide", th), "lisi_dataset", ld)
    add("harmony_theta", sprintf("theta=%g_slide", th), "lisi_slide",   ls)
    rm(o2, he); invisible(gc())
  }
}

# ---- (3) 2-covariate Harmony reproduction (LAPACK check) --------------------
log("2-covariate Harmony reproduction")
o3 <- tryCatch(
  RunHarmony(o, group.by.vars = c("dataset", "slide_id"), theta = c(2, 2),
             reduction.use = "pca", dims.use = 1:NPCS,
             reduction.save = "harmony_2cov", verbose = FALSE),
  error = function(e) conditionMessage(e))
if (is.character(o3)) {
  log(sprintf("  2-cov ERROR (expected if LAPACK-disabled): %s", o3))
  add("harmony_2cov", "dataset+slide", "status", paste("ERROR:", o3))
} else {
  he <- Embeddings(o3, "harmony_2cov")[, 1:NPCS]
  ld <- round(lisi_batch(he, md$dataset), 4)
  ls <- round(lisi_batch(he, md$slide_id), 4)
  log(sprintf("  2-cov OK -> lisi_dataset=%.4f lisi_slide=%.4f", ld, ls))
  add("harmony_2cov", "dataset+slide", "lisi_dataset", ld)
  add("harmony_2cov", "dataset+slide", "lisi_slide",   ls)
}

fwrite(res, OUT_TSV, sep = "\t")
log(sprintf("done -- results in %s", OUT_TSV))
print(res)
