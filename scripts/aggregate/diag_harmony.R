#!/usr/bin/env Rscript
# Harmony-only diagnostic on the embed checkpoint (Leiden question already settled:
# igraph::cluster_leiden over-partitions catastrophically -> switch to reference
# leidenalg in the production script).
#   (1) slide-only Harmony theta sweep  : can a harder diversity penalty (theta 4, 8)
#                                         clear the LISI(dataset) gate (>=1.5) without
#                                         two-covariate Harmony?
#   (2) 2-covariate Harmony reproduction: confirm the LAPACK error is live on this
#                                         install (justifies the version-bump/relink fix).
# Reuses aggregate/embed_checkpoint.rds (PCA + Harmony + meta). One load, both tests.

suppressPackageStartupMessages({
  library(Seurat); library(harmony); library(data.table); library(lisi)
})
set.seed(1)

CKPT    <- "/mnt/data/projects/spatial-rads/aggregate/embed_checkpoint.rds"
LOG     <- "/home/jeszyman/repos/spatial-rads/results/aggregate/diag_harmony.log"
OUT_TSV <- "/home/jeszyman/repos/spatial-rads/results/aggregate/diag_harmony.tsv"
NPCS    <- 30
N_LISI  <- 20000
PERPLEX <- 30

cat("", file = LOG)
log <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg),
                         file = LOG, append = TRUE)
res <- data.table(test = character(), metric = character(), value = character())
add <- function(t, m, v) res <<- rbind(res, data.table(test = t, metric = m,
                                                        value = as.character(v)))
lisi_batch <- function(emb, labels, n_sub = N_LISI, perplexity = PERPLEX) {
  i <- sample(nrow(emb), min(n_sub, nrow(emb)))
  L <- labels[i]; if (length(unique(L)) < 2) return(NA_real_)
  mean(lisi::compute_lisi(emb[i, , drop = FALSE], data.frame(.b = factor(L)), ".b",
                          perplexity = perplexity)[[".b"]])
}

log("loading checkpoint")
o        <- readRDS(CKPT)
md       <- o@meta.data
pca_emb  <- Embeddings(o, "pca")[, 1:NPCS]
harm_emb <- Embeddings(o, "harmony")[, 1:NPCS]
log(sprintf("loaded %d cells | datasets %s", ncol(o), paste(table(md$dataset), collapse = "/")))

# baselines (perfect-mixing ceiling for the 0.29/0.71 imbalance is ~1.70)
add("baseline_pca",        "lisi_dataset", round(lisi_batch(pca_emb,  md$dataset), 4))
add("baseline_theta2_slide", "lisi_dataset", round(lisi_batch(harm_emb, md$dataset), 4))
add("baseline_theta2_slide", "lisi_slide",   round(lisi_batch(harm_emb, md$slide_id), 4))
log(sprintf("baselines: pca lisi_dataset=%s, theta2 lisi_dataset=%s",
            res[test=="baseline_pca" & metric=="lisi_dataset", value],
            res[test=="baseline_theta2_slide" & metric=="lisi_dataset", value]))

# ---- (1) slide-only Harmony at higher theta --------------------------------
for (th in c(4, 8)) {
  log(sprintf("RunHarmony slide-only theta=%g", th))
  o2 <- tryCatch(
    RunHarmony(o, group.by.vars = "slide_id", theta = th, reduction.use = "pca",
               dims.use = 1:NPCS, reduction.save = "harmony_t", verbose = FALSE),
    error = function(e) { log(sprintf("  theta=%g ERROR: %s", th, conditionMessage(e))); NULL })
  if (!is.null(o2)) {
    he <- Embeddings(o2, "harmony_t")[, 1:NPCS]
    ld <- round(lisi_batch(he, md$dataset), 4); ls <- round(lisi_batch(he, md$slide_id), 4)
    log(sprintf("  theta=%g -> lisi_dataset=%.4f lisi_slide=%.4f", th, ld, ls))
    add(sprintf("theta%g_slide", th), "lisi_dataset", ld)
    add(sprintf("theta%g_slide", th), "lisi_slide",   ls)
    rm(o2, he); invisible(gc())
  }
}

# ---- (2) two-covariate Harmony (LAPACK check) ------------------------------
log("RunHarmony 2-covariate dataset+slide_id (LAPACK check)")
o3 <- tryCatch(
  RunHarmony(o, group.by.vars = c("dataset", "slide_id"), theta = c(2, 2),
             reduction.use = "pca", dims.use = 1:NPCS, reduction.save = "harmony_2cov",
             verbose = FALSE),
  error = function(e) conditionMessage(e))
if (is.character(o3)) {
  log(sprintf("  2-cov ERROR (LAPACK still broken): %s", o3))
  add("harmony_2cov", "status", paste("ERROR:", o3))
} else {
  he <- Embeddings(o3, "harmony_2cov")[, 1:NPCS]
  add("harmony_2cov", "lisi_dataset", round(lisi_batch(he, md$dataset), 4))
  add("harmony_2cov", "lisi_slide",   round(lisi_batch(he, md$slide_id), 4))
  log(sprintf("  2-cov OK -> lisi_dataset=%s", res[test=="harmony_2cov" & metric=="lisi_dataset", value]))
}

fwrite(res, OUT_TSV, sep = "\t")
log(sprintf("done -- %s", OUT_TSV)); print(res)
