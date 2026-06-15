#!/usr/bin/env Rscript
# Characterize negprobe (10) + falsecode (184) controls for BOTH datasets, from raw RDS.
# Per-cell metrics + per-(dataset,slide,fov) summary + per-dataset distribution quantiles.
# Report-only calibration step; sets no thresholds. Args: <m01_dir> <m02_dir> <out_tsv> <out_png>
suppressMessages({library(Seurat); library(data.table); library(Matrix); library(ggplot2)})
a <- commandArgs(trailingOnly = TRUE)
M01_DIR <- a[1]; M02_DIR <- a[2]; OUT_TSV <- a[3]; OUT_PNG <- a[4]
NNEG <- 10; NFALSE <- 184

per_fov <- list(); quant <- list()
load_set <- function(dir, pat, ds) {
  files <- sort(list.files(dir, pattern = pat, full.names = TRUE))
  out <- vector("list", length(files))
  for (i in seq_along(files)) {
    o <- UpdateSeuratObject(readRDS(files[i]))     # M02 objects are older-format; harmless no-op for M01
    nc_rna   <- as.numeric(o$nCount_RNA)
    nc_neg   <- if ("nCount_negprobes" %in% colnames(o@meta.data)) as.numeric(o$nCount_negprobes)
               else as.numeric(Matrix::colSums(LayerData(o, assay = "negprobes", layer = "counts")))
    nc_false <- if ("nCount_falsecode" %in% colnames(o@meta.data)) as.numeric(o$nCount_falsecode)
               else as.numeric(Matrix::colSums(LayerData(o, assay = "falsecode", layer = "counts")))
    slide <- unique(as.character(o$Run_Tissue_name))[1]
    dt <- data.table(dataset = ds, slide = slide, fov = as.integer(o$fov),
                     rna = nc_rna, neg = nc_neg, false = nc_false,
                     neg_per_code = nc_neg / NNEG, false_per_code = nc_false / NFALSE,
                     false_frac = nc_false / (nc_rna + nc_neg + nc_false + 1e-9))
    out[[i]] <- dt; rm(o); invisible(gc())
  }
  rbindlist(out)
}
m01 <- load_set(M01_DIR, "Mutter_01_CosMmR\\.RDS$", "Mutter_01")
m02 <- load_set(M02_DIR, "Mutter_02_CosMmR\\.RDS$", "Mutter_02")
all <- rbind(m01, m02)

# per-dataset distribution quantiles (report)
qs <- c(0.5, 0.75, 0.9, 0.95, 0.99)
quant <- all[, as.list(c(
  n_cells = .N,
  setNames(quantile(false_frac, qs), paste0("false_frac_q", qs)),
  setNames(quantile(false_per_code, qs), paste0("false_per_code_q", qs)),
  setNames(quantile(neg_per_code, qs), paste0("neg_per_code_q", qs)))), by = dataset]
dir.create(dirname(OUT_TSV), recursive = TRUE, showWarnings = FALSE)
fwrite(quant, OUT_TSV, sep = "\t")
cat("per-dataset control quantiles:\n"); print(quant)

# cohort-vs-cohort: per-FOV mean false_frac, paired isn't possible (different FOVs), so
# overlay the two per-FOV distributions on a common axis (density), M01 vs M02.
fov <- all[, .(mean_false_frac = mean(false_frac), n = .N), by = .(dataset, slide, fov)]
p <- ggplot(fov, aes(x = mean_false_frac, colour = dataset)) +
  geom_density() +
  labs(title = "Per-FOV falsecode fraction by cohort",
       x = "mean falsecode fraction per FOV", y = "density", colour = NULL) +
  theme_bw()
dir.create(dirname(OUT_PNG), recursive = TRUE, showWarnings = FALSE)
ggsave(OUT_PNG, p, width = 7, height = 4.5, dpi = 150)
cat(sprintf("wrote %s and %s\n", OUT_TSV, OUT_PNG))
