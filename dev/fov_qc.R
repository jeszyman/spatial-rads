#!/usr/bin/env Rscript
# Cells-per-FOV QC: per-sample distribution of cells per CosMx FOV, flagging sparse FOVs
# (tissue edges / segmentation gaps). Caches per-(sample, fov) counts so the plot can be
# re-rendered without re-reading the raw objects.
suppressMessages({library(data.table); library(ggplot2)})
source("/home/jeszyman/repos/science/R/figure_schema.R")

RAW      <- "/mnt/data/projects/spatial-rads/processing/raw"
AGG      <- "results/aggregate"
CACHE    <- file.path(AGG, "qc_cells_per_fov.tsv")
OUT      <- file.path(AGG, "plots", "qc_cells_per_fov")
FLAG_MIN <- 20   # FOVs with fewer cells flagged sparse (edge / segmentation)

if (file.exists(CACHE)) {
  fc <- fread(CACHE)
} else {
  suppressMessages(library(Seurat))
  files <- sort(list.files(RAW, pattern = "\\.raw\\.rds$", full.names = TRUE))
  stopifnot(length(files) > 0)
  fc <- rbindlist(lapply(files, function(f) {
    m <- as.data.table(readRDS(f)@meta.data)
    dt <- m[, .(n_cells = .N), by = fov]
    dt[, `:=`(sample_id = as.character(m$sample_id[1]), dataset = as.character(m$dataset[1]))]
    rm(m); gc(); dt
  }))
  fwrite(fc, CACHE, sep = "\t")
}

fc[, sample_id := factor(sample_id, levels = unique(fc[order(dataset, sample_id), sample_id]))]
n_sparse <- fc[n_cells < FLAG_MIN, .N]; n_fov <- nrow(fc)

key <- paste(strwrap(paste(
  "Cells per field of view (FOV), one box per sample; raw cells before QC filtering.",
  sprintf("Dashed line: sparse-FOV flag at %d cells (%d of %d FOVs below).", FLAG_MIN, n_sparse, n_fov),
  "Mutter_01 + Mutter_02, CosMx."), width = 130), collapse = "\n")

p <- ggplot(fc, aes(sample_id, n_cells, fill = dataset)) +
  geom_boxplot(outlier.size = 0.4, linewidth = 0.3) +
  geom_hline(yintercept = FLAG_MIN, linetype = "dashed", colour = "grey40") +
  facet_grid(~dataset, scales = "free_x", space = "free_x") +
  scale_y_log10(breaks = c(10, 30, 100, 300, 1000, 3000)) +
  labs(x = NULL, y = "Cells per FOV (log10)", caption = key) +
  theme_scifig(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
        legend.position = "none",
        plot.caption = element_text(size = 9, hjust = 0, colour = "gray30", margin = margin(t = 12)),
        plot.caption.position = "plot")

save_plot(p, OUT, w = 9, h = 5)
cat(sprintf("%d FOVs across %d samples; %d sparse (<%d cells)\n",
            n_fov, uniqueN(fc$sample_id), n_sparse, FLAG_MIN))
