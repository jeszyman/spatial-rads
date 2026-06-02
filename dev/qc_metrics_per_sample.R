#!/usr/bin/env Rscript
# Per-sample QC metric distributions: per-cell nCount / nFeature / propNegative / Area by sample,
# raw cells before filtering, with the QC gate lines. Caches a downsampled per-cell metric table
# so the figure can be re-rendered without re-reading the raw objects. Exact cells-removed counts
# come from the full-data removal-attribution table, not this display sample.
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
source("/home/jeszyman/repos/science/R/figure_schema.R")

RAW   <- "/mnt/data/projects/spatial-rads/processing/raw"
AGG   <- "results/aggregate"
CACHE <- file.path(AGG, "qc_metrics_per_cell.tsv")
OUT   <- file.path(AGG, "plots", "qc_metrics_per_sample")
N_SUB <- 4000   # cells/sample sampled for display

# QC gates (config/config.yaml qc:)
MIN_COUNTS <- 20; MIN_FEATURES <- 10; MAX_PROPNEG <- 0.5

if (file.exists(CACHE)) {
  d <- fread(CACHE)
} else {
  suppressMessages(library(Seurat))
  files <- sort(list.files(RAW, pattern = "\\.raw\\.rds$", full.names = TRUE))
  stopifnot(length(files) > 0)
  set.seed(1)
  d <- rbindlist(lapply(files, function(f) {
    m   <- as.data.table(readRDS(f)@meta.data)
    idx <- sample(seq_len(nrow(m)), min(N_SUB, nrow(m)))
    dt  <- m[idx, .(nCount_RNA, nFeature_RNA, propNegative, Area)]
    dt[, `:=`(sample_id = as.character(m$sample_id[1]), dataset = as.character(m$dataset[1]))]
    rm(m); gc(); dt
  }))
  fwrite(d, CACHE, sep = "\t")
}

d[, sample_id := factor(sample_id, levels = unique(d[order(dataset, sample_id), sample_id]))]

mk <- function(col, lab, logy = TRUE, hline = NA, show_x = FALSE) {
  p <- ggplot(d, aes(sample_id, .data[[col]], fill = dataset)) +
    geom_boxplot(outlier.size = 0.15, linewidth = 0.25) +
    labs(x = NULL, y = lab) +
    theme_scifig(base_size = 10)
  if (logy)          p <- p + scale_y_log10()
  if (!is.na(hline)) p <- p + geom_hline(yintercept = hline, linetype = "dashed", colour = "grey40")
  if (show_x) p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7))
  else        p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p
}

pn_max <- max(d$propNegative, na.rm = TRUE)
p1 <- mk("nCount_RNA",   "nCount (log10)",   TRUE, MIN_COUNTS)
p2 <- mk("nFeature_RNA", "nFeature (log10)", TRUE, MIN_FEATURES)
p3 <- mk("propNegative", "propNegative",     FALSE) +
  coord_cartesian(ylim = c(0, pn_max * 1.08)) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.04, vjust = 1.5, size = 2.8,
           label = sprintf("gate < %g (data well below)", MAX_PROPNEG))
p4 <- mk("Area", "Area (log10)", TRUE, show_x = TRUE)

key <- paste(strwrap(paste(
  "Per-cell QC metrics by sample; raw cells before filtering,",
  sprintf("downsampled to %d cells/sample for display.", N_SUB),
  "Dashed lines: count and complexity floors (nCount > 20, nFeature > 10).",
  "Cell area gated per-sample at 3 MAD of log10(area).",
  "Mutter_01 + Mutter_02, CosMx."), width = 150), collapse = "\n")

fig <- (p1 / p2 / p3 / p4) +
  plot_layout(guides = "collect", heights = c(1, 1, 0.55, 1)) +
  plot_annotation(caption = key,
                  theme = theme(plot.caption = element_text(size = 9, hjust = 0, colour = "gray30",
                                                            margin = margin(t = 10)),
                                plot.caption.position = "plot"))
fig <- fig & theme(legend.position = "bottom", legend.title = element_blank())

save_plot(fig, OUT, w = 9, h = 9)
cat(sprintf("%d cells across %d samples (<=%d/sample)\n", nrow(d), uniqueN(d$sample_id), N_SUB))
