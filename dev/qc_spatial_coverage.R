#!/usr/bin/env Rscript
# Spatial tissue-coverage QC: each CosMx FOV drawn at its real slide position
# (x_slide_mm, y_slide_mm) and filled by median transcript counts per cell, one
# small-multiple per sample. Reveals coverage problems a per-FOV *distribution*
# (fov_qc.R) cannot localize -- tissue edges, detachment / necrosis, and
# segmentation dead zones show up as their actual spatial footprint. Caches the
# per-(sample, fov) footprint + median so the plot re-renders without re-reading
# the raw objects.
suppressMessages({library(data.table); library(ggplot2); library(viridisLite)
                  library(scales); library(stringr)})
source("/home/jeszyman/repos/science/R/figure_schema.R")

RAW   <- "/mnt/data/projects/spatial-rads/processing/raw"
AGG   <- "results/aggregate"
CACHE <- file.path(AGG, "qc_spatial_coverage.tsv")
OUT   <- file.path(AGG, "plots", "qc_spatial_coverage")

if (file.exists(CACHE)) {
  fc <- fread(CACHE)
} else {
  suppressMessages(library(Seurat))
  files <- sort(list.files(RAW, pattern = "\\.raw\\.rds$", full.names = TRUE))
  stopifnot(length(files) > 0)
  fc <- rbindlist(lapply(files, function(f) {
    m  <- as.data.table(readRDS(f)@meta.data)
    dt <- m[, .(n_cells     = .N,
                med_counts  = as.double(median(nCount_RNA)),
                xmin = min(x_slide_mm), xmax = max(x_slide_mm),
                ymin = min(y_slide_mm), ymax = max(y_slide_mm)), by = fov]
    dt[, `:=`(sample_id = as.character(m$sample_id[1]),
              dataset   = as.character(m$dataset[1]))]
    rm(m); gc(); dt
  }))
  fwrite(fc, CACHE, sep = "\t")
}

fc[, sample_id := factor(sample_id, levels = unique(fc[order(dataset, sample_id), sample_id]))]
# Slide coords are per-slide (samples share a slide are offset apart); recenter each
# sample on its own tissue so all panels sit at a common origin. Fixed (shared) scales +
# coord_equal then keep true mm aspect AND honest relative tissue sizes across panels.
fc[, `:=`(cx = (min(xmin) + max(xmax)) / 2, cy = (min(ymin) + max(ymax)) / 2), by = sample_id]
fc[, `:=`(xmin = xmin - cx, xmax = xmax - cx, ymin = ymin - cy, ymax = ymax - cy)]
n_samp <- uniqueN(fc$sample_id); n_fov <- nrow(fc)
ds_n   <- fc[, .(s = uniqueN(sample_id)), by = dataset][order(dataset)]
ds_str <- paste(sprintf("%s n=%d", ds_n$dataset, ds_n$s), collapse = ", ")
cap_hi <- as.numeric(quantile(fc$med_counts, 0.99))   # squish top 1% so a few rich FOVs don't wash out the scale

m01 <- round(median(fc[dataset == "Mutter_01", med_counts]))
m02 <- round(median(fc[dataset == "Mutter_02", med_counts]))
key <- paste(
  "Transcript capture is spatially uniform within each tissue: low-count regions track the tissue edge, with no interior dead zones, so coverage and segmentation read clean across all samples.",
  sprintf("Mutter_01 + Mutter_02 (%s), CosMx 950-gene common panel; %d FOVs, raw cells before QC.", ds_str, n_fov),
  "Each panel is one sample; each rectangle is one field of view (FOV) drawn at its true slide position and filled by median counts per cell.",
  sprintf("Capture floors are comparable across cohorts (per-FOV median %d vs %d counts); Mutter_02 adds high-capture hotspots (yellow cores) that widen its range. Fill saturates at the 99th percentile so those few rich FOVs do not wash out the scale.", m01, m02))

p <- ggplot(fc, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = med_counts)) +
  geom_rect() +
  facet_wrap(~ sample_id, ncol = 6) +
  coord_equal() +
  scale_fill_viridis_c(name = "Median counts / cell",
                       limits = c(0, cap_hi), oob = scales::squish,
                       breaks = pretty(c(0, cap_hi), 4)) +
  labs(x = "x (mm)", y = "y (mm)", caption = str_wrap(key, 150)) +
  theme_scifig(base_size = 11) +
  theme(legend.position = "bottom",
        legend.key.height = unit(0.3, "cm"), legend.key.width = unit(1.1, "cm"),
        legend.title = element_text(size = 9), legend.text = element_text(size = 8),
        strip.text = element_text(size = 8),
        axis.text = element_text(size = 6),
        axis.title.x = element_text(margin = margin(t = 6)),
        axis.title.y = element_text(margin = margin(r = 6)),
        plot.caption = element_text(size = 9, hjust = 0, colour = "gray30", margin = margin(t = 12)),
        plot.caption.position = "plot")

save_plot(p, OUT, w = 12, h = 9)
cat(sprintf("spatial coverage: %d FOVs across %d samples | median-count range %.0f-%.0f (cap %.0f) | %s\n",
            n_fov, n_samp, min(fc$med_counts), max(fc$med_counts), cap_hi, ds_str))
