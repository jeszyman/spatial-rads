#!/usr/bin/env Rscript
# scVI twin of umap_integration.R: same two-panel integration-QC UMAP, but on the
# scVI latent embedding + the merged-scale resolved atlas (cluster-then-annotate),
# instead of the Harmony embedding + per-sample transfer labels.
#   (A) colored by dataset -> do the two cohorts mix on the scVI manifold;
#   (B) colored by resolved cell_subtype -> does shared biology co-localize.
# Reads the 200k scVI-UMAP subsample emitted by scripts/aggregate/full_umap.py.
suppressMessages({library(arrow); library(data.table); library(ggplot2); library(scattermore)
                  library(patchwork); library(ggsci); library(stringr)})
source("/home/jeszyman/repos/science/R/figure_schema.R")

PARQ <- "/mnt/data/projects/spatial-rads/aggregate/full/umap_coords.parquet"
OUT  <- "results/aggregate/plots/umap_integration_scvi"
ILISI <- 0.67; GATE <- 0.5   # scVI scaled iLISI from results/aggregate/full_gates.json
stopifnot(file.exists(PARQ))

d <- as.data.table(read_parquet(PARQ))
setnames(d, c("umap_int_1", "umap_int_2", "cell_subtype"), c("UMAP1", "UMAP2", "cell_type"))
n_cell <- nrow(d)
ds_lvls <- sort(unique(d$dataset))
pct <- d[, .(p = round(100 * .N / n_cell)), by = dataset]
pct_str <- paste(sprintf("%s %d%%", pct$dataset, pct$p), collapse = ", ")

# Grey the unassigned catch-all + rare types; color only the commonest RESOLVED
# subtypes on top so panel B reads as populations, not one blob. Curated max-
# distinctness hues (Trubetskoy) -- NO grey/brown/olive (those dissolve into backdrop).
catchall <- c("unassigned")
GREY_LAB <- "unassigned / rare"
res <- d[!cell_type %in% catchall, .N, by = cell_type][order(-N)][seq_len(min(8, .N)), cell_type]
d[, ct := factor(fifelse(cell_type %in% res, cell_type, GREY_LAB), levels = c(res, GREY_LAB))]

distinct8 <- c("#e6194b", "#4363d8", "#3cb44b", "#f58231",
               "#911eb4", "#42d4f4", "#f032e6", "#000075")
ds_cols <- setNames(c("#1b6ca8", "#e07b39")[seq_along(ds_lvls)], ds_lvls)
ct_cols <- setNames(c(distinct8[seq_along(res)], "grey85"), c(res, GREY_LAB))

backdrop <- d[, .(UMAP1, UMAP2)]   # no `dataset` col -> repeats in every facet of A

umap_base <- function(g, sub) {
  g + coord_equal() +
    labs(x = "UMAP 1", y = "UMAP 2", subtitle = sub) +
    theme_scifig(base_size = 11) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.border = element_rect(colour = "grey80", fill = NA),
          legend.position = "bottom", legend.title = element_blank(),
          legend.key.size = unit(0.3, "cm"), legend.text = element_text(size = 8),
          plot.subtitle = element_text(face = "bold", size = 11, hjust = 0),
          axis.title.x = element_text(margin = margin(t = 4)),
          axis.title.y = element_text(margin = margin(r = 4)))
}

g1 <- umap_base(
  ggplot() +
    geom_scattermore(data = backdrop, aes(UMAP1, UMAP2), color = "grey86", pointsize = 1.4) +
    geom_scattermore(data = d, aes(UMAP1, UMAP2, color = dataset), pointsize = 1.9, alpha = 0.6) +
    facet_wrap(~ dataset, nrow = 1) +
    scale_color_manual(values = ds_cols, guide = "none"),
  "by dataset (each cohort vs. all cells, grey)")

g2 <- umap_base(
  ggplot() +
    geom_scattermore(data = d[ct == GREY_LAB], aes(UMAP1, UMAP2),
                     color = "grey85", pointsize = 1.3) +
    geom_scattermore(data = d[ct != GREY_LAB], aes(UMAP1, UMAP2, color = ct),
                     pointsize = 2.7, alpha = 0.8) +
    scale_color_manual(values = ct_cols, breaks = res) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 4)),
  "by cell type (resolved; grey = unassigned / rare)")

key <- paste(
  "Batch-mixing QC for the aggregate scVI embedding (companion to the LISI figure), one shared UMAP shown two ways.",
  "(A) Split by dataset: each cohort (color) over a grey backdrop of all cells, so the regions each fills are directly comparable -- a cohort filling the same shape as the backdrop is well mixed.",
  sprintf("Residual under-mixing is quantified separately (scVI scaled iLISI %.2f, gate >= %.1f).", ILISI, GATE),
  "(B) Colored by cell type: the 8 commonest resolved subtypes are colored over a grey backdrop of unassigned + rare types, so resolved populations stand out rather than a single blob. Shared identities co-locate across cohorts.",
  "Labels are the merged-scale resolved atlas (scVI cluster-then-annotate), which supersedes the per-sample transfer labels.",
  sprintf("%s-cell random subsample; %s.", format(n_cell, big.mark = ","), pct_str))

combo <- (g1 | g2) +
  plot_layout(widths = c(2, 1.2)) +
  plot_annotation(tag_levels = "A", caption = str_wrap(key, 200)) &
  theme(plot.tag = element_text(face = "bold", size = 12),
        plot.caption = element_text(size = 9, hjust = 0, colour = "gray30", margin = margin(t = 12)),
        plot.caption.position = "plot")

save_plot(combo, OUT, w = 13.5, h = 6.2)
cat(sprintf("scVI UMAP: %d cells | datasets %s | %d resolved subtypes colored + grey catch-all\n",
            n_cell, pct_str, length(res)))
