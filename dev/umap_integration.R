#!/usr/bin/env Rscript
# Integration-QC UMAP: the visual companion to lisi_integration.R. Two panels share
# one Harmony embedding -- (A) colored by dataset shows whether the two cohorts mix;
# (B) colored by cell type shows whether shared biology co-localizes across them.
# Reads the slim TSV from umap_integration_extract.R (150k-cell subsample).
suppressMessages({library(data.table); library(ggplot2); library(scattermore)
                  library(patchwork); library(ggsci); library(stringr)})
source("/home/jeszyman/repos/science/R/figure_schema.R")

TSV <- "results/aggregate/umap_integration.tsv"
OUT <- "results/aggregate/plots/umap_integration"
ILISI_POST <- 1.28; GATE <- 1.5   # from results/aggregate/celltype_embed_metrics.tsv
stopifnot(file.exists(TSV))

d <- fread(TSV)
n_cell <- nrow(d)
ds_lvls <- sort(unique(d$dataset))
pct <- d[, .(p = round(100 * .N / n_cell)), by = dataset]
pct_str <- paste(sprintf("%s %d%%", pct$dataset, pct$p), collapse = ", ")

# Grey the heterogeneous catch-all buckets (a/b are de-novo ImmGen, not real types)
# + rare types; color only the commonest RESOLVED types on top, so panel B reads as
# populations not one blob. Palette is a curated max-distinctness set (Trubetskoy)
# of fully saturated hues -- deliberately NO grey/brown/olive, which the category10
# palette includes and which dissolve into the grey backdrop (that was the blob).
catchall <- c("a", "b")
GREY_LAB <- "a / unresolved"
res <- d[!cell_type %in% catchall, .N, by = cell_type][order(-N)][seq_len(min(8, .N)), cell_type]
d[, ct := factor(fifelse(cell_type %in% res, cell_type, GREY_LAB), levels = c(res, GREY_LAB))]

distinct8 <- c("#e6194b", "#4363d8", "#3cb44b", "#f58231",
               "#911eb4", "#42d4f4", "#f032e6", "#000075")
ds_cols <- setNames(c("#1b6ca8", "#e07b39")[seq_along(ds_lvls)], ds_lvls)
ct_cols <- setNames(c(distinct8[seq_along(res)], "grey85"), c(res, GREY_LAB))

# Panel A is split by dataset (small multiples on shared coords): each cohort is
# drawn over a grey backdrop of ALL cells, so the reader compares which regions
# each fills. backdrop has no `dataset` col -> it repeats in every facet.
backdrop <- d[, .(UMAP1, UMAP2)]

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
  "by cell type (resolved; grey = a / unresolved)")

key <- paste(
  "Batch-mixing QC for the aggregate Harmony embedding (companion to the LISI figure), one shared UMAP shown two ways.",
  "(A) Split by dataset: each cohort (color) over a grey backdrop of all cells, so the regions each fills are directly comparable -- a cohort filling the same shape as the backdrop is well mixed.",
  sprintf("Residual under-mixing is quantified separately (post-Harmony iLISI %.2f, gate %.1f).", ILISI_POST, GATE),
  "(B) Colored by cell type: the 8 commonest resolved types are colored over a grey backdrop of the heterogeneous `a`/de-novo bucket and rare types, so resolved populations stand out rather than a single blob. Shared identities co-locate across cohorts.",
  "Labels are per-sample atlas calls (Yi-ImmGen for Mutter_01, anchor transfer for Mutter_02).",
  sprintf("%s-cell random subsample; %s.", format(n_cell, big.mark = ","), pct_str))

combo <- (g1 | g2) +
  plot_layout(widths = c(2, 1.2)) +
  plot_annotation(tag_levels = "A", caption = str_wrap(key, 200)) &
  theme(plot.tag = element_text(face = "bold", size = 12),
        plot.caption = element_text(size = 9, hjust = 0, colour = "gray30", margin = margin(t = 12)),
        plot.caption.position = "plot")

save_plot(combo, OUT, w = 13.5, h = 6.2)
cat(sprintf("UMAP: %d cells | datasets %s | %d resolved types colored + grey catch-all\n",
            n_cell, pct_str, length(res)))
