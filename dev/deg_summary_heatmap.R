#!/usr/bin/env Rscript
# DEG-count heatmap: significant pseudobulk genes per cell type x contrast (Mutter_02 day2).
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
source("/home/jeszyman/repos/science/R/figure_schema.R")

AGG <- "results/aggregate"
OUT <- file.path(AGG, "plots", "deg_summary_heatmap_m02day2")

s   <- fread(file.path(AGG, "deg_summary_m02day2.tsv"))
skp <- fread(file.path(AGG, "degs_pseudobulk_skipped.tsv"))
contr_levels <- c("MBRT_vs_Ctrl", "SBRT_vs_Ctrl", "MBRT_vs_SBRT")

tested <- s[, .(cell_type, contrast, n_sig = n_padj_05, n_strong = n_padj_05_lfc_1,
                status = "tested")]
skipped <- unique(skp[, .(cell_type, contrast)])
skipped[, `:=`(n_sig = NA_integer_, n_strong = NA_integer_, status = "skipped")]

d <- rbind(tested, skipped)
d[, contrast := factor(contrast, levels = contr_levels)]

# Row order: most-responsive cell types toward the top; skip-only types pinned to the bottom.
ord       <- d[status == "tested", .(tot = sum(n_sig, na.rm = TRUE)), by = cell_type][order(tot)]
skip_only <- setdiff(unique(d$cell_type), ord$cell_type)
d[, cell_type := factor(cell_type, levels = c(skip_only, ord$cell_type))]

d[, lab := fifelse(status == "skipped", "skip",
            fifelse(n_strong > 0, sprintf("%d (%d)", n_sig, n_strong), as.character(n_sig)))]
thr <- max(d$n_sig, na.rm = TRUE) / 2
d[, txt_col := fifelse(!is.na(n_sig) & n_sig > thr, "white", "black")]

key <- paste(strwrap(paste(
  "Differentially expressed genes per cell type across three contrasts.",
  "Tile shade and number: genes at padj < 0.05; (n): genes also passing |log2FC| > 1;",
  "grey: cell type below the abundance floor (not tested).",
  "Mutter_02 day2, pseudobulk DESeq2."), width = 98), collapse = "\n")

p <- ggplot(d, aes(contrast, cell_type, fill = n_sig)) +
  geom_tile(colour = "grey90") +
  geom_text(aes(label = lab, colour = txt_col), size = 2.4) +
  scale_fill_gradient(low = "white", high = "firebrick", na.value = "grey85",
                      name = "Genes\n(padj < 0.05)") +
  scale_colour_identity() +
  labs(x = NULL, y = NULL, caption = key) +
  theme_scifig(base_size = 9) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        plot.caption = element_text(size = 8, hjust = 0, colour = "gray30",
                                    margin = margin(t = 12)),
        plot.caption.position = "plot")

save_plot(p, OUT, w = 6.5, h = 9.5)
cat(sprintf("%d tested + %d skipped tiles, %d cell types\n",
            nrow(tested), nrow(skipped), uniqueN(d$cell_type)))
