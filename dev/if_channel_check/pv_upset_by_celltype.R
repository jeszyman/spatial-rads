suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(UpSetR); library(patchwork)
})

OUT <- "/home/jeszyman/repos/spatial-rads/dev/if_channel_check"
dge <- fread("/tmp/dge_pv_by_celltype_fixed.tsv")

# ---- Peak-UP and valley-UP gene sets per cell class ----
# Effect-size threshold (p-values not all usable at small n for some classes)
LFC <- 0.20
PADJ <- 0.05
peak_up <- dge[avg_log2FC >  LFC & p_val_adj < PADJ, .(cell_class, gene)]
valley_up <- dge[avg_log2FC < -LFC & p_val_adj < PADJ, .(cell_class, gene)]

# Fallback: some cell classes have few significant hits at p_adj<0.05 due to small n;
# relax to nominal p<0.01 and keep effect size for exploratory upset
peak_up_relaxed <- dge[avg_log2FC >  LFC & p_val < 0.01, .(cell_class, gene)]
valley_up_relaxed <- dge[avg_log2FC < -LFC & p_val < 0.01, .(cell_class, gene)]

cat("Peak-UP genes per class (p_adj<0.05, |lfc|>0.2):\n")
print(peak_up[, .N, by = cell_class])
cat("\nRelaxed (nominal p<0.01):\n")
print(peak_up_relaxed[, .N, by = cell_class])

# ---- Build list for UpSet ----
classes <- unique(dge$cell_class)
peak_lists <- split(peak_up_relaxed$gene, peak_up_relaxed$cell_class)
valley_lists <- split(valley_up_relaxed$gene, valley_up_relaxed$cell_class)
# Ensure every class has an entry (empty vectors for classes with 0 hits)
for (cls in classes) {
  if (is.null(peak_lists[[cls]]))   peak_lists[[cls]]   <- character(0)
  if (is.null(valley_lists[[cls]])) valley_lists[[cls]] <- character(0)
}

# ---- Save UpSet plots ----
png(file.path(OUT, "upset_peak_up.png"), width = 1400, height = 800, res = 130)
print(upset(fromList(peak_lists), nsets = length(peak_lists),
            order.by = "freq", nintersects = 30,
            sets.bar.color = "#e41a1c",
            main.bar.color = "#333333",
            matrix.color = "#333333",
            mainbar.y.label = "Peak-UP genes in intersection",
            sets.x.label = "Peak-UP genes per class",
            text.scale = c(1.4, 1.2, 1.2, 1.0, 1.3, 1.0)))
dev.off()

png(file.path(OUT, "upset_valley_up.png"), width = 1400, height = 800, res = 130)
print(upset(fromList(valley_lists), nsets = length(valley_lists),
            order.by = "freq", nintersects = 30,
            sets.bar.color = "#377eb8",
            main.bar.color = "#333333",
            matrix.color = "#333333",
            mainbar.y.label = "Valley-UP genes in intersection",
            sets.x.label = "Valley-UP genes per class",
            text.scale = c(1.4, 1.2, 1.2, 1.0, 1.3, 1.0)))
dev.off()

# ---- Write out "cell-type-specific peak state" markers ----
count_in <- function(gene, gene_lists) {
  sum(sapply(gene_lists, function(x) gene %in% x))
}
all_peak_genes <- unique(unlist(peak_lists))
peak_specific <- data.table(gene = all_peak_genes,
                            n_classes = sapply(all_peak_genes, count_in,
                                               gene_lists = peak_lists),
                            classes = sapply(all_peak_genes, function(g)
                              paste(names(peak_lists)[sapply(peak_lists, function(x) g %in% x)], collapse = ",")))
peak_specific <- peak_specific[order(n_classes, gene)]

all_valley_genes <- unique(unlist(valley_lists))
valley_specific <- data.table(gene = all_valley_genes,
                              n_classes = sapply(all_valley_genes, count_in,
                                                 gene_lists = valley_lists),
                              classes = sapply(all_valley_genes, function(g)
                                paste(names(valley_lists)[sapply(valley_lists, function(x) g %in% x)], collapse = ",")))
valley_specific <- valley_specific[order(n_classes, gene)]

fwrite(peak_specific, "/tmp/peak_up_specificity.tsv", sep = "\t")
fwrite(valley_specific, "/tmp/valley_up_specificity.tsv", sep = "\t")

cat("\n=== Peak-UP genes specific to EXACTLY ONE cell class ===\n")
print(peak_specific[n_classes == 1])
cat("\n=== Valley-UP genes specific to EXACTLY ONE cell class ===\n")
print(valley_specific[n_classes == 1])
cat("\n=== Peak-UP genes shared across all 6 classes ===\n")
print(peak_specific[n_classes == 6])
cat("\n=== Valley-UP genes shared across all 6 classes ===\n")
print(valley_specific[n_classes == 6])

cat("\nSaved:\n")
cat("  ", file.path(OUT, "upset_peak_up.png"), "\n")
cat("  ", file.path(OUT, "upset_valley_up.png"), "\n")
cat("  /tmp/peak_up_specificity.tsv\n")
cat("  /tmp/valley_up_specificity.tsv\n")
