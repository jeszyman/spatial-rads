#!/usr/bin/env Rscript
# aggregate.smk M01<->M02 day-2 effect-size concordance on the unified labels.
# Qualitative concordance gate (NOT validation): do the two independently-collected
# cohorts agree on the DIRECTION/magnitude of day-2 treatment effects, now that both
# are jointly typed? Supersedes dev/mbrt_vs_sbrt/08_set2_validation.R (which re-typed
# M02 separately). M02 effect = pseudobulk DESeq2 log2FC (degs_pseudobulk_m02day2.tsv,
# n=4, the formal track). M01 effect = descriptive mean-ratio, n=1: log2 of the ratio
# of mean linear-normalized expression (Seurat AverageExpression on the data layer,
# which de-logs LogNormalize) between arm and Control, per cell_type x gene. Spearman
# rho of M01-vs-M02 log2FC per cell_type x contrast (MBRT/SBRT vs Control). No p-values
# on the M01 side.
# Args: <merged.rds> <full_labels.parquet> <degs_m02.tsv> <out_tsv> <plot_scatter>
suppressPackageStartupMessages({
  library(Seurat); library(arrow); library(data.table); library(ggplot2)
})
a <- commandArgs(trailingOnly = TRUE)
merged_path <- a[1]; labels_path <- a[2]; degs_m02_path <- a[3]
out_tsv <- a[4]; plot_sc <- a[5]
EPS <- 1e-3
CONDS <- c("Control", "MBRT_day2", "SBRT_day2")
CONTRASTS <- c(MBRT_vs_Ctrl = "MBRT_day2", SBRT_vs_Ctrl = "SBRT_day2")

dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(plot_sc), recursive = TRUE, showWarnings = FALSE)

# --- M02 side: pseudobulk DESeq2 log2FC (the two vs-Control contrasts) -----------
m02 <- fread(degs_m02_path)[contrast %in% names(CONTRASTS),
                            .(cell_type, gene, contrast, log2FC_m02 = log2FC)]

# --- M01 side: descriptive mean-ratio log2FC (n=1) -------------------------------
o   <- readRDS(merged_path)
lab <- as.data.table(read_parquet(labels_path))[, .(cell, cell_subtype)]
lv  <- setNames(lab$cell_subtype, lab$cell)
o$cell_subtype <- unname(lv[colnames(o)])
m01 <- subset(o, cells = colnames(o)[o$dataset == "Mutter_01" &
                                     !is.na(o$cell_subtype) &
                                     o$cell_subtype != "unassigned" &
                                     o$condition %in% CONDS])
rm(o); invisible(gc())
m01$grp <- paste(m01$cell_subtype, m01$condition, sep = "@@")  # '@@' survives Seurat ident renaming
Idents(m01) <- "grp"
ae <- AverageExpression(m01, assays = "RNA", layer = "data")$RNA   # genes x grp (linear)
ae <- as.matrix(ae); rm(m01); invisible(gc())

# Seurat maps every '_' -> '-' in ident strings (so MBRT_day2 -> MBRT-day2); the '@@'
# separator and the underscore-free cell_subtype names survive intact.
gl       <- strsplit(colnames(ae), "@@", fixed = TRUE)
col_ct   <- vapply(gl, `[`, "", 1L)
col_cond <- vapply(gl, `[`, "", 2L)
dash     <- function(x) gsub("_", "-", x)
cts <- sort(unique(col_ct))
m01_rows <- rbindlist(lapply(cts, function(ct) {
  cc <- which(col_ct == ct & col_cond == "Control")
  if (length(cc) != 1L) return(NULL)
  rbindlist(lapply(names(CONTRASTS), function(cn) {
    tc <- which(col_ct == ct & col_cond == dash(CONTRASTS[[cn]]))
    if (length(tc) != 1L) return(NULL)
    data.table(cell_type = ct, gene = rownames(ae), contrast = cn,
               log2FC_m01 = log2((ae[, tc] + EPS) / (ae[, cc] + EPS)))
  }))
}))

# --- join + per cell_type x contrast Spearman -----------------------------------
j <- merge(m01_rows, m02, by = c("cell_type", "gene", "contrast"))
conc <- j[, .(spearman_rho = if (.N >= 10L)
                 cor(log2FC_m01, log2FC_m02, method = "spearman") else NA_real_,
              n_genes = .N), by = .(cell_type, contrast)]
setorder(conc, contrast, -spearman_rho)
fwrite(conc, out_tsv, sep = "\t")

# --- scatter: M01 vs M02 log2FC, faceted by cell_type, coloured by contrast ------
lim <- quantile(abs(c(j$log2FC_m01, j$log2FC_m02)), 0.99, na.rm = TRUE)
p <- ggplot(j, aes(log2FC_m01, log2FC_m02, colour = contrast)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_point(size = 0.4, alpha = 0.3) +
  facet_wrap(~ cell_type) +
  coord_cartesian(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  scale_colour_brewer(palette = "Set1") +
  labs(x = "M01 log2FC (mean-ratio, n=1)", y = "M02 log2FC (pseudobulk DESeq2, n=4)",
       title = "M01 vs M02 day-2 effect-size concordance (unified labels)") +
  theme_bw(base_size = 9) + theme(legend.position = "bottom")
ggsave(plot_sc, p, width = 10, height = 8, dpi = 150)

cat(sprintf("concordance: %d genes matched | %d cell_type x contrast | median rho %.3f (range %.3f-%.3f)\n",
            nrow(j), nrow(conc), median(conc$spearman_rho, na.rm = TRUE),
            min(conc$spearman_rho, na.rm = TRUE), max(conc$spearman_rho, na.rm = TRUE)))
print(conc)
