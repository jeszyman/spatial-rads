#!/usr/bin/env Rscript
# aggregate.smk Track 2 inference -- pseudobulk DE (M02 day2, DESeq2).
# Abundance floor (approved 2026-05-31): a cell type is tested only if it has
# >=10 cells in >=3 of 4 samples in EVERY condition -- otherwise the whole type
# is skipped (logged), because a pseudobulk profile of a handful of cells is
# sampling noise. Within a tested type, only sample columns with >=10 cells are
# used. Per type: gene filter (>=10 counts in >=4 samples) -> DESeq2
# (~slide_id+condition) -> 3 contrasts via results() + apeglm lfcShrink (two fits:
# Control ref, then SBRT ref for MBRT-vs-SBRT). padj is DESeq2's per-contrast BH.
# Args: <pseudobulk_se.rds> <out_degs.tsv> <out_summary.tsv> <out_skipped.tsv> <volcano_dir>
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(DESeq2)
  library(apeglm)
  library(data.table)
  library(ggplot2)
})

args        <- commandArgs(trailingOnly = TRUE)
se_path     <- args[1]
out_degs    <- args[2]
out_summary <- args[3]
out_skipped <- args[4]
volcano_dir <- args[5]

MIN_CELLS    <- 10L     # per (sample x cell_type) to count a usable replicate
MIN_SAMPLES  <- 3L      # usable replicates required per condition
GENE_MIN_CT  <- 10L     # gene-level: >= this many counts ...
GENE_MIN_SMP <- 4L      # ... in >= this many samples
CONDS        <- c("Control", "MBRT_day2", "SBRT_day2")
CONTRASTS    <- list(
  MBRT_vs_Ctrl = c("MBRT_day2", "Control"),
  SBRT_vs_Ctrl = c("SBRT_day2", "Control"),
  MBRT_vs_SBRT = c("MBRT_day2", "SBRT_day2"))

dir.create(dirname(out_degs), recursive = TRUE, showWarnings = FALSE)
dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)

se <- readRDS(se_path)
cd <- as.data.table(as.data.frame(colData(se)), keep.rownames = "col")

deg_rows  <- list()
summ_rows <- list()
skip_rows <- list()

sanitize <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

# --- abundance floor: usable replicates per (cell_type, condition) ---
usable <- cd[n_cells >= MIN_CELLS, .(n_pass = .N), by = .(cell_type, condition)]
all_types <- sort(unique(cd$cell_type))

for (ct in all_types) {
  u  <- usable[cell_type == ct]
  np <- setNames(rep(0L, length(CONDS)), CONDS)
  np[u$condition] <- u$n_pass
  failing <- names(np)[np < MIN_SAMPLES]

  if (length(failing) > 0) {                         # floor failure -> skip type
    for (cn in names(CONTRASTS)) {
      grp <- CONTRASTS[[cn]]
      bad <- intersect(grp, failing)
      if (length(bad) == 0) bad <- failing[1]        # type-wide skip, note worst
      ncells_bad <- cd[cell_type == ct & condition == bad[1], sum(n_cells)]
      skip_rows[[length(skip_rows) + 1]] <- data.table(
        cell_type = ct, contrast = cn,
        reason = sprintf("abundance_floor: condition %s has %d/4 samples >=%d cells (need %d)",
                         bad[1], np[[bad[1]]], MIN_CELLS, MIN_SAMPLES),
        n_cells_in_failed_group = ncells_bad)
    }
    next
  }

  # --- subset to usable sample columns for this cell type ---
  cols <- cd[cell_type == ct & n_cells >= MIN_CELLS, col]
  sub  <- se[, cols]
  colData(sub)$condition <- factor(colData(sub)$condition, levels = CONDS)
  colData(sub)$slide_id  <- factor(colData(sub)$slide_id)

  # design estimability guard
  mm <- model.matrix(~ slide_id + condition, data = as.data.frame(colData(sub)))
  if (qr(mm)$rank < ncol(mm)) {
    for (cn in names(CONTRASTS)) skip_rows[[length(skip_rows) + 1]] <- data.table(
      cell_type = ct, contrast = cn,
      reason = "design not full rank after abundance filter (slide_id confounded)",
      n_cells_in_failed_group = NA_integer_)
    next
  }

  # gene filter: >=GENE_MIN_CT counts in >=GENE_MIN_SMP samples
  cm     <- assay(sub)
  keep_g <- rowSums(cm >= GENE_MIN_CT) >= GENE_MIN_SMP
  sub    <- sub[keep_g, ]

  n_used   <- ncol(sub)
  min_libs <- min(colSums(assay(sub)))
  mean_nc  <- round(mean(colData(sub)$n_cells), 1)

  dds <- DESeqDataSetFromMatrix(assay(sub), as.data.frame(colData(sub)),
                                design = ~ slide_id + condition)
  dds <- DESeq(dds, quiet = TRUE)
  dds_sb <- dds
  colData(dds_sb)$condition <- relevel(colData(dds_sb)$condition, ref = "SBRT_day2")
  dds_sb <- DESeq(dds_sb, quiet = TRUE)

  for (cn in names(CONTRASTS)) {
    coef <- if (cn == "MBRT_vs_Ctrl") "condition_MBRT_day2_vs_Control"
            else if (cn == "SBRT_vs_Ctrl") "condition_SBRT_day2_vs_Control"
            else "condition_MBRT_day2_vs_SBRT_day2"
    obj  <- if (cn == "MBRT_vs_SBRT") dds_sb else dds
    res  <- results(obj, name = coef)                      # baseMean, stat, pvalue, padj
    shr  <- lfcShrink(obj, coef = coef, type = "apeglm", quiet = TRUE)  # log2FC, lfcSE

    dt <- data.table(
      cell_type = ct, contrast = cn, gene = rownames(res),
      log2FC = shr$log2FoldChange, lfcSE = shr$lfcSE,
      stat = res$stat, pvalue = res$pvalue, padj = res$padj,
      baseMean = res$baseMean,
      n_samples_used = n_used, min_counts_per_sample = min_libs,
      mean_n_cells_per_sample = mean_nc, dataset = "Mutter_02")
    deg_rows[[length(deg_rows) + 1]] <- dt

    sig  <- dt[!is.na(padj) & padj < 0.05]
    up   <- sig[log2FC > 0][order(-log2FC), head(gene, 5)]
    down <- sig[log2FC < 0][order(log2FC), head(gene, 5)]
    summ_rows[[length(summ_rows) + 1]] <- data.table(
      cell_type = ct, contrast = cn,
      n_genes_tested = dt[!is.na(pvalue), .N],
      n_padj_05 = nrow(sig),
      n_padj_05_lfc_1 = sig[abs(log2FC) > 1, .N],
      top5_up = paste(up, collapse = ","),
      top5_down = paste(down, collapse = ","))

    # volcano (build plot data as a copy -- do not mutate the stored degs rows)
    vd <- dt[!is.na(pvalue), .(log2FC, pvalue, padj)]
    vd[, sig := !is.na(padj) & padj < 0.05]
    p <- ggplot(vd, aes(log2FC, -log10(pvalue), colour = sig)) +
      geom_point(size = 0.6, alpha = 0.6) +
      geom_vline(xintercept = c(-1, 1), linetype = 2, colour = "grey60") +
      scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "firebrick"),
                          name = "padj<0.05") +
      labs(title = sprintf("%s  |  %s", ct, cn), x = "log2FC (apeglm)",
           y = "-log10 p") +
      theme_bw(base_size = 9)
    ggsave(file.path(volcano_dir, sprintf("volcano_%s_%s.png", sanitize(ct), cn)),
           p, width = 5, height = 4, dpi = 130)
  }
}

degs <- rbindlist(deg_rows)
setorder(degs, cell_type, contrast, padj, na.last = TRUE)
fwrite(degs, out_degs, sep = "\t")
fwrite(rbindlist(summ_rows), out_summary, sep = "\t")
skipped <- if (length(skip_rows)) rbindlist(skip_rows) else data.table(
  cell_type = character(), contrast = character(),
  reason = character(), n_cells_in_failed_group = integer())
fwrite(skipped, out_skipped, sep = "\t")

cat(sprintf("deg_pseudobulk: %d cell types tested, %d skipped | %d (gene x ct x contrast) rows | %d padj<0.05\n",
            uniqueN(degs$cell_type), uniqueN(skipped$cell_type),
            nrow(degs), degs[!is.na(padj) & padj < 0.05, .N]))
