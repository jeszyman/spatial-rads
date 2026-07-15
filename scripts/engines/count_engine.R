#!/usr/bin/env Rscript
# Negative-binomial pseudobulk count engine (modularization Task 3). Per-cell-type DE on
# the (sample x cell_type) pseudobulk SummarizedExperiment via DESeq2, reading the arm
# contrasts + their reference levels from the resolved comparison registry -- no embedded
# arm literals. Ports deg_pseudobulk.R's statistics verbatim, with three changes: contrasts
# from the registry, thresholds from engine_params.yaml, and NO padj emitted (correction is
# the results tier's job).
#
# apeglm two-fit pattern, preserved: lfcShrink(type="apeglm") accepts only coef=, not an
# arbitrary contrast vector, so a shrunk MBRT-vs-SBRT estimate needs SBRT as the DESeq2
# reference. The engine therefore builds ONE DESeq fit per distinct reference level named
# in the registry (here Control and SBRT_day2) and pulls each contrast from the fit whose
# reference is that contrast's denominator. The DESeq2 design (~ slide_id + condition, with
# intercept) is a method-internal parameterization, not a per-comparison literal; only the
# arm identities and reference levels come from the registry.
# Args: <pseudobulk_se.rds> <cohort> <comparisons.tsv> <engine_params.yaml>
#       <out_degs.tsv> <out_skipped.tsv>
suppressPackageStartupMessages({
  library(SummarizedExperiment); library(DESeq2); library(apeglm)
  library(data.table); library(yaml)
})
a <- commandArgs(trailingOnly = TRUE)
se_path <- a[1]; coh <- a[2]; comp_path <- a[3]; params_path <- a[4]
out_degs <- a[5]; out_skipped <- a[6]

## ---- thresholds (study-wide, from engine config) ----
P <- read_yaml(params_path)$count_engine
MIN_CELLS <- P$min_cells; MIN_SAMPLES <- P$min_samples
GENE_MIN_CT <- P$gene_min_count; GENE_MIN_SMP <- P$gene_min_samples

## ---- contrasts from the registry: name -> (num, den=ref) ----
comp <- fread(comp_path)
fam  <- comp[cohort == coh & kind == "sample" & !is.na(contrast_num_level)]
fam  <- unique(fam, by = c("name", "contrast_num_level", "contrast_den_level"))
stopifnot(nrow(fam) >= 1)
CONTRASTS <- lapply(seq_len(nrow(fam)), function(i)
  list(name = fam$name[i], num = fam$contrast_num_level[i], den = fam$contrast_den_level[i]))
# condition levels: primary reference (den shared by the vs-baseline contrasts) first, so the
# default fit's reference matches deg_pseudobulk.R (Control), then the remaining arms sorted.
den_counts  <- sort(table(vapply(CONTRASTS, `[[`, "", "den")), decreasing = TRUE)
primary_ref <- names(den_counts)[1]
all_levels  <- unique(c(vapply(CONTRASTS, `[[`, "", "num"), vapply(CONTRASTS, `[[`, "", "den")))
CONDS       <- c(primary_ref, sort(setdiff(all_levels, primary_ref)))
ref_levels  <- unique(vapply(CONTRASTS, `[[`, "", "den"))     # one DESeq fit per distinct ref

dir.create(dirname(out_degs), recursive = TRUE, showWarnings = FALSE)
se <- readRDS(se_path)
cd <- as.data.table(as.data.frame(colData(se)), keep.rownames = "col")

deg_rows <- list(); skip_rows <- list()
usable   <- cd[n_cells >= MIN_CELLS, .(n_pass = .N), by = .(cell_type, condition)]
all_types <- sort(unique(cd$cell_type))

for (ct in all_types) {
  u  <- usable[cell_type == ct]
  np <- setNames(rep(0L, length(CONDS)), CONDS); np[u$condition] <- u$n_pass
  failing <- names(np)[np < MIN_SAMPLES]
  if (length(failing) > 0) {                                   # abundance floor -> skip type
    for (cc in CONTRASTS) {
      bad <- intersect(c(cc$num, cc$den), failing); if (!length(bad)) bad <- failing[1]
      ncells_bad <- cd[cell_type == ct & condition == bad[1], sum(n_cells)]
      skip_rows[[length(skip_rows) + 1]] <- data.table(
        cell_type = ct, contrast = cc$name,
        reason = sprintf("abundance_floor: condition %s has %d/4 samples >=%d cells (need %d)",
                         bad[1], np[[bad[1]]], MIN_CELLS, MIN_SAMPLES),
        n_cells_in_failed_group = ncells_bad)
    }
    next
  }

  cols <- cd[cell_type == ct & n_cells >= MIN_CELLS, col]
  sub  <- se[, cols]
  colData(sub)$condition <- factor(colData(sub)$condition, levels = CONDS)
  colData(sub)$slide_id  <- factor(colData(sub)$slide_id)

  mm <- model.matrix(~ slide_id + condition, data = as.data.frame(colData(sub)))
  if (qr(mm)$rank < ncol(mm)) {                                # design not full rank
    for (cc in CONTRASTS) skip_rows[[length(skip_rows) + 1]] <- data.table(
      cell_type = ct, contrast = cc$name,
      reason = "design not full rank after abundance filter (slide_id confounded)",
      n_cells_in_failed_group = NA_integer_)
    next
  }

  cm     <- assay(sub)
  keep_g <- rowSums(cm >= GENE_MIN_CT) >= GENE_MIN_SMP
  sub    <- sub[keep_g, ]
  n_used <- ncol(sub)
  df_res <- n_used - qr(mm)$rank                               # residual df for the fit

  # one DESeq fit per distinct reference level (generalizes the Control/SBRT two-fit)
  fits <- setNames(lapply(ref_levels, function(rl) {
    d <- DESeqDataSetFromMatrix(assay(sub), as.data.frame(colData(sub)),
                                design = ~ slide_id + condition)
    colData(d)$condition <- relevel(factor(colData(d)$condition, levels = CONDS), ref = rl)
    DESeq(d, quiet = TRUE)
  }), ref_levels)

  for (cc in CONTRASTS) {
    coef <- sprintf("condition_%s_vs_%s", cc$num, cc$den)
    obj  <- fits[[cc$den]]
    res  <- results(obj, name = coef, independentFiltering = TRUE)  # baseMean, stat, pvalue
    shr  <- lfcShrink(obj, coef = coef, type = "apeglm", quiet = TRUE)  # log2FC, lfcSE
    deg_rows[[length(deg_rows) + 1]] <- data.table(
      comparison = coh, contrast = cc$name, unit = ct, feature_type = "gene",
      feature_id = rownames(res),
      estimate = shr$log2FoldChange, se = shr$lfcSE, df = df_res,
      stat = res$stat, p = res$pvalue,
      baseMean = res$baseMean, n_samples_used = n_used)
  }
}

degs <- rbindlist(deg_rows)
setorder(degs, unit, contrast, p, na.last = TRUE)
fwrite(degs, out_degs, sep = "\t")
skipped <- if (length(skip_rows)) rbindlist(skip_rows) else data.table(
  cell_type = character(), contrast = character(),
  reason = character(), n_cells_in_failed_group = integer())
fwrite(skipped, out_skipped, sep = "\t")

cat(sprintf("count_engine[%s]: %d cell types tested, %d skipped | %d (gene x ct x contrast) rows\n",
            coh, uniqueN(degs$unit), uniqueN(skipped$cell_type), nrow(degs)))
