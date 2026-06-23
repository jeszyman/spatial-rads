#!/usr/bin/env Rscript
# aggregate.smk pathway track -- COMPUTE half (plan-differential-robustness Fix/plumbing).
# Per-cell module scoring (UCell + Seurat AddModuleScore) on the merged 3.27M-cell object
# over primary + Hallmark sets, summarized per (sample x cell_type x pathway x score_type),
# with the M02 day2 limma test and the UCell-vs-AMS concordance. This is the ~5h step; it
# is split from plotting (pathway_plots.R) so a plot-layer bug can never roll back and
# force the recompute. Outputs are the three cached TSVs the plots read back.
# Args: <merged.rds> <full_labels.parquet> <pathway_sets.tsv> <out_summary> <out_test> <out_conc>
suppressPackageStartupMessages({
  library(Seurat); library(arrow); library(UCell); library(limma); library(data.table)
})
args        <- commandArgs(trailingOnly = TRUE)
merged_path <- args[1]; labels_path <- args[2]; sets_path <- args[3]
out_summary <- args[4]; out_test <- args[5]; out_conc <- args[6]

MIN_CELLS <- 10L; MIN_SAMPLES <- 3L; CONDS <- c("Control", "MBRT_day2", "SBRT_day2")
UCELL_CORES <- 8L; AMS_NBIN <- 24L; AMS_CTRL <- 20L; SEED <- 42L
dir.create(dirname(out_summary), recursive = TRUE, showWarnings = FALSE)

gs_long  <- fread(sets_path)                            # set, tier, source, gene
all_sets <- lapply(split(gs_long$gene, gs_long$set), unique)
set_meta <- unique(gs_long[, .(pathway = set, pathway_source = source, tier)])
set_meta[, n_set_genes := vapply(all_sets[pathway], length, integer(1))]

o   <- readRDS(merged_path)
lab <- as.data.table(read_parquet(labels_path))[, .(cell, cell_subtype)]
lv  <- setNames(lab$cell_subtype, lab$cell)
o$cell_type <- unname(lv[colnames(o)])
o <- subset(o, cells = colnames(o)[!is.na(o$cell_type) & o$cell_type != "unassigned"])
panel  <- rownames(o)
sets_p <- lapply(all_sets, intersect, panel)
np     <- lengths(sets_p)
dropped <- names(sets_p)[np == 0L]
sets_p <- sets_p[np > 0L]
set_meta[, n_panel_genes := np[pathway]]
set_meta[, panel_coverage_frac := round(n_panel_genes / n_set_genes, 4)]
cat(sprintf("gene sets: %d total, %d scored, %d dropped (0 panel genes)\n",
            length(all_sets), length(sets_p), length(dropped)))

set.seed(SEED)
o <- AddModuleScore_UCell(o, features = sets_p, assay = "RNA", slot = "data",
                          ncores = UCELL_CORES, name = "__UC", force.gc = TRUE)
o <- AddModuleScore(o, features = sets_p, name = "__AMS",
                    nbin = AMS_NBIN, ctrl = AMS_CTRL, seed = SEED, assay = "RNA")
md       <- as.data.table(o@meta.data, keep.rownames = "cell")
uc_cols  <- paste0(names(sets_p), "__UC")
ams_cols <- paste0("__AMS", seq_along(sets_p))
uc_mat   <- as.matrix(md[, ..uc_cols]);  colnames(uc_mat)  <- names(sets_p)
ams_mat  <- as.matrix(md[, ..ams_cols]); colnames(ams_mat) <- names(sets_p)
keys     <- md[, .(cell, sample_id, cell_type, condition, treatment, timepoint_h, dataset, slide_id)]
samp_lkp <- unique(keys[, .(sample_id, condition, treatment, timepoint_h, dataset, slide_id)])
rm(o); invisible(gc())

summarize_one <- function(vec, st) {
  d <- data.table(sample_id = keys$sample_id, cell_type = keys$cell_type, v = vec)
  a <- d[, .(mean = mean(v), sd = sd(v), median = median(v), n_cells = .N),
         by = .(sample_id, cell_type)]
  a[, score_type := st]; a
}
summ <- rbindlist(lapply(seq_len(ncol(uc_mat)), function(j) {
  p <- colnames(uc_mat)[j]
  s <- rbind(summarize_one(uc_mat[, j], "UCell"), summarize_one(ams_mat[, j], "AMS"))
  s[, pathway := p]; s
}))
summ <- merge(summ, samp_lkp, by = "sample_id", sort = FALSE)
summ <- merge(summ, set_meta, by = "pathway", sort = FALSE)

summary_out <- summ[, .(sample_id, cell_type, pathway_name = pathway, pathway_source,
                        tier, score_type, mean, sd, median, n_cells, condition, treatment,
                        timepoint_h, dataset, slide_id, n_set_genes, n_panel_genes,
                        panel_coverage_frac)]   # treatment added so pathway_plots can build the M01 timecourse
setorder(summary_out, dataset, cell_type, pathway_name, score_type, sample_id)
fwrite(summary_out, out_summary, sep = "\t")

m2 <- summ[dataset == "Mutter_02" & condition %in% CONDS]
m2[, condition := factor(condition, levels = CONDS)]
cm_names <- c("MBRT_vs_Ctrl", "SBRT_vs_Ctrl", "MBRT_vs_SBRT")
test_rows <- list()
for (ct in sort(unique(m2$cell_type))) {
  for (st in c("UCell", "AMS")) {
    sub <- m2[cell_type == ct & score_type == st & n_cells >= MIN_CELLS]
    if (nrow(sub) == 0) next
    sc  <- unique(sub[, .(sample_id, condition, slide_id)])
    cnt <- setNames(rep(0L, length(CONDS)), CONDS)
    cc  <- sc[, .N, by = condition]; cnt[as.character(cc$condition)] <- cc$N
    nmin <- min(cnt); if (nmin < MIN_SAMPLES) next
    w   <- dcast(sub, pathway ~ sample_id, value.var = "mean")
    mat <- as.matrix(w[, -1]); rownames(mat) <- w$pathway
    mat <- mat[complete.cases(mat), , drop = FALSE]; if (nrow(mat) == 0) next
    samp <- sc[match(colnames(mat), sample_id)]
    samp[, condition := factor(condition, levels = CONDS)]; samp[, slide_id := factor(slide_id)]
    design <- model.matrix(~ 0 + condition + slide_id, data = samp)
    colnames(design) <- make.names(colnames(design))
    if (qr(design)$rank < ncol(design)) next
    cont <- makeContrasts(
      MBRT_vs_Ctrl = conditionMBRT_day2 - conditionControl,
      SBRT_vs_Ctrl = conditionSBRT_day2 - conditionControl,
      MBRT_vs_SBRT = conditionMBRT_day2 - conditionSBRT_day2, levels = design)
    fit2 <- eBayes(contrasts.fit(lmFit(mat, design), cont), robust = TRUE)
    for (cn in cm_names) {
      est <- fit2$coefficients[, cn]
      se  <- fit2$stdev.unscaled[, cn] * sqrt(fit2$s2.post)
      test_rows[[length(test_rows) + 1]] <- data.table(
        cell_type = ct, pathway_name = rownames(mat), score_type = st, contrast = cn,
        estimate = est, se = se, t_stat = fit2$t[, cn], pvalue = fit2$p.value[, cn],
        n_samples_per_group = nmin)
    }
  }
}
test <- rbindlist(test_rows)
test[, padj_bh := p.adjust(pvalue, method = "BH")]
test <- merge(test, set_meta[, .(pathway_name = pathway, pathway_source, tier,
                                 n_panel_genes, panel_coverage_frac)], by = "pathway_name", sort = FALSE)
test[, dataset := "Mutter_02"]
test_out <- test[, .(cell_type, pathway_name, pathway_source, tier, score_type, contrast,
                     estimate, se, t_stat, pvalue, padj_bh, n_samples_per_group,
                     n_panel_genes, panel_coverage_frac, dataset)]
setorder(test_out, cell_type, score_type, contrast, padj_bh, na.last = TRUE)
fwrite(test_out, out_test, sep = "\t")

wide <- dcast(summ, dataset + cell_type + pathway + sample_id ~ score_type, value.var = "mean")
conc <- wide[, .(pearson_r = if (.N >= 3 && sd(UCell) > 0 && sd(AMS) > 0)
                              round(cor(UCell, AMS), 4) else NA_real_, n_samples = .N),
             by = .(cell_type, pathway, dataset)]
setnames(conc, "pathway", "pathway_name")
setorder(conc, dataset, cell_type, pathway_name)
fwrite(conc, out_conc, sep = "\t")

cat(sprintf("pathway_scores: %d summary rows | test %d rows (%d padj<0.05) | conc %d rows | dropped: %s\n",
            nrow(summary_out), nrow(test_out), test_out[!is.na(padj_bh) & padj_bh < 0.05, .N],
            nrow(conc), if (length(dropped)) paste(dropped, collapse = ",") else "none"))
