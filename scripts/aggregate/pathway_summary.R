#!/usr/bin/env Rscript
# aggregate.smk pathway track -- per-cell module scoring + summaries (flank cohort).
# Recomputes BOTH scores (UCell and Seurat AddModuleScore) on the merged 3.27M-cell
# object for one consistent, auditable pass over 4 project-priority pathways
# (config yaml, tier=primary) + 50 MSigDB Hallmark sets (msigdbr human->mouse
# orthologs, tier=exploratory). No coverage filtering: every set with >=1 panel gene
# is scored and n_set_genes / n_panel_genes / panel_coverage_frac are attached so
# downstream applies thresholds. AddModuleScore ctrl is lowered to 20 because the
# 950-gene panel gives ~40 genes/bin (nbin=24) and Seurat errors if ctrl > bin size.
# UCell runs ncores=4: that is fork (BiocParallel) parallelism, not BLAS, so the
# project's threads:1 BLAS-hygiene intent is preserved (BLAS stays pinned to 1).
# Summaries per (sample x cell_type x pathway x score_type). M02 day2 gets a limma
# test on per-sample means (~slide_id+condition, same abundance floor as the DE
# track, global BH); M01 is descriptive timecourse only. UCell-vs-AMS concordance is
# Pearson r across samples per (cell_type x pathway x dataset).
# Args: <merged.rds> <full_labels.parquet> <pathway_yaml> <out_summary.tsv>
#       <out_test.tsv> <out_concordance.tsv> <plot_heatmap> <plot_timecourse> <plot_scatter>
suppressPackageStartupMessages({
  library(Seurat)
  library(arrow)
  library(UCell)
  library(msigdbr)
  library(yaml)
  library(limma)
  library(data.table)
  library(ggplot2)
})

args        <- commandArgs(trailingOnly = TRUE)
merged_path <- args[1]
labels_path <- args[2]
yaml_path   <- args[3]
out_summary <- args[4]
out_test    <- args[5]
out_conc    <- args[6]
plot_heat   <- args[7]
plot_tc     <- args[8]
plot_scatter<- args[9]

MIN_CELLS   <- 10L
MIN_SAMPLES <- 3L
CONDS       <- c("Control", "MBRT_day2", "SBRT_day2")
UCELL_CORES <- 8L
AMS_NBIN    <- 24L
AMS_CTRL    <- 20L
SEED        <- 42L

dir.create(dirname(out_summary), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(plot_heat),  recursive = TRUE, showWarnings = FALSE)

# --- gene sets: 4 project-priority (primary) + 50 MSigDB Hallmark (exploratory) ---
prim_lists <- lapply(read_yaml(yaml_path), as.character)
prim_meta  <- data.table(pathway = names(prim_lists),
                         pathway_source = "project", tier = "primary")
hm       <- as.data.table(msigdbr(species = "Mus musculus", collection = "H"))
hm_lists <- lapply(split(hm$gene_symbol, hm$gs_name), unique)
hm_meta  <- data.table(pathway = names(hm_lists),
                       pathway_source = "MSigDB_Hallmark", tier = "exploratory")

all_sets    <- c(prim_lists, hm_lists)
set_meta    <- rbind(prim_meta, hm_meta)
set_meta[, n_set_genes := vapply(all_sets[pathway], length, integer(1))]

# --- load merged object, attach canonical aggregate labels, drop unassigned ------
# cell_subtype from the unified tier-1/2 label table (full_labels.parquet) replaces
# the stale per-sample cell_type; barcodes are the globally-unique sample_cell keys.
o      <- readRDS(merged_path)
lab    <- as.data.table(read_parquet(labels_path))[, .(cell, cell_subtype)]
lv     <- setNames(lab$cell_subtype, lab$cell)
o$cell_type <- unname(lv[colnames(o)])
o <- subset(o, cells = colnames(o)[!is.na(o$cell_type) & o$cell_type != "unassigned"])
panel  <- rownames(o)
sets_p <- lapply(all_sets, intersect, panel)            # panel-restricted gene sets
np     <- lengths(sets_p)
dropped <- names(sets_p)[np == 0L]
sets_p <- sets_p[np > 0L]                               # need >=1 panel gene to score
set_meta[, n_panel_genes := np[pathway]]
set_meta[, panel_coverage_frac := round(n_panel_genes / n_set_genes, 4)]
cat(sprintf("gene sets: %d total, %d scored, %d dropped (0 panel genes)\n",
            length(all_sets), length(sets_p), length(dropped)))

# --- score: UCell (rank-based) + Seurat AddModuleScore (binned control) ---
set.seed(SEED)
o <- AddModuleScore_UCell(o, features = sets_p, assay = "RNA", slot = "data",
                          ncores = UCELL_CORES, name = "__UC", force.gc = TRUE)
o <- AddModuleScore(o, features = sets_p, name = "__AMS",
                    nbin = AMS_NBIN, ctrl = AMS_CTRL, seed = SEED, assay = "RNA")

md   <- as.data.table(o@meta.data, keep.rownames = "cell")
uc_cols  <- paste0(names(sets_p), "__UC")
ams_cols <- paste0("__AMS", seq_along(sets_p))           # Seurat names by index, in set order
uc_mat   <- as.matrix(md[, ..uc_cols]);  colnames(uc_mat)  <- names(sets_p)
ams_mat  <- as.matrix(md[, ..ams_cols]); colnames(ams_mat) <- names(sets_p)

keys <- md[, .(cell, sample_id, cell_type, condition, treatment,
               timepoint_h, dataset, slide_id)]
samp_lkp <- unique(keys[, .(sample_id, condition, treatment, timepoint_h,
                            dataset, slide_id)])
rm(o); invisible(gc())

# --- summarize per (sample x cell_type x pathway x score_type): loop pathways ---
# (avoids a 3.27M x 54 x 2 melt; one grouped aggregation per pathway/score_type)
summarize_one <- function(vec, st) {
  d <- data.table(sample_id = keys$sample_id, cell_type = keys$cell_type, v = vec)
  a <- d[, .(mean = mean(v), sd = sd(v), median = median(v), n_cells = .N),
         by = .(sample_id, cell_type)]
  a[, score_type := st]
  a
}
summ <- rbindlist(lapply(seq_len(ncol(uc_mat)), function(j) {
  p  <- colnames(uc_mat)[j]
  s  <- rbind(summarize_one(uc_mat[, j], "UCell"), summarize_one(ams_mat[, j], "AMS"))
  s[, pathway := p]
  s
}))
summ <- merge(summ, samp_lkp, by = "sample_id", sort = FALSE)
summ <- merge(summ, set_meta, by = "pathway", sort = FALSE)

summary_out <- summ[, .(sample_id, cell_type, pathway_name = pathway, pathway_source,
                        tier, score_type, mean, sd, median, n_cells, condition,
                        timepoint_h, dataset, slide_id, n_set_genes, n_panel_genes,
                        panel_coverage_frac)]
setorder(summary_out, dataset, cell_type, pathway_name, score_type, sample_id)
fwrite(summary_out, out_summary, sep = "\t")

# --- M02 day2 limma test on per-sample means (per cell_type x score_type) -------
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
    cc  <- sc[, .N, by = condition]
    cnt[as.character(cc$condition)] <- cc$N
    nmin <- min(cnt)                                     # min usable samples per condition
    if (nmin < MIN_SAMPLES) next                         # abundance floor

    w   <- dcast(sub, pathway ~ sample_id, value.var = "mean")
    mat <- as.matrix(w[, -1]); rownames(mat) <- w$pathway
    mat <- mat[complete.cases(mat), , drop = FALSE]
    if (nrow(mat) == 0) next

    samp <- sc[match(colnames(mat), sample_id)]
    samp[, condition := factor(condition, levels = CONDS)]
    samp[, slide_id  := factor(slide_id)]
    design <- model.matrix(~ 0 + condition + slide_id, data = samp)
    colnames(design) <- make.names(colnames(design))
    if (qr(design)$rank < ncol(design)) next             # slide_id confounded

    cont <- makeContrasts(
      MBRT_vs_Ctrl = conditionMBRT_day2 - conditionControl,
      SBRT_vs_Ctrl = conditionSBRT_day2 - conditionControl,
      MBRT_vs_SBRT = conditionMBRT_day2 - conditionSBRT_day2,
      levels = design)
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
test[, padj_bh := p.adjust(pvalue, method = "BH")]       # global across all rows
test <- merge(test, set_meta[, .(pathway_name = pathway, pathway_source, tier,
                                 n_panel_genes, panel_coverage_frac)],
              by = "pathway_name", sort = FALSE)
test[, dataset := "Mutter_02"]
test_out <- test[, .(cell_type, pathway_name, pathway_source, tier, score_type,
                     contrast, estimate, se, t_stat, pvalue, padj_bh,
                     n_samples_per_group, n_panel_genes, panel_coverage_frac, dataset)]
setorder(test_out, cell_type, score_type, contrast, padj_bh, na.last = TRUE)
fwrite(test_out, out_test, sep = "\t")

# --- UCell-vs-AMS concordance: Pearson r across samples per cell_type x pathway --
wide <- dcast(summ, dataset + cell_type + pathway + sample_id ~ score_type,
              value.var = "mean")
conc <- wide[, .(pearson_r = if (.N >= 3 && sd(UCell) > 0 && sd(AMS) > 0)
                              round(cor(UCell, AMS), 4) else NA_real_,
                 n_samples = .N),
             by = .(cell_type, pathway, dataset)]
setnames(conc, "pathway", "pathway_name")
setorder(conc, dataset, cell_type, pathway_name)
fwrite(conc, out_conc, sep = "\t")

# --- plot 1: M02 primary-pathway test heatmap (UCell), * = padj<0.05 ------------
ph <- test_out[score_type == "UCell" & pathway_source == "project"]
ph[, sig := !is.na(padj_bh) & padj_bh < 0.05]
ph[, contrast := factor(contrast, levels = cm_names)]
ph[, pathway_name := factor(pathway_name, levels = rev(sort(unique(pathway_name))))]
ct_ord <- ph[contrast == "SBRT_vs_Ctrl", .(m = mean(estimate, na.rm = TRUE)),
             by = cell_type][order(m), cell_type]
ph[, cell_type := factor(cell_type, levels = ct_ord)]
p1 <- ggplot(ph, aes(cell_type, pathway_name, fill = estimate)) +
  geom_tile(colour = "grey92") +
  geom_text(data = ph[sig == TRUE], aes(label = "*"), size = 3, vjust = 0.75) +
  facet_wrap(~ contrast, ncol = 1) +
  scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick",
                       midpoint = 0, name = "mean-score\nestimate") +
  labs(x = NULL, y = NULL,
       title = "M02 day2 primary-pathway scores (UCell, limma, global BH; * padj<0.05)") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot_heat, p1, width = 11, height = 6.5, dpi = 150)

# --- plot 2: M01 primary-pathway timecourse (UCell, cell-count-weighted) --------
m1 <- summ[dataset == "Mutter_01" & tier == "primary" & score_type == "UCell"]
m1w <- m1[, .(wmean = sum(mean * n_cells) / sum(n_cells)),
          by = .(sample_id, pathway, timepoint_h, treatment)]
p2 <- ggplot(m1w, aes(timepoint_h, wmean, colour = treatment)) +
  geom_line(aes(group = treatment), linewidth = 0.4) +
  geom_point(size = 1.1) +
  facet_wrap(~ pathway, scales = "free_y") +
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "timepoint (h)", y = "mean UCell score (cell-count weighted)",
       title = "M01 primary-pathway score over time (n=1, descriptive)") +
  theme_bw(base_size = 9)
ggsave(plot_tc, p2, width = 8, height = 6, dpi = 150)

# --- plot 3: UCell vs AMS per-sample-mean concordance ---------------------------
wide[, tier := setNames(set_meta$tier, set_meta$pathway)[pathway]]
r_all <- wide[is.finite(UCell) & is.finite(AMS), round(cor(UCell, AMS), 3)]
p3 <- ggplot(wide, aes(UCell, AMS, colour = tier)) +
  geom_point(size = 0.4, alpha = 0.15) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  scale_colour_manual(values = c(primary = "firebrick", exploratory = "grey50")) +
  labs(x = "per-sample mean UCell", y = "per-sample mean AddModuleScore",
       title = sprintf("UCell vs AddModuleScore concordance (overall Pearson r = %.3f)", r_all)) +
  theme_bw(base_size = 9)
ggsave(plot_scatter, p3, width = 7, height = 6, dpi = 150)

cat(sprintf("pathway_summary: %d summary rows | test %d rows (%d padj<0.05) | conc %d rows | dropped sets: %s\n",
            nrow(summary_out), nrow(test_out), test_out[!is.na(padj_bh) & padj_bh < 0.05, .N],
            nrow(conc), if (length(dropped)) paste(dropped, collapse = ",") else "none"))
