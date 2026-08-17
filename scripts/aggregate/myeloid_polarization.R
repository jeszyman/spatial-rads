#!/usr/bin/env Rscript
# aggregate.smk myeloid M1/M2 polarization. Re-tests the dev "MBRT day-2 M2 skew"
# finding (dev/mbrt_vs_sbrt/11_m1_m2_polarization.R) on the unified cross-dataset
# labels. Restricts to the Macrophages cell_subtype (the canonical M1/M2 substrate),
# UCell-scores the 16+16 canonical M1/M2 marker panels (panel-filtered, kept genes
# logged), and aggregates per-cell scores to per-sample means + M2/M1 ratio. M02
# day-2 gets a limma test on the per-sample metric matrix (M1, M2, M2/M1 ratio;
# ~0+condition+slide_id, 3 contrasts, 95% CIs, BH) - few outcome rows so eBayes
# moderation is negligible. M01 is descriptive (n=1).
# Args: <merged.rds> <full_labels.parquet> <samples.tsv> <out_scores.tsv> <out_test.tsv> <plot_m02>
suppressPackageStartupMessages({
  library(Seurat); library(arrow); library(UCell)
  library(data.table); library(ggplot2)
})
a <- commandArgs(trailingOnly = TRUE)
merged_path <- a[1]; labels_path <- a[2]; samples_path <- a[3]
out_scores <- a[4]; out_lm_input <- a[5]; plot_m02 <- a[6]
SEED <- 42L; UCELL_CORES <- 8L

m1_markers <- c("Nos2","Tnf","Il6","Il1b","Il12b","Cxcl9","Cxcl10","Cxcl11","Cd86",
                "Cd80","Tlr2","Tlr4","Stat1","Irf5","Hif1a","Nfkb1")
m2_markers <- c("Cd163","Mrc1","Arg1","Il4ra","Mgl2","Retnla","Chi3l3","Ym1",
                "Stab1","Vegfa","Ccl22","Mmp9","Tgfb1","Stat6","Klf4","Cd206")

dir.create(dirname(out_scores), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(plot_m02), recursive = TRUE, showWarnings = FALSE)

o   <- readRDS(merged_path)
lab <- as.data.table(read_parquet(labels_path))[, .(cell, cell_subtype)]
lv  <- setNames(lab$cell_subtype, lab$cell)
o$cell_subtype <- unname(lv[colnames(o)])

# samples.tsv is the single source of truth for sample-level metadata (condition,
# timepoint, dataset, slide). merged.rds is a heavy intermediate rebuilt independently of
# it and its baked-in condition/timepoint_h/dataset/slide_id copies can drift out of sync;
# join all four from samples.tsv by sample_id (dataset = its "name" column, values
# "Mutter_01"/"Mutter_02") and overwrite the merged object's copies rather than trust them.
ss <- fread(samples_path)
scol <- names(ss)[grepl("sample", names(ss), ignore.case = TRUE)][1]
smeta <- unique(ss[, .(sample_id = get(scol), condition, timepoint_h, dataset = name, slide_id)])
midx <- match(o$sample_id, smeta$sample_id)
o$condition   <- smeta$condition[midx]
o$timepoint_h <- smeta$timepoint_h[midx]
o$dataset     <- smeta$dataset[midx]
o$slide_id    <- smeta$slide_id[midx]
stopifnot(!anyNA(o$sample_id), !anyNA(o$condition), !anyNA(o$slide_id),
          !anyNA(o$dataset), !anyNA(o$timepoint_h))

mye <- subset(o, cells = colnames(o)[!is.na(o$cell_subtype) &
                                     o$cell_subtype == "Macrophages"])
rm(o); invisible(gc())

panel <- rownames(mye)
m1_in <- intersect(m1_markers, panel); m2_in <- intersect(m2_markers, panel)
cat(sprintf("Macrophages: %d cells | M1 %d/%d on panel (%s) | M2 %d/%d on panel (%s)\n",
            ncol(mye), length(m1_in), length(m1_markers), paste(m1_in, collapse=","),
            length(m2_in), length(m2_markers), paste(m2_in, collapse=",")))

set.seed(SEED)
mye <- AddModuleScore_UCell(mye, features = list(M1 = m1_in, M2 = m2_in),
                            assay = "RNA", slot = "data",
                            ncores = UCELL_CORES, name = "_score", force.gc = TRUE)

md <- as.data.table(mye@meta.data, keep.rownames = "cell")
score <- md[!is.na(M1_score) & !is.na(M2_score),
            .(M1_mean = mean(M1_score), M2_mean = mean(M2_score),
              M2_M1_ratio = mean(M2_score) / mean(M1_score), n_cells = .N),
            by = .(sample_id, condition, slide_id, timepoint_h, dataset)]
setorder(score, dataset, timepoint_h, sample_id)
fwrite(score, out_scores, sep = "\t")

# --- engine input: long (sample_id, feature_id, value, condition, slide_id), all flank
# samples. The lm_engine applies cohort_samples.tsv whitelist filtering. ---
metrics <- c("M1_mean", "M2_mean", "M2_M1_ratio")
flank <- score[!grepl("^Tongue", condition)]
lm_input <- melt(flank[, c("sample_id", "condition", "slide_id", metrics), with = FALSE],
                 id.vars = c("sample_id", "condition", "slide_id"),
                 variable.name = "feature_id", value.name = "value")
fwrite(lm_input, out_lm_input, sep = "\t")

# --- plot: M1 / M2 / ratio by arm (M02 day2) -----------------------------------
m2 <- score[dataset == "Mutter_02" & timepoint_h == 48L]
pd <- melt(m2[, c("sample_id", "condition", metrics), with = FALSE],
           id.vars = c("sample_id", "condition"),
           variable.name = "metric", value.name = "value")
p <- ggplot(pd, aes(condition, value, fill = condition)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_point(size = 1, position = position_jitter(width = 0.12)) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_fill_brewer(palette = "Set1") +
  labs(x = NULL, y = NULL,
       title = "M02 day-2 macrophage M1/M2 polarization by arm (n=2/arm)") +
  theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(plot_m02, p, width = 8, height = 4.5, dpi = 150)

cat(sprintf("myeloid M1/M2: %d samples | lm input %d rows (%d metrics x %d flank samples)\n",
            nrow(score), nrow(lm_input), length(metrics), nrow(flank)))
