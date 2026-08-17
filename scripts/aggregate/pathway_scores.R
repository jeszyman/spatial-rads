#!/usr/bin/env Rscript
# aggregate.smk pathway track -- COMPUTE half (plan-differential-robustness Fix/plumbing;
# arm test split out to pathway_arm_test.R for cross-timepoint parity, task 3).
# Per-cell module scoring (UCell + Seurat AddModuleScore) on the merged 3.27M-cell object
# over primary + Hallmark sets, summarized per (sample x cell_type x pathway x score_type),
# plus the UCell-vs-AMS concordance. This is the ~5h step; it is split from plotting
# (pathway_plots.R) and from the per-cohort arm test (pathway_arm_test.R) so neither a
# plot-layer bug nor a design/contrast change can ever roll back and force the recompute.
# Outputs are the two cached TSVs pathway_arm_test.R and pathway_plots.R read back.
# Args: <merged.rds> <full_labels.parquet> <pathway_sets.tsv> <out_summary> <out_conc>
suppressPackageStartupMessages({
  library(Seurat); library(arrow); library(UCell); library(data.table)
})
args        <- commandArgs(trailingOnly = TRUE)
merged_path <- args[1]; labels_path <- args[2]; sets_path <- args[3]
out_summary <- args[4]; out_conc <- args[5]

UCELL_CORES <- 8L; AMS_NBIN <- 24L; AMS_CTRL <- 20L; SEED <- 42L
# maxRank truncates the zero-tail (genes ranked below it are treated as unexpressed).
# UCell requires maxRank >= the largest signature; Fibrosis_remodeling has 136 on-panel
# genes, so 136 is the floor. The default 1500 is tuned for ~20k-gene scRNA-seq and never
# truncates on this 950-gene panel (washing out sparse scores); UCell's spatial guidance
# is "at most the number of probes in the panel". 136 is the tightest valid truncation
# without altering the frozen gene sets (UCell/pyUCell 2026 CosMx guidance).
UCELL_MAXRANK <- 136L
dir.create(dirname(out_summary), recursive = TRUE, showWarnings = FALSE)

gs_long  <- fread(sets_path)                            # set, tier, source, gene
all_sets <- lapply(split(gs_long$gene, gs_long$set), unique)
# Myeloid M1/M2 polarization marker sets (macrophage state programs) scored alongside
# the curated pathway sets so the myeloid_M2 hypothesis has pathway evidence. Markers
# match scripts/aggregate/myeloid_polarization.R; UCell intersects with the panel.
all_sets[["M1_polarization"]] <- c("Nos2","Tnf","Il6","Il1b","Il12b","Cxcl9","Cxcl10",
  "Cxcl11","Cd86","Cd80","Tlr2","Tlr4","Stat1","Irf5","Hif1a","Nfkb1")
all_sets[["M2_polarization"]] <- c("Cd163","Mrc1","Arg1","Il4ra","Mgl2","Retnla","Chi3l3",
  "Ym1","Stab1","Vegfa","Ccl22","Mmp9","Tgfb1","Stat6","Klf4","Cd206")
set_meta <- unique(gs_long[, .(pathway = set, pathway_source = source, tier)])
set_meta <- rbind(set_meta, data.table(
  pathway = c("M1_polarization","M2_polarization"),
  pathway_source = "myeloid_markers", tier = "primary"))
set_meta[, n_set_genes := vapply(all_sets[pathway], length, integer(1))]

o   <- readRDS(merged_path)
lab <- as.data.table(read_parquet(labels_path))[, .(cell, cell_subtype)]
lv  <- setNames(lab$cell_subtype, lab$cell)
o$cell_type <- unname(lv[colnames(o)])

# samples.tsv is the single source of truth for sample-level metadata (condition,
# timepoint, dataset, slide). merged.rds is a heavy intermediate rebuilt independently of
# it and its baked-in condition/timepoint_h/dataset/slide_id copies can drift out of sync;
# join all four from samples.tsv by sample_id (dataset = its "name" column, values
# "Mutter_01"/"Mutter_02") and overwrite the merged object's copies rather than trust them.
ss <- fread("results/data_model/samples.tsv")
scol <- names(ss)[grepl("sample", names(ss), ignore.case = TRUE)][1]
smeta <- unique(ss[, .(sample_id = get(scol), condition, timepoint_h, dataset = name, slide_id)])
midx <- match(o$sample_id, smeta$sample_id)
o$condition   <- smeta$condition[midx]
o$timepoint_h <- smeta$timepoint_h[midx]
o$dataset     <- smeta$dataset[midx]
o$slide_id    <- smeta$slide_id[midx]
# Fail-fast: assert EVERY column the downstream arm test (pathway_arm_test.R) needs is
# populated, NOW (seconds) not after the ~5h scoring.
stopifnot(!anyNA(o$sample_id), !anyNA(o$condition), !anyNA(o$slide_id),
          !anyNA(o$dataset), !anyNA(o$cell_type))

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
                          maxRank = UCELL_MAXRANK,
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

wide <- dcast(summ, dataset + cell_type + pathway + sample_id ~ score_type, value.var = "mean")
conc <- wide[, .(pearson_r = if (.N >= 3 && sd(UCell) > 0 && sd(AMS) > 0)
                              round(cor(UCell, AMS), 4) else NA_real_, n_samples = .N),
             by = .(cell_type, pathway, dataset)]
setnames(conc, "pathway", "pathway_name")
setorder(conc, dataset, cell_type, pathway_name)
fwrite(conc, out_conc, sep = "\t")

cat(sprintf("pathway_scores: %d summary rows | conc %d rows | dropped: %s\n",
            nrow(summary_out), nrow(conc),
            if (length(dropped)) paste(dropped, collapse = ",") else "none"))
