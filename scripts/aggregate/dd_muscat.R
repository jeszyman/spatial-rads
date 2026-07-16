#!/usr/bin/env Rscript
# muscat differential-detection (runs INSIDE muscat.sif via apptainer). Tests whether
# the FRACTION of cells expressing a gene changes between arms, per cell type. Uses
# muscat's detection aggregation (num.detected) + edgeR pseudobulk DS (pbDS), which is
# the benchmarked differential-detection method (Gilis/Clement 2025). Reported ONLY as
# a detection change, never composition-vs-regulation.
suppressPackageStartupMessages({
  library(muscat); library(SingleCellExperiment); library(data.table)
})
a <- commandArgs(trailingOnly = TRUE)
sce_path <- a[1]; out_p <- a[2]

sce <- readRDS(sce_path)
sce$condition <- factor(sce$condition, levels = c("Control","MBRT_day2","SBRT_day2"))
sce <- prepSCE(sce, kid = "cell_subtype", sid = "sample_id", gid = "condition", drop = TRUE)

# Detection counts: number of cells expressing each gene per sample x cluster.
pb_det <- aggregateData(sce, assay = "counts", fun = "num.detected",
                        by = c("cluster_id", "sample_id"))

# Design keyed to the pseudobulk's OWN sample order (colData(pb_det)); rownames must
# match pb sample columns or pbDS subscripting fails. Cell-means on group_id: this is
# the exploratory detection descriptor; the slide-adjusted arm effect lives in the
# primary DESeq2 DE (count_engine ~0+condition+slide_id), not here.
pbi <- as.data.frame(colData(pb_det))
pbi$group_id <- factor(pbi$group_id, levels = c("Control","MBRT_day2","SBRT_day2"))
design <- model.matrix(~ 0 + group_id, data = pbi)
rownames(design) <- rownames(pbi)
colnames(design) <- levels(pbi$group_id)
CN <- c(MBRT_vs_Ctrl = "MBRT_day2 - Control",
        SBRT_vs_Ctrl = "SBRT_day2 - Control",
        MBRT_vs_SBRT = "MBRT_day2 - SBRT_day2")
cm <- limma::makeContrasts(contrasts = CN, levels = design); colnames(cm) <- names(CN)

res <- pbDS(pb_det, design = design, contrast = cm, method = "edgeR", verbose = FALSE)

rows <- list()
for (cn in names(CN)) {
  tb <- res$table[[cn]]
  for (ct in names(tb)) {
    d <- as.data.table(tb[[ct]])
    if (!nrow(d)) next
    rows[[length(rows)+1]] <- d[, .(cell_type = ct, contrast = cn, gene = gene,
                                    dd_log2fc = logFC, dd_p = p_val)]
  }
}
out <- rbindlist(rows)
out[, dd_padj := p.adjust(dd_p, "BH"), by = .(cell_type, contrast)]
fwrite(out, out_p, sep = "\t")
cat(sprintf("differential_detection: %d rows, %d sig (BH<0.05)\n",
            nrow(out), out[dd_padj < 0.05, .N]))
