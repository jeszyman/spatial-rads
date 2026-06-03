#!/usr/bin/env Rscript
# aggregate.smk Track 2 -- pseudobulk construction (M02 day2 only).
# Sums raw counts within each (sample x cell_type) group into one column via a
# sparse cells-by-groups indicator multiply (counts %*% G), then packages a
# SummarizedExperiment with quality columns. No filtering here -- the abundance
# floor and gene filter are applied downstream in deg_pseudobulk.R so the choice
# is explicit and auditable. Per plan-aggregate.md Track 2 inference.
# Args: <merged_typed.rds> <full_labels.parquet> <out_se.rds> <out_qc.tsv>
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(SummarizedExperiment)
  library(data.table)
  library(arrow)
})

args        <- commandArgs(trailingOnly = TRUE)
merged      <- args[1]
labels_path <- args[2]
out_se      <- args[3]
out_qc      <- args[4]

o   <- readRDS(merged)
md  <- o@meta.data
lab <- as.data.table(read_parquet(labels_path))[, .(cell, cell_subtype)]
md$cell_type <- lab$cell_subtype[match(rownames(md), lab$cell)]   # unified labels (full_labels.parquet)
keep <- which(md$dataset == "Mutter_02" & !is.na(md$cell_type))
stopifnot(length(keep) > 0)

cnt <- LayerData(o, assay = "RNA", layer = "counts")[, keep, drop = FALSE]  # 950 x Nm2
md2 <- md[keep, ]
rm(o); invisible(gc())

# --- sparse group-sum: (950 x Ncells) %*% (Ncells x Ngroups) = 950 x Ngroups ---
grp   <- factor(paste(md2$sample_id, md2$cell_type, sep = "__"))
G     <- sparseMatrix(i = seq_along(grp), j = as.integer(grp), x = 1,
                      dims = c(length(grp), nlevels(grp)),
                      dimnames = list(NULL, levels(grp)))
pb    <- as.matrix(cnt %*% G)                       # 950 x Ngroups, dense counts
storage.mode(pb) <- "integer"

# --- colData: parse group key, join sample-level fields ---
key      <- data.table(group = colnames(pb))
key[, c("sample_id", "cell_type") := tstrsplit(group, "__", fixed = TRUE)]
samp_lkp <- unique(as.data.table(md2)[, .(sample_id, condition, slide_id, timepoint_h)])
cd       <- merge(key, samp_lkp, by = "sample_id", sort = FALSE)
setkey(cd, group); cd <- cd[colnames(pb)]           # align to assay columns
stopifnot(identical(cd$group, colnames(pb)))

n_cells      <- as.integer(colSums(G))
total_counts <- as.numeric(colSums(pb))
cd[, `:=`(n_cells = n_cells,
          total_counts = total_counts,
          mean_libsize = round(total_counts / n_cells, 1))]

se <- SummarizedExperiment(
  assays  = list(counts = pb),
  colData = DataFrame(cd[, .(sample_id, cell_type, condition, slide_id,
                             timepoint_h, n_cells, total_counts, mean_libsize)],
                      row.names = cd$group),
  rowData = DataFrame(gene_symbol = rownames(pb)))

dir.create(dirname(out_se), recursive = TRUE, showWarnings = FALSE)
saveRDS(se, out_se)

# --- QC artifact: per (sample x cell_type) quality ---
qc <- cd[, .(sample_id, cell_type, condition, slide_id, n_cells, total_counts,
             mean_libsize)]
qc[, frac_genes_with_count_gt_5 := round(colMeans(pb > 5), 4)]
setorder(qc, cell_type, sample_id)
dir.create(dirname(out_qc), recursive = TRUE, showWarnings = FALSE)
fwrite(qc, out_qc, sep = "\t")

cat(sprintf("pseudobulk: %d M02 cells -> %d (sample x cell_type) columns x %d genes | %d cell types, %d samples\n",
            length(keep), ncol(pb), nrow(pb),
            uniqueN(cd$cell_type), uniqueN(cd$sample_id)))
