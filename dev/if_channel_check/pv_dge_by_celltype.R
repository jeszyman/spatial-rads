suppressPackageStartupMessages({
  library(arrow); library(data.table); library(Seurat); library(dplyr)
})

peak_fovs   <- unique(c(224, 209, 201, 186, 172, 225, 229, 215, 175,
                        176, 206, 181, 168, 167))
valley_fovs <- c(199, 213, 227, 188, 204, 218, 178, 165, 170, 231, 189,
                 173, 154)

OUT <- "/home/jeszyman/repos/spatial-rads/dev/if_channel_check"

meta <- as.data.table(read_parquet("/tmp/mutter01_meta.parquet"))
counts_df <- as.data.table(read_parquet("/tmp/mutter01_counts.parquet"))

# Subset to Block_21 (MBRT 4h) only
s4 <- meta[Slide == "20250529_214712_S4" & Block == "Block_21"]
s4[, class := fifelse(fov %in% peak_fovs, "peak",
                      fifelse(fov %in% valley_fovs, "valley", NA_character_))]
keep <- s4[!is.na(class)]
cat("Peak cells:", sum(keep$class == "peak"),
    "| Valley cells:", sum(keep$class == "valley"), "\n")

# Map Main_cell_Types to radiobiology classes
setnames(keep, "ImmuneAtlas_ImmGen_Main_cell_Types", "main_type")
classify <- function(x) {
  out <- rep(NA_character_, length(x))
  out[x == "a"] <- "Tumor cells (a)"
  out[x %in% c("B.cell", "Memory.B", "Plasma", "Plasmablast",
               "GC_centroblasts", "GC_centrocyes", "Spleen.CD19")] <- "B cells"
  out[x %in% c("CD8.T.cell", "Spleen.Naive.CD4", "Spleen.Naive.CD8",
               "Spleen.CD4Act.48hrs", "Spleen.LN.Naive.CD4",
               "Spleen.Treg", "Colon.Treg.Nrplo")] <- "T cells"
  out[x %in% c("Macrophage", "PerC.macrophage", "spleen.red.pulp.macs",
               "small.peritoneal.macs", "Microglia")] <- "Macrophages"
  out[x %in% c("Ly6Clo.blood.monocytes", "Ly6Chi.blood.monocytes",
               "Dendritic")] <- "Monocytes/DC"
  out[x %in% c("Fibroblastic.reticular", "Pericyte")] <- "Stroma"
  out
}
keep[, cell_class := classify(main_type)]
keep <- keep[!is.na(cell_class)]
cat("\nCells per class x peak/valley:\n")
print(table(keep$cell_class, keep$class))

# Load counts for these cells
counts_sub <- counts_df[cell_id %in% keep$cell_id]
setorderv(keep, "cell_id"); setorderv(counts_sub, "cell_id")
stopifnot(all(keep$cell_id == counts_sub$cell_id))

gene_cols <- setdiff(colnames(counts_sub), c("Slide", "fov", "cell_id"))
mat <- t(as.matrix(counts_sub[, ..gene_cols]))
colnames(mat) <- counts_sub$cell_id; rownames(mat) <- gene_cols
md <- as.data.frame(keep); rownames(md) <- md$cell_id

obj <- CreateSeuratObject(counts = mat, meta.data = md)
obj <- NormalizeData(obj, verbose = FALSE)

# ---- Per-class DGE ----
classes <- levels(factor(keep$cell_class))
all_res <- list()
for (cls in classes) {
  cat(sprintf("\n=== %s ===\n", cls))
  sub <- subset(obj, subset = cell_class == cls)
  np <- sum(sub$class == "peak"); nv <- sum(sub$class == "valley")
  cat(sprintf("peak n = %d, valley n = %d\n", np, nv))
  if (np < 30 || nv < 30) { cat("Skipped (small n)\n"); next }
  Idents(sub) <- sub$class
  m <- FindMarkers(sub, ident.1 = "peak", ident.2 = "valley",
                   test.use = "wilcox", logfc.threshold = 0.1,
                   min.pct = 0.05, verbose = FALSE)
  m$gene <- rownames(m); m$cell_class <- cls
  m <- as.data.table(m)[, .(cell_class, gene, avg_log2FC, pct.1, pct.2,
                            p_val, p_val_adj)]
  setorder(m, -avg_log2FC)
  cat("Top 10 peak-UP:\n"); print(head(m, 10))
  cat("Top 10 valley-UP:\n"); print(head(m[order(avg_log2FC)], 10))
  all_res[[cls]] <- m
}
combined <- rbindlist(all_res)
fwrite(combined, "/tmp/dge_pv_by_celltype_fixed.tsv", sep = "\t")
cat("\nSaved /tmp/dge_pv_by_celltype_fixed.tsv\n")
