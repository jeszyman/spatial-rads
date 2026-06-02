suppressPackageStartupMessages({
  library(arrow); library(data.table); library(Seurat); library(dplyr)
})

peak_fovs   <- unique(c(224, 209, 201, 186, 172, 225, 229, 215, 175, 176, 206, 181, 168, 167))
valley_fovs <- c(199, 213, 227, 188, 204, 218, 178, 165, 170, 231, 189, 173, 154)

meta <- as.data.table(read_parquet("/tmp/mutter01_meta.parquet"))
s4 <- meta[Slide == "20250529_214712_S4" & Block == "Block_21"]
s4[, class := fifelse(fov %in% peak_fovs, "peak",
                      fifelse(fov %in% valley_fovs, "valley", NA_character_))]
keep <- s4[!is.na(class)]
cat("Peak cells:", sum(keep$class == "peak"),
    "| Valley cells:", sum(keep$class == "valley"), "\n")

counts_df <- read_parquet("/tmp/mutter01_counts.parquet")
counts_dt <- as.data.table(counts_df)
counts_dt <- counts_dt[cell_id %in% keep$cell_id]
setkey(counts_dt, cell_id)
keep <- keep[cell_id %in% counts_dt$cell_id]
setkey(keep, cell_id)

gene_cols <- setdiff(colnames(counts_dt), c("Slide", "fov", "cell_id"))
mat <- t(as.matrix(counts_dt[, ..gene_cols]))
colnames(mat) <- counts_dt$cell_id
rownames(mat) <- gene_cols

md <- as.data.frame(keep[match(colnames(mat), cell_id)])
rownames(md) <- md$cell_id

obj <- CreateSeuratObject(counts = mat, meta.data = md)
obj <- NormalizeData(obj, verbose = FALSE)
Idents(obj) <- obj$class

cat("\n=== Bulk DGE: peak vs valley ===\n")
markers <- FindMarkers(obj, ident.1 = "peak", ident.2 = "valley",
                       test.use = "wilcox", logfc.threshold = 0.1,
                       min.pct = 0.05, verbose = FALSE)
markers$gene <- rownames(markers)
markers <- as.data.table(markers)[, .(gene, avg_log2FC, pct.1, pct.2, p_val, p_val_adj)]
setorder(markers, -avg_log2FC)

cat("\nTop 25 UP in peak:\n")
print(head(markers, 25))
cat("\nTop 25 UP in valley:\n")
print(head(markers[order(avg_log2FC)], 25))

fwrite(markers, "/tmp/dge_peak_vs_valley_bulk.tsv", sep = "\t")

cat("\n=== Stratified DGE by cell family ===\n")
fam_col <- "ImmuneAtlas_ImmGen_cellFamily_Main_cell_Types"
obj$fam <- obj@meta.data[[fam_col]]
fam_counts <- table(obj$fam, obj$class)
print(fam_counts)

# Only run DGE for families with >=50 cells per class
ok_fams <- names(which(rowSums(fam_counts >= 50) == 2))
cat("\nFamilies with >=50 cells in both classes:", paste(ok_fams, collapse = ", "), "\n")

strat_results <- lapply(ok_fams, function(fam) {
  sub <- subset(obj, subset = fam == !!fam)
  Idents(sub) <- sub$class
  m <- FindMarkers(sub, ident.1 = "peak", ident.2 = "valley",
                   test.use = "wilcox", logfc.threshold = 0.1,
                   min.pct = 0.05, verbose = FALSE)
  m$gene <- rownames(m); m$family <- fam
  as.data.table(m)[, .(family, gene, avg_log2FC, pct.1, pct.2, p_val, p_val_adj)]
})
strat <- rbindlist(strat_results)
setorder(strat, family, -avg_log2FC)
fwrite(strat, "/tmp/dge_peak_vs_valley_by_family.tsv", sep = "\t")

cat("\nTop peak-UP per family:\n")
for (fam in ok_fams) {
  cat("---", fam, "---\n")
  print(head(strat[family == fam], 10))
}
