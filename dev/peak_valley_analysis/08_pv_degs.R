if (!exists("PATHWAY_COLS")) source("00_load_data.R")

obj <- readRDS(file.path(ANALYSIS_DIR, "objects", "seurat_clustered.rds"))
pv <- read_tsv(file.path(DATA_DIR, "mbrt4h_peak_valley.tsv"), show_col_types = FALSE)

# Subset to MBRT_4h with zone annotation
mbrt4h_cells <- pv$cell_id
obj_4h <- subset(obj, cells = mbrt4h_cells)
obj_4h$zone <- pv$zone[match(colnames(obj_4h), pv$cell_id)]
Idents(obj_4h) <- "zone"
cat(sprintf("MBRT_4h with zones: %d cells (peak=%d, valley=%d)\n",
            ncol(obj_4h), sum(obj_4h$zone == "peak"), sum(obj_4h$zone == "valley")))

# ---- Bulk DEGs: peak vs valley ----
cat("Running peak vs valley DEGs (all cells)...\n")
bulk_degs <- FindMarkers(obj_4h, ident.1 = "peak", ident.2 = "valley",
                         min.pct = 0.1, logfc.threshold = 0.05)
bulk_degs$gene <- rownames(bulk_degs)
bulk_degs <- bulk_degs %>% arrange(desc(abs(avg_log2FC)))
write_tsv(bulk_degs, file.path(DATA_DIR, "pv_degs_bulk.tsv"))
cat(sprintf("Bulk peak vs valley DEGs: %d\n", nrow(bulk_degs)))

cat("\nTop 20 by effect size:\n")
print(head(bulk_degs %>% select(gene, avg_log2FC, pct.1, pct.2, p_val_adj), 20))

# Volcano retired: -log10(p_val) pseudoreplicates the n=1 design. Mutter_02 is 2d-only and
# H2AX peak/valley ground truth is 4h-only, so NO biological replicates are coming for this
# 4h contrast -- the n=1 spatial-descriptive views ARE final: across-stripe profiles +
# spatial pathway maps (05), zone maps (07), gradient + dot/box plots. bulk_degs (written
# above) is the descriptive effect-size table.

# ---- Per-cell-type DEGs ----
cat("\nRunning per-cell-type peak vs valley DEGs...\n")
ct_counts <- table(obj_4h$cell_type_validated, obj_4h$zone)
major_cts <- rownames(ct_counts)[apply(ct_counts, 1, min) >= 250]
cat(sprintf("Cell types with >=250 cells in both zones: %s\n", paste(major_cts, collapse = ", ")))

ct_degs <- map_dfr(major_cts, function(ct) {
  cells <- colnames(obj_4h)[obj_4h$cell_type_validated == ct]
  sub <- subset(obj_4h, cells = cells)
  Idents(sub) <- "zone"
  tryCatch({
    degs <- FindMarkers(sub, ident.1 = "peak", ident.2 = "valley",
                        min.pct = 0.1, logfc.threshold = 0.05)
    degs$gene <- rownames(degs)
    degs$cell_type <- ct
    degs %>% arrange(desc(abs(avg_log2FC)))
  }, error = function(e) { cat(sprintf("  %s: error — %s\n", ct, e$message)); NULL })
})
write_tsv(ct_degs, file.path(DATA_DIR, "pv_degs_by_celltype.tsv"))
cat(sprintf("Per-cell-type DEGs: %d rows across %d types\n",
            nrow(ct_degs), length(unique(ct_degs$cell_type))))

cat("Peak vs valley DEG analysis complete.\n")
