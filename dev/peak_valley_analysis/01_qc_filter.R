if (!exists("obj")) source("00_load_data.R")

# --- QC filtering ---
pre_qc <- table(obj$Condition)
obj <- subset(obj, subset = qcFlagsCell == "Pass")
obj <- subset(obj, subset = nCount_RNA > 20 & nFeature_RNA > 10 & propNegative < 0.5)
post_qc <- table(obj$Condition)

qc_summary <- data.frame(
  condition = names(pre_qc),
  pre_filter = as.integer(pre_qc),
  post_filter = as.integer(post_qc[names(pre_qc)]),
  pct_retained = round(as.integer(post_qc[names(pre_qc)]) / as.integer(pre_qc) * 100, 1)
)
print(qc_summary)
write_tsv(qc_summary, file.path(DATA_DIR, "qc_summary.tsv"))
cat(sprintf("Post-QC: %d cells retained (%.0f%%)\n", ncol(obj),
            ncol(obj) / sum(as.integer(pre_qc)) * 100))

# --- QC violin plots ---
p_qc <- VlnPlot(obj, features = c("nCount_RNA", "nFeature_RNA", "propNegative"),
                 group.by = "Condition", pt.size = 0, ncol = 3)
ggsave(file.path(PLOT_DIR, "qc_violins.png"), plot = p_qc, width = 16, height = 5, dpi = 150)

# --- Save QC'd object ---
saveRDS(obj, file.path(ANALYSIS_DIR, "objects", "seurat_qc.rds"))
cat("QC'd Seurat object saved.\n")
