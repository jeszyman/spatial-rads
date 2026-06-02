suppressPackageStartupMessages({
  library(arrow); library(data.table); library(Matrix); library(ggplot2)
  library(patchwork)
})

OUT <- "/home/jeszyman/repos/spatial-rads/dev/if_channel_check"
meta <- as.data.table(read_parquet("/tmp/mutter01_meta.parquet"))
counts_df <- as.data.table(read_parquet("/tmp/mutter01_counts.parquet"))
setnames(meta, "ImmuneAtlas_ImmGen_Main_cell_Types", "main_type")

# Focus on Block_21 "a" cells
m <- meta[Slide == "20250529_214712_S4" & Block == "Block_21" & main_type == "a"]
c_df <- counts_df[cell_id %in% m$cell_id]
setkeyv(m, "cell_id"); setkeyv(c_df, "cell_id")
m <- m[cell_id %in% c_df$cell_id]
setorderv(m, "cell_id"); setorderv(c_df, "cell_id")
cat("'a' cells in Block_21:", nrow(m), "\n")

gene_cols <- setdiff(colnames(c_df), c("Slide","fov","cell_id"))
mat <- t(as.matrix(c_df[, ..gene_cols]))
lib <- colSums(mat); lib[lib == 0] <- 1
normed <- log2(t(t(mat) / lib) * 1e4 + 1)

# Marker panels
markers <- list(
  epithelial  = c("Epcam","Krt8","Krt18","Krt19","Cdh1","Krt5","Krt14","Krt1"),
  endothelial = c("Pecam1","Cdh5","Vwf","Esam","Eng","Ramp2","Kdr","Flt1",
                  "Emcn","Cldn5","Tie1","Tek"),
  fibroblast  = c("Col1a1","Col3a1","Col1a2","Pdgfra","Fap","Dcn","Lum","Postn"),
  sm_muscle   = c("Acta2","Myh11","Tagln","Cnn1","Des"),
  pericyte    = c("Pdgfrb","Rgs5","Cspg4","Mcam","Kcnj8"),
  erythroid   = c("Hba-a1","Hba-a2","Hbb-b1","Hbb-bs"),
  immune_res  = c("Ptprc","Cd68","Cd163","Cd14","Lyz2"),
  tumor_mes   = c("Vim","S100a4","Twist1","Zeb1","Snai1"),  # mes/EMT
  proliferation = c("Mki67","Top2a","Mcm2","Pcna","Birc5","Ube2c")
)

# Keep only markers that exist
markers <- lapply(markers, function(g) intersect(g, rownames(normed)))
cat("Marker panel sizes:\n"); print(sapply(markers, length))

# Score each cell for each marker set (mean log-norm expression)
score_list <- list()
for (nm in names(markers)) {
  if (length(markers[[nm]]) == 0) { score_list[[nm]] <- rep(NA, ncol(normed)); next }
  score_list[[nm]] <- colMeans(normed[markers[[nm]], , drop = FALSE])
}
scores <- as.data.table(score_list)
scores[, cell_id := colnames(normed)]
scores[, x := m$x_slide_mm]
scores[, y := m$y_slide_mm]

# Primary contamination classes: epithelial vs endothelial vs fibroblast
scores[, primary := {
  vals <- cbind(epithelial, endothelial, fibroblast, sm_muscle, erythroid, immune_res)
  cls  <- c("epithelial","endothelial","fibroblast","sm_muscle","erythroid","immune_res")
  apply(vals, 1, function(r) {
    if (all(is.na(r) | r < 0.1)) "unassigned" else cls[which.max(r)]
  })
}]
cat("\nPrimary marker-based class:\n"); print(scores[, .N, by = primary][order(-N)])

# Histogram of endothelial score vs epithelial score
p1 <- ggplot(scores, aes(x = epithelial, y = endothelial, color = primary)) +
  geom_point(size = 0.3, alpha = 0.5) +
  labs(title = "Block_21 'a' cells: epithelial vs endothelial marker score",
       subtitle = sprintf("%d cells. Colored by dominant marker class.", nrow(scores))) +
  theme_minimal() + theme(legend.position = "right")

# Spatial plot colored by primary class
p2 <- ggplot(scores, aes(x = x, y = y, color = primary)) +
  geom_point(size = 0.3, alpha = 0.6) +
  coord_fixed() +
  scale_color_manual(values = c(epithelial = "#e41a1c", endothelial = "#377eb8",
                                fibroblast = "#4daf4a", sm_muscle = "#984ea3",
                                erythroid = "#a65628", immune_res = "#ff7f00",
                                unassigned = "grey70")) +
  labs(title = "Spatial distribution of marker-based classes within 'a'",
       x = "x mm", y = "y mm") +
  theme_minimal() + theme(legend.position = "bottom")

# Top genes in "a" vs non-"a" cells (to see what "a" is enriched for)
other_m <- meta[Slide == "20250529_214712_S4" & Block == "Block_21" & main_type != "a"]
other_c <- counts_df[cell_id %in% other_m$cell_id]
other_mat <- t(as.matrix(other_c[, ..gene_cols]))
other_lib <- colSums(other_mat); other_lib[other_lib == 0] <- 1
other_norm <- log2(t(t(other_mat) / other_lib) * 1e4 + 1)

a_mean <- rowMeans(normed)
other_mean <- rowMeans(other_norm)
de <- data.table(gene = rownames(normed),
                 a_mean = a_mean, other_mean = other_mean,
                 logfc = a_mean - other_mean)
setorder(de, -logfc)
cat("\nTop 25 genes enriched in 'a' vs non-'a':\n")
print(de[1:25])

ggsave(file.path(OUT, "tumor_contamination_audit.png"),
       p1 + p2, width = 16, height = 7, dpi = 130)
cat("Saved tumor_contamination_audit.png\n")
fwrite(de, "/tmp/a_vs_other_dge.tsv", sep = "\t")
fwrite(scores, "/tmp/a_cell_marker_scores.tsv", sep = "\t")
