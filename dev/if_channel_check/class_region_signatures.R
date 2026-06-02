suppressPackageStartupMessages({
  library(arrow); library(data.table); library(ggplot2); library(patchwork)
})

OUT <- "/home/jeszyman/repos/spatial-rads/dev/if_channel_check"

peak_fovs   <- unique(c(224, 209, 201, 186, 172, 225, 229, 215, 175,
                        176, 206, 181, 168, 167))
valley_fovs <- c(199, 213, 227, 188, 204, 218, 178, 165, 170, 231, 189,
                 173, 154)

# ---- Load data ----
meta <- as.data.table(read_parquet("/tmp/mutter01_meta.parquet"))
counts <- as.data.table(read_parquet("/tmp/mutter01_counts.parquet"))
design <- as.data.table(readxl::read_excel("/home/jeszyman/repos/spatial-rads/data/metadata.xlsx"))
meta <- merge(meta, design, by.x = c("Slide", "Block"),
              by.y = c("slide", "block_id"))
meta <- meta[is.na(model) | model == "flank"]

setnames(meta, "ImmuneAtlas_ImmGen_Main_cell_Types", "main_type")
classify <- function(x) {
  out <- rep(NA_character_, length(x))
  out[x == "a"] <- "Tumor cells (a)"
  out[x %in% c("B.cell","Memory.B","Plasma","Plasmablast",
               "GC_centroblasts","GC_centrocyes","Spleen.CD19")] <- "B cells"
  out[x %in% c("CD8.T.cell","Spleen.Naive.CD4","Spleen.Naive.CD8",
               "Spleen.CD4Act.48hrs","Spleen.LN.Naive.CD4",
               "Spleen.Treg","Colon.Treg.Nrplo")] <- "T cells"
  out[x %in% c("Macrophage","PerC.macrophage","spleen.red.pulp.macs",
               "small.peritoneal.macs","Microglia")] <- "Macrophages"
  out[x %in% c("Ly6Clo.blood.monocytes","Ly6Chi.blood.monocytes",
               "Dendritic")] <- "Monocytes/DC"
  out[x %in% c("Fibroblastic.reticular","Pericyte")] <- "Stroma"
  out
}
meta[, cell_class := classify(main_type)]
meta <- meta[!is.na(cell_class)]

# ---- Load DGE to extract class × region signatures ----
dge <- fread("/tmp/dge_pv_by_celltype_fixed.tsv")
# For each (class, region), take top 10 DE genes by effect size
signatures <- list()
for (cls in unique(dge$cell_class)) {
  pk <- dge[cell_class == cls & avg_log2FC > 0.1][order(-avg_log2FC)][1:10, gene]
  pk <- pk[!is.na(pk)]
  vl <- dge[cell_class == cls & avg_log2FC < -0.1][order(avg_log2FC)][1:10, gene]
  vl <- vl[!is.na(vl)]
  signatures[[paste(cls, "peak", sep = "__")]]   <- pk
  signatures[[paste(cls, "valley", sep = "__")]] <- vl
}
cat("Signature sizes (genes):\n")
print(sapply(signatures, length))

# ---- Score each cell ----
# For efficiency, build ONE big matrix of all signature genes
all_sig_genes <- unique(unlist(signatures))
all_sig_genes <- all_sig_genes[all_sig_genes %in% colnames(counts)]
cat("Unique signature genes across all:", length(all_sig_genes), "\n")

counts <- counts[cell_id %in% meta$cell_id]
setkeyv(meta, "cell_id"); setkeyv(counts, "cell_id")
meta <- meta[cell_id %in% counts$cell_id]
setorderv(meta, "cell_id"); setorderv(counts, "cell_id")

gene_cols <- setdiff(colnames(counts), c("Slide", "fov", "cell_id"))
lib <- rowSums(counts[, ..gene_cols]); lib[lib == 0] <- 1
# For each signature gene, compute log-normalized expression per cell
sig_expr <- sapply(all_sig_genes, function(g)
  log2((counts[[g]] / lib) * 1e4 + 1))  # cells x genes
rownames(sig_expr) <- counts$cell_id

# Score = mean signature-gene expression per cell
for (sig_name in names(signatures)) {
  genes_in_sig <- intersect(signatures[[sig_name]], colnames(sig_expr))
  if (length(genes_in_sig) == 0) next
  meta[, (sig_name) := rowMeans(sig_expr[, genes_in_sig, drop = FALSE])]
}

# ---- Reshape for plotting ----
sig_cols <- names(signatures)
scores <- meta[, c("cell_id", "cell_class", "condition", "treatment",
                    "timepoint_h", "fov", sig_cols), with = FALSE]

# Assign peak/valley labels for MBRT_4h
scores[, pv := NA_character_]
scores[condition == "MBRT_4h" & fov %in% peak_fovs,   pv := "peak"]
scores[condition == "MBRT_4h" & fov %in% valley_fovs, pv := "valley"]

# Melt into long form for each signature
long_list <- list()
for (sig_name in sig_cols) {
  parts <- strsplit(sig_name, "__")[[1]]
  sig_class <- parts[1]; sig_region <- parts[2]
  sub <- scores[cell_class == sig_class,
                .(cell_id, cell_class, condition, treatment, timepoint_h,
                  fov, pv, score = get(sig_name))]
  sub[, signature := sig_region]
  sub[, signature_class := sig_class]
  long_list[[sig_name]] <- sub
}
long <- rbindlist(long_list)

# ---- Summary: mean score per (class, signature, condition) + peak/valley at 4h ----
tp_map <- c("0" = "0 Gy", "1" = "1h", "4" = "4h", "48" = "2d", "144" = "6d")
long[, tp := factor(tp_map[as.character(timepoint_h)],
                    levels = c("0 Gy", "1h", "4h", "2d", "6d"))]
long[, treat := fifelse(is.na(treatment) | treatment == "NT", "Control", treatment)]

summ <- long[, .(n = .N, mean_score = mean(score, na.rm = TRUE)),
             by = .(signature_class, signature, condition, treat, tp)]
setorder(summ, signature_class, signature, treat, tp)

# Peak/Valley overlay at MBRT_4h
pv_summ <- long[condition == "MBRT_4h" & !is.na(pv),
                .(n = .N, mean_score = mean(score, na.rm = TRUE)),
                by = .(signature_class, signature, pv)]
pv_summ[, condition := paste0("MBRT_4h_", pv)]
pv_summ[, treat := paste0("MBRT ", pv)]
pv_summ[, tp := factor("4h", levels = c("0 Gy", "1h", "4h", "2d", "6d"))]
pv_summ[, pv := NULL]
summ <- rbind(summ, pv_summ, fill = TRUE)

# ---- Plot decay curves for each class × signature ----
treat_levels <- c("Control", "MBRT", "SBRT", "MBRT peak", "MBRT valley")
summ[, treat := factor(treat, levels = treat_levels)]
line_s <- summ[treat %in% c("Control", "MBRT", "SBRT")]
pt_s   <- summ[treat %in% c("MBRT peak", "MBRT valley")]

treat_colors <- c("Control" = "gray40", "MBRT" = "#e41a1c",
                  "SBRT" = "#377eb8", "MBRT peak" = "#b30000",
                  "MBRT valley" = "#fc9272")
treat_shapes <- c("Control" = 16, "MBRT" = 17, "SBRT" = 15,
                  "MBRT peak" = 8, "MBRT valley" = 4)

classes_order <- c("Tumor cells (a)", "Macrophages", "Monocytes/DC",
                   "T cells", "B cells", "Stroma")
summ[, signature_class := factor(signature_class, levels = classes_order)]
line_s[, signature_class := factor(signature_class, levels = classes_order)]
pt_s[,   signature_class := factor(signature_class, levels = classes_order)]

p <- ggplot() +
  geom_line(data = line_s,
            aes(x = tp, y = mean_score, color = treat, group = treat),
            linewidth = 0.7) +
  geom_point(data = line_s,
             aes(x = tp, y = mean_score, color = treat, shape = treat),
             size = 2.5) +
  geom_point(data = pt_s,
             aes(x = tp, y = mean_score, color = treat, shape = treat),
             size = 3.5, stroke = 1.2) +
  facet_grid(signature ~ signature_class, scales = "free_y") +
  scale_color_manual(values = treat_colors) +
  scale_shape_manual(values = treat_shapes) +
  labs(x = "Timepoint",
       y = "Signature score (mean log-norm expression of top 10 DE genes)",
       title = "Class × region signatures across time",
       subtitle = "Rows: signature region. Columns: cell class. Each panel scores cells of that class.") +
  theme_minimal() +
  theme(legend.position = "bottom", strip.text = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                                   size = 8),
        panel.spacing.x = unit(0.8, "lines"),
        plot.margin = margin(8, 14, 8, 20),
        axis.title.y = element_text(margin = margin(r = 12)))

ggsave(file.path(OUT, "class_region_signatures_kinetics.png"), p,
       width = 22, height = 8, dpi = 140)
cat("Saved class_region_signatures_kinetics.png\n")

fwrite(summ, "/tmp/class_region_signature_scores.tsv", sep = "\t")
cat("Saved /tmp/class_region_signature_scores.tsv\n")
