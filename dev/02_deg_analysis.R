library(Seurat)
library(tidyverse)
library(future)

obj <- readRDS("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/objects/seurat_clustered.rds")
cat(sprintf("Loaded: %d genes x %d cells\n", nrow(obj), ncol(obj)))

# --- DEG function (effect-size ranked, NOT p-value based) ---
run_deg <- function(obj, ident.1, ident.2, group.by = "Condition") {
  Idents(obj) <- group.by
  markers <- FindMarkers(obj, ident.1 = ident.1, ident.2 = ident.2,
                         min.pct = 0.1, logfc.threshold = 0.1)
  markers$gene <- rownames(markers)
  markers$comparison <- paste0(ident.1, "_vs_", ident.2)
  markers %>% arrange(desc(abs(avg_log2FC)))
}

# --- MBRT vs Control ---
cat("Running MBRT vs Control DEGs...\n")
mbrt_tp <- c("MBRT_1h", "MBRT_4h", "MBRT_day2", "MBRT_day6")
deg_mbrt_ctrl <- map_dfr(mbrt_tp, function(tp) {
  cat(sprintf("  %s vs Control...\n", tp))
  sub <- subset(obj, cells = colnames(obj)[obj$Condition %in% c(tp, "Control")])
  run_deg(sub, tp, "Control")
})
cat(sprintf("MBRT vs Control: %d DEGs\n", nrow(deg_mbrt_ctrl)))

# --- SBRT vs Control ---
cat("Running SBRT vs Control DEGs...\n")
sbrt_tp <- c("SBRT_4h", "SBRT_day2", "SBRT_day6")
deg_sbrt_ctrl <- map_dfr(sbrt_tp, function(tp) {
  cat(sprintf("  %s vs Control...\n", tp))
  sub <- subset(obj, cells = colnames(obj)[obj$Condition %in% c(tp, "Control")])
  run_deg(sub, tp, "Control")
})
cat(sprintf("SBRT vs Control: %d DEGs\n", nrow(deg_sbrt_ctrl)))

# --- MBRT vs SBRT at matched timepoints ---
cat("Running MBRT vs SBRT DEGs...\n")
matched <- list(c("MBRT_4h", "SBRT_4h"), c("MBRT_day2", "SBRT_day2"), c("MBRT_day6", "SBRT_day6"))
deg_mbrt_sbrt <- map_dfr(matched, function(pair) {
  cat(sprintf("  %s vs %s...\n", pair[1], pair[2]))
  sub <- subset(obj, cells = colnames(obj)[obj$Condition %in% pair])
  run_deg(sub, pair[1], pair[2])
})
cat(sprintf("MBRT vs SBRT: %d DEGs\n", nrow(deg_mbrt_sbrt)))

# --- Combine and save ---
all_degs <- bind_rows(
  mbrt_vs_ctrl = deg_mbrt_ctrl,
  sbrt_vs_ctrl = deg_sbrt_ctrl,
  mbrt_vs_sbrt = deg_mbrt_sbrt,
  .id = "comparison_type"
)
write.csv(all_degs, "/mnt/gcs/jeszyman/projects/spatial-rads/analysis/tables/deg_all_cells.csv", row.names = FALSE)
cat(sprintf("Total DEG rows: %d\n", nrow(all_degs)))

# Top 10 by effect size for each MBRT vs SBRT comparison
cat("\nTop DEGs by effect size (MBRT vs SBRT):\n")
all_degs %>%
  filter(comparison_type == "mbrt_vs_sbrt") %>%
  group_by(comparison) %>%
  slice_max(abs(avg_log2FC), n = 5) %>%
  select(comparison, gene, avg_log2FC, pct.1, pct.2) %>%
  print(n = 30)

# --- Per-cell-type DEGs (MBRT vs SBRT at 4h) ---
cat("\nRunning per-cell-type DEGs (MBRT vs SBRT at 4h)...\n")
ct_counts <- table(obj$cell_type_validated[obj$Condition %in% c("MBRT_4h", "SBRT_4h")])
major_cts <- names(ct_counts[ct_counts > 500])
cat(sprintf("Cell types with >500 cells at 4h: %s\n", paste(major_cts, collapse = ", ")))

deg_by_ct <- map_dfr(major_cts, function(ct) {
  cells <- colnames(obj)[obj$Condition %in% c("MBRT_4h", "SBRT_4h") &
                          obj$cell_type_validated == ct]
  if (length(cells) < 50) return(NULL)
  cat(sprintf("  %s (%d cells)...\n", ct, length(cells)))
  sub <- subset(obj, cells = cells)
  tryCatch(run_deg(sub, "MBRT_4h", "SBRT_4h") %>% mutate(cell_type = ct),
           error = function(e) { cat(sprintf("    Error: %s\n", e$message)); NULL })
})
write.csv(deg_by_ct, "/mnt/gcs/jeszyman/projects/spatial-rads/analysis/tables/deg_by_celltype_4h.csv", row.names = FALSE)
cat(sprintf("Per-cell-type DEGs: %d rows across %d cell types\n",
            nrow(deg_by_ct), length(unique(deg_by_ct$cell_type))))
cat("DEG analysis complete.\n")
