## 06_cellchat_spatial.R
## Spatial CellChat at 4h and day2 (MBRT, SBRT, Control). Uses spatial coordinates
## to filter LR interactions to physically proximal cell pairs.

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(CellChat)
  library(patchwork)
})

DATA_DIR <- "data"
PLOT_DIR <- "plots"
OBJECTS_DIR <- "/mnt/data/projects/spatial-rads/analysis/objects/mbrt_vs_sbrt"
tmp_data <- "/mnt/data/projects/spatial-rads/analysis/tmp"
if (!dir.exists(tmp_data)) dir.create(tmp_data, recursive = TRUE)
Sys.setenv(TMPDIR = tmp_data, TMP = tmp_data, TEMP = tmp_data)

obj <- readRDS(file.path(OBJECTS_DIR, "seurat_filtered.rds"))
cat(sprintf("Loaded: %d cells\n", ncol(obj)))

# Apply default necrosis exclusion at day2 (no exclusion at 4h)
keep_main <- !(obj$timepoint_h %in% c(48) & obj$necrosis_zone %in% TRUE)
obj <- subset(obj, cells = colnames(obj)[keep_main])

run_cellchat_spatial <- function(obj_sub, label) {
  cat(sprintf("\n--- CellChat: %s (n=%d cells) ---\n", label, ncol(obj_sub)))
  if (ncol(obj_sub) < 100) return(NULL)
  data.input <- GetAssayData(obj_sub, layer = "data")
  meta <- obj_sub@meta.data
  meta$cell_type_major <- factor(meta$cell_type_major)
  spatial.locs <- as.matrix(meta[, c("x_slide_mm", "y_slide_mm")])
  rownames(spatial.locs) <- colnames(obj_sub)

  # CellChat spatial: requires spatial.factors with conversion factor (mm -> um)
  spatial.factors <- data.frame(ratio = 1000, tol = 5)  # 1 mm = 1000 um, tol = 5 um
  rownames(spatial.factors) <- "single_slide"

  cc <- tryCatch({
    cc <- createCellChat(object = data.input, meta = meta, group.by = "cell_type_major",
                         datatype = "spatial",
                         coordinates = spatial.locs,
                         spatial.factors = spatial.factors)
    cc@DB <- CellChatDB.mouse
    cc <- subsetData(cc)
    cc <- identifyOverExpressedGenes(cc)
    cc <- identifyOverExpressedInteractions(cc)
    cc <- computeCommunProb(cc, type = "truncatedMean", trim = 0.1,
                            distance.use = TRUE, interaction.range = 250,
                            scale.distance = 3,
                            contact.range = 10,
                            contact.knn.k = 6)
    cc <- filterCommunication(cc, min.cells = 10)
    cc <- computeCommunProbPathway(cc)
    cc <- aggregateNet(cc)
    cc
  }, error = function(e) {
    cat(sprintf("  CellChat error: %s\n", conditionMessage(e)))
    NULL
  })
  cc
}

# Run for 4h and day2 conditions
conditions <- list(
  list(label = "Control", cond = c("Control")),
  list(label = "MBRT_4h",  cond = c("MBRT_4h")),
  list(label = "SBRT_4h",  cond = c("SBRT_4h")),
  list(label = "MBRT_day2", cond = c("MBRT_day2")),
  list(label = "SBRT_day2", cond = c("SBRT_day2"))
)

cc_list <- list()
for (cd in conditions) {
  cells_keep <- colnames(obj)[obj$Condition %in% cd$cond]
  if (length(cells_keep) < 100) next
  # Subsample to keep CellChat tractable: 30000 cells per condition
  if (length(cells_keep) > 30000) {
    set.seed(42)
    cells_keep <- sample(cells_keep, 30000)
  }
  obj_sub <- subset(obj, cells = cells_keep)
  cc <- run_cellchat_spatial(obj_sub, cd$label)
  if (!is.null(cc)) cc_list[[cd$label]] <- cc
}

saveRDS(cc_list, file.path(OBJECTS_DIR, "cellchat_list.rds"))
cat(sprintf("\nSaved %d CellChat objects\n", length(cc_list)))

# Compare MBRT vs SBRT at 4h and day2
compare_pair <- function(cc_a, cc_b, label_a, label_b, suffix) {
  if (is.null(cc_a) || is.null(cc_b)) return(NULL)
  cc_merged <- mergeCellChat(list(cc_a, cc_b), add.names = c(label_a, label_b))

  # Pathway scores summary
  pat_a <- cc_a@netP$pathways
  pat_b <- cc_b@netP$pathways
  all_pat <- union(pat_a, pat_b)
  cat(sprintf("\nPathways in %s: %d\n", label_a, length(pat_a)))
  cat(sprintf("Pathways in %s: %d\n", label_b, length(pat_b)))
  cat(sprintf("Union: %d\n", length(all_pat)))

  # Pathway prob for each
  prob_a <- if (length(pat_a) > 0) sapply(pat_a, function(p) sum(cc_a@netP$prob[, , p])) else c()
  prob_b <- if (length(pat_b) > 0) sapply(pat_b, function(p) sum(cc_b@netP$prob[, , p])) else c()
  pathway_df <- tibble(
    pathway = all_pat,
    prob_a  = prob_a[all_pat],
    prob_b  = prob_b[all_pat]) %>%
    mutate(prob_a = replace_na(prob_a, 0),
           prob_b = replace_na(prob_b, 0),
           delta = prob_a - prob_b,
           ident_a = label_a, ident_b = label_b)
  write_tsv(pathway_df, file.path(DATA_DIR, sprintf("cellchat_pathway_%s.tsv", suffix)))

  # Plot top differential pathways
  top_diff <- pathway_df %>%
    arrange(desc(abs(delta))) %>%
    slice_head(n = 15) %>%
    mutate(pathway = factor(pathway, levels = pathway))
  p <- top_diff %>%
    pivot_longer(c(prob_a, prob_b), names_to = "ident", values_to = "prob") %>%
    mutate(ident = ifelse(ident == "prob_a", label_a, label_b)) %>%
    ggplot(aes(x = pathway, y = prob, fill = ident)) +
    geom_col(position = "dodge") +
    coord_flip() +
    scale_fill_manual(values = c("MBRT_4h"="#E74C3C","SBRT_4h"="#3498DB",
                                  "MBRT_day2"="#E74C3C","SBRT_day2"="#3498DB",
                                  "Control"="#7F8C8D")) +
    labs(title = sprintf("Top differential pathways: %s vs %s", label_a, label_b),
         y = "Sum probability", x = NULL) +
    theme_bw()
  ggsave(file.path(PLOT_DIR, sprintf("cellchat_top_pathways_%s.png", suffix)),
         plot = p, width = 9, height = 6, dpi = 130)

  cc_merged
}

if (!is.null(cc_list[["MBRT_4h"]]) && !is.null(cc_list[["SBRT_4h"]]))
  compare_pair(cc_list[["MBRT_4h"]], cc_list[["SBRT_4h"]], "MBRT_4h", "SBRT_4h", "mbrt_vs_sbrt_4h")
if (!is.null(cc_list[["MBRT_day2"]]) && !is.null(cc_list[["SBRT_day2"]]))
  compare_pair(cc_list[["MBRT_day2"]], cc_list[["SBRT_day2"]], "MBRT_day2", "SBRT_day2", "mbrt_vs_sbrt_day2")

cat("\n=== 06_cellchat_spatial.R complete ===\n")
