#!/usr/bin/env Rscript
# Cell-assignment spatial map: cell centroids colored by cell type -- the centroid-based
# analog of the He-2022 CosMx SMI Fig 1c "Cell assignment" panel. Dataset-agnostic: both
# Mutter_01 and Mutter_02 scored objects share x_slide_mm/y_slide_mm coords and a cell_type
# label. Input is one scored .rds or a directory of per-sample scored .rds (metadata bound,
# faceted). Per-panel subsample keeps large cohorts renderable.
# Standalone:
#   conda run -n spatial-rads Rscript scripts/cell_assignment_map.R <input> <out.png> [color_col] [facet_col] [max_cells]
# Defaults: input=processing/scored dir, color=cell_type, facet=auto (sample_id when >1 sample), max_cells=75000/panel.
suppressMessages({library(Seurat); library(dplyr); library(ggplot2)})

args      <- commandArgs(trailingOnly = TRUE)
INPUT     <- if (length(args) >= 1) args[1] else "/mnt/data/projects/spatial-rads/processing/scored"
OUT       <- if (length(args) >= 2) args[2] else "results/processing/plots/cell_assignment_map.png"
COLOR     <- if (length(args) >= 3) args[3] else "cell_type"
FACET     <- if (length(args) >= 4 && nzchar(args[4])) args[4] else NA_character_  # NA = auto
MAX_CELLS <- if (length(args) >= 5) as.integer(args[5]) else 75000L

files <- if (dir.exists(INPUT)) sort(list.files(INPUT, pattern = "\\.rds$", full.names = TRUE)) else INPUT
stopifnot(length(files) > 0)

keep <- c("x_slide_mm", "y_slide_mm", COLOR, "sample_id", "dataset", "condition", "treatment", "model")
md <- bind_rows(lapply(files, function(f) {
  m <- readRDS(f)@meta.data
  if (!"sample_id" %in% names(m)) m$sample_id <- sub("\\.scored\\.rds$|\\.rds$", "", basename(f))
  m[, intersect(keep, names(m)), drop = FALSE]
}))
stopifnot(all(c("x_slide_mm", "y_slide_mm") %in% names(md)), COLOR %in% names(md))

# auto-facet by sample_id only when more than one sample is present
if (is.na(FACET) && dplyr::n_distinct(md$sample_id) > 1) FACET <- "sample_id"
if (!is.na(FACET) && (!FACET %in% names(md) || dplyr::n_distinct(md[[FACET]]) < 2)) FACET <- NA_character_

# per-panel (or global) subsample for render tractability over the ~3.3M-cell cohort
grp <- if (is.na(FACET)) rep("all", nrow(md)) else md[[FACET]]
set.seed(1)
idx <- unlist(lapply(split(seq_len(nrow(md)), grp), function(ii)
  if (length(ii) > MAX_CELLS) sample(ii, MAX_CELLS) else ii), use.names = FALSE)
n_total <- nrow(md); md <- md[idx, , drop = FALSE]

n_lev   <- dplyr::n_distinct(md[[COLOR]])
leg_col <- min(6L, max(1L, ceiling(n_lev / 8)))
p <- ggplot(md, aes(x = x_slide_mm, y = y_slide_mm, color = .data[[COLOR]])) +
  geom_point(size = 0.1, alpha = 0.3) +
  coord_fixed() +
  labs(title = "Cell-assignment spatial map",
       subtitle = sprintf("%s plotted at slide-mm centroids%s (%d of %d cells shown)",
                           COLOR, if (is.na(FACET)) "" else paste0(", faceted by ", FACET),
                           nrow(md), n_total),
       x = "x (mm)", y = "y (mm)", color = COLOR) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.5, "lines")) +
  guides(color = guide_legend(ncol = leg_col, override.aes = list(size = 2, alpha = 1)))
if (!is.na(FACET)) p <- p + facet_wrap(stats::as.formula(paste0("~", FACET)), scales = "free")

dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
n_panel <- if (is.na(FACET)) 1L else dplyr::n_distinct(md[[FACET]])
ncols   <- ceiling(sqrt(n_panel))
width   <- if (n_panel == 1L) 11 else ncols * 5
height  <- (if (n_panel == 1L) 7 else ceiling(n_panel / ncols) * 5) + ceiling(n_lev / leg_col) * 0.22
ggsave(OUT, p, width = width, height = height, dpi = 150, limitsize = FALSE)
cat(sprintf("wrote %s (%d cells, color=%s, facet=%s, %d panel(s))\n",
            OUT, nrow(md), COLOR, ifelse(is.na(FACET), "none", FACET), n_panel))
