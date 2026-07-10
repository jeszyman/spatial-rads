#!/usr/bin/env Rscript
# How much each arm moves the transcriptome at day 2 -- one clean bar per contrast:
# the count of cell-type DE genes whose |log2FC| clears a modest effect-size threshold
# (NO p-values). Args: <results_master.tsv> <out.png>
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
a <- commandArgs(trailingOnly = TRUE); m <- fread(a[1]); out <- a[2]
IMMUNE <- c("T cells","NK cells","ILC","Plasma cells","Macrophages","DC","Mast cells","Neutrophils")
STROMA <- c("Fibroblast","SmoothMuscle","Adipocyte","Endothelial")
TH <- 0.25
d <- m[readout_class == "DE" & !is.na(effect)]
d[, comp := fifelse(unit %in% IMMUNE, "immune",
            fifelse(unit %in% STROMA, "stroma",
             fifelse(unit == "Tumor", "tumor", "other")))]
d <- d[comp != "other"]
n <- d[abs(effect) >= TH, .N, by = contrast]
n[, contrast := factor(contrast, levels = c("MBRT_vs_Ctrl","MBRT_vs_SBRT","SBRT_vs_Ctrl"))]
print(n)
p <- ggplot(n, aes(contrast, N, fill = contrast)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = N), vjust = -0.4, size = 5) +
  scale_fill_manual(values = c(MBRT_vs_Ctrl = "#1f77b4", MBRT_vs_SBRT = "#2ca02c",
                               SBRT_vs_Ctrl = "#d62728"), guide = "none") +
  expand_limits(y = max(n$N) * 1.1) +
  labs(y = "# cell-type DE genes moved  (|log2FC| >= 0.25)", x = NULL,
       title = "How much each arm moves the transcriptome (day 2)",
       subtitle = "effect-size count, no p-values: SBRT moves ~19x more genes than MBRT vs Control") +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.x = element_blank())
ggsave(out, p, width = 6.5, height = 4.6, dpi = 175)
cat("wrote", out, "\n")
