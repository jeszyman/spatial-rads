#!/usr/bin/env Rscript
# Batch-integration QC: iLISI before vs after Harmony, for the dataset and slide batch axes.
# Higher iLISI = better mixing (floor 1 = batches fully separated; ceiling = number of batches).
# Reads the long-format embed metrics table written by scripts/aggregate/embed_celltype.R.
suppressMessages({library(data.table); library(ggplot2)})
source("/home/jeszyman/repos/science/R/figure_schema.R")

AGG     <- "results/aggregate"
METRICS <- file.path(AGG, "celltype_embed_metrics.tsv")
OUT     <- file.path(AGG, "plots", "lisi_integration")
stopifnot(file.exists(METRICS))

FLOOR <- 1   # iLISI floor: 1 = batches fully separated (no mixing). Stems start here.

m <- fread(METRICS)
d <- m[grepl("^lisi_(dataset|slide)_(pre|post)$", metric)]
d[, axis  := fifelse(grepl("dataset", metric), "Dataset (2 batches)", "Slide")]
d[, stage := factor(fifelse(grepl("pre$", metric), "Pre-Harmony", "Post-Harmony"),
                    levels = c("Pre-Harmony", "Post-Harmony"))]
gate    <- m[metric == "lisi_min_gate", value]
gate_df <- data.table(axis = "Dataset (2 batches)", y = gate)

key <- paste(strwrap(paste(
  "Batch mixing before vs after Harmony integration, by batch axis.",
  "iLISI = mean local inverse-Simpson index; higher = better mixed.",
  "Stems begin at the floor iLISI = 1 (batches fully separated, no mixing);",
  "stem length is mixing gained above that floor (ceiling = number of batches).",
  sprintf("Dashed line: integration gate, LISI(dataset) >= %.2f.", gate),
  "Mutter_01 + Mutter_02, CosMx; 20k-cell subsample."), width = 118), collapse = "\n")

p <- ggplot(d, aes(stage, value, colour = stage)) +
  geom_hline(yintercept = FLOOR, colour = "grey75", linewidth = 0.4) +
  geom_segment(aes(xend = stage, y = FLOOR, yend = value), linewidth = 1.1) +
  geom_point(size = 4.2) +
  geom_text(aes(label = sprintf("%.2f", value)), vjust = -1.0, size = 3.2, colour = "black") +
  geom_hline(data = gate_df, aes(yintercept = y), linetype = "dashed", colour = "grey40") +
  facet_wrap(~axis, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.14))) +
  labs(x = NULL, y = "iLISI (higher = better mixed)", caption = key) +
  theme_scifig(base_size = 11) +
  theme(legend.position = "none",
        plot.caption = element_text(size = 9, hjust = 0, colour = "gray30", margin = margin(t = 12)),
        plot.caption.position = "plot")

save_plot(p, OUT, w = 7, h = 4.5)
cat(sprintf("LISI(dataset): pre %.2f -> post %.2f (gate %.2f)\n",
            d[grepl("dataset_pre", metric), value],
            d[grepl("dataset_post", metric), value], gate))
