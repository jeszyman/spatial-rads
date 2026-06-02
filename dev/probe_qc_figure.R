#!/usr/bin/env Rscript
# Probe signal-to-background QC: per-gene mean count / mean negative-probe count, Mutter_01 vs
# Mutter_02. Genes at/below background (S2B <= 1) in BOTH datasets are candidate failed probes
# (Ozirmak Lermi 2025, 10.1038/s41467-025-63414-1) -- flagged for review, NOT dropped. Reads the
# per-gene report written by scripts/probe_qc.R.
suppressMessages({library(data.table); library(ggplot2); library(ggrepel)})
source("/home/jeszyman/repos/science/R/figure_schema.R")

REP <- "results/processing/probe_qc_report.tsv"
OUT <- "results/aggregate/plots/probe_qc_s2b"
stopifnot(file.exists(REP))

d <- fread(REP)
n_flag <- sum(d$flagged_failed); n_gene <- nrow(d)
d[, status := fifelse(flagged_failed, "Flagged: at background in both", "Above background")]

key <- paste(strwrap(paste(
  sprintf("CosMx panel performance: %d of %d genes sit at or below the negative-probe background in both datasets", n_flag, n_gene),
  "(candidate failed probes, flagged for review, not dropped).",
  "Signal-to-background = per-gene mean count divided by mean negative-probe count, per dataset.",
  "Dashed lines mark the background floor (S2B = 1); genes below both are flagged.",
  "Mutter_01 + Mutter_02, CosMx 950-gene shared panel."), width = 96), collapse = "\n")

p <- ggplot(d, aes(m01_s2b, m02_s2b)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey45") +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey45") +
  geom_point(aes(fill = status), shape = 21, size = 2.3, stroke = 0.4, colour = "grey20", alpha = 0.85) +
  geom_text_repel(data = d[flagged_failed == TRUE], aes(label = gene),
                  size = 2.8, colour = "#b3123b", min.segment.length = 0,
                  box.padding = 0.4, max.overlaps = Inf, seed = 1) +
  scale_x_log10(breaks = c(0.5, 1, 3, 10, 30, 100)) +
  scale_y_log10(breaks = c(0.5, 1, 3, 10, 30, 100, 300)) +
  scale_fill_manual(values = c("Above background" = "grey65",
                               "Flagged: at background in both" = "#d1495b")) +
  labs(x = "Signal-to-background  (Mutter_01)", y = "Signal-to-background  (Mutter_02)",
       fill = NULL, caption = key) +
  theme_scifig(base_size = 11) +
  theme(legend.position = "bottom",
        axis.title.x = element_text(margin = margin(t = 8)),
        axis.title.y = element_text(margin = margin(r = 8)),
        plot.caption = element_text(size = 9, hjust = 0, colour = "gray30", margin = margin(t = 12)),
        plot.caption.position = "plot")

save_plot(p, OUT, w = 6.5, h = 6.2)
cat(sprintf("%d genes; %d flagged at-background in both datasets\n", n_gene, n_flag))
