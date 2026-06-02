suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggrepel); library(patchwork)
})

OUT <- "/home/jeszyman/repos/spatial-rads/dev/if_channel_check"
d <- fread("/tmp/tumor_pv_mbrt_vs_sbrt_null.tsv")
d <- d[pct_expr_mbrt > 0.05 & pct_expr_sbrt > 0.05]
sd_null <- sd(d$peak_minus_valley_sbrt)
mbrt_mean <- mean(d$peak_minus_valley_mbrt)
cat(sprintf("SBRT sd = %.3f | MBRT mean shift = %+.3f | n_genes = %d\n",
            sd_null, mbrt_mean, nrow(d)))

# Label genes: top peak-specific, valley-specific
top_peak <- d[order(-delta)][1:12, gene]
top_valley <- d[order(delta)][1:5, gene]
d[, label := ifelse(gene %in% c(top_peak, top_valley), gene, NA_character_)]
d[, category := "most genes (no strong signal)"]
d[abs(peak_minus_valley_mbrt) > 0.08 &
  abs(delta) > 2 * sd_null & peak_minus_valley_mbrt > 0,
  category := "MBRT-specific peak enrichment"]
d[peak_minus_valley_mbrt < -0.08 & delta < -2 * sd_null,
  category := "MBRT-specific valley enrichment"]

# Main scatter
xlim <- range(d$peak_minus_valley_sbrt); ylim <- range(d$peak_minus_valley_mbrt)

# Layers for annotation:
# a: horizontal band of SBRT null = ±2*sd_null
# b: identity y=x line (shared architecture)
# c: horizontal line at MBRT global mean shift
main <- ggplot(d, aes(x = peak_minus_valley_sbrt,
                      y = peak_minus_valley_mbrt)) +
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = -2 * sd_null, ymax = 2 * sd_null,
           fill = "grey90", alpha = 0.6) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_hline(yintercept = mbrt_mean, linetype = "dotted",
             color = "darkorange", linewidth = 0.5) +
  geom_point(aes(color = category), alpha = 0.7, size = 0.9) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 50,
                  color = "black", min.segment.length = 0,
                  box.padding = 0.3, seed = 1) +
  scale_color_manual(values = c(
    "MBRT-specific peak enrichment" = "#e41a1c",
    "MBRT-specific valley enrichment" = "#377eb8",
    "most genes (no strong signal)" = "grey60")) +
  annotate("text", x = max(xlim), y = mbrt_mean + 0.008,
           label = sprintf("MBRT global peak shift (+%.3f)", mbrt_mean),
           hjust = 1, color = "darkorange", size = 3.3) +
  annotate("text", x = max(xlim) * 0.85, y = 2 * sd_null + 0.008,
           label = sprintf("SBRT null: ±2 sd (±%.3f)", 2 * sd_null),
           hjust = 1, color = "grey40", size = 3.3) +
  labs(x = sprintf("Peak − Valley under SBRT (uniform dose → this is NULL noise, sd=%.3f)", sd_null),
       y = sprintf("Peak − Valley under MBRT (real peaks, sd=%.3f)", sd(d$peak_minus_valley_mbrt)),
       color = NULL,
       title = "Each dot = one gene (tumor cells, 4h). MBRT has real peaks, SBRT does not.",
       subtitle = paste0(
         "HOW TO READ: SBRT is our null. Under uniform dose, peak-valley differences should be near zero — and they are (narrow grey band). ",
         "Under MBRT, many genes jump UP far above the grey band (red dots). Those are real MBRT peak-enriched genes. ",
         "Only the vertical direction (y) matters for signal; horizontal (x) just tells us the null is well-behaved.")) +
  theme_minimal() +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 11),
        plot.subtitle = element_text(size = 9, lineheight = 1.2))

# Right panel: marginal density of MBRT vs SBRT
dens <- rbind(
  data.table(arm = "SBRT (null)",
             val = d$peak_minus_valley_sbrt),
  data.table(arm = "MBRT (real stripes)",
             val = d$peak_minus_valley_mbrt))
dens[, arm := factor(arm, levels = c("SBRT (null)", "MBRT (real stripes)"))]
dplot <- ggplot(dens, aes(x = val, fill = arm, color = arm)) +
  geom_density(alpha = 0.35, linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  scale_fill_manual(values = c("SBRT (null)" = "grey50",
                               "MBRT (real stripes)" = "#e41a1c")) +
  scale_color_manual(values = c("SBRT (null)" = "grey30",
                                "MBRT (real stripes)" = "#b30000")) +
  labs(x = "Peak − Valley (log-norm expression)",
       y = "Density of genes",
       title = "Distribution of per-gene peak−valley across all panel genes",
       subtitle = "SBRT is narrow around 0 (no real stripes). MBRT is wider and shifted right (real peaks exist; some genes strongly up).",
       fill = NULL, color = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 11),
        plot.subtitle = element_text(size = 9, lineheight = 1.2))

combined <- main / dplot + plot_layout(heights = c(2.5, 1))
ggsave(file.path(OUT, "tumor_pv_mbrt_vs_sbrt_annotated.png"),
       combined, width = 12, height = 13, dpi = 130)
cat("Saved tumor_pv_mbrt_vs_sbrt_annotated.png\n")
