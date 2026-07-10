#!/usr/bin/env Rscript
# A priori program panel -- effect size + 95% CI per curated program x cell type x
# contrast. Effect-size-FORWARD: nothing is gated by padj/MDE; points whose |effect|
# is below the magnitude floor (0.5x MDE) are faded, not dropped. Args: <results_master.tsv> <out.png>
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
a <- commandArgs(trailingOnly = TRUE); m <- fread(a[1]); out <- a[2]
PRIMARY <- c("TypeI_interferon","TypeII_interferon","STING","DNA_Damage_Repair",
             "Angiogenesis","Hypoxia","Fibrosis_remodeling","Stromal_stress_senescence")
d <- m[readout_class == "pathway" & feature %in% PRIMARY & !is.na(effect)]
d[, feature := factor(feature, levels = PRIMARY)]
cov <- d[, .(cov = unique(panel_cov_frac)[1]), by = feature]
d <- merge(d, cov, by = "feature")
d[, flab := sprintf("%s  (%.0f%% panel)", feature, 100 * cov)]
d[, flab := factor(flab, levels = unique(d$flab[order(d$feature)]))]
d[, fade := ifelse(trend_call == "below-floor", "below", "ok")]
p <- ggplot(d, aes(effect, unit, color = contrast, alpha = fade)) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_pointrange(aes(xmin = ci_low, xmax = ci_high),
                  position = position_dodge(width = 0.6), size = 0.3, fatten = 1.4) +
  facet_wrap(~ flab, scales = "free", ncol = 2) +
  scale_alpha_manual(values = c(ok = 1, below = 0.22), guide = "none") +
  scale_color_manual(values = c(MBRT_vs_Ctrl = "#1f77b4", SBRT_vs_Ctrl = "#d62728",
                                MBRT_vs_SBRT = "#2ca02c"), name = NULL) +
  labs(x = "program-score effect (limma estimate) +/- 95% CI", y = NULL,
       title = "A priori program activity by cell type -- all three contrasts (effect-size, not padj-gated)",
       subtitle = "faded = |effect| below the magnitude floor (0.5x n=4 MDE); panel-coverage shown per program") +
  theme_bw(base_size = 9) +
  theme(legend.position = "top", strip.text = element_text(size = 7.3),
        panel.grid.minor = element_blank())
ggsave(out, p, width = 11, height = 9.5, dpi = 170)
cat("wrote", out, "|", nrow(d), "program rows\n")
