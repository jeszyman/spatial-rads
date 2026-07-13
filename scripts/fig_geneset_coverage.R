#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# fig_geneset_coverage.R
# Pathway-set panel coverage by provenance: for each study concept, the fraction of
# its genes present on the common 950-gene panel, shown per source (NanoString module,
# MSigDB Hallmark, custom). Concepts with both computational sources get two rows so
# the panel-designed module and the genome-wide Hallmark set are compared directly.
# The confirmatory set per concept (highest coverage) is marked; the rest are scored
# exploratory. Interpretation aid: a null in a low-coverage set is a panel blind spot,
# not necessarily biology. Reads the coverage table from build_gene_sets.R.
# Args (all optional, canonical defaults): <gene_set_panel_coverage.tsv> <out.png>
# -----------------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))
source("/home/jeszyman/repos/science/R/figure_schema.R")

a   <- commandArgs(trailingOnly = TRUE)
cov <- if (length(a) >= 1) a[1] else "results/data_model/gene_set_panel_coverage.tsv"
out <- if (length(a) >= 2) sub("\\.png$", "", a[2]) else "results/data_model/plots/geneset_coverage"

prov_lab <- c(nanostring = "NanoString module", hallmark = "MSigDB Hallmark", custom = "Custom (cited)")
prov_col <- c("NanoString module" = "#3c78d8", "MSigDB Hallmark" = "#e6a817", "Custom (cited)" = "#c0504d")

# concept id -> display label (fix the run-together TypeI/TypeII, strip underscores)
nice_concept <- function(x) x %>%
  str_replace("^TypeI_", "Type I ") %>% str_replace("^TypeII_", "Type II ") %>%
  str_replace_all("_", " ")

d <- read_tsv(cov, show_col_types = FALSE, comment = "#") %>%   # drop msigdbr_version line
  mutate(provenance  = factor(prov_lab[source], levels = prov_lab),
         concept_lab = nice_concept(concept)) %>%
  # order rows: group by concept (confirmatory coverage descending), source within concept
  group_by(concept) %>% mutate(concept_cov = max(coverage)) %>% ungroup() %>%
  arrange(concept_cov, desc(source == "nanostring")) %>%
  mutate(row_lab = paste0(concept_lab, "  [", recode(source, nanostring = "NS", hallmark = "HM", custom = "cu"), "]"),
         band    = as.integer(factor(concept, levels = unique(concept))) %% 2)  # alternating concept shade
d$row_lab <- factor(d$row_lab, levels = d$row_lab)
d$ri <- as.integer(d$row_lab)
# full-width shading spanning every other concept's row pair (baked alpha; no scale to wash it out)
bands <- d %>% filter(band == 1) %>%
  summarise(ymin = min(ri) - 0.5, ymax = max(ri) + 0.5, .by = concept)

legend_text <- str_c(
  "Every confirmatory pathway signature is a panel-designed NanoString module with near-complete ",
  "coverage (83-97%); the genome-wide MSigDB Hallmark sets, shown for comparison, are diluted off ",
  "the 950-gene panel and stay exploratory. A null in a low-coverage set is a panel blind spot, ",
  "not evidence of absent biology. Position gives the on-panel gene fraction per concept and source; dot ",
  "size is the number of on-panel genes (small = statistically thin even at high coverage, e.g. ",
  "STING 2/3). Ring marks the confirmatory set (highest coverage) per concept; STING has no ",
  "computational source and is a cited custom set.")

p <- ggplot(d, aes(coverage, row_lab)) +
  geom_rect(data = bands, inherit.aes = FALSE,
            aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
            fill = "grey90", alpha = 0.6) +
  geom_segment(aes(x = 0, xend = coverage, yend = row_lab), color = "grey75") +
  geom_point(aes(fill = provenance, size = n_panel), shape = 21, stroke = 0.4, color = "white") +
  geom_point(data = filter(d, confirmatory), aes(size = n_panel),
             shape = 21, stroke = 1.6, color = "grey20", fill = NA) +
  geom_text(aes(label = sprintf("%d/%d", n_panel, n_total)), hjust = -0.5, size = 3.2) +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1.1), breaks = seq(0, 1, 0.25)) +
  scale_size_area(max_size = 8, limits = c(0, max(d$n_panel)),
                  breaks = c(5, 50, 100), name = "On-panel genes (n)") +
  scale_fill_manual(values = prov_col, name = NULL, drop = FALSE) +
  guides(fill = guide_legend(override.aes = list(size = 4.5), order = 1),
         size = guide_legend(order = 2)) +
  labs(x = "Fraction of Set Genes on the Common 950-Gene Panel", y = NULL,
       caption = str_wrap(legend_text, width = round(9 * 15))) +
  theme_scifig(base_size = 12) +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.box.just = "left",
        legend.key.size = unit(0.4, "cm"),
        axis.title.x = element_text(margin = margin(t = 8)),
        plot.caption = element_text(size = 10, hjust = 0, color = "gray30",
                                    lineheight = 1.3, margin = margin(t = 20)),
        plot.caption.position = "plot")

save_plot(p, out, w = 9, h = 6.5)
