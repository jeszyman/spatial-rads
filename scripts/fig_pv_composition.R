#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# fig_pv_composition.R
# No cell-type composition shift between beam peak and valley zones at 4h.
# Paired per-FOV compartment fractions (core bands), peak vs valley, 39 FOVs.
# M01 4h block (sam0003), current QC'd counts + merged-scale atlas labels.
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({library(Seurat);library(arrow);library(dplyr)
  library(readr);library(tidyr);library(ggplot2)})
source("/home/jeszyman/repos/science/R/figure_schema.R")

o <- readRDS("/mnt/data/projects/spatial-rads/processing/norm/sam0003.norm.rds")
md <- o@meta.data; md$cell <- rownames(md)
a <- read_parquet("results/aggregate/full_labels.parquet")
a <- a[grepl("^sam0003_",a$cell),c("cell","compartment")]; a$cell<-sub("^sam0003_","",a$cell)
z <- read_tsv("dev/peak_valley_analysis/data/mbrt4h_peak_valley.tsv",show_col_types=FALSE)[,c("cell_id","dist_to_peak")]
md <- md %>% left_join(a,by="cell") %>% left_join(z,by=c("cell"="cell_id"))
md <- md[!is.na(md$dist_to_peak) & !is.na(md$compartment),]
md$zone <- ifelse(md$dist_to_peak<0.10,"peak",ifelse(md$dist_to_peak>0.40,"valley",NA))
md <- md[!is.na(md$zone),]
keepf <- md %>% count(fov,zone) %>% filter(n>=30) %>% count(fov) %>% filter(n==2) %>% pull(fov)
md <- md[md$fov %in% keepf,]

comp <- md %>% count(fov,zone,compartment) %>% group_by(fov,zone) %>%
  mutate(frac=n/sum(n)) %>% ungroup()
comp$compartment <- factor(comp$compartment, levels=c("tumor","stroma","immune"))
comp$zone <- factor(comp$zone, levels=c("valley","peak"))

## paired t on logit for annotation
lg <- function(p) log((p+.005)/(1-p+.005))
ann <- comp %>% select(fov,zone,compartment,frac) %>%
  pivot_wider(names_from=zone,values_from=frac,values_fill=0) %>%
  group_by(compartment) %>%
  summarise(p=t.test(lg(peak),lg(valley),paired=TRUE)$p.value,
            ymax=max(c(peak,valley)), .groups="drop") %>%
  mutate(lab=sprintf("paired p=%.2f", p))

p <- ggplot(comp, aes(x=zone, y=frac, fill=zone)) +
  geom_line(aes(group=fov), color="grey75", linewidth=0.3, alpha=0.6) +
  geom_point(aes(group=fov), color="grey55", size=0.7, alpha=0.6) +
  geom_boxplot(width=0.5, alpha=0.55, outlier.shape=NA) +
  facet_wrap(~compartment, scales="free_y") +
  geom_text(data=ann, aes(x=1.5, y=ymax*1.05, label=lab), inherit.aes=FALSE, size=4) +
  scale_fill_manual(values=c("valley"="#85C1E9","peak"="#E74C3C"), guide="none") +
  labs(x=NULL, y="fraction of cells in FOV",
       title="No peak/valley cell-type composition shift at 4h (39 paired FOVs)") +
  theme_scifig() + theme(plot.title=element_text(size=12))

save_plot(p, "results/aggregate/plots/pv_composition", w=10, h=4.8)
