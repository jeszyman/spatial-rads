# Programme B redo: rebuild cGAS-STING as a panel-available signaling-OUTPUT module
# (no sensors on panel; this is IRF3 + NF-kB arms + STING-induced chemokines).
# Type I IFN reused (existing 9-ISG UCell). MBRT peak/valley vs SBRT, 4h, tumor cells.
suppressPackageStartupMessages({library(Seurat); library(UCell); library(data.table); library(arrow); library(ggplot2)})

sting_v2 <- c("Irf3","Nfkb1","Rela","Nfkbia","Cxcl9","Cxcl10","Ccl5","Tnf","Il6")
scored <- "/mnt/data/projects/spatial-rads/processing/scored"
fl <- as.data.table(read_parquet("results/aggregate/full_labels.parquet", col_select=c("cell","compartment")))

get_tumor <- function(s) {
  o <- readRDS(file.path(scored, paste0(s, ".scored.rds")))
  o <- AddModuleScore_UCell(o, features = list(cGAS_STING = sting_v2), name = "")
  md <- o@meta.data
  d <- data.table(cell_id = as.character(md$cell_id),
                  IFN = md$TypeI_interferon_UCell, STING = md$cGAS_STING)
  rm(o); gc()
  comp <- fl[startsWith(cell, paste0(s, "_"))][, cell_id := sub(paste0("^", s, "_"), "", cell)]
  merge(d, comp[, .(cell_id, compartment)], by = "cell_id")[compartment == "tumor"]
}

mbrt <- get_tumor("sam0003")
zone <- fread("dev/peak_valley_analysis/data/mbrt4h_peak_valley.tsv", select = c("cell_id","zone"))
mbrt <- merge(mbrt, zone, by = "cell_id")
mbrt[, group := ifelse(zone == "peak", "MBRT\npeak", "MBRT\nvalley")]
sbrt <- get_tumor("sam0006")[, group := "SBRT"]

d <- rbind(mbrt[, .(group, IFN, STING)], sbrt[, .(group, IFN, STING)])
d[, group := factor(group, levels = c("MBRT\npeak","MBRT\nvalley","SBRT"))]
long <- melt(d, id.vars="group", measure.vars=c("IFN","STING"), variable.name="pathway", value.name="score")
long[, pathway := factor(pathway, c("IFN","STING"),
   labels = c("Type I interferon (9 ISGs)", "cGAS-STING output (IRF3+NF-kB+chemokines, 9 genes)"))]

summ <- long[, .(n=.N, mean=round(mean(score),4), median=round(median(score),4)), by=.(pathway, group)]
print(summ); fwrite(summ, "dev/peak_valley_analysis/data/maldi_programB_pv_sbrt.tsv", sep="\t")

cols <- c("MBRT\npeak"="#B2182B","MBRT\nvalley"="#2166AC","SBRT"="grey50")
p <- ggplot(long, aes(group, score, fill=group)) +
  geom_violin(scale="width", trim=TRUE, width=0.85, alpha=0.85, draw_quantiles=c(.25,.5,.75), linewidth=0.3) +
  stat_summary(fun=mean, geom="point", shape=23, size=2.6, fill="white", colour="black") +
  facet_wrap(~pathway, nrow=1, scales="free_y") +
  scale_fill_manual(values=cols, guide="none") +
  labs(title="Programme B readouts at 4 h (tumor compartment)",
       subtitle="MBRT peak vs valley vs SBRT  -  n=1 sample/group, descriptive; lines=quartiles, diamond=mean",
       x=NULL, y="UCell score") +
  theme_bw(base_size=12) +
  theme(panel.grid.major.x=element_blank(), strip.background=element_rect(fill="grey92"))

ggsave("dev/peak_valley_analysis/maldi_programB_pv_sbrt.png", p, width=10, height=4.6, dpi=150)
cat("\nsaved (rebuilt STING)\n")
