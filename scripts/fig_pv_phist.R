#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# fig_pv_phist.R
# Peak/valley group-mean DE is null (flat p-histogram) in every compartment AND
# all-cells, while the same FOV-pseudobulk machinery detects tumor-vs-stroma
# emphatically (spike at 0). M01 4h block (sam0003), current QC'd + atlas labels.
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({library(Seurat);library(arrow);library(edgeR)
  library(limma);library(dplyr);library(readr);library(Matrix);library(ggplot2)})
source("/home/jeszyman/repos/science/R/figure_schema.R")
DDR <- c("Cdkn1a","Gadd45a","Gadd45b","Mdm2","Bax","Pmaip1","Bbc3","Ddit3","Atm","Atr","Chek1","Chek2","Rad51","Brca1","Brca2","Xrcc5","Xrcc6","Parp1")

o <- readRDS("/mnt/data/projects/spatial-rads/processing/norm/sam0003.norm.rds")
cts <- GetAssayData(o,assay="RNA",layer="counts"); md <- o@meta.data; md$cell <- rownames(md)
a <- read_parquet("results/aggregate/full_labels.parquet")
a <- a[grepl("^sam0003_",a$cell),c("cell","compartment")]; a$cell<-sub("^sam0003_","",a$cell)
z <- read_tsv("dev/peak_valley_analysis/data/mbrt4h_peak_valley.tsv",show_col_types=FALSE)[,c("cell_id","dist_to_peak")]
md <- md %>% left_join(a,by="cell") %>% left_join(z,by=c("cell"="cell_id"))
md <- md[!is.na(md$compartment),]

## generic paired FOV-pseudobulk voom DE; returns p-values + n sig
pv_de <- function(sub, split_col, lvls, floor=30) {
  sub$grp <- paste0("f",sub$fov,"_",sub[[split_col]])
  ag <- t(fac2sparse(sub$grp) %*% t(cts[,sub$cell]))
  g <- data.frame(grp=colnames(ag)); g$fov<-sub("_.*","",g$grp); g$k<-sub(".*_","",g$grp)
  nc<-table(sub$grp); g$ncell<-as.integer(nc[g$grp]); g<-g[g$ncell>=floor,]
  both<-names(which(tapply(g$k,g$fov,function(x)length(unique(x)))==2)); g<-g[g$fov%in%both,]
  ag<-ag[,g$grp]; g$fov<-factor(g$fov); g$k<-factor(g$k,levels=lvls)
  y<-DGEList(ag); des<-model.matrix(~fov+k,g); kg<-filterByExpr(y,des); y<-y[kg,,keep.lib.sizes=FALSE]
  y<-calcNormFactors(y); v<-voom(y,des); fit<-eBayes(lmFit(v,des))
  tt<-topTable(fit,coef=paste0("k",lvls[2]),number=Inf,sort.by="none"); tt<-tt[!rownames(tt)%in%DDR,]
  list(p=tt$P.Value, nsig=sum(tt$adj.P.Val<0.05), nfov=length(both), ngene=nrow(tt))
}

## peak/valley all-cells + per compartment (core bands)
md$zone <- ifelse(md$dist_to_peak<0.10,"peak",ifelse(md$dist_to_peak>0.40,"valley",NA))
res <- list()
res[["Peak vs valley: all cells"]] <- pv_de(md[!is.na(md$zone),],"zone",c("valley","peak"))
for(cp in c("tumor","stroma","immune")){
  s <- md[md$compartment==cp & !is.na(md$zone),]
  res[[paste0("Peak vs valley: ",cp)]] <- pv_de(s,"zone",c("valley","peak"))
}
## positive control: tumor vs stroma (same machinery)
s2 <- md[md$compartment %in% c("tumor","stroma"),]
res[["Positive control: tumor vs stroma"]] <- pv_de(s2,"compartment",c("stroma","tumor"))

df <- do.call(rbind, lapply(names(res), function(nm){
  r<-res[[nm]]
  lab <- sprintf("%s\n(%d FOVs, %d genes; %d at FDR<0.05)", nm, r$nfov, r$ngene, r$nsig)
  data.frame(panel=lab, p=r$p, pos=grepl("Positive",nm))
}))
ord <- unique(df$panel); ord <- c(ord[!grepl("Positive",ord)], ord[grepl("Positive",ord)])
df$panel <- factor(df$panel, levels=ord)

p <- ggplot(df, aes(x=p, fill=pos)) +
  geom_histogram(breaks=seq(0,1,0.05), color="grey30", linewidth=0.2) +
  geom_hline(yintercept=0, linewidth=0.3) +
  facet_wrap(~panel, scales="free_y", ncol=3) +
  scale_fill_manual(values=c("FALSE"="grey65","TRUE"="#0072B2"), guide="none") +
  labs(x="p-value", y="genes",
       title="Peak/valley DE is null in every compartment and all-cells; same machinery detects tumor vs stroma") +
  theme_scifig() +
  theme(strip.text=element_text(size=10), plot.title=element_text(size=12))

save_plot(p, "results/aggregate/plots/pv_phist_null", w=13, h=7)
