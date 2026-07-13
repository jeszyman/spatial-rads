# Raw p-value histograms for the Step 2 FOV-pseudobulk nulls, to test the true-null
# claim (flat = well-calibrated null; left-spike = real signal; right-slope/U = mis-cal).
suppressPackageStartupMessages({
  library(arrow); library(data.table); library(Matrix)
  library(DESeq2); library(limma); library(edgeR); library(ggplot2); library(patchwork)
})
OUT   <- "/home/jeszyman/repos/spatial-rads/dev/if_channel_check"
META  <- "/mnt/data/projects/spatial-rads/inputs/mutter01/Analysis_Mutter_01_CosMmR_Mutter_updated_metadata.parquet"
CNT   <- "/mnt/data/projects/spatial-rads/inputs/mutter01/projects_Mutter_01_CosMmR_app_data_raw_tx_counts_matrix.parquet"
ATLAS <- "/home/jeszyman/repos/spatial-rads/results/aggregate/full_labels.parquet"
SLIDE <- "20250529_214712_S4"; BLOCK <- "Block_21"; MIN_CELLS_UNIT <- 20
set.seed(1)
peak_fovs <- c(224,209,201,186,172,225,229,215,175,176,206,181,168,167)

meta <- as.data.table(read_parquet(META)); setnames(meta,"ImmuneAtlas_ImmGen_Main_cell_Types","main_type")
blk  <- meta[Slide==SLIDE & Block==BLOCK]
pc <- blk[fov %in% peak_fovs, .(cx=mean(x_slide_mm), cy=mean(y_slide_mm)), by=fov]
search_theta <- function(cx,cy,n=4,grid=seq(-pi/2,pi/2,length.out=181)){best<-list(within=Inf)
  for(th in grid){d<--sin(th)*cx+cos(th)*cy;km<-kmeans(d,n,nstart=10)
    if(km$tot.withinss<best$within)best<-list(theta=th,within=km$tot.withinss,centers=sort(km$centers[,1]),cluster=km$cluster)};best}
fit<-search_theta(pc$cx,pc$cy);theta<-fit$theta;spacing<-median(diff(fit$centers))
pc[,d_perp:=-sin(theta)*cx+cos(theta)*cy][,stripe:=fit$cluster]
peak_half<-pc[,.(sd=sd(d_perp-fit$centers[stripe])),by=stripe][,mean(sd)]+0.15;c0<-fit$centers[1]
label_fast<-function(d,off,sp,hw){r<-(d-off)%%sp;dtp<-pmin(r,sp-r)
  list(zone=ifelse(dtp<hw,"peak",ifelse(dtp>sp/2-hw,"valley","transition")),dtp=dtp)}

atl<-as.data.table(read_parquet(ATLAS))[dataset=="Mutter_01"];atl[,mcid:=sub("^sam[0-9]+_","",cell)]
tumor_ids<-intersect(atl[compartment=="tumor",mcid],blk$cell_id)
cnt<-as.data.table(read_parquet(CNT))[cell_id %in% tumor_ids]
gene_cols<-setdiff(colnames(cnt),c("Slide","fov","cell_id"))
pos<-blk[cell_id %in% cnt$cell_id,.(cell_id,fov,x_slide_mm,y_slide_mm)]
setkey(cnt,cell_id);setkey(pos,cell_id);pos<-pos[cnt$cell_id]
d<--sin(theta)*pos$x_slide_mm+cos(theta)*pos$y_slide_mm;lab<-label_fast(d,c0,spacing,peak_half)
pos[,`:=`(zone=lab$zone,dtp=lab$dtp)];cmat<-as.matrix(cnt[,..gene_cols]);rownames(cmat)<-cnt$cell_id

## Design A raw p
keep<-pos$zone %in% c("peak","valley");agg<-paste(pos$fov[keep],pos$zone[keep],sep="__")
units<-rowsum(cmat[keep,,drop=FALSE],agg);nc<-as.integer(table(agg)[rownames(units)])
udt<-data.table(unit=rownames(units),ncell=nc);udt[,c("fov","zone"):=tstrsplit(unit,"__")]
udt<-udt[ncell>=MIN_CELLS_UNIT];pf<-udt[,.N,by=fov][N==2,fov];udt<-udt[fov %in% pf];units<-units[udt$unit,]
cd<-t(units);gk<-rowSums(cd>=1)>=0.1*ncol(cd)&rowSums(cd)>=10;cd<-cd[gk,]
coldata<-data.frame(fov=factor(udt$fov),zone=factor(udt$zone,levels=c("valley","peak")),row.names=udt$unit)
dds<-DESeq(DESeqDataSetFromMatrix(round(cd),coldata,~fov+zone),quiet=TRUE)
pA<-results(dds,name="zone_peak_vs_valley")$pvalue;pA<-pA[!is.na(pA)]

## Design B raw p
fu<-rowsum(cmat,as.character(pos$fov))
fm<-pos[,.(fov=as.character(fov),dtp)][,.(mean_dtp=mean(dtp),ncell=.N),by=fov][match(rownames(fu),fov)]
kU<-fm$ncell>=MIN_CELLS_UNIT;fu<-fu[kU,];fm<-fm[kU]
cB<-t(fu);gkB<-rowSums(cB>=1)>=0.1*ncol(cB)&rowSums(cB)>=10;cB<-cB[gkB,]
dge<-calcNormFactors(DGEList(cB));des<-model.matrix(~fm$mean_dtp)
fitB<-eBayes(lmFit(voom(dge,des),des));pB<-topTable(fitB,coef=2,number=Inf)$P.Value

cat(sprintf("Design A: %d genes, KS-vs-uniform p=%.3f, frac p<0.05=%.3f\n",
            length(pA),ks.test(pA,"punif")$p.value,mean(pA<0.05)))
cat(sprintf("Design B: %d genes, KS-vs-uniform p=%.3f, frac p<0.05=%.3f\n",
            length(pB),ks.test(pB,"punif")$p.value,mean(pB<0.05)))

mk<-function(p,title){dt<-data.table(p=p)
  ggplot(dt,aes(p))+geom_histogram(breaks=seq(0,1,0.05),fill="#4a7fb5",color="white",linewidth=0.2)+
    geom_hline(yintercept=length(p)/20,linetype="dashed",color="grey40")+
    labs(x="raw p-value",y="genes",title=title,
         subtitle=sprintf("%d genes | frac p<0.05 = %.3f (uniform: 0.050) | KS p = %.2f",
                          length(p),mean(p<0.05),ks.test(p,"punif")$p.value))+
    theme_minimal(base_size=11)+theme(plot.subtitle=element_text(size=9))}
g<-mk(pA,"Design A: within-FOV peak vs valley (~ FOV + zone)")/
   mk(pB,"Design B: continuous phase regression (FOV unit)")
ggsave(file.path(OUT,"pv_step2_phist.png"),g,width=7,height=7,dpi=140)
cat("saved pv_step2_phist.png\n")
