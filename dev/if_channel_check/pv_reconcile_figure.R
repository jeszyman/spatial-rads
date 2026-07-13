# Reconciliation figure: peak/valley is null under the honest within-FOV test.
# A: Block_21 atlas-tumor cells on the real fixed beam geometry (peak/valley/transition).
# B: raw p-value histogram from the paired within-FOV test (~ FOV + zone), flat = null.
suppressPackageStartupMessages({
  library(arrow); library(data.table); library(Matrix)
  library(DESeq2); library(ggplot2); library(patchwork)
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

atl<-as.data.table(read_parquet(ATLAS))[dataset=="Mutter_01"];atl[,mcid:=sub("^sam[0-9]+_","",cell)]
tumor_ids<-intersect(atl[compartment=="tumor",mcid],blk$cell_id)
cnt<-as.data.table(read_parquet(CNT))[cell_id %in% tumor_ids]
gene_cols<-setdiff(colnames(cnt),c("Slide","fov","cell_id"))
pos<-blk[cell_id %in% cnt$cell_id,.(cell_id,fov,x_slide_mm,y_slide_mm)]
setkey(cnt,cell_id);setkey(pos,cell_id);pos<-pos[cnt$cell_id]
d<--sin(theta)*pos$x_slide_mm+cos(theta)*pos$y_slide_mm
r<-(d-c0)%%spacing;dtp<-pmin(r,spacing-r)
pos[,zone:=ifelse(dtp<peak_half,"peak",ifelse(dtp>spacing/2-peak_half,"valley","transition"))]

## ---- panel B: within-FOV paired test raw p ----
cmat<-as.matrix(cnt[,..gene_cols]);rownames(cmat)<-cnt$cell_id
keep<-pos$zone %in% c("peak","valley");agg<-paste(pos$fov[keep],pos$zone[keep],sep="__")
units<-rowsum(cmat[keep,,drop=FALSE],agg);nc<-as.integer(table(agg)[rownames(units)])
udt<-data.table(unit=rownames(units),ncell=nc);udt[,c("fv","zn"):=tstrsplit(unit,"__")]
udt<-udt[ncell>=MIN_CELLS_UNIT];pf<-udt[,.N,by=fv][N==2,fv];udt<-udt[fv %in% pf];units<-units[udt$unit,]
cd<-t(units);gk<-rowSums(cd>=1)>=0.1*ncol(cd)&rowSums(cd)>=10;cd<-cd[gk,]
coldata<-data.frame(fov=factor(udt$fv),zone=factor(udt$zn,levels=c("valley","peak")),row.names=udt$unit)
dds<-DESeq(DESeqDataSetFromMatrix(round(cd),coldata,~fov+zone),quiet=TRUE)
pv<-results(dds,name="zone_peak_vs_valley")$pvalue;pv<-pv[!is.na(pv)]
frac<-mean(pv<0.05);n_pf<-length(pf)

## ---- palette ----
cols<-c(peak="#c0392b",valley="#2c6fbb",transition="grey82")
th_deg<-theta*180/pi

## ---- rotation arc + label (centroid) ----
xr<-range(pos$x_slide_mm); yr<-range(pos$y_slide_mm); cx<-mean(xr); cy<-mean(yr)
arc_r<-1.05
arc_ang<-seq(theta+pi/2, theta+pi/2+70*pi/180, length.out=60)
arc<-data.table(x=cx+arc_r*cos(arc_ang), y=cy+arc_r*sin(arc_ang))

## ---- panel A ----
lines_dt<-data.table(int=fit$centers/cos(theta))   # y = tan(theta) x + c/cos(theta)
pA<-ggplot(pos,aes(x_slide_mm,y_slide_mm,color=zone))+
  geom_point(size=0.15,alpha=0.5,shape=16)+
  geom_abline(data=lines_dt,aes(slope=tan(theta),intercept=int),
              linetype="dashed",linewidth=0.35,color="grey25")+
  geom_path(data=arc,aes(x,y),inherit.aes=FALSE,color="grey10",linewidth=1.3,
            arrow=arrow(length=unit(0.45,"cm"),type="closed"))+
  annotate("text",x=cx,y=cy+arc_r+0.42,label="randomize\nbeam angle",
           size=6.5,fontface="bold.italic",color="grey10",lineheight=0.9)+
  scale_color_manual(values=cols,breaks=c("peak","valley","transition"),
                     labels=c("peak","valley","transition"))+
  coord_equal()+
  guides(color=guide_legend(override.aes=list(size=2.5,alpha=1)))+
  labs(x="slide x (mm)",y="slide y (mm)",color=NULL,
       title="A  Dose peak/valley zones on the irradiated tumor",
       subtitle=sprintf("MBRT 4h, atlas-defined tumor cells (n=%s); fitted beam stripes %.0f°, %.2f mm spacing",
                        format(nrow(pos),big.mark=","),th_deg,spacing))+
  theme_classic(base_size=12)+
  theme(legend.position=c(0.99,0.02),legend.justification=c(1,0),
        legend.background=element_rect(fill="white",color="grey80",linewidth=0.3),
        plot.title=element_text(face="bold",size=12),plot.subtitle=element_text(size=9.5))

## ---- panel B ----
pB<-ggplot(data.table(p=pv),aes(p))+
  geom_histogram(breaks=seq(0,1,0.05),fill="#4a7fb5",color="white",linewidth=0.25)+
  geom_hline(yintercept=length(pv)/20,linetype="dashed",color="grey30")+
  annotate("text",x=0.97,y=length(pv)/20,label="uniform (true null)",
           hjust=1,vjust=-0.6,size=3,color="grey30")+
  labs(x="raw p-value (peak vs valley)",y="genes",
       title="B  No gene distinguishes peak from valley",
       subtitle=sprintf("Paired within-field test (~ field + zone), %d fields; fraction p<0.05 = %.3f (null expects 0.05)",
                        n_pf,frac))+
  theme_classic(base_size=12)+
  theme(plot.title=element_text(face="bold",size=12),plot.subtitle=element_text(size=9.5))

fig<-pA/pB+plot_layout(heights=c(1.35,1))
ggsave(file.path(OUT,"pv_reconcile_figure.png"),fig,width=8,height=10,dpi=150)
cat(sprintf("saved. panel A n=%d cells (peak %d/valley %d/transition %d); panel B %d genes, %d fields, frac p<0.05=%.3f\n",
            nrow(pos),sum(pos$zone=="peak"),sum(pos$zone=="valley"),sum(pos$zone=="transition"),
            length(pv),n_pf,frac))
