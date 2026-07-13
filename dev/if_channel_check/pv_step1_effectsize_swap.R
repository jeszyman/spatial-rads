# Step 1: does the 37-gene peak signature survive swapping the retired 'a' bucket +
# old counts for current atlas-tumor + current-QC cells? Isolates confounds #1 (label)
# and #2 (QC), holding the statistic (pooled per-cell log-CPM delta) and the wide real
# bands fixed. Primary readout = effect-size stability (collapse-toward-zero vs stable),
# per the reframe; the rotation-null p is retained only to regenerate the original
# survivor set and to check whether the swapped set produces any survivors at all.
suppressPackageStartupMessages({ library(arrow); library(data.table); library(Matrix) })

OUT   <- "/home/jeszyman/repos/spatial-rads/dev/if_channel_check"
META  <- "/mnt/data/projects/spatial-rads/inputs/mutter01/Analysis_Mutter_01_CosMmR_Mutter_updated_metadata.parquet"
CNT   <- "/mnt/data/projects/spatial-rads/inputs/mutter01/projects_Mutter_01_CosMmR_app_data_raw_tx_counts_matrix.parquet"
ATLAS <- "/home/jeszyman/repos/spatial-rads/results/aggregate/full_labels.parquet"
SLIDE <- "20250529_214712_S4"; BLOCK <- "Block_21"
set.seed(1)

peak_fovs <- c(224,209,201,186,172,225,229,215,175,176,206,181,168,167)

meta <- as.data.table(read_parquet(META)); setnames(meta, "ImmuneAtlas_ImmGen_Main_cell_Types", "main_type")
blk  <- meta[Slide == SLIDE & Block == BLOCK]

## ---- shared stripe geometry (cell-set independent) ----
pc <- blk[fov %in% peak_fovs, .(cx = mean(x_slide_mm), cy = mean(y_slide_mm)), by = fov]
search_theta <- function(cx, cy, n = 4, grid = seq(-pi/2, pi/2, length.out = 181)) {
  best <- list(within = Inf)
  for (th in grid) { d <- -sin(th)*cx + cos(th)*cy; km <- kmeans(d, n, nstart = 10)
    if (km$tot.withinss < best$within) best <- list(theta=th, within=km$tot.withinss, centers=sort(km$centers[,1]), cluster=km$cluster) }
  best }
fit <- search_theta(pc$cx, pc$cy); theta <- fit$theta; spacing <- median(diff(fit$centers))
pc[, d_perp := -sin(theta)*cx + cos(theta)*cy][, stripe := fit$cluster]
peak_half <- pc[, .(sd = sd(d_perp - fit$centers[stripe])), by=stripe][, mean(sd)] + 0.15
c0 <- fit$centers[1]
cat(sprintf("Geometry: theta=%.1f, spacing=%.3f, peak_half=%.3f\n", theta*180/pi, spacing, peak_half))

## ---- vectorized nearest-peak labeling on a regular grid (semantics match original) ----
label_fast <- function(d, offset, spacing, half_w) {
  r <- (d - offset) %% spacing; dtp <- pmin(r, spacing - r)
  ifelse(dtp < half_w, "peak", ifelse(dtp > spacing/2 - half_w, "valley", "transition")) }

## ---- atlas tumor cell ids on this block (current QC + tumor label) ----
atl <- as.data.table(read_parquet(ATLAS))[dataset == "Mutter_01"]
atl[, mcid := sub("^sam[0-9]+_", "", cell)]
atlas_tumor_ids <- atl[compartment == "tumor", mcid]

## ---- load counts once, subset per cell set ----
cnt <- as.data.table(read_parquet(CNT))
cnt <- cnt[cell_id %in% blk$cell_id]
gene_cols <- setdiff(colnames(cnt), c("Slide","fov","cell_id"))
posmap <- blk[, .(cell_id, x_slide_mm, y_slide_mm)]

# per-cell-set: returns real delta per gene + rotation-null survivor calls
run_set <- function(ids, tag, n_null = 200) {
  cd <- cnt[cell_id %in% ids]
  m  <- posmap[cell_id %in% cd$cell_id]
  setkey(cd, cell_id); setkey(m, cell_id); m <- m[cd$cell_id]
  mat <- t(as.matrix(cd[, ..gene_cols])); lib <- colSums(mat); lib[lib==0] <- 1
  normed <- log2(t(t(mat)/lib) * 1e4 + 1)
  keep <- rowSums(mat > 0)/ncol(mat) > 0.05; normed <- normed[keep, ]
  d <- -sin(theta)*m$x_slide_mm + cos(theta)*m$y_slide_mm
  delta <- function(lab) { pi_ <- which(lab=="peak"); vi <- which(lab=="valley")
    if (length(pi_)<50 || length(vi)<50) return(NULL)
    rowMeans(normed[,pi_,drop=FALSE]) - rowMeans(normed[,vi,drop=FALSE]) }
  real_lab <- label_fast(d, c0, spacing, peak_half)
  rd <- delta(real_lab)
  # rotation null: random angle (>=30deg from real) + random phase
  raw <- runif(n_null*3, -pi/2, pi/2)
  tn  <- head(raw[pmin(abs(raw-theta), abs(raw-(theta+pi)), abs(raw-(theta-pi))) > pi/6], n_null)
  nd  <- matrix(NA_real_, nrow=length(rd), ncol=length(tn), dimnames=list(names(rd), NULL))
  for (i in seq_along(tn)) { di <- -sin(tn[i])*m$x_slide_mm + cos(tn[i])*m$y_slide_mm
    dd <- delta(label_fast(di, runif(1,0,spacing), spacing, peak_half)); if(!is.null(dd)) nd[,i] <- dd }
  osd <- sd(as.vector(nd), na.rm=TRUE)
  empp <- sapply(seq_along(rd), function(g){ x<-nd[g,]; x<-x[!is.na(x)]; if(length(x)<10) NA_real_ else sum(abs(x)>=abs(rd[g]))/length(x) })
  surv <- !is.na(empp) & empp < 0.05 & abs(rd) > 2*osd
  cat(sprintf("[%s] %d cells, %d genes, null sd=%.3f, survivors=%d\n", tag, ncol(normed), nrow(normed), osd, sum(surv)))
  data.table(gene=names(rd), delta=as.numeric(rd), emp_p=empp, survives=surv, cells=ncol(normed))
}

old <- run_set(blk[main_type=="a", cell_id],   "OLD a-bucket")
new <- run_set(atlas_tumor_ids,                "NEW atlas-tumor")

## ---- the survivor set is defined by the OLD run (reproduces the original ~37) ----
sig_genes <- old[survives==TRUE, gene]
cat(sprintf("\nOriginal-method survivors on OLD set: %d genes\n", length(sig_genes)))

cmp <- merge(old[, .(gene, old_delta=delta, old_p=emp_p, old_surv=survives)],
             new[, .(gene, new_delta=delta, new_p=emp_p, new_surv=survives)], by="gene")
cmp[, is_sig := gene %in% sig_genes]

## ---- bootstrap CI on the delta for the survivor genes, both cell sets ----
boot_ci <- function(ids, genes, nboot=500) {
  cd <- cnt[cell_id %in% ids]; m <- posmap[cell_id %in% cd$cell_id]
  setkey(cd, cell_id); setkey(m, cell_id); m <- m[cd$cell_id]
  mat <- t(as.matrix(cd[, ..genes])); lib <- colSums(t(as.matrix(cd[, ..gene_cols]))); lib[lib==0]<-1
  normed <- log2(t(t(mat)/lib)*1e4 + 1)
  d <- -sin(theta)*m$x_slide_mm + cos(theta)*m$y_slide_mm
  lab <- label_fast(d, c0, spacing, peak_half); pi_<-which(lab=="peak"); vi<-which(lab=="valley")
  bs <- matrix(NA_real_, nrow=length(genes), ncol=nboot)
  for (b in 1:nboot) { pb<-sample(pi_,replace=TRUE); vb<-sample(vi,replace=TRUE)
    bs[,b] <- rowMeans(normed[,pb,drop=FALSE]) - rowMeans(normed[,vb,drop=FALSE]) }
  data.table(gene=genes, lo=apply(bs,1,quantile,0.025), hi=apply(bs,1,quantile,0.975)) }

if (length(sig_genes) > 0) {
  ob <- boot_ci(blk[main_type=="a", cell_id], sig_genes); setnames(ob,c("gene","old_lo","old_hi"))
  nb <- boot_ci(atlas_tumor_ids, sig_genes);              setnames(nb,c("gene","new_lo","new_hi"))
  cmp <- merge(merge(cmp, ob, by="gene", all.x=TRUE), nb, by="gene", all.x=TRUE)
  # collapse = new CI includes 0 while old CI excluded 0
  cmp[, collapsed := is_sig & (old_lo>0 | old_hi<0) & (new_lo<=0 & new_hi>=0)]
}

setorder(cmp, -old_delta)
fwrite(cmp, file.path(OUT, "pv_step1_effectsize_swap.tsv"), sep="\t")

sig <- cmp[is_sig==TRUE]
cat(sprintf("\n=== Survivor genes (%d): old vs new peak-valley delta ===\n", nrow(sig)))
print(sig[, .(gene, old_delta=round(old_delta,3), new_delta=round(new_delta,3),
              new_ci=sprintf("[%.2f,%.2f]", new_lo, new_hi),
              new_surv, collapsed)][order(-old_delta)])
cat(sprintf("\nSurvivors retained by new-set rotation null: %d / %d\n", sum(sig$new_surv), nrow(sig)))
cat(sprintf("Survivors whose effect collapsed to a 0-spanning CI: %d / %d\n", sum(sig$collapsed, na.rm=TRUE), nrow(sig)))
cat(sprintf("Genome-wide delta correlation old vs new: r=%.3f (Spearman rho=%.3f)\n",
            cor(cmp$old_delta, cmp$new_delta, use="complete.obs"),
            cor(cmp$old_delta, cmp$new_delta, use="complete.obs", method="spearman")))
cat(sprintf("Median |delta| shrinkage on survivors: old=%.3f new=%.3f (ratio %.2f)\n",
            median(abs(sig$old_delta)), median(abs(sig$new_delta)),
            median(abs(sig$new_delta))/median(abs(sig$old_delta))))
