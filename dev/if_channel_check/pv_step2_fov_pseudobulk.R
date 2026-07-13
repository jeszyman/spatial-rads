# Step 2: swap the statistic from pooled-per-cell to FOV-pseudobulk, holding cell set
# (current atlas-tumor + current QC, as in Step 1) and geometry fixed. Isolates confound
# #4 (pseudoreplication). Two aggregation designs:
#   (A) paired within-FOV peak-vs-valley, ~ FOV + zone (FOV as blocking factor) -- the
#       honest within-section contrast that removes FOV-level architecture.
#   (B) continuous phase regression at the FOV unit: FOV-pseudobulk expression ~ mean
#       distance-to-nearest-peak, limma-voom. Judged on slope effect size + CI, no phase-null p.
# Readout: do the Step-1 residual survivors hold up when the replicate unit is the FOV?
suppressPackageStartupMessages({
  library(arrow); library(data.table); library(Matrix)
  library(DESeq2); library(limma); library(edgeR)
})

OUT   <- "/home/jeszyman/repos/spatial-rads/dev/if_channel_check"
META  <- "/mnt/data/projects/spatial-rads/inputs/mutter01/Analysis_Mutter_01_CosMmR_Mutter_updated_metadata.parquet"
CNT   <- "/mnt/data/projects/spatial-rads/inputs/mutter01/projects_Mutter_01_CosMmR_app_data_raw_tx_counts_matrix.parquet"
ATLAS <- "/home/jeszyman/repos/spatial-rads/results/aggregate/full_labels.parquet"
STEP1 <- file.path(OUT, "pv_step1_effectsize_swap.tsv")
SLIDE <- "20250529_214712_S4"; BLOCK <- "Block_21"
MIN_CELLS_UNIT <- 20
set.seed(1)

peak_fovs <- c(224,209,201,186,172,225,229,215,175,176,206,181,168,167)

meta <- as.data.table(read_parquet(META)); setnames(meta, "ImmuneAtlas_ImmGen_Main_cell_Types", "main_type")
blk  <- meta[Slide == SLIDE & Block == BLOCK]

## ---- geometry (identical to Steps 0/1) ----
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
label_fast <- function(d, offset, spacing, half_w) {
  r <- (d - offset) %% spacing; dtp <- pmin(r, spacing - r)
  list(zone = ifelse(dtp < half_w, "peak", ifelse(dtp > spacing/2 - half_w, "valley", "transition")), dtp = dtp) }
cat(sprintf("Geometry: theta=%.1f, spacing=%.3f, peak_half=%.3f\n", theta*180/pi, spacing, peak_half))

## ---- atlas-tumor cells on this block ----
atl <- as.data.table(read_parquet(ATLAS))[dataset == "Mutter_01"]
atl[, mcid := sub("^sam[0-9]+_", "", cell)]
tumor_ids <- intersect(atl[compartment == "tumor", mcid], blk$cell_id)   # Block_21 only
cat(sprintf("atlas-tumor cells on Block_21: %d\n", length(tumor_ids)))

## ---- counts for atlas-tumor cells; attach fov, zone, dtp ----
cnt <- as.data.table(read_parquet(CNT))[cell_id %in% tumor_ids]
gene_cols <- setdiff(colnames(cnt), c("Slide","fov","cell_id"))
pos <- blk[cell_id %in% cnt$cell_id, .(cell_id, fov, x_slide_mm, y_slide_mm)]
setkey(cnt, cell_id); setkey(pos, cell_id); pos <- pos[cnt$cell_id]
d <- -sin(theta)*pos$x_slide_mm + cos(theta)*pos$y_slide_mm
lab <- label_fast(d, c0, spacing, peak_half)
pos[, `:=`(zone = lab$zone, dtp = lab$dtp)]
cat(sprintf("atlas-tumor cells with counts: %d\n", nrow(pos)))
cat("cell zone table:\n"); print(table(pos$zone))

cmat <- as.matrix(cnt[, ..gene_cols])          # cells x genes
rownames(cmat) <- cnt$cell_id

## ================= Design A: paired within-FOV, ~ FOV + zone =================
keep <- pos$zone %in% c("peak","valley")
agg_key <- paste(pos$fov[keep], pos$zone[keep], sep="__")
units <- rowsum(cmat[keep, , drop=FALSE], group = agg_key)     # units x genes
ncell <- as.integer(table(agg_key)[rownames(units)])
udt <- data.table(unit = rownames(units), ncell = ncell)
udt[, c("fov","zone") := tstrsplit(unit, "__")]
udt <- udt[ncell >= MIN_CELLS_UNIT]
# keep FOVs with BOTH zones (paired)
paired_fovs <- udt[, .N, by=fov][N == 2, fov]
udt <- udt[fov %in% paired_fovs]
units <- units[udt$unit, , drop=FALSE]
cat(sprintf("\nDesign A: %d FOV-zone units across %d paired FOVs (%d peak / %d valley)\n",
            nrow(udt), length(paired_fovs), sum(udt$zone=="peak"), sum(udt$zone=="valley")))

cd <- t(units)                                  # genes x units
gkeep <- rowSums(cd >= 1) >= (0.1 * ncol(cd)) & rowSums(cd) >= 10
cd <- cd[gkeep, ]
coldata <- data.frame(fov = factor(udt$fov), zone = factor(udt$zone, levels=c("valley","peak")),
                      row.names = udt$unit)
dds <- DESeqDataSetFromMatrix(round(cd), coldata, design = ~ fov + zone)
dds <- DESeq(dds, quiet = TRUE)
resA <- as.data.table(as.data.frame(results(dds, name = "zone_peak_vs_valley")), keep.rownames = "gene")
resA[, `:=`(ci_lo = log2FoldChange - 1.96*lfcSE, ci_hi = log2FoldChange + 1.96*lfcSE)]
setorder(resA, padj)
nsurv_A <- resA[padj < 0.05, .N]
cat(sprintf("Design A survivors (padj<0.05, ~FOV+zone): %d\n", nsurv_A))

## ================= Design B: continuous phase regression, FOV unit =================
# all Block_21 atlas-tumor cells contribute to FOV-level pseudobulk
fov_units <- rowsum(cmat, group = as.character(pos$fov))
fmeta <- pos[, .(fov = as.character(fov), dtp)][, .(mean_dtp = mean(dtp), ncell = .N), by = fov]
fmeta <- fmeta[match(rownames(fov_units), fov)]   # align to matrix row order (character keys)
keepU <- fmeta$ncell >= MIN_CELLS_UNIT
fov_units <- fov_units[keepU, ]; fmeta <- fmeta[keepU]
cat(sprintf("\nDesign B: %d FOV pseudobulk units, mean_dtp range [%.3f, %.3f]\n",
            nrow(fmeta), min(fmeta$mean_dtp), max(fmeta$mean_dtp)))
cB <- t(fov_units); gkeepB <- rowSums(cB >= 1) >= (0.1*ncol(cB)) & rowSums(cB) >= 10
cB <- cB[gkeepB, ]
dge <- DGEList(cB); dge <- calcNormFactors(dge)
# phase covariate: mean distance-to-nearest-peak; NEGATIVE slope = peak-up
design <- model.matrix(~ fmeta$mean_dtp)
v <- voom(dge, design); fitB <- eBayes(lmFit(v, design))
resB <- as.data.table(topTable(fitB, coef=2, number=Inf), keep.rownames="gene")
resB[, peak_up_slope := -logFC]                   # sign so positive = peak-up
setnames(resB, c("logFC","adj.P.Val"), c("phase_slope","phase_padj"))

## ================= interest gene set from Step 1 =================
s1 <- fread(STEP1)
old_surv <- s1[old_surv==TRUE, gene]; new_surv <- s1[new_surv==TRUE, gene]
ifn_mhc  <- c("B2m","H2-K1","H2-D1","H2-Q7","Ifitm2","Ifitm3","Ifi27","Oas1a/g","Stat1","Irf1","Tap1")
interest <- unique(c(old_surv, new_surv, ifn_mhc, "Clu","Tpt1","Itm2b","Cdkn1a"))

tab <- merge(resA[, .(gene, A_lfc=log2FoldChange, A_lo=ci_lo, A_hi=ci_hi, A_padj=padj)],
             resB[, .(gene, B_slope=peak_up_slope, B_padj=phase_padj)], by="gene", all=TRUE)
tab[, in_step1_new_surv := gene %in% new_surv]
tab[, is_ifn_mhc := gene %in% ifn_mhc]
fwrite(tab, file.path(OUT, "pv_step2_fov_pseudobulk.tsv"), sep="\t")

show <- tab[gene %in% interest][order(-A_lfc)]
cat("\n=== Step-1 survivors + IFN/MHC under FOV-pseudobulk (Design A) + phase slope (Design B) ===\n")
cat("A_lfc = within-FOV peak-vs-valley log2FC (paired); positive=peak-up. B_slope: positive=peak-up.\n")
print(show[, .(gene,
               A_lfc = round(A_lfc,3), A_ci = sprintf("[%.2f,%.2f]", A_lo, A_hi),
               A_padj = sprintf("%.2g", A_padj),
               B_slope = round(B_slope,3), B_padj = sprintf("%.2g", B_padj),
               new_surv = in_step1_new_surv)])

cat(sprintf("\n=== SUMMARY ===\n"))
cat(sprintf("Design A (~FOV+zone) genome-wide survivors padj<0.05: %d\n", nsurv_A))
cat(sprintf("Design B (phase regression) genome-wide padj<0.05: %d\n", resB[phase_padj<0.05, .N]))
cat(sprintf("Step-1 residual new-survivors also A-significant: %d / %d\n",
            tab[in_step1_new_surv==TRUE & A_padj<0.05, .N], length(new_surv)))
cat(sprintf("IFN/MHC members with A CI excluding 0 (peak-up): %d / %d\n",
            tab[is_ifn_mhc==TRUE & A_lo>0, .N], length(ifn_mhc)))
cat("Key single genes:\n")
for (g in c("Clu","H2-D1","B2m","Cdkn1a","Tpt1","Itm2b")) {
  r <- tab[gene==g]; if (nrow(r)) cat(sprintf("  %-8s A_lfc=%.3f CI[%.2f,%.2f] padj=%.2g | B_slope=%.3f padj=%.2g\n",
      g, r$A_lfc, r$A_lo, r$A_hi, r$A_padj, r$B_slope, r$B_padj)) }
