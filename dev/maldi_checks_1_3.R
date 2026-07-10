# Checks 1 & 3 for the MALDI dual-programme hypothesis (Baranwal lipidomics deck).
# Check 1: tumor-restricted Type-I IFN (Programme B) by treatment x timepoint (MBRT vs SBRT, 4h vs day2).
# Check 3: day-2 tumor per-cell correlation of DDR (lingering damage) vs Type-I IFN (immune signal).
# Tumor = per-sample v1 labels a + tumor_epithelial (documented convention); flank only.
suppressPackageStartupMessages({library(Seurat); library(data.table)})

scored_dir <- "/mnt/data/projects/spatial-rads/processing/scored"
flank <- sprintf("sam%04d", c(1:8, 12:23))   # exclude tongue 9,10,11
tumor_set <- c("a", "tumor_epithelial")
score_cols <- c("TypeI_interferon_UCell", "DNA_Damage_Repair_UCell", "STING_UCell")
meta_cols  <- c("sample_id","dataset","treatment","condition","timepoint_h","model","cell_type", score_cols)

per_sample <- list(); day2_tumor <- list()
for (s in flank) {
  o <- readRDS(file.path(scored_dir, paste0(s, ".scored.rds")))
  md <- as.data.table(o@meta.data)[, ..meta_cols]
  rm(o); gc()
  tu <- md[cell_type %in% tumor_set]
  per_sample[[s]] <- data.table(
    sample_id = s, dataset = md$dataset[1], treatment = md$treatment[1],
    timepoint_h = md$timepoint_h[1], model = md$model[1],
    n_cells = nrow(md), n_tumor = nrow(tu),
    ifn = mean(tu$TypeI_interferon_UCell, na.rm=TRUE),
    ddr = mean(tu$DNA_Damage_Repair_UCell, na.rm=TRUE),
    sting = mean(tu$STING_UCell, na.rm=TRUE))
  if (md$timepoint_h[1] == 48)
    day2_tumor[[s]] <- tu[, .(sample_id=s, dataset=md$dataset[1], treatment=md$treatment[1],
                              ddr=DNA_Damage_Repair_UCell, ifn=TypeI_interferon_UCell)]
  cat("done", s, "\n")
}
ps <- rbindlist(per_sample)
fwrite(ps, "dev/maldi_check1_per_sample.tsv", sep="\t")

cat("\n================ CHECK 1: tumor Type-I IFN (UCell) by sample ================\n")
print(ps[order(dataset, timepoint_h, treatment),
         .(sample_id, dataset, treatment, timepoint_h, n_tumor,
           ifn=round(ifn,4), ddr=round(ddr,4), sting=round(sting,4))])

cat("\n---- M01 flank 4h vs day2 (n=1, descriptive) ----\n")
print(ps[dataset=="Mutter_01" & timepoint_h %in% c(4,48),
         .(treatment, timepoint_h, ifn=round(ifn,4), sting=round(sting,4))][order(timepoint_h,treatment)])

cat("\n---- M02 day-2 tumor IFN, mean across n=4 replicates (inferential cohort) ----\n")
print(ps[dataset=="Mutter_02",
         .(n=.N, ifn_mean=round(mean(ifn),4), ifn_sd=round(sd(ifn),4),
           sting_mean=round(mean(sting),4), ddr_mean=round(mean(ddr),4)), by=treatment])

d2 <- rbindlist(day2_tumor)
cat("\n================ CHECK 3: day-2 tumor per-cell Spearman(DDR, IFN) ================\n")
cat("---- pooled by treatment (all day-2 flank tumor cells) ----\n")
print(d2[, .(n_cells=.N, rho=round(cor(ddr, ifn, method="spearman"),3)), by=treatment])
cat("\n---- per-sample, M02 day-2 ----\n")
print(d2[dataset=="Mutter_02",
         .(n=.N, rho=round(cor(ddr, ifn, method="spearman"),3)), by=.(treatment, sample_id)][order(treatment, sample_id)])
cat("\nALL DONE\n")
