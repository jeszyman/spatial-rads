suppressPackageStartupMessages({
  library(arrow); library(data.table)
})

meta <- as.data.table(read_parquet("/tmp/mutter01_meta.parquet"))
counts_df <- as.data.table(read_parquet("/tmp/mutter01_counts.parquet"))
design <- as.data.table(readxl::read_excel("/home/jeszyman/repos/spatial-rads/data/metadata.xlsx"))

meta <- merge(meta, design, by.x = c("Slide", "Block"), by.y = c("slide", "block_id"))
fam <- "ImmuneAtlas_ImmGen_cellFamily_Main_cell_Types"
setnames(meta, fam, "family")

# Focus on tumor "a" cells across conditions
a <- meta[family == "a" & (is.na(model) | model == "flank"), .(cell_id, condition, Slide, Block, treatment, timepoint_h)]
cat("Tumor 'a' cells per condition:\n")
print(a[, .N, by = condition][order(condition)])

# Match counts for these cells
ac <- counts_df[cell_id %in% a$cell_id]
stopifnot(nrow(ac) == nrow(a))
setkeyv(a, "cell_id"); setkeyv(ac, "cell_id")

# Library size and raw Cdkn1a per cell
a[, lib := rowSums(ac[, !c("Slide", "fov", "cell_id"), with = FALSE])[match(cell_id, ac$cell_id)]]
a[, cdkn1a_raw := ac$Cdkn1a[match(cell_id, ac$cell_id)]]
a[, cdkn1a_norm := log2((cdkn1a_raw / pmax(lib, 1)) * 1e4 + 1)]

cat("\nTumor 'a' cell summary per condition (raw + normalized):\n")
print(a[, .(n = .N,
            mean_lib = round(mean(lib)),
            mean_cdkn1a_raw = round(mean(cdkn1a_raw), 3),
            median_cdkn1a_raw = median(cdkn1a_raw),
            sum_cdkn1a = sum(cdkn1a_raw),
            pct_expr = round(100 * mean(cdkn1a_raw > 0), 2),
            mean_cdkn1a_norm = round(mean(cdkn1a_norm), 3)), by = condition])

# For comparison: Cdkn1a in OTHER families in SBRT_4h
cat("\nSBRT_4h: Cdkn1a by family (top 10 families by n):\n")
s4h_all <- meta[condition == "SBRT_4h", .(cell_id, family)]
s4h_c <- counts_df[cell_id %in% s4h_all$cell_id]
s4h_all <- merge(s4h_all, s4h_c[, .(cell_id, Cdkn1a)], by = "cell_id")
s4h_all[, lib := rowSums(s4h_c[match(cell_id, s4h_c$cell_id), !c("Slide", "fov", "cell_id"), with=FALSE])]
tab <- s4h_all[, .(n = .N, pct_expr = 100*mean(Cdkn1a > 0), mean_raw = mean(Cdkn1a)), by = family][order(-n)][1:10]
print(tab)

# Global: how many cells have Cdkn1a > 0 in SBRT_4h block at all?
cat("\nSBRT_4h (Block_17): cells with Cdkn1a > 0 across ALL families:\n")
b17 <- meta[Slide == "20250529_214712_S4" & Block == "Block_17"]
b17_c <- counts_df[cell_id %in% b17$cell_id]
cat("Total cells:", nrow(b17), "| with Cdkn1a > 0:", sum(b17_c$Cdkn1a > 0),
    "(", round(100*mean(b17_c$Cdkn1a > 0), 2), "%)\n")

# Same for MBRT_4h (Block_21)
cat("\nMBRT_4h (Block_21): cells with Cdkn1a > 0 across ALL families:\n")
b21 <- meta[Slide == "20250529_214712_S4" & Block == "Block_21"]
b21_c <- counts_df[cell_id %in% b21$cell_id]
cat("Total cells:", nrow(b21), "| with Cdkn1a > 0:", sum(b21_c$Cdkn1a > 0),
    "(", round(100*mean(b21_c$Cdkn1a > 0), 2), "%)\n")

# Distribution of ALL gene counts per block: is SBRT_4h generally lower-count?
cat("\nLibrary size distribution per block (Slide S4):\n")
for (b in c("Block_12", "Block_17", "Block_21")) {
  bcells <- meta[Slide == "20250529_214712_S4" & Block == b, cell_id]
  bc <- counts_df[cell_id %in% bcells]
  libs <- rowSums(bc[, !c("Slide", "fov", "cell_id"), with = FALSE])
  cat(sprintf("%s (n=%d): median lib = %.0f, mean = %.0f\n",
              b, length(libs), median(libs), mean(libs)))
}
