# Programme A (ferroptosis) readable proxies, MBRT 4h tumor, peak vs valley (locked labels).
# Execution arm (ACSL4/ALOX15/GPX4) is OFF-panel; readable = buffer (Gpx1, Fasn) + PUFA-handling (Fabp5,Cd36,Hilpda,Ldlr).
suppressPackageStartupMessages({library(Seurat); library(data.table); library(arrow)})

o <- readRDS("/mnt/data/projects/spatial-rads/processing/scored/sam0003.scored.rds")
genes <- intersect(c("Gpx1","Fasn","Fabp5","Cd36","Hilpda","Ldlr"), rownames(o))
expr <- FetchData(o, vars=genes, layer="data")          # LogNormalized, rows = Cells(o)
md <- o@meta.data
key <- if ("cell_id" %in% colnames(md)) as.character(md$cell_id) else rownames(md)
dt <- data.table(cell_id=key, as.data.table(expr))
cat("scored key head:", paste(head(key,2), collapse=" | "), "\n")

pv <- fread("dev/peak_valley_analysis/data/mbrt4h_peak_valley.tsv", select=c("cell_id","zone"))
cat("pv key head:   ", paste(head(pv$cell_id,2), collapse=" | "), "\n")

fl <- as.data.table(read_parquet("results/aggregate/full_labels.parquet", col_select=c("cell","compartment")))
fl <- fl[startsWith(cell,"sam0003_")][, cell_id := sub("^sam0003_","",cell)]

m <- merge(dt, pv, by="cell_id")
m <- merge(m, fl[,.(cell_id,compartment)], by="cell_id", all.x=TRUE)
cat("matched to zone:", nrow(m), "| with compartment:", sum(!is.na(m$compartment)), "\n")

tum <- m[compartment=="tumor"]
cat("\n=== MBRT 4h TUMOR (locked labels): mean LogNorm expression by zone ===\n")
print(tum[, c(.(n=.N), lapply(.SD, function(x) round(mean(x),4))), by=zone, .SDcols=genes])

if (all(c("Gpx1","Fasn") %in% genes)) {
  tum[, defense := as.numeric(scale(Gpx1)) + as.numeric(scale(Fasn))]
  cat("\ndefense composite z(Gpx1)+z(Fasn) by zone (higher = more shielded / ferroptosis buffered):\n")
  print(tum[, .(n=.N, defense=round(mean(defense),4)), by=zone])
}
cat("\nDONE\n")
