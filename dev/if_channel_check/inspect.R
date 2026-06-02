library(arrow)
library(data.table)

meta <- as.data.table(read_parquet("/tmp/mutter01_meta.parquet"))
cat("All columns:\n")
print(colnames(meta))

cat("\nSample per-slide:\n")
print(meta[, .N, by = Slide])

cat("\nCell type counts (Main_cell_Types):\n")
print(meta[, .N, by = ImmuneAtlas_ImmGen_Main_cell_Types][order(-N)])

cat("\nSample FOV counts per slide (first 3):\n")
print(meta[, .N, by = .(Slide, fov)][order(Slide, fov)][, head(.SD, 3), by = Slide])
