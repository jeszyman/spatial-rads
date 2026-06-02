library(arrow)
library(data.table)

meta <- as.data.table(read_parquet("/tmp/mutter01_meta.parquet"))

s4 <- meta[Slide == "20250529_214712_S4"]
cat("\n=== Block distribution on S4 ===\n")
print(s4[, .(n_cells = .N, fov_min = min(fov), fov_max = max(fov), n_fovs = length(unique(fov))), by = Block])

# Focus on MBRT_4h: Block_21 per metadata.xlsx
mbrt4h_fovs <- s4[Block == "Block_21"]
cat("\nMBRT_4h (Block_21) FOVs:", length(unique(mbrt4h_fovs$fov)),
    "range:", range(mbrt4h_fovs$fov), "\n")

# Per-FOV bounding box in slide mm
fov_boxes <- mbrt4h_fovs[, .(
  xmin = min(x_slide_mm),
  xmax = max(x_slide_mm),
  ymin = min(y_slide_mm),
  ymax = max(y_slide_mm),
  xmid = mean(range(x_slide_mm)),
  ymid = mean(range(y_slide_mm)),
  width_um  = 1000 * (max(x_slide_mm) - min(x_slide_mm)),
  height_um = 1000 * (max(y_slide_mm) - min(y_slide_mm)),
  n_cells = .N
), by = fov][order(fov)]

cat("\n=== MBRT_4h FOV summary ===\n")
print(summary(fov_boxes$width_um))
print(summary(fov_boxes$height_um))
print(head(fov_boxes, 3))

fwrite(fov_boxes, "/tmp/mutter01_s4_mbrt4h_fov_boxes.tsv", sep = "\t")
cat("\nWrote /tmp/mutter01_s4_mbrt4h_fov_boxes.tsv with", nrow(fov_boxes), "FOVs\n")
