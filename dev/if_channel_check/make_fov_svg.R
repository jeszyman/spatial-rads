library(data.table)

args <- commandArgs(trailingOnly = TRUE)
vflip <- "--vflip" %in% args
out_path <- if (vflip) "/tmp/mutter01_s4_mbrt4h_fovs_vflip.svg" else "/tmp/mutter01_s4_mbrt4h_fovs.svg"

fov <- fread("/tmp/mutter01_s4_mbrt4h_fov_boxes.tsv")
if (vflip) {
  ymax_global <- max(fov$ymax)
  ymin_global <- min(fov$ymin)
  fov[, `:=`(
    ymid = ymax_global + ymin_global - ymid,
    ymin_new = ymax_global + ymin_global - ymax,
    ymax_new = ymax_global + ymin_global - ymin
  )]
  fov[, `:=`(ymin = ymin_new, ymax = ymax_new, ymin_new = NULL, ymax_new = NULL)]
}

# Use cell-extent centroid; box is fixed 510 x 510 um (CosMx spec)
box_um <- 510
fov[, `:=`(
  box_xmin_mm = xmid - box_um / 2000,
  box_xmax_mm = xmid + box_um / 2000,
  box_ymin_mm = ymid - box_um / 2000,
  box_ymax_mm = ymid + box_um / 2000
)]

# SVG viewport: overall bbox of all FOVs, with small margin
margin <- 0.5
vb_xmin <- min(fov$box_xmin_mm) - margin
vb_xmax <- max(fov$box_xmax_mm) + margin
vb_ymin <- min(fov$box_ymin_mm) - margin
vb_ymax <- max(fov$box_ymax_mm) + margin
w <- vb_xmax - vb_xmin
h <- vb_ymax - vb_ymin

# Write SVG in mm
svg <- c(
  sprintf('<?xml version="1.0" encoding="UTF-8"?>'),
  sprintf('<svg xmlns="http://www.w3.org/2000/svg" version="1.1"'),
  sprintf('  width="%.3fmm" height="%.3fmm"', w, h),
  sprintf('  viewBox="%.4f %.4f %.4f %.4f">', vb_xmin, vb_ymin, w, h),
  '  <style>',
  '    .fov { fill: none; stroke: #e91e63; stroke-width: 0.012; stroke-opacity: 0.7; }',
  '    .lbl { font-family: sans-serif; font-size: 0.12px; fill: #e91e63; text-anchor: middle; dominant-baseline: central; font-weight: bold; }',
  '  </style>',
  sprintf('  <g id="fovs">')
)

for (i in seq_len(nrow(fov))) {
  svg <- c(svg, sprintf(
    '    <rect class="fov" x="%.4f" y="%.4f" width="%.4f" height="%.4f" data-fov="%d" />',
    fov$box_xmin_mm[i], fov$box_ymin_mm[i],
    box_um / 1000, box_um / 1000, fov$fov[i]
  ))
  svg <- c(svg, sprintf(
    '    <text class="lbl" x="%.4f" y="%.4f">%d</text>',
    fov$xmid[i], fov$ymid[i], fov$fov[i]
  ))
}

svg <- c(svg, "  </g>", "</svg>")

writeLines(svg, out_path)
cat("Wrote", out_path, "with", nrow(fov), "FOVs\n")
cat("SVG size (mm):", round(w, 2), "x", round(h, 2), "\n")
cat("FOV size: 510 x 510 um (fixed per CosMx spec)\n")
