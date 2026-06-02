library(officer)
library(tidyverse)
library(png)

PLOT_DIR <- "plots"
DATA_DIR <- "data"

SLIDE_W <- 13.333
SLIDE_H <- 7.5
MARGIN  <- 0.5
TITLE_H <- 1.2

img_fit <- function(path, max_w = SLIDE_W - 2 * MARGIN,
                    max_h = SLIDE_H - TITLE_H - MARGIN) {
  info <- readPNG(path, info = TRUE)
  img_w <- attr(info, "dim")[2]
  img_h <- attr(info, "dim")[1]
  ar <- img_w / img_h
  w <- max_w
  h <- w / ar
  if (h > max_h) {
    h <- max_h
    w <- h * ar
  }
  left <- (SLIDE_W - w) / 2
  top  <- TITLE_H + (max_h - h) / 2
  list(path = path, w = w, h = h, left = left, top = top)
}

add_image_slide <- function(pptx, title, img_path) {
  pptx <- add_slide(pptx, layout = "Title Only", master = "Office Theme")
  pptx <- ph_with(pptx, value = title,
                   location = ph_location(left = MARGIN, top = 0.2,
                                          width = SLIDE_W - 2 * MARGIN, height = 0.9))
  fit <- img_fit(img_path)
  pptx <- ph_with(pptx, value = external_img(fit$path, width = fit$w, height = fit$h),
                   location = ph_location(left = fit$left, top = fit$top,
                                          width = fit$w, height = fit$h))
  pptx
}

pptx <- read_pptx("/tmp/widescreen_template.pptx")

title_fp <- fp_text(font.size = 32, bold = TRUE, font.family = "Calibri", color = "#2C3E50")
subtitle_fp <- fp_text(font.size = 18, font.family = "Calibri", color = "#7F8C8D")
body_fp <- fp_text(font.size = 16, font.family = "Calibri", color = "#2C3E50")
bold_fp <- fp_text(font.size = 16, font.family = "Calibri", color = "#2C3E50", bold = TRUE)
small_fp <- fp_text(font.size = 13, font.family = "Calibri", color = "#555555")

# --- Slide 1: Title ---
pptx <- add_slide(pptx, layout = "Blank", master = "Office Theme")
pptx <- ph_with(pptx,
  value = block_list(
    fpar(ftext("Microbeam Peak-Valley", title_fp)),
    fpar(ftext("Spatial Transcriptomics", title_fp)),
    fpar(ftext("")),
    fpar(ftext("CosMx single-cell analysis of MBRT dose heterogeneity", subtitle_fp)),
    fpar(ftext("4T1 flank tumors | Mutter_01 dataset", subtitle_fp)),
    fpar(ftext("")),
    fpar(ftext("Jeffrey Szymanski | April 2026", small_fp))
  ),
  location = ph_location(left = 1.5, top = 1.5, width = 10, height = 5))

# --- Slide 2: Experimental Design ---
pptx <- add_slide(pptx, layout = "Blank", master = "Office Theme")
pptx <- ph_with(pptx,
  value = fpar(ftext("Experimental Design", title_fp)),
  location = ph_location(left = MARGIN, top = 0.2, width = 12, height = 0.9))

pptx <- ph_with(pptx,
  value = block_list(
    fpar(ftext("Tumor model", bold_fp), ftext("   4T1 mouse mammary carcinoma (Balb/c flank)", body_fp)),
    fpar(ftext("")),
    fpar(ftext("Platform", bold_fp), ftext("   NanoString CosMx (~950 genes + 21 custom targets)", body_fp)),
    fpar(ftext("")),
    fpar(ftext("Conditions", bold_fp), ftext("   11 total \u2014 MBRT / SBRT / Control \u00d7 5 timepoints", body_fp)),
    fpar(ftext("")),
    fpar(ftext("Models", bold_fp), ftext("   8 flank + 3 tongue conditions", body_fp)),
    fpar(ftext("")),
    fpar(ftext("Replication", bold_fp), ftext("   n = 1 per condition (exploratory; Mutter_02 pending)", body_fp)),
    fpar(ftext("")),
    fpar(ftext("")),
    fpar(ftext("Key comparisons:", bold_fp)),
    fpar(ftext("  \u2022  MBRT vs SBRT at matched timepoints (4h, day 2, day 6)", body_fp)),
    fpar(ftext("  \u2022  MBRT peak vs valley within 4h (H2AX-validated)", body_fp)),
    fpar(ftext("  \u2022  Kinetic trajectories: 1h \u2192 4h \u2192 day 2 \u2192 day 6", body_fp))
  ),
  location = ph_location(left = 1.0, top = 1.2, width = 11, height = 5.5))

# --- Slide 3: Sample Table ---
pptx <- add_slide(pptx, layout = "Blank", master = "Office Theme")
pptx <- ph_with(pptx,
  value = fpar(ftext("Samples and Cell Counts", title_fp)),
  location = ph_location(left = MARGIN, top = 0.2, width = 12, height = 0.9))

qc <- read_tsv(file.path(DATA_DIR, "qc_summary.tsv"), show_col_types = FALSE)
treatments <- c("Control", "MBRT", "MBRT", "MBRT", "MBRT",
                 "SBRT", "SBRT", "SBRT", "MBRT", "MBRT", "SBRT")
timepoints <- c("0h", "1h", "4h", "48h", "144h",
                 "4h", "48h", "144h", "240h", "192h", "240h")
models <- c("Flank", rep("Flank", 7), rep("Tongue", 3))

mono_fp <- fp_text(font.size = 14, font.family = "Courier New", color = "#2C3E50")
mono_bold <- fp_text(font.size = 14, font.family = "Courier New", color = "#2C3E50", bold = TRUE)

header_line <- sprintf("%-24s %-8s %-7s %-7s %10s", "Condition", "Tx", "Time", "Model", "Cells")
table_lines <- list(fpar(ftext(header_line, mono_bold)))

for (i in seq_len(nrow(qc))) {
  line <- sprintf("%-24s %-8s %-7s %-7s %10s",
                   qc$condition[i], treatments[i], timepoints[i],
                   models[i], format(qc$post_filter[i], big.mark = ","))
  table_lines <- c(table_lines, list(fpar(ftext(line, mono_fp))))
}

table_lines <- c(table_lines,
  list(fpar(ftext(""))),
  list(fpar(ftext(sprintf("Total post-QC: %s cells  (from 1.43M raw)",
                          format(sum(qc$post_filter), big.mark = ",")),
                  fp_text(font.size = 15, bold = TRUE, font.family = "Calibri", color = "#E74C3C")))))

pptx <- ph_with(pptx,
  value = do.call(block_list, table_lines),
  location = ph_location(left = 1.5, top = 1.3, width = 10.5, height = 5.5))

# --- Slide 4: QC Violins ---
pptx <- add_image_slide(pptx, "Quality Control: Metric Distributions",
                          file.path(PLOT_DIR, "qc_violins.png"))

# --- Slide 5: UMAP ---
pptx <- add_image_slide(pptx, "Clustering and Cell Type Landscape",
                          file.path(PLOT_DIR, "umap_landscape.png"))

# --- Slide 6: Cell type validation table ---
pptx <- add_slide(pptx, layout = "Title Only", master = "Office Theme")
pptx <- ph_with(pptx,
  value = fpar(ftext("Cell Type Validation Summary", title_fp)),
  location = ph_location(left = MARGIN, top = 0.2, width = 12, height = 0.9))

val_header <- sprintf("%-22s %-30s %-10s", "Cell Type", "Expected Marker(s)", "Detected")
val_rows <- list(
  sprintf("%-22s %-30s %-10s", "4T1 Tumor",   "Krt8, Krt18, Epcam, Cdh1",  "Yes"),
  sprintf("%-22s %-30s %-10s", "Macrophage",   "Adgre1 (F4/80), Cd68",      "Yes"),
  sprintf("%-22s %-30s %-10s", "CD8 T cell",   "Cd8a, Cd3e",                "Yes"),
  sprintf("%-22s %-30s %-10s", "B cell",        "Cd19, Ms4a1 (CD20)",       "Yes"),
  sprintf("%-22s %-30s %-10s", "Dendritic",     "Itgax (CD11c), H2-Aa",     "Yes"),
  sprintf("%-22s %-30s %-10s", "Neutrophil",    "Ly6g, S100a8",             "Yes"),
  sprintf("%-22s %-30s %-10s", "Fibroblast",    "Col1a1, Dcn",              "Yes"),
  sprintf("%-22s %-30s %-10s", "Endothelial",   "Pecam1, Cdh5, Vwf",       "Yes"),
  sprintf("%-22s %-30s %-10s", "Pericyte",      "Pdgfrb, Acta2",            "Yes")
)

val_lines <- list(fpar(ftext(val_header, mono_bold)))
for (row in val_rows)
  val_lines <- c(val_lines, list(fpar(ftext(row, mono_fp))))
val_lines <- c(val_lines,
  list(fpar(ftext(""))),
  list(fpar(ftext("All major cell types confirmed by canonical marker expression in CosMx panel",
    fp_text(font.size = 14, font.family = "Calibri", color = "#27AE60", italic = TRUE)))))

pptx <- ph_with(pptx,
  value = do.call(block_list, val_lines),
  location = ph_location(left = 1.5, top = 1.3, width = 10.5, height = 5.5))

# --- Slide 7: Spatial cell class maps ---
pptx <- add_image_slide(pptx, "Spatial Cell Class Maps (Flank)",
                          file.path(PLOT_DIR, "spatial_cellclass_flank.png"))

# --- Slide 8: Peak/Valley Classification ---
pptx <- add_slide(pptx, layout = "Blank", master = "Office Theme")
pptx <- ph_with(pptx,
  value = fpar(ftext("Peak/Valley Classification: Manual FOV Annotation", title_fp)),
  location = ph_location(left = MARGIN, top = 0.2, width = 12, height = 0.9))

fit_overlay <- img_fit(file.path(PLOT_DIR, "h2ax_fov_overlay.png"),
                       max_w = 7.5, max_h = 5.5)
pptx <- ph_with(pptx,
  value = external_img(fit_overlay$path, width = fit_overlay$w, height = fit_overlay$h),
  location = ph_location(left = 0.5, top = fit_overlay$top,
                          width = fit_overlay$w, height = fit_overlay$h))

pptx <- ph_with(pptx,
  value = block_list(
    fpar(ftext("Method", bold_fp)),
    fpar(ftext("")),
    fpar(ftext("1. Gamma-H2AX IHC of serial section", body_fp)),
    fpar(ftext("   reveals MBRT beam stripe pattern", body_fp)),
    fpar(ftext("")),
    fpar(ftext("2. CosMx FOV grid overlaid on IHC", body_fp)),
    fpar(ftext("   image to identify FOV positions", body_fp)),
    fpar(ftext("")),
    fpar(ftext("3. FOVs classified as peak or valley", body_fp)),
    fpar(ftext("   based on H2AX signal intensity", body_fp)),
    fpar(ftext("")),
    fpar(ftext("4. Cells inherit zone from their FOV", body_fp)),
    fpar(ftext("")),
    fpar(ftext("MBRT 4h only — H2AX clears by day 2", small_fp))
  ),
  location = ph_location(left = 8.5, top = 1.3, width = 4.3, height = 5.5))

# --- Save ---
out_path <- "presentation.pptx"
print(pptx, target = out_path)
cat("Saved:", out_path, "\n")
