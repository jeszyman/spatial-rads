#!/usr/bin/env Rscript
# One-time migration: flat metadata.xlsx (Sheet1) -> relational 4-sheet metadata.xlsx
# conforming to the sci_proj_repo Data layout standard (science.org 08049f58).
# Preview:  conda run -n basecamp Rscript scripts/build_metadata_xlsx.R
# Write:    conda run -n basecamp Rscript scripts/build_metadata_xlsx.R --write
suppressMessages({library(readxl); library(dplyr); library(tibble); library(tidyr); library(writexl)})

WRITE <- "--write" %in% commandArgs(trailingOnly = TRUE)
SRC   <- "data/sources/2026-05-29-legacy-flat-metadata.xlsx"  # legacy flat metadata (write-protected source)
XLSX  <- "data/metadata.xlsx"                                  # relational output, built from SRC
INPUT <- "/mnt/data/projects/spatial-rads/inputs"

flat <- read_excel(SRC)  # 11 Mutter_01 rows (legacy flat)

## datasets
datasets <- tibble(
  dataset_id    = c("dat0001", "dat0002"),
  name          = c("Mutter_01", "Mutter_02"),
  platform      = "CosMx",
  panel_n       = c(1000L, 972L),
  provider      = c("Yi Liu / Mutter Lab", "Mutter Lab"),
  format        = c("parquet", "rds"),
  counts_path   = c(file.path(INPUT, "mutter01/projects_Mutter_01_CosMmR_app_data_raw_tx_counts_matrix.parquet"), NA),
  metadata_path = c(file.path(INPUT, "mutter01/Analysis_Mutter_01_CosMmR_Mutter_updated_metadata.parquet"), NA)
)

## slides (4 Mutter_01 physical slides + 4 Mutter_02)
m01_barcodes <- sort(unique(flat$slide))
slides_m01 <- tibble(
  slide_id         = sprintf("sld%04d", seq_along(m01_barcodes)),
  dataset_id       = "dat0001",
  physical_barcode = m01_barcodes,
  input_path       = NA_character_  # monolithic parquet (see datasets.counts_path)
)
slides_m02 <- tibble(
  slide_id         = sprintf("sld%04d", 4 + 1:4),
  dataset_id       = "dat0002",
  physical_barcode = sprintf("Mutter_02_slide%02d", 1:4),
  input_path       = file.path(INPUT, sprintf("mutter02/seuratObject_%02d_Mutter_02_CosMmR.RDS", 1:4))
)
slides <- bind_rows(slides_m01, slides_m02)
barcode2sld <- setNames(slides_m01$slide_id, slides_m01$physical_barcode)

## Mutter_02 sample layout (Jenn Fazzari 2026-04-09; from 00_load_data_mutter02.R:40-62)
## rank by descending y_slide_mm: 1=top=Control, 2=mid=MBRT, 3=bottom=SBRT
m02 <- tribble(
  ~slide_num, ~rank, ~region_id, ~treatment,
  1, 1, 1,  "Control",
  1, 2, 34, "MBRT",
  1, 3, 28, "SBRT",
  2, 1, 2,  "Control",
  2, 2, 35, "MBRT",
  2, 3, 30, "SBRT",
  3, 1, 3,  "Control",
  3, 2, 19, "MBRT",
  3, 3, 16, "SBRT",
  4, 1, 1,  "Control",
  4, 2, 20, "MBRT",
  4, 3, 18, "SBRT"
)
cond_map_m02 <- c(Control = "Control", MBRT = "MBRT_day2", SBRT = "SBRT_day2")
dose_m02     <- c(Control = 0, MBRT = NA_real_, SBRT = 20)  # SBRT 20 Gy per layout note; MBRT peak dose TBD

## samples
samples_m01 <- flat %>% transmute(
  block_label = block_id, condition, treatment, timepoint_h,
  model, tumor_cell_line, tumor_type, strain = host_strain,
  slide_id = unname(barcode2sld[slide]),
  dose_gy = NA_real_,            # Mutter_01 doses not yet recorded
  extract_key = condition        # Condition field that splits the Mutter_01 parquet
) %>% mutate(  # flat-source Control row has blank biology; it is the flank 4T1 baseline
  model           = coalesce(model,           if_else(treatment == "NT", "flank", NA_character_)),
  tumor_cell_line = coalesce(tumor_cell_line, if_else(treatment == "NT", "4T1", NA_character_)),
  tumor_type      = coalesce(tumor_type,      if_else(treatment == "NT", "mammary carcinoma, TNBC/claudin-low", NA_character_)),
  strain          = coalesce(strain,          if_else(treatment == "NT", "Balb/c", NA_character_))
)
samples_m02 <- m02 %>%
  mutate(  # derive dose/condition from raw label, THEN harmonize control to NT (matches Mutter_01)
    dose_gy   = unname(dose_m02[treatment]),
    condition = unname(cond_map_m02[treatment]),
    treatment = if_else(treatment == "Control", "NT", treatment)
  ) %>% transmute(
    block_label = paste0("s", slide_num, "_region", region_id),
    condition, treatment, timepoint_h = 48,
    model = "flank", tumor_cell_line = "4T1",
    tumor_type = "mammary carcinoma, TNBC/claudin-low", strain = "Balb/c",
    slide_id = slides_m02$slide_id[slide_num],
    dose_gy,
    extract_key = as.character(rank)  # y-band rank within the per-slide RDS
  )
samples_all <- bind_rows(samples_m01, samples_m02) %>%
  mutate(sample_id = sprintf("sam%04d", row_number()),
         mouse_id  = sprintf("mou%04d", row_number()))

## mice (one per sample; n=1 per condition)
mice <- samples_all %>% transmute(mouse_id, strain, sex = NA_character_, notes = NA_character_)

## final samples sheet (strain lives in mice)
samples <- samples_all %>% select(
  sample_id, mouse_id, slide_id, block_label, condition, treatment,
  dose_gy, timepoint_h, model, tumor_cell_line, tumor_type, extract_key
)

cat("=== datasets ===\n"); print(as.data.frame(datasets))
cat("\n=== slides ===\n");  print(as.data.frame(slides))
cat("\n=== mice ===\n");    print(as.data.frame(mice))
cat("\n=== samples ===\n"); print(as.data.frame(samples))

if (WRITE) {
  bak <- "(none)"
  if (file.exists(XLSX)) {
    bak <- sprintf("data/metadata-prev-%s.xlsx.bak", format(Sys.time(), "%Y%m%d-%H%M%S"))
    file.copy(XLSX, bak, overwrite = TRUE)
  }
  write_xlsx(list(mice = mice, datasets = datasets, slides = slides, samples = samples), XLSX)
  cat(sprintf("\nWROTE %s (prev backup: %s)\n", XLSX, bak))
} else {
  cat("\n[preview only -- rerun with --write to write metadata.xlsx]\n")
}
