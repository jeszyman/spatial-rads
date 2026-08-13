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
  panel_base    = "Mouse UCC",                                       # Mouse RNA Universal Cell Characterization
  panel_custom  = c("Alan Fields 21-gene", "1K_Mayo_Thomp-Fields1"), # same Fields add-on; M02 label is the in-object Panel build name
  assay_type    = "RNA",                                             # RNA-only; no protein-expression assay on either run
  n_negprobe    = c(10L, 10L),           # both datasets carry the 10 negprobe assays in the raw per-slide RDS
  n_falsecode   = c(184L, 184L),         # both carry the 184 falsecode assays in the raw per-slide RDS
  provider      = c("Yi Liu / Mutter Lab", "Mutter Lab"),
  format        = c("rds", "rds"),       # M01 re-based onto per-slide raw RDS (Yi 2026-06-10); counts come from slides.input_path
  counts_path   = c(NA, NA),             # counts now sourced from the per-slide RDS, not a monolithic counts parquet
  metadata_path = c(file.path(INPUT, "mutter01/Analysis_Mutter_01_CosMmR_Mutter_updated_metadata.parquet"), NA)  # M01: sole source of Condition + Yi labels
)

## slides (4 Mutter_01 physical slides + 4 Mutter_02)
m01_barcodes <- sort(unique(flat$slide))
slides_m01 <- tibble(
  slide_id         = sprintf("sld%04d", seq_along(m01_barcodes)),
  dataset_id       = "dat0001",
  physical_barcode = m01_barcodes,
  input_path       = file.path(INPUT, sprintf("mutter01/seuratObject_%02d_Mutter_01_CosMmR.RDS", seq_along(m01_barcodes)))  # per-slide raw RDS; sorted barcode _S<i> -> seuratObject_0<i>
)
stopifnot(grepl("_S\\d+$", m01_barcodes))  # positional barcode->RDS mapping holds only if barcodes carry the _S<i> suffix
slides_m02 <- tibble(
  slide_id         = sprintf("sld%04d", 4 + 1:4),
  dataset_id       = "dat0002",
  physical_barcode = sprintf("Mutter_02_slide%02d", 1:4),
  input_path       = file.path(INPUT, sprintf("mutter02/seuratObject_%02d_Mutter_02_CosMmR.RDS", 1:4))
)
slides <- bind_rows(slides_m01, slides_m02)
barcode2sld <- setNames(slides_m01$slide_id, slides_m01$physical_barcode)

## Mutter_02 sample layout: block-to-slide positions from Jenn Fazzari (2026-04-09),
## timepoints from the yH2AX layout doc (slides 1-2 = 2d, slides 3-4 = 4h).
## rank by descending y_slide_mm: 1=top=Control, 2=mid=MBRT, 3=bottom=SBRT
LAYOUT_CSV <- "data/sources/mutter02_block_layout.csv"
layout     <- read.csv(LAYOUT_CSV)
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
## Join timepoint from the layout CSV (keyed on block + slide)
m02 <- m02 %>%
  left_join(layout %>% select(block, slide, timepoint_h),
            by = c("region_id" = "block", "slide_num" = "slide"))
stopifnot(!anyNA(m02$timepoint_h),
          length(unique(m02$timepoint_h)) > 1)
stopifnot(all(layout$block %in% m02$region_id))
## Condition suffix encodes timepoint: MBRT_4h / MBRT_day2 etc; Control stays Control
cond_suffix  <- ifelse(m02$timepoint_h == 4, "_4h", "_day2")
cond_map_m02 <- ifelse(m02$treatment == "Control", "Control",
                       paste0(m02$treatment, cond_suffix))
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
  mutate(
    dose_gy   = unname(dose_m02[treatment]),
    condition = cond_map_m02,
    treatment = if_else(treatment == "Control", "NT", treatment)
  ) %>% transmute(
    block_label = paste0("s", slide_num, "_region", region_id),
    condition, treatment, timepoint_h,
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

## if_channels: CosMx morphology/IF stains used for segmentation + lineage (NOT a protein-expression assay).
## Both datasets verified to carry the identical 5 channels (M01 parquet IF cols; M02 RDS meta.data),
## so spec x datasets via crossing(); append divergent rows by hand if a future run adds a channel (e.g. H2AX).
if_spec <- tribble(
  ~channel,    ~target,           ~role,                   ~meta_cols,
  "DAPI",      "DNA",             "segmentation_nuclear",  "Mean.DAPI;Max.DAPI",
  "CD298.B2M", "B2M/CD298",       "segmentation_membrane", "Mean.CD298.B2M;Max.CD298.B2M",
  "PanCK",     "pan-cytokeratin", "lineage_epithelial",    "Mean.PanCK;Max.PanCK",
  "CD45",      "PTPRC",           "lineage_immune",        "Mean.CD45;Max.CD45",
  "G",         "unspecified",     "redundant_epithelial",  "Mean.G;Max.G"
)
if_channels <- tidyr::crossing(dataset_id = datasets$dataset_id, if_spec) %>%
  mutate(notes = if_else(channel == "G", "r~0.92 with Mean.PanCK; not a literal duplicate", NA_character_)) %>%
  select(dataset_id, channel, target, role, meta_cols, notes)

cat("=== datasets ===\n");    print(as.data.frame(datasets))
cat("\n=== slides ===\n");     print(as.data.frame(slides))
cat("\n=== mice ===\n");       print(as.data.frame(mice))
cat("\n=== samples ===\n");    print(as.data.frame(samples))
cat("\n=== if_channels ===\n"); print(as.data.frame(if_channels))

if (WRITE) {
  bak <- "(none)"
  if (file.exists(XLSX)) {
    bak <- sprintf("data/metadata-prev-%s.xlsx.bak", format(Sys.time(), "%Y%m%d-%H%M%S"))
    file.copy(XLSX, bak, overwrite = TRUE)
  }
  write_xlsx(list(mice = mice, datasets = datasets, slides = slides, samples = samples, if_channels = if_channels), XLSX)
  cat(sprintf("\nWROTE %s (prev backup: %s)\n", XLSX, bak))
} else {
  cat("\n[preview only -- rerun with --write to write metadata.xlsx]\n")
}
