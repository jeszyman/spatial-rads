#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# panel_provenance.R
# Partition each dataset's delivered CosMx gene panel into UCC-standard vs custom
# by diffing the actual object gene lists against the Bruker Mouse UCC vendor gene
# list, and mark the cross-dataset common set. Ground truth = the data objects +
# the vendor panel doc, NOT the hand-entered metadata panel fields.
# Args: <vendor_ucc.xlsx> <m01_rds_dir> <m02_rds_dir> <out_membership.tsv> <out_summary.tsv>
# -----------------------------------------------------------------------------

packages <- c("Seurat", "readxl", "data.table")
suppressPackageStartupMessages(invisible(lapply(packages, require, character.only = TRUE)))

args    <- commandArgs(trailingOnly = TRUE)
vendor  <- args[1]
m01_dir <- args[2]
m02_dir <- args[3]
out_mem <- args[4]
out_sum <- args[5]

#' Genes present in a dataset's delivered per-slide object (RNA assay rownames).
dataset_genes <- function(dir, pat) {
  f <- list.files(dir, pattern = pat, full.names = TRUE)[1]
  rownames(UpdateSeuratObject(readRDS(f))[["RNA"]])
}

g01 <- dataset_genes(m01_dir, "Mutter_01.*RDS$")
g02 <- dataset_genes(m02_dir, "Mutter_02.*RDS$")
obj_genes <- union(g01, g02)

# --- vendor UCC standard set: auto-detect the gene-symbol column (max overlap) ---
gp   <- suppressMessages(read_excel(vendor, sheet = "Gene and Probe Details", skip = 1))
ov   <- vapply(gp, function(col) sum(as.character(col) %in% obj_genes), integer(1))
ucc  <- unique(na.omit(as.character(gp[[which.max(ov)]])))
ucc  <- ucc[ucc %in% obj_genes | grepl("^[A-Za-z0-9.-]+$", ucc)]  # gene-like tokens
stopifnot(max(ov) > 500)  # sanity: the detected column really is genes

# Drop platform control probes (negprobes / system controls) -- these are not genes and
# otherwise inflate the standard count (the Mouse UCC is a 1000-plex + 10 negprobes) and
# leak into the "custom" set (e.g. NegativeAdd).
is_control <- function(x) {
  grepl("^Negative[0-9]+$|^NegativeAdd$|^SystemControl|^NegPrb|^FalseCode", x, ignore.case = TRUE)
}
g01 <- g01[!is_control(g01)]
g02 <- g02[!is_control(g02)]
ucc <- ucc[!is_control(ucc)]

# --- per-gene membership across UCC standard + both datasets ---
all_g <- sort(unique(c(ucc, g01, g02)))
mem <- data.table(
  gene         = all_g,
  ucc_standard = all_g %in% ucc,
  mutter_01    = all_g %in% g01,
  mutter_02    = all_g %in% g02)
mem[, class := fifelse(mutter_01 & mutter_02, "common",
                fifelse(mutter_01 | mutter_02, "dataset_specific", "ucc_only"))]
mem[, custom := (mutter_01 | mutter_02) & !ucc_standard]   # in a dataset but not UCC standard
fwrite(mem, out_mem, sep = "\t")

# --- summary counts ---
sm <- data.table(
  metric = c("ucc_standard", "mutter_01_total", "mutter_02_total", "common",
             "mutter_01_custom", "mutter_02_custom",
             "mutter_01_missing_from_ucc", "mutter_02_missing_from_ucc"),
  n = c(length(ucc), length(g01), length(g02), length(intersect(g01, g02)),
        sum(mem$mutter_01 & mem$custom), sum(mem$mutter_02 & mem$custom),
        length(setdiff(g01, ucc)), length(setdiff(g02, ucc))))
fwrite(sm, out_sum, sep = "\t")
print(sm)
cat(sprintf("panel_provenance: UCC standard %d | M01 %d (custom %d) | M02 %d (custom %d) | common %d\n",
            length(ucc), length(g01), sum(mem$mutter_01 & mem$custom),
            length(g02), sum(mem$mutter_02 & mem$custom), length(intersect(g01, g02))))
