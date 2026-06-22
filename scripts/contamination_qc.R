#!/usr/bin/env Rscript
# Contamination QC (report-only): SpatialQM MECR marker-bleed metric at sample + FOV grain.
# MECR = mean Jaccard co-expression of cross-lineage marker pairs (LOW = clean, HIGH = bleed /
# segmentation contamination). FastReseg / per-cell resegmentation is infeasible here (no transcript
# coords or cell polygons delivered), so MECR -- the feasible peer-reviewed metric (SpatialQM,
# Plummer/Segato-Dezem et al., Nat Biotechnol 2025; vendored in scripts/aggregate/spatialqm_metrics.R)
# -- is computed per sample and per FOV. The per-FOV table carries orthogonal triage columns so a
# high-MECR flag separates segmentation contamination (high MECR, normal counts/background) from
# necrosis / low-quality tissue (high MECR + LOW counts + HIGH propNeg) from genuine mixed biology
# (high MECR, normal QC). REPORT-ONLY -- never excludes cells/FOVs (necrosis is itself radiation signal).
# Args: <qc_dir> <spatialqm_metrics.R> <lineage_markers.yaml> <out_sample.tsv> <out_fov.tsv>
suppressMessages({library(Seurat); library(data.table); library(yaml); library(Matrix)})
a <- commandArgs(trailingOnly = TRUE)
QC_DIR <- a[1]; SQM <- a[2]; MARKERS <- a[3]; OUT_SAMPLE <- a[4]; OUT_FOV <- a[5]
MIN_FOV_CELLS <- 50L                                   # skip tiny FOVs -- the Jaccard rate is unstable
source(SQM)                                            # provides sqm_mecr(counts, marker_df, max_markers, seed)

## marker_df: flatten lineage_markers.yaml (Lineage: [genes]) to a long (gene, cell_type) data.frame
ly <- yaml::read_yaml(MARKERS)
marker_df <- rbindlist(lapply(names(ly), function(ct) data.table(gene = as.character(ly[[ct]]), cell_type = ct)))
marker_df <- as.data.frame(unique(marker_df[!is.na(gene) & nzchar(gene)][!duplicated(gene)]))

files <- sort(list.files(QC_DIR, pattern = "\\.qc\\.rds$", full.names = TRUE)); stopifnot(length(files) > 0)
samp <- list(); fov <- list()
for (f in files) {
  o  <- readRDS(f); md <- o@meta.data
  cm <- LayerData(o, assay = "RNA", layer = "counts")
  sid <- as.character(md$sample_id[1]); ds <- as.character(md$dataset[1])
  samp[[f]] <- data.table(sample_id = sid, dataset = ds, n_cells = ncol(cm), mecr = sqm_mecr(cm, marker_df))
  idx_by_fov <- split(seq_len(ncol(cm)), md$fov)
  rows <- lapply(names(idx_by_fov), function(fv) {
    idx <- idx_by_fov[[fv]]
    if (length(idx) < MIN_FOV_CELLS) return(NULL)
    data.table(sample_id = sid, dataset = ds, fov = as.integer(fv), n_cells = length(idx),
               mecr         = sqm_mecr(cm[, idx, drop = FALSE], marker_df),   # off-label per-FOV: may be NaN if a sampled pair has zero union
               med_nCount   = as.numeric(median(md$nCount_RNA[idx])),
               med_nFeature = as.numeric(median(md$nFeature_RNA[idx])),
               mean_propNeg = mean(md$propNegative[idx]),
               med_area     = as.numeric(median(md$Area[idx])))
  })
  fov[[f]] <- rbindlist(rows)
  cat(sprintf("%-28s %-10s n=%7d  mecr=%.3f  fovs=%d\n", sid, ds, ncol(cm), samp[[f]]$mecr, length(rows)))
  rm(o, cm); invisible(gc())
}
sample_tab <- rbindlist(samp); fov_tab <- rbindlist(fov)

## per-FOV 3-MAD outlier flag on MECR, within dataset (mirrors the falsecode FOV flag in control_sidecar.R)
fov_tab[, `:=`(med = median(mecr, na.rm = TRUE), mad = mad(mecr, na.rm = TRUE)), by = dataset]
fov_tab[, mecr_flag := !is.na(mecr) & mecr > med + 3 * mad][, c("med", "mad") := NULL]

dir.create(dirname(OUT_SAMPLE), recursive = TRUE, showWarnings = FALSE)
fwrite(sample_tab, OUT_SAMPLE, sep = "\t"); fwrite(fov_tab, OUT_FOV, sep = "\t")
cat(sprintf("contamination QC: %d samples; %d FOVs (%d valid MECR), %d flagged (report-only)\n",
            nrow(sample_tab), nrow(fov_tab), sum(!is.na(fov_tab$mecr)), sum(fov_tab$mecr_flag)))
