#!/usr/bin/env Rscript
# Peak-vs-valley differential expression at 4h post-MBRT (sam0003, the
# MBRT_4h arm), plus two registration-sensitivity validation reruns and a
# continuous distance-to-peak regression. Companion to
# m01_4h_sample_level.R (whole-slide arm contrasts); this script instead
# splits the single MBRT_4h slide into peak/valley zones using the fitted
# stripe model and tests peak vs valley within-slide, paired on FOV. n=1
# slide (descriptive); FOV is the pseudo-replicate unit for the paired
# design, matching the project's established `~ FOV + zone` DESeq2
# precedent (CLAUDE.md, rotation-null reconcile, 2026-07-13).
suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(Seurat)
  library(DESeq2)
  library(limma)
  library(edgeR)
  library(qvalue)
  library(yaml)
  library(IHW)
})

DATADIR    <- "/mnt/data/projects/spatial-rads/processing/norm"
LABELS     <- "results/aggregate/full_labels.parquet"
STRIPE     <- "dev/peak_valley_analysis/data/stripe_model.rds"
OUTDIR     <- "results/aggregate/m01_4h"
SAMPLE_ID  <- "sam0003"  # MBRT_4h
PEAK_THR   <- 0.10   # mm, core band peak
VALLEY_THR <- 0.40   # mm, core band valley
MIN_CELLS  <- 20     # minimum cells per FOV-zone pseudobulk unit

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# --- Hypothesis-restricted FDR ---
# Map coarse compartments to the cell subtypes in hypotheses.yaml so we can
# tag genes and recompute BH within the confirmatory set per compartment.
HYPO_FILE <- "config/hypotheses.yaml"
hypo_raw  <- read_yaml(HYPO_FILE)

COMP_TO_SUBTYPES <- list(
  tumor  = c("Tumor"),
  stroma = c("Fibroblast", "SmoothMuscle", "Adipocyte", "Endothelial"),
  immune = c("T cells", "NK cells", "ILC", "Plasma cells", "Macrophages",
             "DC", "Mast cells", "Neutrophils")
)

build_confirmatory_genes <- function(hypo_list, subtypes) {
  genes <- character(0)
  for (h in hypo_list$hypotheses) {
    if (any(h$cell_types %in% subtypes)) {
      for (prog in h$programs) {
        genes <- c(genes, prog$genes)
      }
    }
  }
  unique(genes)
}

tag_confirmatory <- function(dt, comp) {
  if (comp == "all") {
    subtypes <- unlist(COMP_TO_SUBTYPES, use.names = FALSE)
  } else {
    subtypes <- COMP_TO_SUBTYPES[[comp]]
  }
  if (is.null(subtypes)) {
    dt[, `:=`(confirmatory = FALSE, padj_confirmatory = NA_real_)]
    return(dt)
  }
  conf_genes <- build_confirmatory_genes(hypo_raw, subtypes)
  dt[, confirmatory := gene %in% conf_genes]
  conf_p <- dt[confirmatory == TRUE & !is.na(p), p]
  dt[, padj_confirmatory := NA_real_]
  if (length(conf_p) > 0) {
    dt[confirmatory == TRUE & !is.na(p),
       padj_confirmatory := p.adjust(p, method = "BH")]
  }
  dt
}

DDR_EXCLUDE <- character(0)

compartments <- c("all", "tumor", "stroma", "immune")

# Perpendicular distance from a cell to the nearest of the stripe model's
# beam centers, for an arbitrary (possibly offset) set of centers. Shared by
# the primary zone call and both registration-sensitivity validations.
dist_to_nearest_beam <- function(d_perp, centers) {
  do.call(pmin, lapply(centers, function(c) abs(d_perp - c)))
}

# --- Load data (once; reused across all parts below) ---
cat("Loading data...\n")
obj <- readRDS(file.path(DATADIR, paste0(SAMPLE_ID, ".norm.rds")))
lab <- as.data.table(read_parquet(LABELS))
lab <- lab[, .(cell, compartment)]
stripe <- readRDS(STRIPE)

# --- Compute distance to peak ---
# full_labels.parquet keys cells as "<sample_id>_<barcode>"; the per-sample
# norm.rds Seurat object's colnames are the bare barcode. Reconstruct the
# prefixed id only to join labels -- keep "cell" as the bare barcode (==
# the Seurat object's colnames) for indexing the count/data matrices
# everywhere else in this script.
meta <- as.data.table(obj@meta.data, keep.rownames = "cell")
meta[, cell_full := paste0(SAMPLE_ID, "_", cell)]
meta <- merge(meta, lab, by.x = "cell_full", by.y = "cell", all.x = TRUE)
stopifnot(sum(is.na(meta$compartment)) == 0)

theta <- stripe$tilt_deg * pi / 180
d_perp <- -sin(theta) * meta$x_slide_mm + cos(theta) * meta$y_slide_mm
meta[, dist_to_peak := dist_to_nearest_beam(d_perp, stripe$beam_centers)]

# Zone assignment (core bands)
meta[, zone := fifelse(dist_to_peak < PEAK_THR, "peak",
                fifelse(dist_to_peak > VALLEY_THR, "valley", "transition"))]

# Report zone counts
cat("\nZone counts by compartment:\n")
print(meta[, .N, by = .(compartment, zone)][order(compartment, zone)])

# ============================================================
# PART 2: FOV-pseudobulk paired DE
# ============================================================
cat("\n--- Part 2: FOV-pseudobulk paired DE ---\n")

raw_counts <- GetAssayData(obj, layer = "counts")
gene_universe <- setdiff(rownames(raw_counts), DDR_EXCLUDE)
cat(sprintf("Gene universe: %d genes (%d panel genes - %d on-panel DDR exclusions)\n",
            length(gene_universe), nrow(raw_counts),
            nrow(raw_counts) - length(gene_universe)))

for (comp in compartments) {
  cat(sprintf("  Compartment: %s\n", comp))

  if (comp == "all") {
    cells_use <- meta[zone %in% c("peak", "valley"), cell]
  } else {
    cells_use <- meta[zone %in% c("peak", "valley") & compartment == comp, cell]
  }
  dt <- meta[cell %in% cells_use, .(cell, fov, zone)]

  # Drop FOV-zones with too few cells for reliable pseudobulk
  fov_zone_n <- dt[, .N, by = .(fov, zone)]
  fov_zone_pass <- fov_zone_n[N >= MIN_CELLS]
  dt <- dt[paste(fov, zone) %in% paste(fov_zone_pass$fov, fov_zone_pass$zone)]

  # FOVs with cells in both zones (paired design; unpaired FOVs are dropped)
  fov_both <- dt[, .(has_peak = any(zone == "peak"),
                     has_valley = any(zone == "valley")), by = fov]
  fov_use <- fov_both[has_peak == TRUE & has_valley == TRUE, fov]
  dt <- dt[fov %in% fov_use]

  if (length(fov_use) < 3) {
    cat(sprintf("    Skipping: only %d paired FOVs\n", length(fov_use)))
    next
  }

  # Pseudobulk per FOV x zone
  pb_list <- dt[, {
    cc <- cell
    list(counts = list(rowSums(raw_counts[gene_universe, cc, drop = FALSE])))
  }, by = .(fov, zone)]

  count_mat <- do.call(cbind, pb_list$counts)
  colnames(count_mat) <- paste(pb_list$fov, pb_list$zone, sep = "_")
  rownames(count_mat) <- gene_universe

  col_data <- data.frame(
    fov  = factor(pb_list$fov),
    zone = factor(pb_list$zone, levels = c("valley", "peak")),
    row.names = colnames(count_mat)
  )

  # --- DESeq2 ---
  dds <- DESeqDataSetFromMatrix(count_mat, col_data, design = ~ fov + zone)
  dds <- tryCatch({
    DESeq(dds, quiet = TRUE)
  }, error = function(e) {
    if (grepl("cannot compute log geometric means", e$message)) {
      dds <- estimateSizeFactors(dds, type = "poscounts")
      dds <- estimateDispersions(dds, quiet = TRUE)
      dds <- nbinomWaldTest(dds, quiet = TRUE)
      dds
    } else stop(e)
  })
  # Derive the coefficient name from resultsNames() rather than hardcoding
  # "zone_peak_vs_valley" -- guards against any surprise in DESeq2's
  # level-ordering/naming while producing the identical contrast.
  zone_coef <- grep("^zone_", resultsNames(dds), value = TRUE)
  res <- as.data.table(as.data.frame(results(dds, name = zone_coef)),
                       keep.rownames = "gene")
  res[, qvalue_storey := NA_real_]

  # Storey q-values
  pvals <- res$pvalue[!is.na(res$pvalue)]
  qobj <- NULL
  if (length(pvals) > 10) {
    qobj <- tryCatch(qvalue(pvals), error = function(e) NULL)
    if (!is.null(qobj)) {
      res[!is.na(pvalue), qvalue_storey := qobj$qvalues]
    }
  }

  # KS uniformity test (diagnostic only, printed not written)
  ks_res <- ks.test(res$pvalue[!is.na(res$pvalue)], "punif")

  out_deseq <- data.table(
    gene          = res$gene,
    baseMean      = res$baseMean,
    log2FC        = res$log2FoldChange,
    SE            = res$lfcSE,
    p             = res$pvalue,
    padj_BH       = res$padj,
    qvalue_storey = res$qvalue_storey,
    compartment   = comp
  )[order(p)]

  # IHW: baseMean as covariate, BH alpha
  ihw_res <- tryCatch({
    ihw(p ~ baseMean, data = out_deseq[!is.na(p)], alpha = 0.05)
  }, error = function(e) { cat(sprintf("    IHW failed: %s\n", e$message)); NULL })
  out_deseq[, padj_ihw := NA_real_]
  if (!is.null(ihw_res)) {
    out_deseq[!is.na(p), padj_ihw := adj_pvalues(ihw_res)]
    n_ihw_sig <- sum(adj_pvalues(ihw_res) < 0.05, na.rm = TRUE)
    cat(sprintf("    IHW: %d rejections at alpha=0.05\n", n_ihw_sig))
  }

  out_deseq <- tag_confirmatory(out_deseq, comp)

  fwrite(out_deseq, file.path(OUTDIR, sprintf("pv_deseq2_%s.tsv", comp)), sep = "\t")
  cat(sprintf("    DESeq2: %d genes, %d paired FOVs, KS p=%.3f (pi0=%.3f)\n",
              nrow(out_deseq), length(fov_use), ks_res$p.value,
              if (!is.null(qobj)) qobj$pi0 else NA_real_))

  # --- limma-voom ---
  y <- DGEList(counts = count_mat, group = col_data$zone)
  y <- calcNormFactors(y)
  design <- model.matrix(~ fov + zone, data = col_data)
  zone_col <- grep("^zone", colnames(design), value = TRUE)
  v <- voom(y, design)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  tt <- as.data.table(topTable(fit, coef = zone_col, number = Inf, sort.by = "none"),
                      keep.rownames = "gene")
  tt[, qvalue_storey := NA_real_]

  # Storey q-values for voom
  pvals_v <- tt$P.Value[!is.na(tt$P.Value)]
  if (length(pvals_v) > 10) {
    qobj_v <- tryCatch(qvalue(pvals_v), error = function(e) NULL)
    if (!is.null(qobj_v)) {
      tt[!is.na(P.Value), qvalue_storey := qobj_v$qvalues]
    }
  }

  out_voom <- data.table(
    gene          = tt$gene,
    log2FC        = tt$logFC,
    SE            = sqrt(fit$s2.post) * fit$stdev.unscaled[, zone_col],
    p             = tt$P.Value,
    padj_BH       = tt$adj.P.Val,
    qvalue_storey = tt$qvalue_storey,
    compartment   = comp
  )[order(p)]
  out_voom <- tag_confirmatory(out_voom, comp)

  fwrite(out_voom, file.path(OUTDIR, sprintf("pv_voom_%s.tsv", comp)), sep = "\t")
  cat(sprintf("    Voom: %d genes\n", nrow(out_voom)))
}

# ============================================================
# PART 3: Validation reruns
# ============================================================
cat("\n--- Part 3: Validation reruns ---\n")

# Module gene sets for the offset scan (p21 kept separate from the rest of
# the DDR module since it's tracked individually elsewhere in the project;
# DDR here is the remaining 17-gene module).
modules <- list(
  p21 = "Cdkn1a",
  DDR = c("Gadd45a", "Gadd45b", "Mdm2", "Bax", "Pmaip1", "Bbc3",
          "Ddit3", "Atm", "Atr", "Chek1", "Chek2", "Rad51",
          "Brca1", "Brca2", "Xrcc5", "Xrcc6", "Parp1"),
  IFN = c("Bst2", "Ifi27", "Ifit1", "Ifitm1", "Ifitm3", "Irf3",
          "Irf4", "Jak1", "Mx1", "Rigi", "Tyk2")
)

# Normalized expression for module scoring
norm_data <- GetAssayData(obj, layer = "data")

compute_module_mean <- function(expr_mat, genes, cells) {
  genes_use <- intersect(genes, rownames(expr_mat))
  if (length(genes_use) == 0) return(NA_real_)
  mean(colMeans(expr_mat[genes_use, cells, drop = FALSE]))
}

# --- Test A: Offset scan ---
# Shift the assumed beam-center grid along the perpendicular axis and
# re-derive peak/valley membership at each offset; a real biological signal
# anchored to the fitted registration should peak in magnitude near
# offset=0 and decay away from it.
offsets <- seq(-0.50, 0.50, by = 0.05)
scan_results <- rbindlist(lapply(offsets, function(off) {
  shifted_centers <- stripe$beam_centers + off
  dtp <- dist_to_nearest_beam(d_perp, shifted_centers)
  peak_cells <- meta$cell[dtp < PEAK_THR]
  valley_cells <- meta$cell[dtp > VALLEY_THR]

  row <- data.table(
    offset_mm = off,
    n_peak = length(peak_cells),
    n_valley = length(valley_cells)
  )
  for (mname in names(modules)) {
    m_peak <- compute_module_mean(norm_data, modules[[mname]], peak_cells)
    m_valley <- compute_module_mean(norm_data, modules[[mname]], valley_cells)
    set(row, j = paste0(mname, "_diff"), value = m_peak - m_valley)
  }
  row
}))

fwrite(scan_results, file.path(OUTDIR, "validation_offset_scan.tsv"), sep = "\t")
cat(sprintf("  Offset scan: %d offsets, peak p21 diff at offset %.2f mm\n",
            nrow(scan_results),
            scan_results[which.max(abs(p21_diff)), offset_mm]))

# --- Test B: Pure-FOV design ---
# Independent aggregation strategy: classify whole FOVs by centroid
# distance to the nearest beam rather than individual cells, avoiding
# within-FOV zone-boundary ambiguity. Swept over two offsets and five
# band-width thresholds (fraction of beam spacing).
offsets_b <- c(0.00, 0.10)
thresholds <- c(0.05, 0.10, 0.15, 0.20, 0.25)

fov_centroids_base <- meta[, .(cx = mean(x_slide_mm), cy = mean(y_slide_mm)), by = fov]

pure_fov_results <- rbindlist(lapply(offsets_b, function(off) {
  shifted_centers <- stripe$beam_centers + off
  fov_centroids <- copy(fov_centroids_base)
  fov_centroids[, d_centroid := {
    dp <- -sin(theta) * cx + cos(theta) * cy
    dist_to_nearest_beam(dp, shifted_centers)
  }]

  rbindlist(lapply(thresholds, function(thr) {
    peak_fovs <- fov_centroids[d_centroid < thr * stripe$spacing_mm, fov]
    valley_fovs <- fov_centroids[d_centroid > (0.5 - thr) * stripe$spacing_mm, fov]

    row <- data.table(offset = off, threshold = thr,
                      n_peak_fovs = length(peak_fovs),
                      n_valley_fovs = length(valley_fovs))

    for (mname in names(modules)) {
      pc <- meta[fov %in% peak_fovs, cell]
      vc <- meta[fov %in% valley_fovs, cell]
      if (length(pc) > 0 && length(vc) > 0) {
        m_p <- compute_module_mean(norm_data, modules[[mname]], pc)
        m_v <- compute_module_mean(norm_data, modules[[mname]], vc)

        # Per-FOV means for a t-test with FOV (not cell) as the unit
        fov_means_p <- sapply(peak_fovs, function(f) {
          fc <- meta[fov == f, cell]
          compute_module_mean(norm_data, modules[[mname]], fc)
        })
        fov_means_v <- sapply(valley_fovs, function(f) {
          fc <- meta[fov == f, cell]
          compute_module_mean(norm_data, modules[[mname]], fc)
        })
        tt <- tryCatch(t.test(fov_means_p, fov_means_v), error = function(e) NULL)

        set(row, j = paste0(mname, "_diff"), value = m_p - m_v)
        set(row, j = paste0(mname, "_t"), value = if (!is.null(tt)) tt$statistic else NA_real_)
      } else {
        set(row, j = paste0(mname, "_diff"), value = NA_real_)
        set(row, j = paste0(mname, "_t"), value = NA_real_)
      }
    }
    row
  }))
}))

fwrite(pure_fov_results, file.path(OUTDIR, "validation_pure_fov.tsv"), sep = "\t")
cat(sprintf("  Pure-FOV: %d configurations\n", nrow(pure_fov_results)))

cat("\nScript 2 complete.\n")
