#!/usr/bin/env Rscript
# Sample-level descriptive log2FC tables for the three M01 4h treatment-arm
# contrasts (MBRT_vs_Ctrl, SBRT_vs_Ctrl, MBRT_vs_SBRT) x 4 compartment cuts
# (all/tumor/stroma/immune). n=1 per arm (no replicates at 4h), so this
# pseudobulk-sums the raw counts per arm and runs edgeR's fixed-dispersion
# exactTest -- the project's standard no-replicate descriptive approach
# (report effect sizes, not p-values; only log2FC/logCPM/pct_expressing are
# written out).
suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(Seurat)
  library(edgeR)
  library(yaml)
})

DATADIR <- "/mnt/data/projects/spatial-rads/processing/norm"
LABELS  <- "results/aggregate/full_labels.parquet"
HYPO    <- "config/hypotheses.yaml"
PANEL   <- "results/data_model/common_genes.tsv"
OUTDIR  <- "results/aggregate/m01_4h"
BCV     <- 0.2

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# sample_id is needed to reconstruct the merged-cohort cell id (labels are
# keyed as "<sample_id>_<barcode>"; the per-sample norm.rds colnames are the
# bare barcode).
samples <- list(
  Control = list(sample_id = "sam0001", path = file.path(DATADIR, "sam0001.norm.rds")),
  MBRT_4h = list(sample_id = "sam0003", path = file.path(DATADIR, "sam0003.norm.rds")),
  SBRT_4h = list(sample_id = "sam0006", path = file.path(DATADIR, "sam0006.norm.rds"))
)

contrasts <- list(
  MBRT_vs_Ctrl = c("MBRT_4h", "Control"),
  SBRT_vs_Ctrl = c("SBRT_4h", "Control"),
  MBRT_vs_SBRT = c("MBRT_4h", "SBRT_4h")
)

compartments <- c("all", "tumor", "stroma", "immune")

# --- Cell-type labels ---
lab <- as.data.table(read_parquet(LABELS))[, .(cell, compartment)]

# --- Hypothesis gene sets for flagging ---
gene_universe <- readLines(PANEL)

# A handful of panel probes cover multiple paralogs and are stored as one
# slash-joined symbol (e.g. "Oas1a/g"); hypotheses.yaml lists them the same
# way, so a literal match against the panel handles them directly. This is a
# fallback for any alias that does NOT match the panel verbatim: split on "/"
# and test the base symbol with each suffix swapped in against the actual
# gene universe, rather than assuming a naming convention.
expand_gene_alias <- function(gene, universe) {
  if (gene %in% universe) return(gene)
  if (!grepl("/", gene, fixed = TRUE)) return(character(0))
  parts <- strsplit(gene, "/", fixed = TRUE)[[1]]
  base  <- parts[1]
  cand  <- vapply(parts[-1], function(suf)
    paste0(substr(base, 1, nchar(base) - nchar(suf)), suf), character(1))
  cand  <- c(base, cand)
  cand[cand %in% universe]
}

hypo_raw  <- read_yaml(HYPO)
hypo_long <- rbindlist(lapply(hypo_raw$hypotheses, function(h) {
  genes <- unique(unlist(lapply(h$programs, function(p) p$genes)))
  data.table(gene_raw = genes, hypothesis = h$name)
}))
hypo_long[, gene := lapply(gene_raw, expand_gene_alias, universe = gene_universe)]
n_unresolved <- hypo_long[lengths(gene) == 0, uniqueN(gene_raw)]
hypo_long <- hypo_long[lengths(gene) > 0, .(gene = unlist(gene)), by = .(hypothesis, gene_raw)]
# Genes claimed by multiple hypotheses: collapse to comma-separated.
hypo_map <- hypo_long[, .(hypothesis = paste(sort(unique(hypothesis)), collapse = ",")), by = gene]
cat(sprintf("Hypothesis genes: %d flagged on-panel, %d off-panel/unresolved\n",
            uniqueN(hypo_long$gene_raw), n_unresolved))

# --- Load Seurat objects; tag cells with the merged-cohort id + compartment ---
cat("Loading Seurat objects...\n")
objs <- lapply(samples, function(s) {
  obj <- readRDS(s$path)
  bc   <- colnames(obj)
  dt   <- data.table(barcode = bc, cell = paste0(s$sample_id, "_", bc))
  dt   <- merge(dt, lab, by = "cell", all.x = TRUE)
  list(obj = obj, cells = dt)
})

# --- Run contrasts ---
for (cname in names(contrasts)) {
  pair    <- contrasts[[cname]]
  s1_name <- pair[1]
  s2_name <- pair[2]
  s1 <- objs[[s1_name]]
  s2 <- objs[[s2_name]]

  for (comp in compartments) {
    cat(sprintf("  %s / %s\n", cname, comp))

    if (comp == "all") {
      keep1 <- s1$cells$barcode
      keep2 <- s2$cells$barcode
    } else {
      keep1 <- s1$cells[compartment == comp, barcode]
      keep2 <- s2$cells[compartment == comp, barcode]
    }

    if (length(keep1) < 10 || length(keep2) < 10) {
      cat(sprintf("    Skipping: too few cells (%d, %d)\n", length(keep1), length(keep2)))
      next
    }

    # Pseudobulk: sum raw counts
    raw1 <- GetAssayData(s1$obj, layer = "counts")[, keep1, drop = FALSE]
    raw2 <- GetAssayData(s2$obj, layer = "counts")[, keep2, drop = FALSE]
    pb1  <- rowSums(raw1)
    pb2  <- rowSums(raw2)

    # Percent expressing
    pct1 <- rowMeans(raw1 > 0)
    pct2 <- rowMeans(raw2 > 0)

    # edgeR with fixed BCV (n=1 per arm; no replicates to estimate dispersion)
    counts_mat <- cbind(pb1, pb2)
    colnames(counts_mat) <- c(s1_name, s2_name)
    group <- factor(c(s1_name, s2_name), levels = c(s2_name, s1_name))
    y  <- DGEList(counts = counts_mat, group = group)
    y  <- calcNormFactors(y)
    et <- exactTest(y, dispersion = BCV^2)
    tt <- as.data.table(topTags(et, n = Inf, sort.by = "none")$table, keep.rownames = "gene")

    # Assemble output
    out <- data.table(
      gene   = tt$gene,
      log2FC = tt$logFC,
      logCPM = tt$logCPM,
      pct_expressing_group1 = pct1[tt$gene],
      pct_expressing_group2 = pct2[tt$gene]
    )

    # Flag hypothesis genes
    out <- merge(out, hypo_map, by = "gene", all.x = TRUE)

    # Sort by absolute log2FC descending
    out <- out[order(-abs(log2FC))]

    fname <- sprintf("arm_%s_%s.tsv", cname, comp)
    fwrite(out, file.path(OUTDIR, fname), sep = "\t")
    cat(sprintf("    Wrote %s (%d genes, %d/%d cells)\n",
                fname, nrow(out), length(keep1), length(keep2)))
  }
}

cat("Script 1 complete.\n")
