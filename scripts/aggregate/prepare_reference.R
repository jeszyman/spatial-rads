#!/usr/bin/env Rscript
# Build the de-anchored external cell-type reference (panel space) for aggregate atlas typing.
# Per the agreed plan this is a MAMMARY-CENTERED reference: NanoString CellProfileLibrary
# Mouse/Adult MammaryGland_Virgin (tissue-matched -- 4T1 is a mammary carcinoma: luminal
# epithelial + mammary stroma) + ImmGen cellFamily (immune depth: T/B/Mac/DC/NK/Plasma/Neutrophil
# plus thin endothelial/pericyte/fibroblast). Lung + Muscle were DROPPED 2026-06-01: they diluted
# the tissue match (lung AT1/AT2 epithelium polluting the Epithelial anchor) and only existed to
# anchor SmoothMuscle, which neither mammary nor ImmGen provides -- so SmoothMuscle is intentionally
# left UNANCHORED (smooth-muscle cells fall to the nearest stromal type, Pericyte/Fibroblast).
# Relabels every profile COLUMN to a coarse lineage via an ordered keyword map, drops junk-drawer
# attractors (Stem/Prog, Erythroblast, Baso/Eos, organ-specific tubule epithelium), and subsets to
# the 950-gene CosMx panel. No M01 anything: identical external knowledge applied to M01 and M02.
# Args: <cpl_adult_dir> <common_genes.tsv> <out_ref.rds> <out_coverage.tsv>
suppressMessages({library(Matrix)})
args <- commandArgs(trailingOnly = TRUE)
CPL   <- args[1]; PANELF <- args[2]; OUT <- args[3]; COV <- args[4]
panel <- readLines(PANELF)

files <- c("ImmuneAtlas_ImmGen_cellFamily.RData", "MammaryGland_Virgin_MCA.RData")

# Ordered keyword map: FIRST match wins (so Acta2-stroma -> SmoothMuscle before Fibroblast,
# endothelial before Dendritic, etc.). Unmatched columns are dropped (returns NA).
map_lineage <- function(col) {
  rules <- list(
    Endothelial  = "endothelial",
    Pericyte     = "pericyte",
    SmoothMuscle = "Acta2|Muscle\\.cell|Muscle\\.progenitor|myocyte|smooth.?muscle",
    Fibroblast   = "stromal|fibroblast|Col3a1|Col1a|Dcn|Pi16|Inmt|Mgp|reticular",
    Epithelial   = "Luminal|AT1|AT2|Clara|Ciliated|epithel|Alveolar\\.bipotent",
    Plasma       = "Plasma",
    B            = "B\\.cell|^B$|Memory\\.B|GC\\.centro|CD19|IgA",
    T            = "T\\.cell|CD8\\.T|CD4|gdT|NKT|Treg|Thymic|DN[0-9]|ISP|preT|Naive\\.CD|CD4Act|^T$|alpha.?beta|gamma.?delta",
    NK           = "^NK|NK\\.cell|ILC|Nuocyte",
    DC           = "Dendritic|pDC|plasmacytoid",
    Macrophage   = "Macrophage|macs|Monocyte|monocyte|Microglia|Kupffer",
    Neutrophil   = "Neutrophil|Granulocyte")
  for (lin in names(rules)) if (grepl(rules[[lin]], col, ignore.case = TRUE)) return(lin)
  NA_character_
}

cols <- list()   # accumulate per-profile vectors (named by panel genes)
meta <- list()   # profile -> c(lineage, source, label)
for (f in files) {
  e <- new.env(); load(file.path(CPL, f), envir = e); pm <- e$profile_matrix
  pm <- as.matrix(pm)
  for (j in seq_len(ncol(pm))) {
    raw <- colnames(pm)[j]; lin <- map_lineage(raw)
    if (is.na(lin)) next
    v <- setNames(rep(NA_real_, length(panel)), panel)
    common <- intersect(rownames(pm), panel)
    v[common] <- pm[common, j]
    id <- paste0(sub("_MCA|\\.RData|ImmuneAtlas_", "", f), ":", raw)
    cols[[id]] <- v
    meta[[id]] <- c(lineage = lin, source = f, label = raw)
  }
}
mat <- do.call(cbind, cols)                       # 950 x P, NA where profile lacks gene
lineage <- vapply(meta, `[[`, "", "lineage"); names(lineage) <- names(meta)
source  <- vapply(meta, `[[`, "", "source");  names(source)  <- names(meta)

ref <- list(mat = mat, lineage = lineage, source = source, panel = panel)
saveRDS(ref, OUT)

# coverage report
covdf <- data.frame(profile = colnames(mat), lineage = lineage[colnames(mat)],
                    source = source[colnames(mat)],
                    n_panel_genes = colSums(!is.na(mat)), row.names = NULL)
agg <- aggregate(n_panel_genes ~ lineage, covdf, function(x) round(mean(x)))
agg$n_profiles <- as.integer(table(covdf$lineage)[agg$lineage])
agg <- agg[order(-agg$n_profiles), c("lineage","n_profiles","n_panel_genes")]
write.table(covdf, COV, sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("ref_profiles: %d profiles x %d panel genes across %d lineages\n",
            ncol(mat), nrow(mat), length(unique(lineage))))
cat("lineage | n_profiles | mean_panel_genes\n")
for (i in seq_len(nrow(agg)))
  cat(sprintf("  %-13s %3d   %4d\n", agg$lineage[i], agg$n_profiles[i], agg$n_panel_genes[i]))
