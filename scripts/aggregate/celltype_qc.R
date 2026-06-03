#!/usr/bin/env Rscript
# aggregate.smk cell-type-label QC. Go/no-go gate before any differential result
# rests on the unified labels: does each canonical lineage marker enrich in its
# expected cell_subtype? Subsamples the merged 3.27M-cell object to <=CAP cells per
# cell_subtype (rare lineages kept whole so Neutrophil/Mast/Plasma calls are vetted),
# sets Idents to the unified cell_subtype, and draws a DotPlot of panel-filtered
# canonical markers. Provenance: dev/peak_valley_analysis/03_cell_type_validation.R.
# Args: <merged.rds> <full_labels.parquet> <out_markers.tsv> <plot_dotplot>
suppressPackageStartupMessages({
  library(Seurat); library(arrow); library(data.table); library(ggplot2)
})
a <- commandArgs(trailingOnly = TRUE)
merged_path <- a[1]; labels_path <- a[2]; out_tsv <- a[3]; plot_dot <- a[4]
CAP <- 8000L; SEED <- 42L

dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(plot_dot), recursive = TRUE, showWarnings = FALSE)

# canonical lineage markers, ordered by lineage (tumor -> stroma -> immune)
markers <- c(
  Tumor          = "Epcam", Tumor2        = "Krt8",   Tumor3 = "Krt18",
  Fibroblast     = "Col1a1", Fibroblast2  = "Pdgfra",
  SmoothMuscle   = "Acta2",  SmoothMuscle2= "Myh11",  SmoothMuscle3 = "Tagln",
  Endothelial    = "Pecam1", Endothelial2 = "Cdh5",
  Adipocyte      = "Cidea",  Adipocyte2   = "Fabp4",
  Tcell          = "Cd3e",   Tcell2       = "Cd8a",   Tcell3 = "Cd4",
  NK             = "Ncr1",   NK2          = "Klrb1c",
  Bcell          = "Cd19",   Bcell2       = "Ms4a1",
  Plasma         = "Mzb1",   Plasma2      = "Jchain",
  Macrophage     = "Lyz2",   Macrophage2  = "Itgam",  Macrophage3 = "Adgre1",
  DC             = "Itgax",  DC2          = "Flt3",
  Mast           = "Cpa3",   Mast2        = "Kit",
  Neutrophil     = "S100a8", Neutrophil2  = "S100a9", Neutrophil3 = "Elane")

# preferred y-axis (cell_subtype) order; intersected with present levels
subtype_order <- c("Tumor", "Epithelial cells", "Fibroblast", "SmoothMuscle",
                   "Endothelial", "Adipocyte", "T cells", "NK cells", "ILC",
                   "B cells", "Plasma cells", "Macrophages", "DC", "Mast cells",
                   "Neutrophils", "unassigned")

o   <- readRDS(merged_path)
lab <- as.data.table(read_parquet(labels_path))[, .(cell, cell_subtype)]
lv  <- setNames(lab$cell_subtype, lab$cell)
o$cell_subtype <- unname(lv[colnames(o)])
o <- subset(o, cells = colnames(o)[!is.na(o$cell_subtype)])

# stratified subsample: <=CAP cells per subtype (rare lineages kept whole)
set.seed(SEED)
idx_by <- split(seq_len(ncol(o)), o$cell_subtype)
keep   <- unlist(lapply(idx_by, function(ix)
  if (length(ix) > CAP) sample(ix, CAP) else ix), use.names = FALSE)
oq <- o[, keep]; rm(o); invisible(gc())

# panel-filter markers, order subtypes
panel <- rownames(oq)
mk    <- unique(markers[markers %in% panel])
dropped_mk <- setdiff(unique(markers), panel)
lvl   <- intersect(subtype_order, unique(oq$cell_subtype))
lvl   <- c(lvl, setdiff(unique(oq$cell_subtype), lvl))     # any extras at the end
oq$cell_subtype <- factor(oq$cell_subtype, levels = rev(lvl))
Idents(oq) <- "cell_subtype"

dp <- DotPlot(oq, features = mk, assay = "RNA") +
  scale_colour_gradient2(low = "#2166ac", mid = "grey90", high = "#b2182b") +
  labs(x = NULL, y = NULL,
       title = sprintf("Cell-subtype marker QC (<=%d cells/subtype, n=%d)",
                       CAP, ncol(oq))) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot_dot, dp, width = 11, height = 5.5, dpi = 150)

qc <- as.data.table(dp$data)
setnames(qc, c("avg.exp", "pct.exp", "features.plot", "id", "avg.exp.scaled"),
         c("avg_exp", "pct_exp", "marker", "cell_subtype", "avg_exp_scaled"),
         skip_absent = TRUE)
setorder(qc, cell_subtype, marker)
fwrite(qc, out_tsv, sep = "\t")

cat(sprintf("celltype QC: %d cells, %d subtypes, %d/%d markers on panel (dropped: %s)\n",
            ncol(oq), length(lvl), length(mk), length(unique(markers)),
            if (length(dropped_mk)) paste(dropped_mk, collapse = ",") else "none"))
