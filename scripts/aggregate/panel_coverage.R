#!/usr/bin/env Rscript
# Per-(set x gene x cell_group) detection report for the pathway readout genes,
# M02 day2 cells only. For each gene in every config/pathway_gene_lists.yaml set:
# fraction of cells with count>0 and mean count, overall ("all"), within the
# tumor/stroma/immune compartments, and within the Endothelial subtype (Angiogenesis
# is an endothelial program; its markers are endothelial-specific, so against the
# broad stroma pool -- endothelial is only ~3% of stroma -- they read as near-floor
# by dilution, not biology). Each set maps to the cell group its hypothesis targets
# (immune for IFN/STING, tumor for DDR, Endothelial for Angiogenesis, stroma for the
# remaining vascular/fibrosis sets); near_floor flags a (set,gene) whose detection in
# that target group is <5%, so a null on it is read as panel-limited, not biological.
# Args: <merged_typed.rds> <full_labels.parquet> <pathway_yaml> <out.tsv>
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(arrow); library(yaml); library(data.table)
})

a  <- commandArgs(trailingOnly = TRUE)
o  <- readRDS(a[1])
md <- o@meta.data
keep <- which(md$dataset == "Mutter_02")
stopifnot(length(keep) > 0)
cnt   <- LayerData(o, assay = "RNA", layer = "counts")[, keep, drop = FALSE]
cells <- colnames(cnt)
rm(o); invisible(gc())

lab  <- as.data.table(read_parquet(a[2]))[, .(cell, compartment, cell_subtype)]
m    <- match(cells, lab$cell)
comp <- lab$compartment[m]
subt <- lab$cell_subtype[m]

sets  <- yaml::read_yaml(a[3])
genes <- sort(unique(unlist(sets)))
genes <- genes[genes %in% rownames(cnt)]          # panel-present only
cm    <- cnt[genes, , drop = FALSE]
det   <- cm > 0

frac_mean <- function(idx, label) data.table(
  gene = genes, cell_group = label,
  detect_frac = as.numeric(rowSums(det[, idx, drop = FALSE]) / length(idx)),
  mean_count  = as.numeric(rowSums(cm[, idx, drop = FALSE]) / length(idx)))

res <- frac_mean(seq_len(ncol(cm)), "all")
for (cc in c("tumor", "stroma", "immune")) {
  j <- which(comp == cc & !is.na(comp))
  if (length(j) > 0) res <- rbind(res, frac_mean(j, cc))
}
je <- which(subt == "Endothelial" & !is.na(subt))
if (length(je) > 0) res <- rbind(res, frac_mean(je, "Endothelial"))

target <- c(TypeI_interferon = "immune", TypeII_interferon = "immune", STING = "immune",
            DNA_Damage_Repair = "tumor",
            Angiogenesis = "Endothelial", Hypoxia = "stroma",
            Fibrosis_remodeling = "stroma", Stromal_stress_senescence = "stroma")
map <- rbindlist(lapply(names(sets), function(s)
  data.table(set = s, gene = sets[[s]])))[gene %in% genes]

out <- merge(map, res, by = "gene", allow.cartesian = TRUE)
out[, target_group := target[set]]
out[, is_target := cell_group == target_group]
tgt <- out[is_target == TRUE, .(set, gene, target_detect = detect_frac)]
out <- merge(out, tgt, by = c("set", "gene"), all.x = TRUE)
out[, near_floor := !is.na(target_detect) & target_detect < 0.05]
setcolorder(out, c("set", "gene", "cell_group", "detect_frac", "mean_count",
                   "target_group", "is_target", "target_detect", "near_floor"))
setorder(out, set, gene, cell_group)
fwrite(out, a[4], sep = "\t")

nf <- unique(out[near_floor == TRUE, .(set, gene)])
cat(sprintf("readout detection: %d panel genes x %d sets, %d M02 cells | %d (set,gene) near-floor (<5%% in target group)\n",
            length(genes), length(sets), ncol(cm), nrow(nf)))
