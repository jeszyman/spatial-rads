#!/usr/bin/env Rscript
# Fibroblast resting/activated sub-state split (plan-differential-robustness Fix 3a).
# UCell argmax assigns resting vs activated on the pre-registered (config) panels; the
# split is GATED at the GROUP level: the UCell-activated group must show signal-to-
# background >=2 (vs the resting group) on the activated markers AND on >=1 fibroblast-
# specific anchor (Col5a1/Col5a3, absent from the SmoothMuscle panel). This blocks an
# Acta2/Tagln-driven smooth-muscle-bleed false split. Gate fail or too-few activated
# cells => no split (all resting), reported transparently.
# Args: <merged.rds> <full_labels.parquet> <substate_markers.yaml> <out.parquet> <out_gate.tsv>
suppressPackageStartupMessages({
  library(Seurat); library(UCell); library(arrow); library(data.table); library(yaml); library(Matrix)
})
a   <- commandArgs(trailingOnly = TRUE)
o   <- readRDS(a[1])
lab <- as.data.table(read_parquet(a[2]))
mk  <- read_yaml(a[3])$fibroblast
fib <- intersect(lab[cell_subtype == "Fibroblast", cell], colnames(o))
o   <- subset(o, cells = fib)
o   <- NormalizeData(o, verbose = FALSE)                       # ensure a log-norm data layer for s2b
o   <- AddModuleScore_UCell(o, features = list(resting = mk$resting, activated = mk$activated),
                            assay = "RNA", ncores = 4)
md  <- o@meta.data
assigned <- ifelse(md$activated_UCell > md$resting_UCell, "activated", "resting")
data <- LayerData(o, assay = "RNA", layer = "data")

# group-level signal-to-background: mean(log-norm expr) in activated vs resting cells
s2b <- function(genes) {
  genes <- intersect(genes, rownames(data))
  if (!length(genes) || sum(assigned == "activated") < 1 || sum(assigned == "resting") < 1)
    return(setNames(rep(NA_real_, length(genes)), genes))
  act <- Matrix::rowMeans(data[genes, assigned == "activated", drop = FALSE])
  res <- Matrix::rowMeans(data[genes, assigned == "resting",   drop = FALSE])
  act / pmax(res, 1e-6)
}
act_s2b    <- s2b(mk$activated)
anchor_s2b <- s2b(mk$specificity_anchor)
gate_ok <- sum(assigned == "activated") >= mk$gate$min_activated_cells &&
           mean(act_s2b, na.rm = TRUE) >= mk$gate$min_signal_to_bg &&
           any(anchor_s2b >= mk$gate$min_signal_to_bg, na.rm = TRUE)
substate <- if (isTRUE(gate_ok)) assigned else rep("resting", length(assigned))   # gate fail -> no split

out <- data.table(cell = rownames(md), substate = substate)
write_parquet(out, a[4])
gate <- data.table(
  n_fibroblast = nrow(out),
  n_activated  = sum(substate == "activated"),
  n_resting    = sum(substate == "resting"),
  activated_marker_s2b_mean = round(mean(act_s2b, na.rm = TRUE), 2),
  anchor_s2b   = paste(sprintf("%s=%.2f", names(anchor_s2b), anchor_s2b), collapse = ";"),
  best_anchor_s2b = round(max(anchor_s2b, na.rm = TRUE), 2),
  gate_passed  = isTRUE(gate_ok))
fwrite(gate, a[5], sep = "\t")
cat(sprintf("substate: %d fibroblasts | gate %s | %d activated / %d resting | act_s2b=%.2f | anchors %s\n",
            nrow(out), ifelse(isTRUE(gate_ok), "PASS", "FAIL"),
            sum(substate == "activated"), sum(substate == "resting"),
            mean(act_s2b, na.rm = TRUE), gate$anchor_s2b))
