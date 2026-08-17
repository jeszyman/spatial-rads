#!/usr/bin/env Rscript
# InSituCor module prioritization -- the documented post-discovery workflow (Danaher et al.
# Genome Biology 2025; CosMx Analysis Scratch Space vignette section 4): rank modules for
# analyst review by (a) which cell types drive them (module-level attribution), (b) which
# genes drive them within those cell types (gene-level attribution), and (c) overlap with
# the curated pathway sets. No statistics here; prioritization, not testing.
# Scope: modules of >= 4 genes. Diffuse background catch-alls (per-gene weights an order
# of magnitude below every real module) are flagged and excluded from candidate status;
# exceeding the max_module_size cap alone is informational, not disqualifying.
# Args: <modules.tsv> <celltype_attribution.tsv> <gene_attribution.tsv> <pathway_sets.tsv>
#       <out.tsv>
suppressPackageStartupMessages({
  library(data.table)
})
a <- commandArgs(trailingOnly = TRUE)
modules_path <- a[1]; ct_attr_path <- a[2]; gene_attr_path <- a[3]
sets_path <- a[4]; out_path <- a[5]

MIN_GENES <- 4L
# max_module_size in the discovery call; anything larger escaped the subclustering cap
OVERSIZE  <- 25L
# diffuse catch-all rule: median |gene_weight| below this fraction of the across-module
# median of medians (the 300-gene background module sits 5x under this cut, the smallest
# real module 1.5x over it)
CATCHALL_WEIGHT_FRAC <- 0.1
# candidate rule: a clear driving cell type and no curated set covering half the module
ATTR_MIN    <- 0.25
OVERLAP_MAX <- 0.5

mod <- fread(modules_path)
ct  <- fread(ct_attr_path)
ga  <- fread(gene_attr_path)
gs  <- fread(sets_path)[on_panel == TRUE, .(set, gene)]

keep <- mod[, .N, by = module][N >= MIN_GENES, module]
mod  <- mod[module %in% keep]

## ---- per-module: top cell types by module-level attribution ----
setorder(ct, module, -attribution)
top_ct <- ct[module %in% keep,
             .(top_cell_types = paste(sprintf("%s(%.2f)", head(cell_type, 3),
                                              head(attribution, 3)), collapse = "; "),
               max_attribution = attribution[1],
               top_cell_type   = cell_type[1]),
             by = module]

## ---- per-module: top genes within the driving cell type ----
setorder(ga, module, -attribution)
top_genes <- ga[top_ct, on = c("module", cell_type = "top_cell_type")][
  , .(driver_genes = paste(sprintf("%s(%.2f)", head(gene, 5), head(attribution, 5)),
                           collapse = "; ")),
  by = module]

## ---- per-module: best containment overlap with any curated set ----
overlap <- rbindlist(lapply(keep, function(m) {
  mg <- mod[module == m, gene]
  ov <- gs[gene %in% mg, .(n_overlap = uniqueN(gene)), by = set]
  if (nrow(ov) == 0) {
    data.table(module = m, best_set = NA_character_, n_overlap = 0L,
               best_overlap_frac = 0)
  } else {
    setorder(ov, -n_overlap)
    data.table(module = m, best_set = ov$set[1], n_overlap = ov$n_overlap[1],
               best_overlap_frac = round(ov$n_overlap[1] / length(mg), 3))
  }
}))

out <- mod[, .(n_genes = .N,
               median_weight = median(abs(gene_weight)),
               top_weight_genes = paste(head(gene[order(-abs(gene_weight))], 5),
                                        collapse = ", ")),
           by = module]
out[, diffuse_catchall := median_weight < CATCHALL_WEIGHT_FRAC * median(median_weight)]
out <- Reduce(function(x, y) merge(x, y, by = "module", all.x = TRUE),
              list(out, top_ct, top_genes, overlap))
out[, exceeds_size_cap := n_genes > OVERSIZE]
out[, discovery_candidate := !diffuse_catchall & max_attribution >= ATTR_MIN &
      best_overlap_frac < OVERLAP_MAX]
out[, top_cell_type := NULL]
setorder(out, diffuse_catchall, -discovery_candidate, best_overlap_frac, -max_attribution)
fwrite(out, out_path, sep = "\t")

cat(sprintf("insitucor_prioritize: %d modules >= %d genes | %d candidates | %d diffuse catch-all\n",
            nrow(out), MIN_GENES, sum(out$discovery_candidate), sum(out$diffuse_catchall)))
