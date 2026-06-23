#!/usr/bin/env Rscript
# Pairwise gene-set overlap diagnostic (plan-differential-robustness Fix 5). Curated
# (tier=="primary") sets only. Exposes sets that are correlated by construction --
# e.g. STING (3 on-panel genes, 2 shared with the IFN sets) is not independent H1 evidence.
# Args: <pathway_sets.tsv> <out.tsv>
suppressPackageStartupMessages(library(data.table))
a  <- commandArgs(trailingOnly = TRUE)
gs <- fread(a[1])[tier == "primary"]
sets  <- split(gs$gene, gs$set)
nm    <- names(sets)
pairs <- CJ(a = nm, b = nm)[a < b]
pairs[, `:=`(
  n_a     = lengths(sets[a]),
  n_b     = lengths(sets[b]),
  shared  = mapply(function(x, y) length(intersect(sets[[x]], sets[[y]])), a, b),
  jaccard = mapply(function(x, y)
              length(intersect(sets[[x]], sets[[y]])) / length(union(sets[[x]], sets[[y]])), a, b),
  shared_genes = mapply(function(x, y) paste(intersect(sets[[x]], sets[[y]]), collapse = ";"), a, b))]
setorder(pairs, -jaccard)
fwrite(pairs, a[2], sep = "\t")
cat(sprintf("geneset_overlap: %d curated-set pairs; max jaccard %.2f; %d pairs share >=2 genes\n",
            nrow(pairs), max(pairs$jaccard), pairs[shared >= 2, .N]))
