#!/usr/bin/env Rscript
# Audit every gene set in config/pathway_gene_lists.yaml against the common 950-gene
# panel. The YAML keeps full curated lists (off-panel genes drop at scoring time, the
# shipped convention); this script does NOT rewrite it, it only logs each gene's panel
# membership so coverage is auditable and panel-thin sets are visible. Reads the panel
# gene list directly (no merged RDS load needed for a name lookup).
# Args: <common_genes.tsv> <pathway_gene_lists.yaml> <out_coverage.tsv>
suppressPackageStartupMessages({library(yaml); library(data.table)})

a     <- commandArgs(trailingOnly = TRUE)
panel <- readLines(a[1])
sets  <- yaml::read_yaml(a[2])

cov <- rbindlist(lapply(names(sets), function(s) data.table(
  set = s, gene = sets[[s]], in_panel = sets[[s]] %in% panel)))
fwrite(cov, a[3], sep = "\t")

per_set <- cov[, .(n_total = .N, n_panel = sum(in_panel)), by = set]
cat(sprintf("gene sets: %d sets; %d/%d genes on panel\n",
            nrow(per_set), sum(cov$in_panel), nrow(cov)))
invisible(per_set[, cat(sprintf("  %-26s %d/%d on panel\n", set, n_panel, n_total)),
                  by = set])
