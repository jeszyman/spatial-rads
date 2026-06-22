#!/usr/bin/env Rscript
# Generate the canonical tier-structured pathway gene-set artifact + panel-coverage
# table. Tier-1 (primary) = curated sets from pathway_gene_lists.yaml (frozen, provenance
# from pathway_sets_provenance.tsv); Tier-2 (exploratory) = MSigDB mouse Hallmark via
# msigdbr (version recorded). Single source of truth for every downstream gene-set reader.
# Args: <pathway_gene_lists.yaml> <provenance.tsv> <common_genes.tsv> <min_panel_genes>
#       <out_pathway_sets.tsv> <out_coverage.tsv>
suppressPackageStartupMessages({library(yaml); library(data.table); library(msigdbr)})

`%||%` <- function(a, b) if (is.null(a)) b else a

a        <- commandArgs(trailingOnly = TRUE)
yaml_p   <- a[1]; prov_p <- a[2]; panel_p <- a[3]
min_pg   <- as.integer(a[4]); out_sets <- a[5]; out_cov <- a[6]

panel <- readLines(panel_p)
prov  <- fread(prov_p)                                  # set, source

## ---- tier-1: curated (frozen) ----
prim   <- lapply(read_yaml(yaml_p), as.character)
miss_prov <- setdiff(names(prim), prov$set)
if (length(miss_prov)) stop("no provenance for primary set(s): ", paste(miss_prov, collapse=", "))
prim_dt <- rbindlist(lapply(names(prim), function(s) data.table(
  set = s, tier = "primary",
  source = prov[set == s, source], gene = prim[[s]])))

## ---- tier-2: MSigDB mouse Hallmark (pinned, namespaced) ----
mver <- as.character(packageVersion("msigdbr"))
hm   <- as.data.table(msigdbr(species = "Mus musculus", collection = "H"))
hm_dt <- hm[, .(set = gs_name, tier = "exploratory",
                source = sprintf("MSigDB_Hallmark (msigdbr %s)", mver),
                gene = gene_symbol)]
hm_dt <- unique(hm_dt)
stopifnot(all(grepl("^HALLMARK_", hm_dt$set)))         # namespaced -> no collision with primary

sets_long <- rbind(prim_dt, hm_dt)

## ---- coverage + usability gate ----
cov <- sets_long[, .(n_total = uniqueN(gene),
                     n_panel = uniqueN(gene[gene %in% panel])),
                 by = .(set, tier, source)]
cov[, usable := tier == "primary" | n_panel >= min_pg] # tier-1 always usable
cov[, thin   := n_panel < min_pg]

## ---- completeness guard: every primary set present ----
miss_prim <- setdiff(names(prim), unique(sets_long[tier == "primary", set]))
if (length(miss_prim)) stop("primary set(s) missing from artifact: ", paste(miss_prim, collapse=", "))

## ---- tier-1 freeze: generated primary must match the git-committed yaml ----
# Pre-registration integrity: the curated membership is frozen; a divergence between the
# working-tree yaml (what we just read) and HEAD's yaml is a silent change to the
# confirmatory family and must hard-error. system2(stdout=TRUE) yields a per-line vector,
# so collapse to one string before yaml.load.
committed <- tryCatch(
  yaml::yaml.load(paste(system2("git", c("show", "HEAD:config/pathway_gene_lists.yaml"),
                                stdout = TRUE), collapse = "\n")),
  error = function(e) NULL)
if (!is.null(committed)) {
  cg <- lapply(committed, as.character)
  for (s in names(prim))
    if (!setequal(prim[[s]], cg[[s]] %||% character()))
      stop("tier-1 freeze: primary set '", s, "' differs from committed config/pathway_gene_lists.yaml")
}

## ---- write usable sets only; coverage logs ALL (incl. excluded) ----
keep <- cov[usable == TRUE, set]
out  <- sets_long[set %in% keep]
dir.create(dirname(out_sets), recursive = TRUE, showWarnings = FALSE)
fwrite(out, out_sets, sep = "\t")
writeLines(sprintf("# msigdbr_version: %s", mver), out_cov)
fwrite(cov, out_cov, sep = "\t", append = TRUE, col.names = TRUE)
cat(sprintf("pathway_sets: %d primary + %d hallmark; %d usable, %d excluded (logged)\n",
            uniqueN(prim_dt$set), uniqueN(hm_dt$set), length(keep), cov[usable==FALSE, .N]))
