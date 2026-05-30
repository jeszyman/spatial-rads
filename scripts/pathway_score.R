#!/usr/bin/env Rscript
# Pathway/module scoring on a typed object. UCell::AddModuleScore_UCell (primary -- rank-based,
# no random background, dataset-independent) + Seurat AddModuleScore (seed=42, secondary, for
# comparison). Gene lists from config/pathway_gene_lists.yaml; genes absent from the panel dropped.
# Adds <pathway>_UCell and <pathway>_AMS columns. Args: <in.typed.rds> <gene_lists.yaml> <out.scored.rds>
suppressMessages({library(Seurat); library(UCell); library(yaml)})
args <- commandArgs(trailingOnly = TRUE)
IN <- args[1]; YAML <- args[2]; OUT <- args[3]
set.seed(42)

obj  <- readRDS(IN)
sets <- yaml::read_yaml(YAML)
present <- rownames(obj)
sets <- lapply(sets, function(g) intersect(g, present))
sets <- sets[vapply(sets, length, integer(1)) >= 2]    # need >=2 detectable genes
stopifnot(length(sets) > 0)
cat(sprintf("scoring %d pathways: %s\n", length(sets),
            paste(sprintf("%s(%d)", names(sets), lengths(sets)), collapse = ", ")))

## UCell (primary)
obj <- AddModuleScore_UCell(obj, features = sets, name = "_UCell")
## AddModuleScore (secondary, seeded) -> rename AMS1..n to <pathway>_AMS
## panel-aware ctrl/nbin: default ctrl=100 over-samples the thin ~870-gene background pool
obj <- AddModuleScore(obj, features = sets, name = "AMS", seed = 42,
                      nbin = 12, ctrl = min(50L, floor(nrow(obj) / 10)))
for (i in seq_along(sets)) {
  obj[[paste0(names(sets)[i], "_AMS")]] <- obj[[paste0("AMS", i)]]
  obj[[paste0("AMS", i)]] <- NULL
}

dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
saveRDS(obj, OUT)
cat(sprintf("scored %s: %d cells\n", sub("\\.typed\\.rds$", "", basename(IN)), ncol(obj)))
