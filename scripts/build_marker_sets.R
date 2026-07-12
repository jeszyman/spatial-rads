#!/usr/bin/env Rscript
# Panel-coverage audit for the cell-lineage + substate marker sets, mirroring
# build_gene_sets.R. Emits per-set n_total / n_panel / frac_panel / thin, carrying the
# provenance source and known weak-separation flags. Enforces two integrity claims the
# marker YAMLs currently only assert in comments: (1) coarse-lineage on-panel completeness
# (every coarse marker present on the common panel) and (2) a git-HEAD freeze on the
# pre-registered marker membership. Stage-A metadata only.
# Args: <lineage.yaml> <substate.yaml> <provenance.tsv> <common_genes.tsv> <min_panel_genes> <out.tsv>
suppressPackageStartupMessages({library(yaml); library(data.table)})

`%||%` <- function(a, b) if (is.null(a)) b else a

a       <- commandArgs(trailingOnly = TRUE)
lin_p   <- a[1]; sub_p <- a[2]; prov_p <- a[3]; panel_p <- a[4]
min_pg  <- as.integer(a[5]); out_p <- a[6]

panel <- readLines(panel_p)
prov  <- fread(prov_p)

## ---- flatten lineage (flat) + substate (nested) into (set, gene) ----
lin <- read_yaml(lin_p)
lin_dt <- rbindlist(lapply(names(lin), function(s)
  data.table(set = s, gene = as.character(lin[[s]]))))

sub <- read_yaml(sub_p)$fibroblast
sub_dt <- rbindlist(list(
  data.table(set = "fibroblast_resting",   gene = as.character(sub$resting)),
  data.table(set = "fibroblast_activated", gene = as.character(sub$activated))))

sets_long <- rbind(lin_dt, sub_dt)

## ---- provenance completeness (mirrors build_gene_sets.R:21-22) ----
miss_prov <- setdiff(unique(sets_long$set), prov$set)
if (length(miss_prov)) stop("no provenance for marker set(s): ", paste(miss_prov, collapse = ", "))

## ---- coverage ----
cov <- sets_long[, .(n_total = uniqueN(gene),
                     n_panel = uniqueN(gene[gene %in% panel])),
                 by = set]
cov <- merge(cov, prov[, .(set, tier, source, weak_separation)], by = "set", all.x = TRUE)
cov[, frac_panel := n_panel / n_total]
cov[, thin := n_panel < min_pg]
setcolorder(cov, c("set","tier","source","n_total","n_panel","frac_panel","thin","weak_separation"))

## ---- on-panel completeness check: coarse lineages claim ALL markers on-panel ----
coarse_off <- lin_dt[!gene %in% panel]
if (nrow(coarse_off))
  warning("coarse-lineage markers OFF the common panel (weakens that lineage): ",
          paste(sprintf("%s:%s", coarse_off$set, coarse_off$gene), collapse = ", "))

## ---- freeze: lineage marker membership must match git HEAD ----
lin_head <- tryCatch(yaml::yaml.load(paste(system2("git",
  c("show", "HEAD:config/lineage_markers.yaml"), stdout = TRUE), collapse = "\n")),
  error = function(e) NULL)
if (!is.null(lin_head)) for (s in names(lin))
  if (!setequal(as.character(lin[[s]]), as.character(lin_head[[s]] %||% character())))
    stop("freeze: lineage set '", s, "' differs from committed lineage_markers.yaml")

dir.create(dirname(out_p), recursive = TRUE, showWarnings = FALSE)
fwrite(cov, out_p, sep = "\t")
cat(sprintf("marker coverage: %d sets (%d coarse off-panel)\n",
            nrow(cov), nrow(coarse_off)))
print(cov[, .(set, tier, n_panel, n_total, thin, weak_separation)])
