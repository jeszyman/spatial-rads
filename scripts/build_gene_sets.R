#!/usr/bin/env Rscript
# Generate the canonical pathway gene-set artifact + panel-coverage table from the
# concept registry (config/pathway_concept_registry.tsv). Every set is one of three
# provenances, resolved to genes at build time so members are never hand-duplicated:
#   nanostring -> genes tagged '+' in the named vendor Annotations module column(s)
#   hallmark   -> msigdbr(Mus musculus, collection H) filtered to the named HALLMARK_* set(s)
#   custom     -> config/pathway_gene_lists.yaml (bespoke concepts; citation required in registry)
# source_id may name several modules/sets joined by ';' -> union. Per concept the
# (concept,source) with the most on-panel genes is flagged confirmatory; ties -> nanostring.
# Freeze: the two hand-authored surfaces (registry + custom yaml) are checked against
# git HEAD; machine-derived members legitimately change with the vendor file / msigdbr.
# Args: <registry.tsv> <custom_yaml> <vendor_xlsx> <common_genes.tsv> <min_panel_genes>
#       <out_pathway_sets.tsv> <out_coverage.tsv>
suppressPackageStartupMessages({
  library(yaml); library(data.table); library(msigdbr); library(readxl)
})

`%||%` <- function(a, b) if (is.null(a)) b else a
split_ids <- function(x) trimws(strsplit(x, ";", fixed = TRUE)[[1]])

a         <- commandArgs(trailingOnly = TRUE)
reg_p     <- a[1]; yaml_p <- a[2]; vendor_p <- a[3]; panel_p <- a[4]
min_pg    <- as.integer(a[5]); out_sets <- a[6]; out_cov <- a[7]

panel <- readLines(panel_p)
reg   <- fread(reg_p, colClasses = "character")            # concept, source, source_id, tier, citation
stopifnot(all(c("concept","source","source_id","tier","citation") %in% names(reg)))
stopifnot(all(reg$source %in% c("nanostring","hallmark","custom")))
stopifnot(all(reg$tier %in% c("primary","exploratory")))

## ---- custom members: bespoke concepts must be present with a citation ----
custom <- lapply(read_yaml(yaml_p), as.character)
reg_custom <- reg[source == "custom"]
miss_cit <- reg_custom[is.na(citation) | citation == "" | citation == "NA", concept]
if (length(miss_cit)) stop("custom set(s) missing a citation in registry: ", paste(miss_cit, collapse=", "))
miss_yaml <- setdiff(reg_custom$concept, names(custom))
if (length(miss_yaml)) stop("custom set(s) in registry but not in yaml: ", paste(miss_yaml, collapse=", "))

## ---- nanostring source: vendor Annotations sheet, module columns hold '+' flags ----
raw <- as.data.table(read_excel(vendor_p, sheet = "Annotations"))
hdr <- as.character(unlist(raw[1, ]))                      # module names live in row 1
ann <- raw[-1, ]; setnames(ann, hdr); setnames(ann, 1, "Gene")
ns_module <- function(mod) {
  if (!mod %in% names(ann)) stop("vendor module not found: '", mod, "'")
  g <- ann$Gene[ann[[mod]] == "+"]; g[!is.na(g)]
}

## ---- hallmark source: pinned msigdbr mouse Hallmark collection ----
mver <- as.character(packageVersion("msigdbr"))
hm   <- as.data.table(msigdbr(species = "Mus musculus", db_species = "MM", collection = "MH"))
hm_set <- function(nm) {
  if (!nm %in% hm$gs_name) stop("hallmark set not found in msigdbr ", mver, ": '", nm, "'")
  unique(hm[gs_name == nm, gene_symbol])
}

## ---- resolve every (concept, source) row to its gene set ----
members <- function(r) {
  ids <- if (r$source == "custom") NA else split_ids(r$source_id)
  g <- switch(r$source,
    nanostring = unique(unlist(lapply(ids, ns_module))),
    hallmark   = unique(unlist(lapply(ids, hm_set))),
    custom     = custom[[r$concept]])
  unique(as.character(g))
}
# set = the scored-set key: concept name for a primary set, else the concept (the
# hallmark twin of a named concept keeps the HALLMARK_* id so it is distinct from its
# primary sibling). Downstream (pathway_scores.R, assemble_results.R) splits genes by `set`.
set_name <- function(r) if (r$source == "hallmark") r$source_id else r$concept
sets_long <- rbindlist(lapply(seq_len(nrow(reg)), function(i) {
  r <- reg[i]; g <- members(r)
  data.table(set = set_name(r), concept = r$concept, source = r$source, source_id = r$source_id,
             provenance = r$source, tier = r$tier, citation = r$citation, gene = g)
}))

## ---- full Hallmark exploratory sweep: every MH set clearing the panel-gene floor ----
## (concept-twin Hallmark sets already emitted above are excluded to avoid duplicate rows)
named_hm <- unique(unlist(lapply(reg[source == "hallmark", source_id], split_ids)))
hm_all   <- hm[, .(gene = unique(gene_symbol)), by = .(set = gs_name)]
hm_all[, n_panel := sum(gene %in% panel), by = set]
hm_sweep <- hm_all[n_panel >= min_pg & !set %in% named_hm]
hm_sweep <- hm_sweep[, .(set, concept = set, source = "hallmark", source_id = set,
                         provenance = "hallmark", tier = "exploratory", citation = NA_character_, gene)]
sets_long <- rbind(sets_long, hm_sweep, fill = TRUE)

## ---- coverage per scored set + confirmatory selection ----
cov <- sets_long[, .(concept = concept[1], source = source[1], source_id = source_id[1],
                     tier = tier[1], citation = citation[1],
                     n_total = uniqueN(gene),
                     n_panel = uniqueN(gene[gene %in% panel])),
                 by = set]
cov[, coverage := n_panel / n_total]
cov[, thin := n_panel < min_pg]
# one confirmatory set per PRIMARY concept = highest panel COVERAGE (n_panel/n_total); the
# panel-designed module beats a genome-wide set diluted off-panel. exploratory sets are
# never confirmatory. ties -> nanostring, hallmark, custom.
src_rank <- c(nanostring = 1L, hallmark = 2L, custom = 3L)
cov[, confirmatory := FALSE]
conf_sets <- cov[tier == "primary"][order(concept, -coverage, src_rank[source]),
                                    .(set = set[1]), by = concept]$set
cov[set %in% conf_sets, confirmatory := TRUE]

## ---- attach confirmatory + coverage back onto the long table ----
sets_long <- merge(sets_long,
                   cov[, .(set, confirmatory, n_panel, n_total, coverage, thin)],
                   by = "set", all.x = TRUE)
sets_long[, on_panel := gene %in% panel]

## ---- freeze: hand-authored surfaces (registry + custom yaml) vs git HEAD ----
# SPATIALRADS_FREEZE_SKIP=1 bypasses the guard for the intended "building the next
# frozen state" step: after a deliberate registry/yaml revision, build once with the
# skip to regenerate artifacts, then COMMIT the new registry+yaml (which becomes HEAD
# and the new frozen baseline). Never set it in the committed workflow.
freeze_skip <- nzchar(Sys.getenv("SPATIALRADS_FREEZE_SKIP"))
freeze_check <- function(path, parse) {
  if (freeze_skip) return(NULL)
  txt <- tryCatch(paste(system2("git", c("show", paste0("HEAD:", path)), stdout = TRUE),
                        collapse = "\n"), error = function(e) NULL)
  if (is.null(txt) || !nzchar(txt)) return(invisible())   # not yet committed -> nothing to guard
  parse(txt)
}
committed_reg <- freeze_check("config/pathway_concept_registry.tsv",
  function(t) fread(text = t, colClasses = "character"))
if (!is.null(committed_reg) && !isTRUE(all.equal(
      committed_reg[order(concept, source)], reg[order(concept, source)], check.attributes = FALSE)))
  stop("freeze: config/pathway_concept_registry.tsv differs from committed HEAD")
committed_yaml <- freeze_check("config/pathway_gene_lists.yaml",
  function(t) lapply(yaml::yaml.load(t), as.character))
if (!is.null(committed_yaml)) for (s in names(custom))
  if (!setequal(custom[[s]], committed_yaml[[s]] %||% character()))
    stop("freeze: custom set '", s, "' differs from committed config/pathway_gene_lists.yaml")

## ---- write ----
# `set` + `tier` lead so downstream consumers (pathway_scores.R splits genes by `set`;
# assemble_results.R filters tier == "primary" and matches cov by `set`) read a stable key.
dir.create(dirname(out_sets), recursive = TRUE, showWarnings = FALSE)
setcolorder(sets_long, c("set","tier","concept","source","source_id","provenance","citation",
                         "confirmatory","gene","on_panel","n_panel","n_total","coverage","thin"))
fwrite(sets_long, out_sets, sep = "\t")
writeLines(sprintf("# msigdbr_version: %s", mver), out_cov)
cov[, usable := tier == "primary" | n_panel >= min_pg]      # legacy column assemble_results reads
fwrite(cov[, .(set, tier, concept, source, source_id, citation,
               n_total, n_panel, coverage, usable, thin, confirmatory)],
       out_cov, sep = "\t", append = TRUE, col.names = TRUE)
cat(sprintf("pathway_sets: %d sets (%d primary, %d exploratory); %d confirmatory | msigdbr %s\n",
            nrow(cov), cov[tier=="primary", .N], cov[tier=="exploratory", .N],
            cov[confirmatory == TRUE, .N], mver))
