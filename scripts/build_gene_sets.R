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
reg   <- fread(reg_p, colClasses = "character")            # concept, source, source_id, citation
stopifnot(all(c("concept","source","source_id","citation") %in% names(reg)))
stopifnot(all(reg$source %in% c("nanostring","hallmark","custom")))

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
sets_long <- rbindlist(lapply(seq_len(nrow(reg)), function(i) {
  r <- reg[i]; g <- members(r)
  data.table(concept = r$concept, source = r$source, source_id = r$source_id,
             provenance = r$source, citation = r$citation, gene = g)
}))

## ---- coverage per (concept, source) + confirmatory selection ----
cov <- sets_long[, .(source_id = source_id[1], citation = citation[1],
                     n_total = uniqueN(gene),
                     n_panel = uniqueN(gene[gene %in% panel])),
                 by = .(concept, source)]
cov[, coverage := n_panel / n_total]
cov[, thin := n_panel < min_pg]
# one confirmatory set per concept = highest panel COVERAGE (n_panel/n_total); the
# panel-designed module beats a genome-wide set diluted off-panel. ties -> nanostring, hallmark
src_rank <- c(nanostring = 1L, hallmark = 2L, custom = 3L)
cov[, confirmatory := FALSE]
sel <- cov[, .I[order(-coverage, src_rank[source])[1]], by = concept]$V1
cov[sel, confirmatory := TRUE]

## ---- attach confirmatory + coverage back onto the long table ----
sets_long <- merge(sets_long,
                   cov[, .(concept, source, confirmatory, n_panel, n_total, coverage, thin)],
                   by = c("concept", "source"), all.x = TRUE)
sets_long[, on_panel := gene %in% panel]

## ---- freeze: hand-authored surfaces (registry + custom yaml) vs git HEAD ----
freeze_check <- function(path, parse) {
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
dir.create(dirname(out_sets), recursive = TRUE, showWarnings = FALSE)
setcolorder(sets_long, c("concept","source","source_id","provenance","citation",
                         "confirmatory","gene","on_panel","n_panel","n_total","coverage","thin"))
fwrite(sets_long, out_sets, sep = "\t")
writeLines(sprintf("# msigdbr_version: %s", mver), out_cov)
fwrite(cov[, .(concept, source, source_id, citation, n_total, n_panel, coverage, thin, confirmatory)],
       out_cov, sep = "\t", append = TRUE, col.names = TRUE)
cat(sprintf("pathway_sets: %d concepts, %d (concept,source) sets; %d confirmatory | msigdbr %s\n",
            uniqueN(cov$concept), nrow(cov), cov[confirmatory == TRUE, .N], mver))
