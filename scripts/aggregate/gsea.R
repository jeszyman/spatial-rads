#!/usr/bin/env Rscript
# aggregate.smk Track 2 pathway -- GSEA on pseudobulk DESeq2 stat-ranked genes.
# Per (cell_type x contrast), rank panel genes by DESeq2 Wald `stat` (descending)
# and run fgseaMultilevel against 4 project-priority pathways (config yaml,
# tier=primary) plus the 50 MSigDB Hallmark sets (msigdbr human->mouse ortholog
# symbols, tier=exploratory). On a ~950-gene targeted panel most Hallmark sets
# overlap the panel only partially, so n_set_genes (full set) and n_panel_genes
# (overlap with the ranked list = the genes fgsea actually scored) are reported per
# row and downstream applies coverage thresholds. minSize=3 retains the small
# curated primary sets (STING has 3 genes). padj_bh is fgsea's BH across the sets
# tested within each (cell_type x contrast) ranking; pvalue is kept so a tier-scoped
# re-adjustment is possible downstream. Args: <degs.tsv> <pathway_yaml> <out_gsea.tsv>
suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(msigdbr)
  library(yaml)
})

args      <- commandArgs(trailingOnly = TRUE)
degs_path <- args[1]
yaml_path <- args[2]
out_gsea  <- args[3]

MIN_SIZE <- 3L     # keep small curated primary sets (STING = 3 genes)
MAX_SIZE <- 500L

# --- gene sets: 4 project-priority (primary) + 50 MSigDB Hallmark (exploratory) ---
prim_lists <- lapply(read_yaml(yaml_path), as.character)
prim_meta  <- data.table(pathway_name = names(prim_lists),
                         pathway_source = "project", tier = "primary")

hm       <- as.data.table(msigdbr(species = "Mus musculus", collection = "H"))
hm_lists <- lapply(split(hm$gene_symbol, hm$gs_name), unique)
hm_meta  <- data.table(pathway_name = names(hm_lists),
                       pathway_source = "MSigDB_Hallmark", tier = "exploratory")

all_sets      <- c(prim_lists, hm_lists)
set_meta      <- rbind(prim_meta, hm_meta)
set_size_full <- vapply(all_sets, length, integer(1))   # full set size (panel-independent)

# --- DE stats: one ranking per (cell_type x contrast) over tested genes ---
degs   <- fread(degs_path)
strata <- unique(degs[!is.na(stat), .(cell_type, contrast)])

run_one <- function(ct, cn) {
  d     <- degs[cell_type == ct & contrast == cn & !is.na(stat)]
  ranks <- d$stat; names(ranks) <- d$gene
  ranks <- ranks[!duplicated(names(ranks))]
  ranks <- sort(ranks, decreasing = TRUE)
  fg <- fgsea(pathways = all_sets, stats = ranks,
              minSize = MIN_SIZE, maxSize = MAX_SIZE, eps = 0, nproc = 1)
  if (nrow(fg) == 0L) return(NULL)
  fg <- as.data.table(fg)
  fg[, `:=`(cell_type = ct, contrast = cn,
            leading_edge_genes = vapply(leadingEdge, paste, character(1), collapse = ","),
            leading_edge_size  = lengths(leadingEdge))]
  fg
}

set.seed(1)
res <- rbindlist(lapply(seq_len(nrow(strata)),
                        function(i) run_one(strata$cell_type[i], strata$contrast[i])),
                 fill = TRUE)

# --- assemble to schema ---
res <- merge(res, set_meta, by.x = "pathway", by.y = "pathway_name", all.x = TRUE)
res[, n_set_genes := set_size_full[pathway]]
setnames(res, c("pathway", "pval", "padj", "size"),
         c("pathway_name", "pvalue", "padj_bh", "n_panel_genes"))
res[, dataset := "Mutter_02"]
out <- res[, .(cell_type, contrast, pathway_name, pathway_source, tier,
               NES, pvalue, padj_bh, leading_edge_genes, leading_edge_size,
               n_set_genes, n_panel_genes, dataset)]
setorder(out, cell_type, contrast, padj_bh, na.last = TRUE)

dir.create(dirname(out_gsea), recursive = TRUE, showWarnings = FALSE)
fwrite(out, out_gsea, sep = "\t")

cat(sprintf("gsea: %d (cell_type x contrast) rankings | %d pathway rows | %d padj_bh<0.05 (%d primary, %d hallmark) | %d primary-tier rows\n",
            nrow(strata), nrow(out),
            out[!is.na(padj_bh) & padj_bh < 0.05, .N],
            out[!is.na(padj_bh) & padj_bh < 0.05 & tier == "primary", .N],
            out[!is.na(padj_bh) & padj_bh < 0.05 & tier == "exploratory", .N],
            out[tier == "primary", .N]))
