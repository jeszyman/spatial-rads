# FDR helper functions for the results tier. Sourced by assemble_results.R and by
# the results-tier tests. Pure functions only (no CLI/pipeline side effects), so
# tests can source this without triggering the assemble_results.R run body.
suppressPackageStartupMessages({library(data.table)})

# Gatekeeping (secondary FDR view). Gates only programs that carry a DE claim; the
# parent is the matching GSEA enrichment row, the child is the claimed DE genes.
# padj_confirmatory (pooled, primary) is computed elsewhere and untouched here.
add_gatekeeping <- function(master, claims, hyp) {
  master <- copy(master)
  master[, c("gate_program","gate_padj","padj_gated") := .(NA_character_, NA_real_, NA_real_)]
  de_claims <- claims[readout_class=="DE", .(unit, contrast, gene=feature)]
  gate_set  <- unique(hyp[de_claims, on=.(cell_type=unit, contrast, gene),
                          .(cell_type, contrast, program), nomatch=NULL])
  if (!nrow(gate_set)) return(master[])
  parents <- master[readout_class=="gsea" & !is.na(pvalue)][
    gate_set, on=.(unit=cell_type, feature=program, contrast), nomatch=NULL]
  if (!nrow(parents)) return(master[])
  parents[, gate_padj := p.adjust(pvalue, "BH")]
  opened <- parents[gate_padj < 0.05, .(unit, feature, contrast, gate_padj)]
  for (i in seq_len(nrow(opened))) {
    genes <- hyp[program==opened$feature[i] & cell_type==opened$unit[i] &
                 contrast==opened$contrast[i], unique(gene)]
    idx <- master[, which(readout_class=="DE" & unit==opened$unit[i] &
                          contrast==opened$contrast[i] & feature %in% genes & !is.na(pvalue))]
    if (length(idx))
      master[idx, `:=`(gate_program = opened$feature[i],
                       gate_padj    = opened$gate_padj[i],
                       padj_gated   = p.adjust(pvalue, "BH"))]
  }
  master[]
}
