# scripts/aggregate/load_claims.R
suppressPackageStartupMessages({library(yaml); library(data.table)})
source("scripts/aggregate/load_hypotheses.R")

load_claims <- function(claims_path, hyp_path) {
  cfg <- yaml::read_yaml(claims_path)
  hyp <- load_hypotheses(hyp_path)                 # hypothesis, cell_type, contrast, program, gene
  prog_genes <- unique(hyp[, .(program, gene)])
  hyp_contrasts <- unique(hyp[, .(hypothesis, contrast)])
  rows <- list()
  for (hname in names(cfg$claims)) {
    for (rc in names(cfg$claims[[hname]])) {
      spec <- cfg$claims[[hname]][[rc]]
      cts  <- unlist(spec$cell_types %||% spec$parent_type)
      cons <- unlist(spec$contrasts %||% hyp_contrasts[hypothesis==hname, contrast])
      if (rc == "DE") {
        genes <- prog_genes[program %in% unlist(spec$programs), unique(gene)]
        rows[[length(rows)+1]] <- CJ(hypothesis=hname, readout_class="DE",
                                     unit=cts, contrast=cons, feature=genes)
      } else if (rc == "pathway") {
        rows[[length(rows)+1]] <- CJ(hypothesis=hname, readout_class="pathway",
                                     unit=cts, contrast=cons, feature=unlist(spec$programs))
      } else {  # composition / substate_composition: feature = NA
        rows[[length(rows)+1]] <- CJ(hypothesis=hname, readout_class=rc,
                                     unit=cts, contrast=cons, feature=NA_character_)
      }
    }
  }
  rbindlist(rows, use.names=TRUE)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

assert_no_multiclaim <- function(claims) {
  dup <- claims[, .N, by=.(readout_class, unit, contrast, feature)][N > 1]
  if (nrow(dup)) stop("claims map a result to >1 hypothesis:\n",
                      paste(capture.output(print(dup)), collapse="\n"))
  invisible(TRUE)
}

apply_claims <- function(master, claims) {
  master <- copy(master)
  if (!"hypothesis" %in% names(master)) master[, hypothesis := NA_character_]
  key <- unique(claims[, .(readout_class, unit, contrast, feature, hypothesis)])
  master[key, on=.(readout_class, unit, contrast, feature), hypothesis := i.hypothesis]
  master[]
}
