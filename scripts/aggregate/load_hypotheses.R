# scripts/aggregate/load_hypotheses.R
suppressPackageStartupMessages({library(yaml); library(data.table)})

load_hypotheses <- function(path) {
  y <- yaml::read_yaml(path)$hypotheses
  rows <- rbindlist(lapply(y, function(h) {
    prog_rows <- rbindlist(lapply(h$programs, function(p) {
      CJ_genes <- data.table(program = p$name, gene = unlist(p$genes))
      CJ_genes
    }))
    # cross cell_types x contrasts x (program,gene)
    grid <- CJ(cell_type = unlist(h$cell_types),
               contrast  = unlist(h$contrasts),
               idx = seq_len(nrow(prog_rows)))
    out <- cbind(hypothesis = h$name,
                 grid[, .(cell_type, contrast)],
                 prog_rows[grid$idx])
    out
  }))
  shas <- setNames(vapply(y, function(h) as.character(h$provenance$frozen_at), ""),
                   vapply(y, function(h) h$name, ""))
  setattr(rows, "frozen_shas", shas)
  rows[]
}

hypotheses_programs <- function(path) {
  dt <- load_hypotheses(path)
  unique(dt[, .(hypothesis, cell_type, contrast, program)])
}
