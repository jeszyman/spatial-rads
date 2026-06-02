#!/usr/bin/env Rscript
# Build the data model from the relational metadata.xlsx, validating against
# metadata_schema.yaml. Emits data/data_model.rda (named list of tibbles) and the
# master denormalized sample sheet results/data_model/samples.tsv.
# Standalone:  conda run -n spatial-rads Rscript scripts/make_data_model.R
# Snakemake:   via the data_model rule's script: directive
suppressMessages({library(readxl); library(dplyr); library(readr); library(yaml)})
`%||%` <- function(a, b) if (is.null(a)) b else a

# positional args: xlsx schema rda samplesheet (fall back to standard paths for standalone runs)
args   <- commandArgs(trailingOnly = TRUE)
XLSX   <- if (length(args) >= 1) args[1] else "data/metadata.xlsx"
SCHEMA <- if (length(args) >= 2) args[2] else "data/metadata_schema.yaml"
RDA    <- if (length(args) >= 3) args[3] else "data/data_model.rda"
TSV    <- if (length(args) >= 4) args[4] else "results/data_model/samples.tsv"

schema <- yaml::read_yaml(SCHEMA)
dm <- lapply(schema$resources, function(r) readxl::read_excel(XLSX, sheet = r$control$sheet))
names(dm) <- vapply(schema$resources, function(r) r$name, character(1))

## ---- schema-driven validation ----
err <- character(0)
norm_na <- function(x, mv) { x <- as.character(x); x[x %in% mv] <- NA; x }
for (r in schema$resources) {
  d <- dm[[r$name]]; mv <- r$schema$missingValues %||% character(0)
  for (f in r$schema$fields) {
    if (!f$name %in% names(d)) { err <- c(err, sprintf("%s: missing column '%s'", r$name, f$name)); next }
    v <- norm_na(d[[f$name]], mv); cons <- f$constraints
    if (isTRUE(cons$required) && anyNA(v))
      err <- c(err, sprintf("%s.%s: %d missing required value(s)", r$name, f$name, sum(is.na(v))))
    if (isTRUE(cons$unique) && any(duplicated(v[!is.na(v)])))
      err <- c(err, sprintf("%s.%s: duplicate value(s)", r$name, f$name))
    if (!is.null(cons$pattern)) {
      bad <- v[!is.na(v) & !grepl(cons$pattern, v)]
      if (length(bad)) err <- c(err, sprintf("%s.%s: %d value(s) fail pattern %s", r$name, f$name, length(bad), cons$pattern))
    }
    if (!is.null(cons$enum)) {
      bad <- setdiff(na.omit(unique(v)), unlist(cons$enum))
      if (length(bad)) err <- c(err, sprintf("%s.%s: value(s) outside enum: %s", r$name, f$name, paste(bad, collapse = ", ")))
    }
  }
  for (fk in r$foreignKeys %||% list()) {
    child  <- norm_na(dm[[r$name]][[fk$fields]], mv)
    parent <- dm[[fk$reference$resource]][[fk$reference$fields]]
    miss   <- setdiff(na.omit(unique(child)), parent)
    if (length(miss)) err <- c(err, sprintf("%s.%s -> %s.%s unresolved: %s",
      r$name, fk$fields, fk$reference$resource, fk$reference$fields, paste(miss, collapse = ", ")))
  }
}
if (length(err)) stop("metadata_schema validation failed:\n  - ", paste(err, collapse = "\n  - "))

## ---- denormalized master sample sheet ----
## sample-grain only: join the dims that resolve 1:1 to a sample (mouse, slide, dataset).
## dataset-grain dims (e.g. if_channels, 5 rows/dataset) are deliberately NOT joined -- they
## would row-multiply the sheet. They live in data_model.rda; query them there by dataset_id.
master <- dm$samples %>%
  left_join(dm$mice,     by = "mouse_id") %>%
  left_join(dm$slides,   by = "slide_id") %>%
  left_join(dm$datasets, by = "dataset_id") %>%
  mutate(raw_input_path = if_else(format == "rds", input_path, counts_path)) %>%
  relocate(sample_id, dataset_id, name, format, treatment, condition, timepoint_h, model)

## ---- hard checks + summary ----
stopifnot(!anyNA(master$raw_input_path[master$format == "rds"]))         # every rds sample has a slide path
stopifnot(!anyNA(master$raw_input_path[master$format == "parquet"]))     # every parquet sample has counts path
cat(sprintf("data model OK: %d samples, %d mice, %d slides, %d datasets\n",
            nrow(dm$samples), nrow(dm$mice), nrow(dm$slides), nrow(dm$datasets)))
print(as.data.frame(count(master, name, treatment)))

data_model <- dm
dir.create(dirname(RDA), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(TSV), recursive = TRUE, showWarnings = FALSE)
save(data_model, file = RDA)
write_tsv(master, TSV)
cat(sprintf("wrote %s and %s\n", RDA, TSV))
