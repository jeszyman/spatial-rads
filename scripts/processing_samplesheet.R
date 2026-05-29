#!/usr/bin/env Rscript
# Workflow-linked sample sheet: scoped view of data_model.rda for the processing workflow.
# Args: <data_model.rda> <out.tsv>
suppressMessages({library(dplyr); library(readr)})
args <- commandArgs(trailingOnly = TRUE)
RDA <- args[1]; OUT <- args[2]
load(RDA)  # -> data_model (named list of tibbles)
dm <- data_model
scoped <- dm$samples %>%
  left_join(dm$slides,   by = "slide_id") %>%
  left_join(dm$datasets, by = "dataset_id") %>%
  mutate(raw_input_path = if_else(format == "rds", input_path, counts_path)) %>%
  transmute(sample_id, dataset = name, format, treatment, condition, timepoint_h, model,
            slide_id, physical_barcode, block_label, extract_key,
            counts_path, metadata_path, raw_input_path)
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
write_tsv(scoped, OUT)
cat(sprintf("wrote %s (%d samples: %d Mutter_01, %d Mutter_02)\n",
            OUT, nrow(scoped), sum(scoped$dataset=="Mutter_01"), sum(scoped$dataset=="Mutter_02")))
