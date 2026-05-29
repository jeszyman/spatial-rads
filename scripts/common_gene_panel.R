#!/usr/bin/env Rscript
# Common gene panel = Mutter_01 panel (counts parquet columns) intersect Mutter_02 panel
# (rownames of one Mutter_02 Seurat RDS). Args: <samplesheet.tsv> <out.tsv>
suppressMessages({library(readr); library(dplyr); library(arrow)})
args <- commandArgs(trailingOnly = TRUE)
SS <- args[1]; OUT <- args[2]
ss <- read_tsv(SS, show_col_types = FALSE)

m01_counts <- ss %>% filter(dataset == "Mutter_01") %>% pull(counts_path) %>% unique()
m02_rds    <- ss %>% filter(dataset == "Mutter_02") %>% pull(raw_input_path) %>% unique()
stopifnot(length(m01_counts) >= 1, length(m02_rds) >= 1)

m01_genes <- setdiff(arrow::open_dataset(m01_counts[1])$schema$names, c("Slide", "fov", "cell_id"))
suppressMessages(library(Seurat))
o <- readRDS(m02_rds[1]); m02_genes <- rownames(o); rm(o); gc()

common <- sort(intersect(m01_genes, m02_genes))
cat(sprintf("Mutter_01 panel: %d | Mutter_02 panel: %d | common: %d\n",
            length(m01_genes), length(m02_genes), length(common)))
stopifnot(length(common) > 500)  # sanity: panels should overlap heavily
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
writeLines(common, OUT)
