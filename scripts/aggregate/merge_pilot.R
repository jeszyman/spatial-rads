#!/usr/bin/env Rscript
# Stage 0 memory pilot for aggregate.smk. Loads a representative subset of flank
# scored.rds, performs the Seurat v5 sparse merge, records peak RSS (VmHWM), and
# projects to the full flank cohort to choose single-pass vs per-dataset merge.
# Args: <out_tsv> <flank_total_cells> <scored.rds...>
suppressPackageStartupMessages(library(Seurat))

args        <- commandArgs(trailingOnly = TRUE)
out_tsv     <- args[1]
flank_cells <- as.numeric(args[2])           # post-QC cells across all 20 flank samples
rds_paths   <- args[-(1:2)]

# Peak resident set size of this process (Linux high-water mark), kB -> GB.
peak_rss_gb <- function() {
  hwm <- grep("^VmHWM:", readLines("/proc/self/status"), value = TRUE)
  as.numeric(gsub("[^0-9]", "", hwm)) / 1024 / 1024
}

sids <- sub("\\.norm\\.rds$", "", basename(rds_paths))
objs <- lapply(rds_paths, readRDS)
names(objs) <- sids
n_cells <- sum(vapply(objs, ncol, numeric(1)))

invisible(gc(reset = TRUE))
merged <- merge(objs[[1]], y = objs[-1], add.cell.ids = sids)
rm(objs); invisible(gc())

peak    <- peak_rss_gb()
obj_gb  <- as.numeric(object.size(merged)) / 1e9
proj_gb <- peak / n_cells * flank_cells       # cell-based projection to full cohort

res <- data.frame(
  n_samples             = length(rds_paths),
  n_cells               = n_cells,
  peak_rss_gb           = round(peak, 2),
  final_object_gb       = round(obj_gb, 2),
  flank_total_cells     = flank_cells,
  projected_full_rss_gb = round(proj_gb, 2),
  decision              = ifelse(proj_gb < 100, "single_pass", "per_dataset_fallback"),
  stringsAsFactors      = FALSE
)
dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)
write.table(res, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("pilot: %d samples, %d cells | peak RSS %.2f GB | projected full %.2f GB -> %s\n",
            length(rds_paths), n_cells, peak, proj_gb, res$decision))
