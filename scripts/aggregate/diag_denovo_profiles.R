#!/usr/bin/env Rscript
# Cheap diagnostic: load the InSituType EM checkpoint only (no 3.27M-cell counts)
# and print marker expression per cluster from the UPDATED profiles, to decide
# whether the dominant de-novo clusters are real tumor (Epcam/Krt8 high) or
# misclassified immune (Cd3/Cd45 high). Read-only.
suppressPackageStartupMessages({ library(data.table) })
res <- readRDS("/mnt/data/projects/spatial-rads/aggregate/merged_typed.insitutype_res.rds")
prof <- res$profiles                       # genes x clusters (updated)
clust <- res$clust
sz <- sort(table(clust), decreasing = TRUE)

markers <- c(Epcam="Epcam", Krt8="Krt8", Krt18="Krt18",
             Cd3e="Cd3e", Cd3d="Cd3d", Cd45="Ptprc",
             Cd68="Cd68", Csf1r="Csf1r", Pecam1="Pecam1",
             Col1a1="Col1a1", Rgs5="Rgs5", Nkg7="Nkg7")
markers <- markers[markers %in% rownames(prof)]
cat("markers present:", paste(names(markers), collapse=" "), "\n\n")

mt <- t(prof[markers, , drop=FALSE])       # clusters x markers
dt <- data.table(cluster = rownames(mt), n = as.integer(sz[rownames(mt)]))
for (m in names(markers)) dt[[m]] <- round(mt[, markers[m]], 3)
setorder(dt, -n)
print(dt, nrows = 100)
cat("\ntotal cells:", length(clust), " | n clusters:", ncol(prof), "\n")
