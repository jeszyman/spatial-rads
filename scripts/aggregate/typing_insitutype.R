#!/usr/bin/env Rscript
# aggregate.smk typing -- field-standard CosMx cell typing via InSituType
# (Danaher et al. 2022), the platform-native semi-supervised classifier. Replaces
# the hand-rolled atlas-correlation annotate.R: InSituType models the CosMx noise
# structure (per-cell negative-probe background + Poisson counts) that flat
# correlation ignores, and runs ONCE on the merged M01+M02 counts so both cohorts
# get one identical rule -- structural comparability by construction.
#
# Method:
#   reference = the external atlas collapsed to 12 LINEAGE-mean profiles (the 93
#     fine profiles are badly collinear on a 950-gene panel -> they thrash the EM;
#     12 clean lineage anchors separate, and n_clusts de novo clusters absorb
#     within-lineage structure + the 4T1 tumor, which matches no normal profile).
#   neg  = per-cell mean negative-probe count (recover_negprobes.R), the background.
#   update_reference_profiles=TRUE adapts the (linear-scale, cross-platform) atlas
#     means to the CosMx data -- fixes the scale mismatch that hurt correlation.
#   semi-supervised: 12 fixed lineage types + n_clusts auto-selected de novo
#     clusters. De novo clusters expressing Epcam/Krt8 -> tumor_epithelial; the
#     rest -> mapped to their nearest lineage by the updated profile.
#
# Validation: gold-marker recall per lineage, split M01/M02 (honest yardstick).
# Args: <merged.rds> <ref_profiles.rds> <neg.tsv> <out_typed.rds> <out_summary.tsv>
#       <out_validation.tsv> <out_labels.tsv> [n_clusts_lo=6] [n_clusts_hi=12]
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(data.table); library(InSituType)
})
set.seed(1)
a <- commandArgs(trailingOnly = TRUE)
in_rds   <- a[1]; ref_rds <- a[2]; neg_tsv <- a[3]
out_typed<- a[4]; out_summary <- a[5]; out_valid <- a[6]; out_labels <- a[7]
nlo <- if (length(a) >= 8) as.integer(a[8]) else 6L
nhi <- if (length(a) >= 9) as.integer(a[9]) else 12L
EPI_GATE <- 0.25
for (d in unique(c(dirname(out_typed), dirname(out_summary), dirname(out_labels))))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- load ----
o <- readRDS(in_rds); DefaultAssay(o) <- "RNA"
o <- tryCatch(JoinLayers(o), error = function(e) o)
cnt <- LayerData(o, layer = "counts")                 # genes x cells, sparse
ref <- readRDS(ref_rds)
neg_dt <- fread(neg_tsv); negv <- setNames(neg_dt$neg, neg_dt$cell)
cat(sprintf("loaded %d cells x %d genes\n", ncol(cnt), nrow(cnt)))

# ---- collapse atlas to 12 lineage-mean profiles, align to panel, fill NA ----
lin <- ref$lineage; lin_levels <- sort(unique(lin))
prof <- sapply(lin_levels, function(L) rowMeans(ref$mat[, lin == L, drop = FALSE], na.rm = TRUE))
prof <- prof[rownames(cnt), , drop = FALSE]           # genes x 12, panel order
gmean <- rowMeans(prof, na.rm = TRUE); gmean[is.nan(gmean)] <- 0
na_ix <- which(is.na(prof), arr.ind = TRUE)           # uninformative fill: gene mean
if (nrow(na_ix)) prof[na_ix] <- gmean[na_ix[, 1]]
prof[!is.finite(prof)] <- 0
cat(sprintf("reference: %d genes x %d lineages (%s)\n", nrow(prof), ncol(prof),
            paste(lin_levels, collapse = ",")))

# ---- align neg to cells ----
neg <- negv[colnames(o)]
if (anyNA(neg)) { med <- median(neg, na.rm = TRUE); neg[is.na(neg)] <- med
  cat(sprintf("WARN %d cells missing neg -> median %.3f\n", sum(is.na(negv[colnames(o)])), med)) }

# ---- InSituType semi-supervised (cells x genes input) ----
x <- Matrix::t(cnt)                                   # cells x genes
res <- insitutype(
  x                          = x,
  neg                        = as.numeric(neg),
  reference_profiles         = prof,
  n_clusts                   = nlo:nhi,
  update_reference_profiles  = TRUE,
  n_starts                   = 5,
  align_genes                = TRUE)
saveRDS(res, sub("\\.rds$", ".insitutype_res.rds", out_typed))   # checkpoint: mapping is cheap, EM is not
clust <- res$clust; names(clust) <- colnames(o)
cat(sprintf("insitutype done: %d cells assigned across %d clusters\n",
            length(clust), length(unique(clust))))

# ---- map clusters -> cell_type_atlas ----
ref_types <- colnames(prof)
denovo    <- setdiff(unique(clust), ref_types)        # InSituType's new clusters
# de novo -> nearest lineage by cosine of the *updated* profile to the lineage anchors
updated   <- res$profiles                             # genes x (ref + denovo); FEWER genes than prof
gg        <- intersect(rownames(updated), rownames(prof))
profA     <- prof[gg, , drop = FALSE]
cos1 <- function(u, V) { u <- u / sqrt(sum(u^2))
  as.numeric(crossprod(sweep(V, 2, sqrt(colSums(V^2)), "/"), u)) }
denovo_lin <- setNames(rep(NA_character_, length(denovo)), denovo)
for (d in denovo) if (d %in% colnames(updated))
  denovo_lin[d] <- colnames(profA)[which.max(cos1(updated[gg, d], profA))]
lab <- ifelse(clust %in% ref_types, clust, denovo_lin[clust])

# ---- 4T1 overlay: de novo (or Epithelial) cluster expressing Epcam/Krt8 -> tumor ----
epi_genes <- intersect(c("Epcam", "Krt8"), rownames(cnt))
epi_det <- if (length(epi_genes)) {
  Matrix::colSums(cnt[epi_genes, , drop = FALSE] > 0) > 0
} else rep(FALSE, ncol(cnt))
epi_frac <- tapply(epi_det, clust, mean)
hot <- names(which(epi_frac >= EPI_GATE))
is_tumor <- clust %in% hot & lab %in% c(denovo_lin[denovo], "Epithelial")
lab[is_tumor] <- "tumor_epithelial"

o$cell_type_atlas <- unname(lab)
o$insitutype_clust <- unname(clust)
o$insitutype_prob  <- as.numeric(res$prob)            # InSituType preserves input cell order
atl <- o$cell_type_atlas; ds <- o$dataset

# ---- de novo cluster diagnostic (audit the force-mapping) ----
dn <- data.table(denovo = denovo,
  n               = as.integer(table(factor(clust, levels = denovo))[denovo]),
  nearest_lineage = denovo_lin[denovo],
  epcam_krt8_frac = round(as.numeric(epi_frac[denovo]), 3),
  final           = ifelse(denovo %in% hot, "tumor_epithelial", denovo_lin[denovo]))
cat("\n=== de novo cluster mapping ===\n"); print(dn)

# ---- composition summary: cell_type_atlas x dataset ----
tab <- table(atl, ds)
m01 <- if ("Mutter_01" %in% colnames(tab)) tab[, "Mutter_01"] else 0
m02 <- if ("Mutter_02" %in% colnames(tab)) tab[, "Mutter_02"] else 0
summ <- data.table(cell_type = rownames(tab),
  n_M01 = as.integer(m01), n_M02 = as.integer(m02),
  frac_M01 = round(as.numeric(m01) / sum(m01), 4),
  frac_M02 = round(as.numeric(m02) / sum(m02), 4),
  frac_overall = round(rowSums(tab) / sum(tab), 4))
setorder(summ, -frac_overall); fwrite(summ, out_summary, sep = "\t")

# ---- gold-marker recall validation (markers as truth, split M01/M02) ----
det <- function(g) { g <- intersect(g, rownames(cnt)); if (!length(g)) return(rep(0L, ncol(cnt)))
  Matrix::colSums(cnt[g, , drop = FALSE] > 0) }
gold <- list(
  tumor_epithelial = det(c("Epcam","Krt8"))>=2, T = det(c("Cd3e","Cd3d"))>=1,
  B = det(c("Cd79a","Ms4a1"))>=2, Plasma = det(c("Jchain","Sdc1","Xbp1"))>=2,
  NK = det("Ncr1")>=1, Macrophage = det(c("Cd68","Csf1r"))>=1,
  DC = det(c("Itgax","Flt3","Cd209a"))>=2, Neutrophil = det(c("S100a8","S100a9"))>=2,
  Endothelial = det(c("Pecam1","Cdh5"))>=2, Fibroblast = det(c("Col1a1","Dcn"))>=2,
  SmoothMuscle = det(c("Myh11","Acta2"))>=2, Pericyte = det(c("Rgs5","Pdgfrb","Kcnj8"))>=2)
m01c <- ds == "Mutter_01"
valid <- rbindlist(lapply(names(gold), function(k) { sub <- gold[[k]]
  if (!sum(sub, na.rm = TRUE)) return(NULL)
  data.table(lineage = k, gold_n = sum(sub),
    recall_pct = round(100*mean(atl[sub]==k),1),
    recall_M01_pct = round(100*mean(atl[sub & m01c]==k),1),
    recall_M02_pct = round(100*mean(atl[sub & !m01c]==k),1)) }))
fwrite(valid, out_valid, sep = "\t")

# ---- slim labels + typed object ----
fwrite(data.table(cell_id = colnames(o), dataset = ds, insitutype_clust = clust,
                  cell_type_atlas = atl, prob = round(o$insitutype_prob, 3)),
       out_labels, sep = "\t")
saveRDS(o, out_typed)

cat("\n=== composition (M01 vs M02) ===\n"); print(summ)
cat("\n=== gold-marker recall ===\n"); print(valid)
cat(sprintf("\ntyping: %d cells, %d types (%d de novo); typed.rds written.\n",
            ncol(o), uniqueN(atl), length(denovo)))
