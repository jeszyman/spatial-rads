#!/usr/bin/env Rscript
# Spatially proximal ligand-receptor communication (CellChat v2, Jin/Plikus/Nie Nat Protoc
# 2025) inferred separately per treatment arm within a cohort, then compared across arms.
#
# Coordinates are millimetres (x_slide_mm/y_slide_mm) while CellChat works in microns, so
# spatial.factors$ratio = 1000 um/mm. The CellChat FAQ's CosMx recipe (ratio = 0.12028)
# applies to raw AtoMx pixel coordinates; ours are already converted to mm upstream.
# Verified on this cohort: the median 1-NN centroid distance is 0.0072 mm = 7.2 um, a cell
# diameter, which is only self-consistent under the mm reading.
#
# Cells are thinned by whole spatial TILES, never at random. CellChat's spatial constraint
# is a nearest-neighbour statistic (computeRegionDistance -> queryKNN k=1), so random
# thinning silently inflates cell-type-to-cell-type distances: measured on sam0012, random
# 20% inflated them up to 5.05x while tile 20% held them at 0.95-1.00x. Random thinning
# also erases contact-dependent signalling (MHC-I) by pushing every pair past
# contact.range. Whole tiles preserve local density, so distances stay physical.
#
# The 950-gene panel covers only a small fraction of CellChatDB, so per-pathway on-panel
# interaction counts are carried on every result row: thin pathways stay visible rather
# than being silently reported as absent signalling.
#
# Args: <merged.rds> <full_labels.parquet> <coords.parquet> <obs.parquet> <samples.tsv>
#       <cohort_samples.tsv> <cohort> <out_prefix>
suppressPackageStartupMessages({
  library(data.table); library(arrow); library(Matrix)
  library(CellChat); library(BiocNeighbors)
})

## Parameters -----------------------------------------------------------------
# Coordinates are mm; every CellChat distance argument below is microns.
UM_PER_COORD_UNIT <- 1000    # x_slide_mm/y_slide_mm -> microns; spatial.factors$ratio
INTERACTION_RANGE <- as.numeric(Sys.getenv("CELLCHAT_INTERACTION_RANGE", "250"))
MIN_CELLS         <- 10L     # CellChat filterCommunication / identifyOverExpressedGenes default
K_MIN             <- 10L     # min interacting cell pairs defining a spatially proximal pair
NBOOT             <- as.integer(Sys.getenv("CELLCHAT_NBOOT", "100"))
TRIM              <- 0.1     # truncatedMean trim, per the spatial tutorials
SEED              <- 1L
TILE_MM           <- as.numeric(Sys.getenv("CELLCHAT_TILE_MM", "0.5"))   # 2x interaction range
MAX_CELLS_PER_ARM <- as.integer(Sys.getenv("CELLCHAT_MAX_CELLS_PER_ARM", "150000"))
MIN_PANEL_LR      <- 3L      # on-panel interactions below which a pathway is flagged thin
WORKERS           <- as.integer(Sys.getenv("CELLCHAT_WORKERS", "1"))

## Helpers --------------------------------------------------------------------
req_cols <- function(dt, cols, what) {
  miss <- setdiff(cols, names(dt))
  if (length(miss))
    stop(sprintf("%s is missing required column(s): %s", what, paste(miss, collapse = ", ")),
         call. = FALSE)
  invisible(TRUE)
}

nonempty <- function(x, what) {
  n <- if (is.null(x)) 0L else nrow(x)
  if (is.null(n) || n == 0L)
    stop(sprintf("%s resolved to 0 rows; refusing to write an empty table", what), call. = FALSE)
  invisible(TRUE)
}

say <- function(...) cat(sprintf(...), sep = "")

# Expand a CellChatDB ligand/receptor name into its constituent genes; a multi-subunit
# complex only counts as on-panel when every subunit is measured.
db_units <- function(x, complex_tbl) {
  if (x %in% rownames(complex_tbl)) {
    s <- unlist(complex_tbl[x, ], use.names = FALSE)
    s <- s[!is.na(s) & nzchar(s)]
    if (length(s)) return(as.character(s))
  }
  x
}

## Input assembly -------------------------------------------------------------
read_inputs <- function(labels_pq, coords_pq, obs_pq, samples_tsv, cohort_tsv, coh) {
  cs <- fread(cohort_tsv)
  req_cols(cs, c("cohort", "sample_id"), basename(cohort_tsv))
  sids <- unique(cs$sample_id[cs$cohort == coh])
  if (!length(sids))
    stop(sprintf("cohort '%s' has no samples in %s; known cohorts: %s", coh,
                 basename(cohort_tsv), paste(sort(unique(cs$cohort)), collapse = ", ")),
         call. = FALSE)

  smp <- fread(samples_tsv)
  req_cols(smp, c("sample_id", "condition", "slide_id", "dataset_id"), basename(samples_tsv))
  smp <- smp[sample_id %in% sids, .(sample_id, arm = condition, slide_id, dataset_id)]
  miss <- setdiff(sids, smp$sample_id)
  if (length(miss))
    stop(sprintf("cohort '%s' names sample(s) absent from %s: %s", coh,
                 basename(samples_tsv), paste(miss, collapse = ", ")), call. = FALSE)

  obs <- as.data.table(read_parquet(obs_pq, col_select = c("cell", "sample_id", "condition")))
  req_cols(obs, c("cell", "sample_id", "condition"), basename(obs_pq))
  obs <- obs[sample_id %in% sids]
  nonempty(obs, sprintf("cohort '%s' rows in %s", coh, basename(obs_pq)))
  # samples.tsv is the authoritative arm assignment; obs must agree or the metadata is stale
  chk <- merge(unique(obs[, .(sample_id, obs_condition = condition)]),
               smp[, .(sample_id, arm)], by = "sample_id")
  bad <- chk[obs_condition != arm]
  if (nrow(bad))
    stop(sprintf("condition disagrees between %s and %s for: %s", basename(samples_tsv),
                 basename(obs_pq), paste(bad$sample_id, collapse = ", ")), call. = FALSE)

  lab <- as.data.table(read_parquet(labels_pq,
                                    col_select = c("cell", "cell_subtype", "compartment")))
  req_cols(lab, c("cell", "cell_subtype", "compartment"), basename(labels_pq))
  # Epithelial cells (1,180 cohort-wide) is a small marker-rescue alias of Tumor in
  # full_labels.parquet; collapse to the single-Tumor roster used elsewhere (verify_arm_tables.R
  # retired list) so it is not analysed as a near-empty separate cell group.
  lab[cell_subtype == "Epithelial cells", cell_subtype := "Tumor"]

  co <- as.data.table(read_parquet(coords_pq,
                                   col_select = c("cell", "x_slide_mm", "y_slide_mm", "necrosis_zone")))
  req_cols(co, c("cell", "x_slide_mm", "y_slide_mm", "necrosis_zone"), basename(coords_pq))

  d <- merge(obs[, .(cell, sample_id)], lab, by = "cell")
  d <- merge(d, co, by = "cell")
  d <- merge(d, smp, by = "sample_id")
  # Inferred signaling from necrotic tissue is not interpretable; drop before any other filter.
  n_necrotic <- sum(d$necrosis_zone, na.rm = TRUE)
  d <- d[!is.na(cell_subtype) & !is.na(x_slide_mm) & !is.na(y_slide_mm) & !necrosis_zone]
  say("cohort '%s': dropped %d necrotic cell(s) (necrosis_zone == TRUE)\n", coh, n_necrotic)
  nonempty(d, sprintf("cohort '%s' cell table", coh))
  if (uniqueN(d$arm) < 2L)
    stop(sprintf("cohort '%s' resolves to %d arm(s) (%s); the cross-arm comparison needs >= 2",
                 coh, uniqueN(d$arm), paste(sort(unique(d$arm)), collapse = ", ")), call. = FALSE)
  d[]
}

## Spatial geometry -----------------------------------------------------------
# Cell diameter estimated as the median centroid-to-centroid 1-NN distance at FULL density.
# CellChat FAQ: for platforms without a fixed spot size compute the centroid distance and
# use half of it as `tol`; contact.range is "approximately the estimated cell diameter".
estimate_cell_geometry <- function(d) {
  g <- d[, {
    xy <- as.matrix(.SD[, .(x_slide_mm, y_slide_mm)])
    nn <- suppressWarnings(findKNN(xy, k = 1, BNPARAM = AnnoyParam(), get.index = FALSE))
    .(n_cells = .N, diam_um = median(nn$distance[, 1]) * UM_PER_COORD_UNIT)
  }, by = sample_id]
  if (any(!is.finite(g$diam_um) | g$diam_um <= 0))
    stop("cell-diameter estimation produced a non-positive distance; check the coordinates",
         call. = FALSE)
  g[]
}

## Tile thinning --------------------------------------------------------------
# Keep every cell inside a seeded selection of whole square tiles until the per-sample
# budget is met, so local density and therefore nearest-neighbour distances stay physical.
tile_subsample <- function(d, max_per_arm, tile_mm, seed = SEED) {
  set.seed(seed)
  d <- copy(d)
  d[, tile := paste(floor(x_slide_mm / tile_mm), floor(y_slide_mm / tile_mm), sep = ":")]
  keep <- list()
  for (a in sort(unique(d$arm))) {
    da <- d[arm == a]
    budget <- max(1L, as.integer(floor(max_per_arm / uniqueN(da$sample_id))))
    for (s in sort(unique(da$sample_id))) {
      ds <- da[sample_id == s]
      if (nrow(ds) <= budget) { keep[[length(keep) + 1L]] <- ds; next }
      tn <- ds[, .N, by = tile][sample(.N)]      # seeded tile order
      tn[, cum := cumsum(N)]
      ntile <- max(1L, sum(tn$cum <= budget))
      keep[[length(keep) + 1L]] <- ds[tile %in% tn$tile[seq_len(ntile)]]
    }
  }
  out <- rbindlist(keep)
  nonempty(out, "tile-thinned cell table")
  out[]
}

## Ligand-receptor panel coverage ---------------------------------------------
lr_panel_coverage <- function(db_full, db_used, panel, coh) {
  cx <- db_full$complex
  ints <- as.data.table(db_full$interaction)
  req_cols(ints, c("interaction_name", "pathway_name", "ligand", "receptor", "annotation"),
           "CellChatDB$interaction")
  cache <- new.env(hash = TRUE, parent = emptyenv())
  units_on_panel <- function(nm) {
    hit <- cache[[nm]]
    if (!is.null(hit)) return(hit)
    v <- all(db_units(nm, cx) %in% panel)
    assign(nm, v, envir = cache)
    v
  }
  lig <- as.character(ints$ligand); rec <- as.character(ints$receptor)
  ints[, on_panel := vapply(seq_len(.N),
                            function(k) units_on_panel(lig[k]) && units_on_panel(rec[k]),
                            logical(1))]
  ints[, in_db_subset := interaction_name %in% unique(as.data.table(db_used$interaction)$interaction_name)]

  cov <- ints[, .(n_interactions_db         = .N,
                  n_interactions_on_panel   = sum(on_panel),
                  n_interactions_db_subset  = sum(in_db_subset),
                  n_interactions_analysable = sum(on_panel & in_db_subset)),
              by = .(pathway_name, annotation)]
  cov[, frac_on_panel := round(n_interactions_on_panel / n_interactions_db, 4)]
  cov[, adequate_panel_support := n_interactions_analysable >= MIN_PANEL_LR]
  cov[, cohort := coh]
  setorder(cov, -n_interactions_analysable, -n_interactions_on_panel, pathway_name)
  nonempty(cov, "ligand-receptor coverage table")
  list(coverage = cov[], total_db = nrow(ints), total_on_panel = sum(ints$on_panel),
       total_analysable = sum(ints$on_panel & ints$in_db_subset))
}

## Per-arm CellChat object ----------------------------------------------------
arm_meta <- function(da, levels_use) {
  data.frame(labels  = factor(da$cell_subtype, levels = levels_use),
             samples = factor(da$sample_id, levels = sort(unique(da$sample_id))),
             row.names = da$cell)
}

# ratio and tol are indexed positionally over levels(samples) inside computeRegionDistance,
# so the rows must follow that level order.
arm_spatial_factors <- function(da, geom) {
  lv <- sort(unique(da$sample_id))
  i <- match(lv, geom$sample_id)
  if (anyNA(i))
    stop("no cell-diameter estimate for sample(s): ", paste(lv[is.na(i)], collapse = ", "),
         call. = FALSE)
  sf <- data.frame(ratio = rep(UM_PER_COORD_UNIT, length(lv)), tol = geom$diam_um[i] / 2)
  rownames(sf) <- lv
  sf
}

# The group-by-group distance matrix CellChat computes internally, precomputed so one
# scale.distance can be shared: the docs require the same scale factor across objects
# that are later compared.
arm_region_distance <- function(da, geom, levels_use, contact_range) {
  m <- arm_meta(da, levels_use)
  sf <- arm_spatial_factors(da, geom)
  computeRegionDistance(coordinates = as.matrix(da[, .(x_slide_mm, y_slide_mm)]),
                        meta = data.frame(group = m$labels, samples = m$samples,
                                          row.names = rownames(m)),
                        interaction.range = INTERACTION_RANGE,
                        ratio = sf$ratio, tol = sf$tol, k.min = K_MIN,
                        contact.dependent = TRUE, contact.range = contact_range,
                        contact.knn.k = NULL)
}

# CellChat stops unless the minimum off-diagonal scaled distance is >= 1 and recommends
# it land in [1,2].
pick_scale_distance <- function(d_min_um) {
  # geometric ladder stepping by <=1.67, so some value always lands inside the 2x-wide
  # acceptance window for any minimum distance between 1 um and 1 mm
  ladder <- c(1, 0.7, 0.5, 0.3, 0.2, 0.15, 0.1, 0.07, 0.05, 0.03, 0.02, 0.015, 0.01,
              0.007, 0.005, 0.003, 0.002, 0.0015, 0.001)
  ok <- ladder[d_min_um * ladder >= 1 & d_min_um * ladder <= 2]
  if (!length(ok)) {
    s <- 1.5 / d_min_um
    warning(sprintf("no ladder value puts the scaled minimum distance in [1,2] (min = %.3f um); using %.5g",
                    d_min_um, s), call. = FALSE)
    return(s)
  }
  max(ok)
}

run_arm <- function(da, expr, geom, levels_use, db_used, contact_range, scale_distance, arm_id) {
  cells <- da$cell
  m <- arm_meta(da, levels_use)
  coords <- as.matrix(da[, .(x_slide_mm, y_slide_mm)])
  rownames(coords) <- cells

  cc <- createCellChat(object = expr[, cells, drop = FALSE], meta = m, group.by = "labels",
                       datatype = "spatial", coordinates = coords,
                       spatial.factors = arm_spatial_factors(da, geom))
  cc@DB <- db_used
  cc <- subsetData(cc)
  # FAQ, "Application to a dataset with a small panel of genes": skip the DE-based gene
  # selection, which would discard most of an already thin ligand-receptor set.
  cc <- identifyOverExpressedGenes(cc, do.DE = FALSE, min.cells = MIN_CELLS)
  cc <- identifyOverExpressedInteractions(cc)
  if (is.null(cc@LR$LRsig) || nrow(cc@LR$LRsig) == 0L)
    stop(sprintf("arm '%s': no ligand-receptor pair survived the panel intersection", arm_id),
         call. = FALSE)
  say("  arm %s: %d cells, %d L-R pairs usable\n", arm_id, length(cells), nrow(cc@LR$LRsig))

  cc <- computeCommunProb(cc, type = "truncatedMean", trim = TRIM,
                          distance.use = TRUE, interaction.range = INTERACTION_RANGE,
                          scale.distance = scale_distance,
                          contact.dependent = TRUE, contact.range = contact_range,
                          k.min = K_MIN, nboot = NBOOT, seed.use = SEED)
  cc <- filterCommunication(cc, min.cells = MIN_CELLS)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  # Centrality feeds the plotting helpers only; a failure there must not sink the run.
  cc <- tryCatch(netAnalysis_computeCentrality(cc, slot.name = "netP"),
                 error = function(e) { warning(sprintf("arm '%s': centrality skipped (%s)",
                                                       arm_id, conditionMessage(e)),
                                               call. = FALSE); cc })
  cc
}

## Result extraction ----------------------------------------------------------
# computeCommunProbPathway zeroes any L-R pair with permutation pval > 0.05 before summing
# to the pathway, so netP$prob is already significance-filtered and CellChat stores no
# pathway-level pval. n_lr_sig / min_lr_pval report how deep that evidence actually runs.
extract_netP <- function(cc, arm_id, coh, census) {
  prob <- cc@netP$prob
  if (is.null(prob) || length(dim(prob)) != 3L || dim(prob)[3] == 0L)
    stop(sprintf("arm '%s': no pathway-level communication was inferred", arm_id), call. = FALSE)
  dn <- dimnames(prob)
  grid <- expand.grid(source = dn[[1]], target = dn[[2]], pathway_name = dn[[3]],
                      stringsAsFactors = FALSE)
  dt <- data.table(cohort = coh, arm = arm_id, grid, prob = as.vector(prob))
  dt <- dt[prob > 0]
  nonempty(dt, sprintf("arm '%s' netP table", arm_id))

  lp <- cc@net$pval; lpr <- cc@net$prob
  lrn <- dimnames(lp)[[3]]
  ints <- as.data.table(cc@DB$interaction)[, .(interaction_name, pathway_name)]
  lr_pw <- ints$pathway_name[match(lrn, ints$interaction_name)]
  if (anyNA(lr_pw))
    stop(sprintf("arm '%s': %d L-R pair(s) could not be mapped to a pathway", arm_id,
                 sum(is.na(lr_pw))), call. = FALSE)
  ev <- rbindlist(lapply(dn[[3]], function(p) {
    k <- which(lr_pw == p)
    sig <- (lp[, , k, drop = FALSE] < 0.05) & (lpr[, , k, drop = FALSE] > 0)
    pmin_ <- apply(lp[, , k, drop = FALSE], c(1, 2), min)
    data.table(expand.grid(source = dn[[1]], target = dn[[2]], stringsAsFactors = FALSE),
               pathway_name = p, n_lr_sig = as.vector(apply(sig, c(1, 2), sum)),
               min_lr_pval = as.vector(pmin_))
  }))
  dt <- merge(dt, ev, by = c("source", "target", "pathway_name"), all.x = TRUE)

  n <- census[census$arm == arm_id, .(cell_subtype, n_cells)]
  dt[n, n_cells_source := i.n_cells, on = c(source = "cell_subtype")]
  dt[n, n_cells_target := i.n_cells, on = c(target = "cell_subtype")]
  dt[]
}

# Count of L-R pairs actually carried into inference for this arm, by pathway.
arm_lr_used <- function(cc, arm_id) {
  lr <- as.data.table(cc@LR$LRsig)
  out <- lr[, .(n_interactions_used = .N), by = pathway_name]
  out[, arm := arm_id][]
}

# CellChat's cross-condition comparison: information flow per pathway per arm, with a
# Wilcoxon test whose observations are the cell-type PAIRS, not biological replicates.
compare_arms <- function(obj_list, coh) {
  arms <- names(obj_list)
  out <- list()
  for (i in seq_along(arms)) for (j in seq_along(arms)) {
    if (j <= i) next
    merged <- mergeCellChat(obj_list[c(i, j)], add.names = arms[c(i, j)])
    r <- rankNet(merged, slot.name = "netP", measure = "weight", mode = "comparison",
                 comparison = c(1, 2), do.stat = TRUE, paired.test = TRUE,
                 return.data = TRUE, stacked = FALSE)
    df <- as.data.table(r$signaling.contribution)
    req_cols(df, c("name", "contribution", "group"), "rankNet signaling.contribution")
    df[, name := as.character(name)]
    df[, group := as.character(group)]
    w <- dcast(df, name ~ group, value.var = "contribution", fun.aggregate = sum)
    req_cols(w, arms[c(i, j)], "rankNet contribution matrix")
    setnames(w, c("name", arms[i], arms[j]), c("pathway_name", "flow_arm1", "flow_arm2"))
    pv <- if ("pvalues" %in% names(df)) {
      unique(df[, .(pathway_name = name, pval_wilcox = as.numeric(pvalues))])
    } else {
      unique(df[, .(pathway_name = name, pval_wilcox = NA_real_)])
    }
    w <- merge(w, unique(pv, by = "pathway_name"), by = "pathway_name")
    w[, `:=`(cohort = coh, arm1 = arms[i], arm2 = arms[j],
             log2FC_flow = log2((flow_arm2 + 1e-12) / (flow_arm1 + 1e-12)))]
    # rankNet writes pvalues = 0 when fewer than four informative cell-type pairs exist;
    # that is "untestable", not "significant".
    w[, pval_untestable := !is.na(pval_wilcox) & pval_wilcox == 0]
    out[[length(out) + 1L]] <- w
  }
  res <- rbindlist(out, fill = TRUE)
  nonempty(res, "cross-arm comparison table")
  res[]
}

## Main -----------------------------------------------------------------------
main <- function(a) {
  if (length(a) < 8L)
    stop("usage: cellchat_communication.R <merged.rds> <full_labels.parquet> <coords.parquet> ",
         "<obs.parquet> <samples.tsv> <cohort_samples.tsv> <cohort> <out_prefix>", call. = FALSE)
  merged_rds <- a[1]; labels_pq <- a[2]; coords_pq <- a[3]; obs_pq <- a[4]
  samples_tsv <- a[5]; cohort_tsv <- a[6]; coh <- a[7]; out_prefix <- a[8]

  if (WORKERS > 1L) {
    future::plan("multisession", workers = WORKERS)
    options(future.globals.maxSize = 8 * 1024^3)
  }

  # Everything cheap and fail-fast runs before the merged object is read.
  d <- read_inputs(labels_pq, coords_pq, obs_pq, samples_tsv, cohort_tsv, coh)
  say("cohort %s: %d cells, %d samples, arms: %s\n", coh, nrow(d), uniqueN(d$sample_id),
      paste(sort(unique(d$arm)), collapse = ", "))

  geom <- estimate_cell_geometry(d)
  contact_range <- median(geom$diam_um)
  say("cell diameter (median 1-NN centroid distance) per sample: %s um; contact.range = %.2f um\n",
      paste(sprintf("%.2f", geom$diam_um), collapse = ", "), contact_range)

  full_n <- d[, .(n_full = .N), by = .(arm, cell_subtype)]
  d <- tile_subsample(d, MAX_CELLS_PER_ARM, TILE_MM)
  kept <- d[, .N, by = arm]
  say("tile thinning (%.2f mm tiles, budget %d cells/arm): %s\n", TILE_MM, MAX_CELLS_PER_ARM,
      paste(sprintf("%s=%d", kept$arm, kept$N), collapse = " "))

  census <- merge(d[, .(n_cells = .N), by = .(arm, cell_subtype)], full_n,
                  by = c("arm", "cell_subtype"), all.y = TRUE)
  census[is.na(n_cells), n_cells := 0L]
  census[, retention := round(n_cells / n_full, 4)]
  census[, cohort := coh]

  # A cell type must clear CellChat's documented minimum in EVERY arm: otherwise the arms
  # carry different factor levels and computeCommunProb/rankNet cannot be compared.
  ok <- census[, .(ok = all(n_cells >= MIN_CELLS)), by = cell_subtype]
  keep <- ok[ok == TRUE, cell_subtype]
  dropped <- setdiff(sort(unique(census$cell_subtype)), keep)
  if (length(dropped)) {
    say("DROPPED %d cell type(s) below the CellChat minimum of %d cells in >=1 arm: %s\n",
        length(dropped), MIN_CELLS, paste(dropped, collapse = ", "))
    print(census[cell_subtype %in% dropped][order(cell_subtype, arm)])
  }
  if (length(keep) < 2L)
    stop(sprintf("only %d cell type(s) clear the %d-cell minimum in every arm; nothing to compare",
                 length(keep), MIN_CELLS), call. = FALSE)
  if (length(dropped) > length(keep))
    stop(sprintf("more cell types dropped (%d) than retained (%d); the budget of %d cells per arm is too small",
                 length(dropped), length(keep), MAX_CELLS_PER_ARM), call. = FALSE)
  d <- d[cell_subtype %in% keep]
  census <- census[cell_subtype %in% keep]
  levels_use <- sort(keep)

  # Read the merged object once, take the log-normalised layer, subset to the cohort cells
  # and CellChatDB genes immediately, then drop it. Densifying the full 950 x 3.28M matrix
  # would allocate ~25 GB, so every step here stays sparse.
  say("reading %s ...\n", merged_rds)
  obj <- readRDS(merged_rds)
  expr <- SeuratObject::LayerData(obj, assay = "RNA", layer = "data")
  if (is.null(expr) || !nrow(expr)) stop("merged object has no RNA 'data' layer", call. = FALSE)
  panel <- rownames(expr)
  miss <- setdiff(d$cell, colnames(expr))
  if (length(miss))
    stop(sprintf("%d cohort cell(s) absent from the merged object, e.g. %s", length(miss),
                 paste(head(miss, 3), collapse = ", ")), call. = FALSE)
  db_full <- CellChatDB.mouse
  db_used <- subsetDB(db_full)   # Secreted Signaling + ECM-Receptor + Cell-Cell Contact
  db_genes <- unique(c(
    unlist(lapply(c(as.character(db_full$interaction$ligand),
                    as.character(db_full$interaction$receptor)),
                  db_units, complex_tbl = db_full$complex)),
    unlist(db_full$cofactor, use.names = FALSE)))
  expr <- expr[intersect(panel, db_genes), d$cell, drop = FALSE]
  rm(obj); gc(FALSE)
  say("expression: %d CellChatDB genes x %d cells retained from a %d-gene panel\n",
      nrow(expr), ncol(expr), length(panel))

  cov <- lr_panel_coverage(db_full, db_used, panel, coh)
  say("panel coverage: %d/%d CellChatDB interactions fully on panel (%.1f%%); %d analysable after DB subsetting\n",
      cov$total_on_panel, cov$total_db, 100 * cov$total_on_panel / cov$total_db,
      cov$total_analysable)

  # One scale.distance for every arm, as the docs require for cross-object comparison.
  say("calibrating the spatial distance scale ...\n")
  dmins <- numeric(0)
  for (a_ in sort(unique(d$arm))) {
    ds <- arm_region_distance(d[arm == a_], geom, levels_use, contact_range)$d.spatial
    diag(ds) <- NaN
    dmins <- c(dmins, min(ds, na.rm = TRUE))
  }
  d_min_um <- min(dmins)
  if (!is.finite(d_min_um) || d_min_um <= 0)
    stop("no spatially proximal cell-type pair was found; check the coordinate units",
         call. = FALSE)
  scale_distance <- pick_scale_distance(d_min_um)
  say("minimum off-diagonal cell-type distance %.2f um -> scale.distance = %g (scaled min %.2f)\n",
      d_min_um, scale_distance, d_min_um * scale_distance)

  objs <- list(); netp <- list(); used <- list()
  for (a_ in sort(unique(d$arm))) {
    cc <- run_arm(d[arm == a_], expr, geom, levels_use, db_used, contact_range,
                  scale_distance, a_)
    objs[[a_]] <- cc
    netp[[length(netp) + 1L]] <- extract_netP(cc, a_, coh, census)
    used[[length(used) + 1L]] <- arm_lr_used(cc, a_)
  }

  cvj <- cov$coverage[, .(pathway_name, annotation, n_interactions_db, n_interactions_on_panel,
                          n_interactions_analysable, adequate_panel_support)]
  netp <- merge(rbindlist(netp), cvj, by = "pathway_name", all.x = TRUE)
  netp <- merge(netp, rbindlist(used), by = c("pathway_name", "arm"), all.x = TRUE)
  netp[is.na(n_interactions_used), n_interactions_used := 0L]
  setcolorder(netp, c("cohort", "arm", "pathway_name", "source", "target", "prob",
                      "n_lr_sig", "min_lr_pval", "n_cells_source", "n_cells_target"))
  setorder(netp, cohort, arm, pathway_name, source, target)

  diff <- merge(compare_arms(objs, coh), cvj, by = "pathway_name", all.x = TRUE)
  setcolorder(diff, c("cohort", "arm1", "arm2", "pathway_name", "flow_arm1", "flow_arm2",
                      "log2FC_flow", "pval_wilcox", "pval_untestable"))
  setorder(diff, cohort, arm1, arm2, -flow_arm2)

  f_cov  <- paste0(out_prefix, "_lr_coverage.tsv")
  f_netp <- paste0(out_prefix, "_netP_by_arm.tsv")
  f_diff <- paste0(out_prefix, "_diff_arms.tsv")
  f_rds  <- paste0(out_prefix, "_objects.rds")
  dir.create(dirname(f_cov), recursive = TRUE, showWarnings = FALSE)
  nonempty(cov$coverage, f_cov); nonempty(netp, f_netp); nonempty(diff, f_diff)
  fwrite(cov$coverage, f_cov, sep = "\t")
  fwrite(netp, f_netp, sep = "\t")
  fwrite(diff, f_diff, sep = "\t")
  saveRDS(list(objects = objs, census = census, coverage = cov$coverage,
               params = list(cohort = coh, ratio_um_per_coord_unit = UM_PER_COORD_UNIT,
                             interaction_range_um = INTERACTION_RANGE,
                             contact_range_um = contact_range, scale_distance = scale_distance,
                             min_cells = MIN_CELLS, k_min = K_MIN, nboot = NBOOT, trim = TRIM,
                             tile_mm = TILE_MM, max_cells_per_arm = MAX_CELLS_PER_ARM,
                             seed = SEED, cell_diameter_um = geom,
                             dropped_cell_types = dropped,
                             cellchat_version = as.character(utils::packageVersion("CellChat")))),
          f_rds)
  say("wrote %s (%d rows), %s (%d), %s (%d), %s\n", f_cov, nrow(cov$coverage), f_netp,
      nrow(netp), f_diff, nrow(diff), f_rds)
  invisible(TRUE)
}

a <- commandArgs(trailingOnly = TRUE)
if (length(a) >= 8L) main(a)
