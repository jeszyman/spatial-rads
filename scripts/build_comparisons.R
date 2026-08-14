#!/usr/bin/env Rscript
# Resolve the curated comparison registry against the master samplesheet into a
# validated design table. Each registry entry DECLARES intent (groups, unit, formula,
# tier); this builder COMPUTES realized facts (resolved n, distinct mice, model rank,
# residual df, inference-capability) so the experimental design is machine-checkable.
# Pre-registration: tier + hypotheses are frozen against git HEAD (silent change to the
# confirmatory family hard-errors). Stage-A metadata only -- never touches per-cell data.
# Args: <comparisons.yaml> <samples.tsv> <out_comparisons.tsv>
suppressPackageStartupMessages({library(yaml); library(data.table)})

`%||%` <- function(a, b) if (is.null(a)) b else a

a       <- commandArgs(trailingOnly = TRUE)
yaml_p  <- a[1]; ss_p <- a[2]; out_p <- a[3]

reg <- read_yaml(yaml_p)$comparisons
ss  <- fread(ss_p)
ss[, timepoint_h := suppressWarnings(as.integer(timepoint_h))]

# cohort label -> row filter on the samplesheet
# The two *_pooledctrl cohorts are pinned to an explicit sample_id whitelist (not a
# predicate) because their membership is a deliberate UNION of two normally-disjoint
# groups (a treated arm's usual cohort plus the OTHER timepoint's control samples) --
# a control-pooling sensitivity check, not a naturally-growing predicate.
cohort_rows <- function(cohort) {
  switch(cohort,
    mutter02_day2   = ss[name == "Mutter_02" & timepoint_h == 48],
    mutter02_4h     = ss[name == "Mutter_02" & timepoint_h == 4],
    combined_4h     = ss[timepoint_h == 4 & model == "flank"],
    combined_4h_treated = ss[timepoint_h == 4 & model == "flank" & treatment != "NT"],
    mutter01_flank  = ss[name == "Mutter_01" & model == "flank"],
    mutter01_tongue = ss[name == "Mutter_01" & model == "tongue"],
    mutter01_mbrt_4h= ss[name == "Mutter_01" & treatment == "MBRT" & timepoint_h == 4],
    # mutter02_day2 (sam0012-sam0017) + the two 4h-slide controls (sam0018, sam0021)
    mutter02_day2_pooledctrl = ss[sample_id %in% c(
      "sam0012", "sam0013", "sam0014", "sam0015", "sam0016", "sam0017",
      "sam0018", "sam0021")],
    # combined_4h whitelist (sam0003, sam0006, sam0018-sam0023) + the two day-2 controls
    # (sam0012, sam0015)
    combined_4h_pooledctrl = ss[sample_id %in% c(
      "sam0003", "sam0006", "sam0018", "sam0019", "sam0020", "sam0021", "sam0022", "sam0023",
      "sam0012", "sam0015")],
    stop("unknown cohort: ", cohort))
}

# apply a group predicate map (col -> value / list-of-values) to a cohort subset
apply_pred <- function(d, pred) {
  if (is.null(pred)) return(d[0])
  keep <- rep(TRUE, nrow(d))
  for (col in names(pred)) {
    if (!col %in% names(d)) stop("predicate column not in samplesheet: ", col)
    vals <- as.character(unlist(pred[[col]]))
    keep <- keep & (as.character(d[[col]]) %in% vals)
  }
  d[keep]
}

# resid df for `~0+condition+slide_id` on a resolved 2-group subset (else NA)
design_facts <- function(d, formula) {
  if (is.null(formula) || nrow(d) == 0) return(list(rank = NA_integer_, resid = NA_integer_))
  d <- droplevels(as.data.frame(d))
  d$condition <- factor(d$condition); d$slide_id <- factor(d$slide_id)
  mm <- tryCatch(model.matrix(as.formula(formula), d), error = function(e) NULL)
  if (is.null(mm)) return(list(rank = NA_integer_, resid = NA_integer_))
  r <- qr(mm)$rank
  list(rank = r, resid = nrow(mm) - r)
}

rows <- list()
for (e in reg) {
  cohorts <- e$cohorts %||% (e$cohort %||% NA_character_)
  resterms <- e$resolution %||% NA_character_
  for (ch in cohorts) for (res in resterms) {
    d  <- if (is.na(ch)) ss[0] else cohort_rows(ch)
    g1 <- apply_pred(d, e$group1); g2 <- apply_pred(d, e$group2)
    # Design facts reflect the model as ACTUALLY fit: one ~0+condition+slide_id fit over
    # the whole cohort (all arms), with each pairwise contrast extracted from it -- a
    # single-model contrast, not a per-pair 2-arm refit. So resid_df is the shared
    # full-cohort residual df (day-2 RCB: 12 samples, rank 6 -> 6), not the 2-arm subset.
    df <- design_facts(d, e$formula)
    unit <- e$unit %||% "sample"
    nu1  <- if (unit == "mouse") uniqueN(g1$mouse_id) else nrow(g1)
    nu2  <- if (unit == "mouse") uniqueN(g2$mouse_id) else nrow(g2)
    min_n <- e$min_samples %||% 3L
    inf  <- (e$tier %||% "") %in% c("confirmatory","exploratory") &&
            e$kind %in% c("sample","region") && nu1 >= min_n && nu2 >= min_n
    gate <- if (!is.null(e$requires)) "unchecked" else NA_character_
    rows[[length(rows)+1]] <- data.table(
      name = e$name, kind = e$kind, cohort = as.character(ch),
      resolution = as.character(res), unit = unit, tier = e$tier %||% NA_character_,
      formula = e$formula %||% NA_character_,
      hypotheses = paste(e$hypotheses %||% character(), collapse = ";"),
      dose_confounded = isTRUE(e$dose_confounded),
      requires = e$requires %||% NA_character_,
      # Resolved contrast levels for the test engines: engine builds condition<num> -
      # condition<den> with den as the DESeq2 reference. transform tags the proportion
      # (logit) vs continuous (identity) path; NA where the entry names no single contrast.
      contrast_num_level = if (!is.null(e$contrast_levels)) e$contrast_levels$num else NA_character_,
      contrast_den_level = if (!is.null(e$contrast_levels)) e$contrast_levels$den else NA_character_,
      ref_level          = if (!is.null(e$contrast_levels)) e$contrast_levels$den else NA_character_,
      transform          = e$transform %||% NA_character_,
      n_group1 = nrow(g1), n_group2 = nrow(g2),
      n_mouse_group1 = nu1, n_mouse_group2 = nu2,
      model_rank = df$rank, resid_df = df$resid,
      inference_capable = inf, gate = gate)
  }
}
out <- rbindlist(rows, fill = TRUE)

## ---- validation: sample/region entries with predicates must resolve ----
named_with_g1 <- sapply(Filter(function(e) !is.null(e$group1), reg), `[[`, "name")
bad <- out[kind %in% c("sample","region") & name %in% named_with_g1 & n_group1 == 0]
if (nrow(bad)) stop("comparison(s) resolved to zero group1 samples: ",
                    paste(unique(bad$name), collapse = ", "))

## ---- freeze: tier + hypotheses must match git HEAD ----
committed <- tryCatch(
  yaml::yaml.load(paste(system2("git", c("show", "HEAD:config/comparisons.yaml"),
                                stdout = TRUE), collapse = "\n")),
  error = function(e) NULL)
if (!is.null(committed)) {
  ch_tier <- setNames(lapply(committed$comparisons, function(x)
    list(tier = x$tier %||% NA, hyp = x$hypotheses %||% character())),
    sapply(committed$comparisons, `[[`, "name"))
  for (e in reg) {
    c0 <- ch_tier[[e$name]]
    if (!is.null(c0)) {
      if (!identical(e$tier %||% NA, c0$tier))
        stop("freeze: tier changed for '", e$name, "' vs committed comparisons.yaml")
      if (!setequal(e$hypotheses %||% character(), c0$hyp))
        stop("freeze: hypotheses changed for '", e$name, "' vs committed comparisons.yaml")
    }
  }
}

dir.create(dirname(out_p), recursive = TRUE, showWarnings = FALSE)
fwrite(out, out_p, sep = "\t")

## ---- emit per-cohort sample lists (consumed by engines for SE filtering) ----
cohort_samples <- rbindlist(lapply(unique(na.omit(out$cohort)), function(ch) {
  d <- cohort_rows(ch)
  data.table(cohort = ch, sample_id = d$sample_id)
}))
cohort_samples_path <- file.path(dirname(out_p), "cohort_samples.tsv")
fwrite(cohort_samples, cohort_samples_path, sep = "\t")

cat(sprintf("comparisons: %d rows from %d registry entries\n", nrow(out), length(reg)))
print(out[, .(name, cohort, resolution, n_mouse_group1, n_mouse_group2, resid_df, inference_capable, gate)])
