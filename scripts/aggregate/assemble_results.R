#!/usr/bin/env Rscript
# aggregate.smk T13: tier-tagged master results table with tiered FDR.
# Collates every cohort's differential readout (composition, pseudobulk DE, GSEA,
# pathway program scores, niche composition, spatial mixing, myeloid polarization,
# fibroblast substate composition, smiDE per-cell DE) into one long table, one row per
# readout x unit x feature x contrast x comparison. Every readout family except the
# day-2-only muscat differential-detection step is parameterized over the full cohort
# roster (mutter02_day2, combined_4h_treated, combined_4h, mutter02_day2_pooledctrl,
# combined_4h_pooledctrl): each producer writes one file per cohort (e.g.
# de_engine_<cohort>.tsv), and this script globs + rbinds them per readout family,
# deriving -- or, for the two producers that don't emit one themselves, injecting --
# a `comparison` column from the filename.
#
# Each row is tagged tier = confirmatory|exploratory via a frozen a-priori claims-join
# (the results-tier config config/confirmatory_claims.yaml declares each hypothesis's
# pre-registered evidence; producers emit hypothesis=NA and the join fills it, so
# nothing is retro pattern-matched). The claims join keys on
# (readout_class, unit, contrast, feature) only, not comparison -- this stays correct
# because every confirmatory claim's contrasts come from hypotheses.yaml
# (SBRT_vs_Ctrl / MBRT_vs_Ctrl / MBRT_vs_SBRT), and the comparison registry gives those
# three contrast-name strings only to the confirmatory day-2 cohort; every other
# cohort's contrast names carry a _4h / _pooled / _d2pooled suffix. A hard assertion
# below re-verifies that invariant at run time rather than trusting it silently.
#
# Primary confirmatory FDR = pooled BH over the confirmatory family -> padj_confirmatory.
# Two auxiliary FDR views ride alongside: gatekeeping (gate_padj/padj_gated) and IHW,
# fit per comparison (padj_ihw). Effect size + 95% CI carried through; the MDE floor from
# power_mde.tsv is design-matched per (comparison, contrast) and joined on
# (comparison, contrast, unit), so each cohort is measured against its own n/arm and
# formula; rows with no matching floor get NA mde/abs_effect_lt_mde/clears_mde/trend_call.
# n_per_arm and dataset are likewise derived per (comparison, contrast) from the
# registry / cohort_samples.tsv rather than assumed uniform across the cohort roster.
#
# The six pre-registered hypotheses and their evidence live in config/hypotheses.yaml
# (biological definition) + config/confirmatory_claims.yaml (readout->confirmatory map).
# Args: <agg_dir> <comparisons.tsv> <power_mde.tsv> <pathway_sets.tsv>
#       <gene_set_panel_coverage.tsv> <differential_detection.tsv> <overlap_ratio_qc.tsv>
#       <out_master.tsv>
# agg_dir must contain engine/{composition,de,niche,mixing,myeloid,substate}_<cohort>.tsv
# plus engine/voom_engine_<cohort>.tsv for the cohorts whose reported DE inference is
# the moderated t (see the DE block below)
# plus engine/smide_de_<cohort>_{screen,spatial}.tsv, and gsea_pseudobulk_<cohort>.tsv
# and pathway_test_<cohort>.tsv at its top level.
suppressPackageStartupMessages({ library(data.table) })

# --- gatekeeping (secondary FDR view) -------------------------------------------
# Auxiliary to the pooled-primary padj_confirmatory. A gate exists only where a
# program is claimed as gene-level confirmatory evidence (a DE claim). Parent = the
# matching GSEA enrichment row (one BH family over gate-eligible programs); child =
# the claimed DE genes, BH'd within each opened program's cell_type x contrast.
# Programs with no DE claim create no gate. Defined at top so `source(<this file>)`
# exposes it to the gatekeeping test without running the pipeline body.
source("scripts/aggregate/fdr_helpers.R")   # add_gatekeeping()

a <- commandArgs(trailingOnly = TRUE)
agg_p <- a[1]; reg_p <- a[2]; mde_p <- a[3]; sets_p <- a[4]
cov_p <- a[5]; ddet_p <- a[6]; overlap_p <- a[7]; out_p <- a[8]
engine_p <- file.path(agg_p, "engine")

# Panel coverage table for the symmetric per-program coverage columns below.
cov_dt <- fread(cov_p)                                   # set, tier, source, n_total, n_panel, usable, thin

dir.create(dirname(out_p), recursive = TRUE, showWarnings = FALSE)
COLS <- c("comparison","readout_class","unit","feature","contrast","effect","effect_type",
          "ci_low","ci_high","se","stat","pvalue","padj_own","hypothesis","tier",
          "n_per_arm","n_samples_used","dataset")

# --- cohort roster, read from the resolved registry rather than a copy of
# Snakemake's DE_COHORTS/LM_COHORTS lists -- single source of truth for both the glob
# whitelist below and the confirmatory-leak assertion. ---
REG <- fread(reg_p)
KNOWN_COHORTS <- unique(REG[kind == "sample", cohort])
CONFIRMATORY_COHORTS <- unique(REG[kind == "sample" & tier == "confirmatory", cohort])
stopifnot(length(CONFIRMATORY_COHORTS) > 0)

# --- glob-and-rbind one readout family across the cohort roster --------------------
# Producers name their per-cohort file "<prefix><cohort>.tsv". A file whose extracted
# cohort token isn't in KNOWN_COHORTS is a stray/legacy artifact (a pre-rename file
# under an old cohort abbreviation, or a count_engine "*_skipped.tsv" sibling) and is
# silently skipped -- tolerant of the pre-flight reality that not every cohort has run
# yet, but loud if a whole family matches nothing real. Injects `comparison` from the
# filename for the two producers (gsea.R, pathway_arm_test.R) that don't emit the
# column themselves; for producers that do (every engine script), asserts the file's
# own value agrees with the filename.
# `suffix` covers producers that emit more than one file per cohort (smide_de
# writes one per model mode, "<prefix><cohort>_<mode>.tsv"); it is stripped along
# with the prefix before the cohort token is checked.
# `required = FALSE` lets a family that legitimately covers only part of the roster
# (the moderated-t engine runs on the low-df 4h cohorts only) return an empty table
# instead of erroring.
read_cohort_family <- function(dir, prefix, label, suffix = "", required = TRUE) {
  files <- list.files(dir, pattern = paste0("^", prefix, ".*", suffix, "\\.tsv$"),
                      full.names = TRUE)
  hits <- list()
  for (f in files) {
    coh <- sub("\\.tsv$", "", sub(paste0("^", prefix), "", basename(f)))
    if (nzchar(suffix)) coh <- sub(paste0(suffix, "$"), "", coh)
    if (!coh %in% KNOWN_COHORTS) next
    dt <- fread(f)
    if ("comparison" %in% names(dt)) {
      stopifnot(all(dt$comparison == coh))
    } else {
      dt[, comparison := coh]
    }
    hits[[coh]] <- dt
  }
  if (length(hits) == 0) {
    if (!required) {
      cat(sprintf("assemble_results: %-22s no cohort files (optional family)\n", label))
      return(data.table())
    }
    stop(sprintf("assemble_results: readout family '%s' matched zero cohort files (pattern '%s*' in %s)",
                 label, prefix, dir))
  }
  cat(sprintf("assemble_results: %-22s %d cohort file(s): %s\n",
              label, length(hits), paste(names(hits), collapse = ", ")))
  rbindlist(hits, use.names = TRUE)
}

# Engine outputs carry sufficient statistics only (estimate/se/df/stat/p) -- correction is
# recomputed HERE, reproducing each source's BH grouping (extended with `comparison` so
# BH is never pooled across cohorts). Engine `unit`/`feature_id` semantics differ by
# engine: lm_engine puts the feature in feature_id (unit="mouse"), while count_engine
# puts the cell-type stratum in `unit` and the gene in feature_id. The adapter below
# reads the correct column per readout. (Schema-harmonization deferred to daylight.)

# --- composition (propeller logit log2FC; lm_engine proportion path) -------------
comp <- read_cohort_family(engine_p, "composition_engine_", "composition")
comp[, padj := p.adjust(p, "BH"), by = comparison]       # composition.R used global BH, within cohort
comp_m <- comp[, .(comparison, readout_class="composition", unit=feature_id, feature=NA_character_,
  contrast, effect=estimate, effect_type="log2FC_logit",
  ci_low=estimate - qt(0.975, df)*se, ci_high=estimate + qt(0.975, df)*se,
  se, stat, pvalue=p, padj_own=padj,
  hypothesis=NA_character_, n_samples_used=NA_integer_, baseMean=NA_real_)]

# --- substate composition (Fibroblast resting/activated; lm_engine proportion path) --
# Same schema as composition; unit is the sub-state label (e.g. Fibroblast_activated),
# which the myofibroblast_expansion claim targets.
sub <- read_cohort_family(engine_p, "substate_engine_", "substate_composition")
sub[, padj := p.adjust(p, "BH"), by = comparison]
sub_m <- sub[, .(comparison, readout_class="substate_composition", unit=feature_id, feature=NA_character_,
  contrast, effect=estimate, effect_type="log2FC_logit",
  ci_low=estimate - qt(0.975, df)*se, ci_high=estimate + qt(0.975, df)*se,
  se, stat, pvalue=p, padj_own=padj,
  hypothesis=NA_character_, n_samples_used=NA_integer_, baseMean=NA_real_)]

# --- pseudobulk DE (unit=cell_type, feature_id=gene) ------------------------------
# Two engines fit the SAME registry design on the SAME pseudobulk SE with the same
# filters: count_engine (DESeq2, Wald test) and voom_engine (limma-voom, moderated t).
# Which one is REPORTED is a per-cohort property of the realized residual df, and the
# workflow declares it by running voom_engine for those cohorts only. Cohorts with a
# voom_engine file report the moderated t; every other cohort reports the Wald test.
# The DESeq2 Wald estimate and p ride along on the swapped rows as
# effect_deseq2_wald / pvalue_deseq2_wald so the two are never confused, and baseMean
# (the IHW covariate, a DESeq2 quantity) is carried over from the same rows.
de_wald <- read_cohort_family(engine_p, "de_engine_", "DE (DESeq2 Wald)")
de_mod  <- read_cohort_family(engine_p, "voom_engine_", "DE (moderated t)", required = FALSE)
MOD_T_COHORTS <- if (nrow(de_mod)) unique(de_mod$comparison) else character()
if (length(MOD_T_COHORTS))
  cat(sprintf("assemble_results: reported DE inference = limma-voom moderated t for: %s\n",
              paste(sort(MOD_T_COHORTS), collapse = ", ")))
DEKEY <- c("comparison", "contrast", "unit", "feature_id")
wald_aux <- de_wald[comparison %in% MOD_T_COHORTS,
                    .(comparison, contrast, unit, feature_id,
                      effect_deseq2_wald = estimate, pvalue_deseq2_wald = p,
                      baseMean_deseq2 = baseMean)]
de <- rbindlist(list(de_wald[!comparison %in% MOD_T_COHORTS], de_mod),
                use.names = TRUE, fill = TRUE)
de[, inference_method := fifelse(comparison %in% MOD_T_COHORTS,
                                 "limma_voom_moderated_t", "DESeq2_Wald")]
de <- wald_aux[de, on = DEKEY]
de[is.na(baseMean), baseMean := baseMean_deseq2]   # voom rows carry no DESeq2 baseMean
de[, baseMean_deseq2 := NULL]
de[, padj := p.adjust(p, "BH"), by = .(comparison, unit, contrast)]  # per cohort x cell_type x contrast
de_m <- de[, .(comparison, readout_class="DE", unit=unit, feature=feature_id, contrast,
  effect=estimate, effect_type="log2FC", ci_low=estimate-1.96*se,
  ci_high=estimate+1.96*se, se, stat, pvalue=p, padj_own=padj,
  hypothesis=NA_character_, n_samples_used=n_samples_used, baseMean=baseMean,
  inference_method, effect_deseq2_wald, pvalue_deseq2_wald)]

# --- pathway program scores (UCell primary; limma estimate, normal CI) -----------
# padj_bh is already correctly scoped: pathway_arm_test.R computes it once per cohort
# invocation, so each per-cohort file's own BH never saw another cohort's rows -- no
# recompute needed here, just pass through.
pw <- read_cohort_family(agg_p, "pathway_test_", "pathway")[score_type == "UCell"]
pw_m <- pw[, .(comparison, readout_class="pathway", unit=cell_type, feature=pathway_name, contrast,
  effect=estimate, effect_type="limma_score_estimate", ci_low=estimate-1.96*se,
  ci_high=estimate+1.96*se, se, stat=t_stat, pvalue, padj_own=padj_bh,
  hypothesis=NA_character_, n_samples_used=n_samples_per_group, baseMean=NA_real_)]

# --- GSEA (Hallmark sweep, always exploratory; NES, no CI) -----------------------
# padj_bh is already correctly scoped: fgsea corrects within each (cell_type x
# contrast) ranking inside a single cohort's gsea.R invocation -- no recompute needed.
gsea <- read_cohort_family(agg_p, "gsea_pseudobulk_", "gsea")
gsea_m <- gsea[, .(comparison, readout_class="gsea", unit=cell_type, feature=pathway_name, contrast,
  effect=NES, effect_type="NES", ci_low=NA_real_, ci_high=NA_real_, se=NA_real_,
  stat=NES, pvalue, padj_own=padj_bh, hypothesis=NA_character_,
  n_samples_used=NA_integer_, baseMean=NA_real_)]

# --- niche composition (lm_engine proportion path; exploratory) -----------------
ni <- read_cohort_family(engine_p, "niche_engine_", "niche_composition")
ni[, padj := p.adjust(p, "BH"), by = comparison]          # niches.R used global BH, within cohort
ni_m <- ni[, .(comparison, readout_class="niche_composition", unit=as.character(feature_id),
  feature=NA_character_, contrast, effect=estimate, effect_type="log2FC_logit",
  ci_low=estimate - qt(0.975, df)*se, ci_high=estimate + qt(0.975, df)*se,
  se, stat, pvalue=p, padj_own=padj, hypothesis=NA_character_,
  n_samples_used=NA_integer_, baseMean=NA_real_)]

# --- spatial mixing + myeloid polarization (lm_engine identity path; exploratory) --
mx <- read_cohort_family(engine_p, "mixing_engine_", "spatial_mixing")
mx[, padj := p.adjust(p, "BH"), by = comparison]
mx_m <- mx[, .(comparison, readout_class="spatial_mixing", unit="global", feature=feature_id, contrast,
  effect=estimate, effect_type="limma_estimate",
  ci_low=estimate - qt(0.975, df)*se, ci_high=estimate + qt(0.975, df)*se, se,
  stat, pvalue=p, padj_own=padj, hypothesis=NA_character_,
  n_samples_used=NA_integer_, baseMean=NA_real_)]
my <- read_cohort_family(engine_p, "myeloid_engine_", "myeloid_polarization")
my[, padj := p.adjust(p, "BH"), by = comparison]
my_m <- my[, .(comparison, readout_class="myeloid_polarization", unit="Macrophages", feature=feature_id,
  contrast, effect=estimate, effect_type="limma_estimate",
  ci_low=estimate - qt(0.975, df)*se, ci_high=estimate + qt(0.975, df)*se, se,
  stat, pvalue=p, padj_own=padj, hypothesis=NA_character_,
  n_samples_used=NA_integer_, baseMean=NA_real_)]

# --- differential detection (muscat: fraction of cells expressing gene changed) ------
# Detection-only descriptor (never composition-vs-regulation); always exploratory.
# Stays a single day-2-only file (no per-cohort DD in this build).
dd <- fread(ddet_p)
dd_m <- dd[, .(comparison="mutter02_day2", readout_class="detection", unit=cell_type, feature=gene, contrast,
  effect=dd_log2fc, effect_type="detection_log2fc", ci_low=NA_real_, ci_high=NA_real_,
  se=NA_real_, stat=NA_real_, pvalue=dd_p, padj_own=dd_padj,
  hypothesis=NA_character_, n_samples_used=NA_integer_, baseMean=NA_real_)]

# --- smiDE genome-wide per-cell NB mixed-model DE ---------------------------------
# smide_de.R writes one file per (cohort, model mode). The two modes are not
# interchangeable evidence and enter the master under different readout_classes:
#   screen  -> "smiDE_screen". Sample random intercept, no spatial random effect;
#              anti-conservative for a contrast that varies between samples, so it
#              is a discovery screen only. No confirmatory claim keys on it.
#   spatial -> "smiDE_spatial". Per-spatial-unit Gaussian-process fits combined
#              by inverse-variance meta-analysis (mode == "spatial_meta"); the
#              per-unit rows (mode == "spatial") are the meta-analysis inputs and
#              are not cohort-level results, so they are dropped here. The
#              default backend (GP_INLA) reports posterior summaries, so a
#              spatial row's pvalue is a normal approximation to its credible
#              interval, not a frequentist test.
# The spatial fitting geometry (fit_unit, n_fits, n_cells_fit) rides through to
# the master: a spatial effect is only readable next to how many units combined
# into it and how many cells each fit actually saw after the subsample cap.
smide_screen  <- read_cohort_family(engine_p, "smide_de_", "smiDE (screen)",  suffix = "_screen")
smide_spatial <- read_cohort_family(engine_p, "smide_de_", "smiDE (spatial)", suffix = "_spatial")
if (nrow(smide_spatial) > 0 && "mode" %in% names(smide_spatial)) {
  no_meta <- setdiff(unique(smide_spatial$comparison),
                     unique(smide_spatial[mode == "spatial_meta", comparison]))
  if (length(no_meta) > 0)
    stop(sprintf(paste0("assemble_results: smiDE spatial cohort(s) %s contributed a file with no ",
                        "spatial_meta row -- the per-unit fits never combined into a cohort result"),
                 paste(no_meta, collapse = ", ")))
}
smide <- rbindlist(list(smide_screen, smide_spatial), use.names = TRUE, fill = TRUE)
NEED_SMIDE_COLS <- c("comparison","contrast","unit","feature_id","estimate","se","stat","p","mode")
if (nrow(smide) == 0 || !all(NEED_SMIDE_COLS %in% names(smide))) {
  smide_m <- data.table(comparison=character(), readout_class=character(), unit=character(),
    feature=character(), contrast=character(), effect=numeric(), effect_type=character(),
    ci_low=numeric(), ci_high=numeric(), se=numeric(), stat=numeric(), pvalue=numeric(),
    padj_own=numeric(), hypothesis=character(), n_samples_used=integer(), baseMean=numeric(),
    fit_unit=character(), n_fits=integer(), n_cells_fit=integer())
} else {
  smide <- smide[mode != "spatial"]
  smide[, readout_class := fifelse(mode == "screen", "smiDE_screen", "smiDE_spatial")]
  smide[, padj := p.adjust(p, "BH"), by = .(comparison, readout_class, unit, contrast)]
  smide_m <- smide[, .(comparison, readout_class, unit=unit, feature=feature_id, contrast,
    effect=estimate, effect_type="log2FC", ci_low=estimate-1.96*se,
    ci_high=estimate+1.96*se, se, stat, pvalue=p, padj_own=padj,
    hypothesis=NA_character_, n_samples_used=NA_integer_, baseMean=NA_real_,
    fit_unit, n_fits, n_cells_fit)]
}

# fill = TRUE because only the smiDE rows carry the three spatial-geometry
# columns; every other readout family leaves them NA.
master <- rbindlist(list(comp_m, sub_m, de_m, pw_m, gsea_m, ni_m, mx_m, my_m, dd_m, smide_m),
                    use.names = TRUE, fill = TRUE)

# --- n_per_arm + dataset: derived per (comparison, contrast) / cohort, not hardcoded --
# n_per_arm was a blanket 4L (the day-2 n before the 2026-08-13 M02 relabel correction)
# on every readout family; each cohort now has its own registry-declared arm size
# (day-2 n=2/arm, combined_4h n=2, combined_4h_treated n=3, pooled variants n=2/n=3).
# Joined from the registry's `sample`-kind rows, deduped across the whole/compartment
# resolution duplicates (which always carry identical n_group1/n_group2).
reg_n <- unique(REG[kind == "sample", .(comparison = cohort, contrast = name, n_group1, n_group2)])
stopifnot(!anyDuplicated(reg_n, by = c("comparison", "contrast")))
master <- reg_n[master, on = .(comparison, contrast)]
master[, n_per_arm := pmin(n_group1, n_group2)]
master[, c("n_group1", "n_group2") := NULL]
n_no_reg_n <- master[is.na(n_per_arm), .N]
if (n_no_reg_n > 0)
  cat(sprintf("assemble_results: %d row(s) had no registry (comparison,contrast) match -- n_per_arm left NA\n",
              n_no_reg_n))

# dataset was a blanket "Mutter_02"; the combined_4h* cohorts actually pool samples
# from both raw datasets. Derived from which raw dataset(s) contributed samples to
# each cohort (cohort_samples.tsv -> samples.tsv$name), so it reflects cohort
# composition directly rather than a per-readout-family assertion.
cohort_samples <- fread("results/data_model/cohort_samples.tsv")
samp_names <- fread("results/data_model/samples.tsv")[, .(sample_id, name)]
cohort_dataset <- unique(samp_names[cohort_samples, on = "sample_id"][
  , .(dataset = paste(sort(unique(name)), collapse = "+")), by = cohort])
setnames(cohort_dataset, "cohort", "comparison")
master <- cohort_dataset[master, on = "comparison"]
n_no_dataset <- master[is.na(dataset), .N]
if (n_no_dataset > 0)
  cat(sprintf("assemble_results: %d row(s) had no cohort_samples.tsv entry -- dataset left NA\n",
              n_no_dataset))

# --- segmentation-contamination QC (smiDE overlap ratio) -------------------------
# Per-gene, per-cell_subtype: does neighboring-cell-type expression exceed self
# expression (ratio >= 1 = contamination-dominated in that cell type)? Only
# applies to per-gene per-cell-type readouts (DE rows); all other readout_classes
# (composition, niche, mixing, pathway, GSEA) get contamination_ratio = NA. The QC
# table is gene x cell_subtype only (cohort-independent), so it joins onto every
# cohort's DE rows alike.
orm <- fread(overlap_p)
setnames(orm, old = c("gene", "cell_subtype", "ratio"),
         new = c("feature", "unit", "contamination_ratio"),
         skip_absent = TRUE)
master <- orm[, .(feature, unit, contamination_ratio)][master, on = .(feature, unit)]

# --- confirmatory tagging by frozen a-priori claims (not retro pattern-match) -----
# Producers emit hypothesis=NA; the frozen results-tier config declares which
# (readout_class, unit, contrast, feature) rows are each hypothesis's pre-registered
# evidence, and the claims-join fills `hypothesis` for exactly those rows.
source("scripts/aggregate/load_claims.R")
CLAIMS <- load_claims("config/confirmatory_claims.yaml", "config/hypotheses.yaml")
HYP    <- load_hypotheses("config/hypotheses.yaml")   # program<->gene<->cell_type<->contrast, for gatekeeping
assert_no_multiclaim(CLAIMS)
master <- apply_claims(master, CLAIMS)

# Collapse the tumor-compartment unit alias to the single locked roster name.
master[unit == "Epithelial cells", unit := "Tumor"]

master[, tier := fifelse(!is.na(hypothesis), "confirmatory", "exploratory")]

# --- confirmatory-tag scope guard (cross-cohort claims-join leak check) -----------
# The claims join above keys on (readout_class, unit, contrast, feature) only, not
# comparison. That is safe only because every confirmatory hypothesis's contrasts
# (hypotheses.yaml: SBRT_vs_Ctrl/MBRT_vs_Ctrl/MBRT_vs_SBRT) are contrast-name strings
# unique to CONFIRMATORY_COHORTS in the registry -- every other cohort's contrasts carry
# a _4h/_pooled/_d2pooled suffix (config/comparisons.yaml). Re-verify that invariant at
# run time: a confirmatory tag on any other cohort means a contrast-name collision
# slipped into the registry and the claims-join needs a comparison key added.
leaked <- master[tier == "confirmatory" & !comparison %in% CONFIRMATORY_COHORTS]
if (nrow(leaked) > 0)
  stop(sprintf(
    "assemble_results: %d confirmatory-tagged row(s) belong to a non-confirmatory comparison %s -- claims-join leaked across a contrast-name collision",
    nrow(leaked), paste(unique(leaked$comparison), collapse = ", ")))

# --- tiered FDR -----------------------------------------------------------------
# Primary confirmatory FDR = pooled BH over the confirmatory family (guarded above to
# stay within CONFIRMATORY_COHORTS). Exploratory rows keep their within-analysis BH in
# padj_own (no separate padj_exploratory: it was a verbatim copy of padj_own).
master[, padj_confirmatory := NA_real_]
master[tier == "confirmatory", padj_confirmatory := p.adjust(pvalue, method = "BH")]

# --- gatekeeping (secondary FDR view; padj_confirmatory above stays primary) -----
master <- add_gatekeeping(master, CLAIMS, HYP)

# --- IHW auxiliary FDR (keyed on DE baseMean; padj_confirmatory stays primary) ----
master <- add_ihw(master)

# --- join each row's own design-matched MDE from power_mde.tsv --------------------
# power_mde.tsv now carries one row per (comparison, contrast) x readout x unit, each
# built from that cohort's registry n/arm, residual df and formula. The join key is
# therefore (comparison, contrast, unit) -- keying on unit alone would hand a cohort a
# floor computed under a design that isn't its own (day-2 is blocked n=2/arm; the 4h
# cohorts are unblocked n=2-3/arm). Rows with no matching MDE row keep NA for all four
# MDE-derived columns rather than a spurious floor.
mde <- fread(mde_p)
mde_comp <- mde[readout_class == "composition",
                .(key = paste(comparison, contrast, cell_type), mde, mde_scale)]
mde_de   <- mde[readout_class == "pseudobulk_DE",
                .(key = paste(comparison, contrast, cell_type), mde, mde_scale)]
mde_prog <- mde[readout_class == "program_score",
                .(key = paste(comparison, contrast, cell_type, detail), mde, mde_scale)]
master[, `:=`(mde = NA_real_, mde_scale = NA_character_)]
mkey3 <- paste(master$comparison, master$contrast, master$unit)
mkey4 <- paste(mkey3, master$feature)
ic <- master$readout_class == "composition"
mm <- mde_comp[match(mkey3[ic], mde_comp$key)]
master[ic, `:=`(mde = mm$mde, mde_scale = mm$mde_scale)]
id <- master$readout_class == "DE"
mm <- mde_de[match(mkey3[id], mde_de$key)]
master[id, `:=`(mde = mm$mde, mde_scale = mm$mde_scale)]
ip <- master$readout_class == "pathway"
mm <- mde_prog[match(mkey4[ip], mde_prog$key)]
master[ip, `:=`(mde = mm$mde, mde_scale = mm$mde_scale)]
master[, abs_effect_lt_mde := NA]
master[is.finite(mde), abs_effect_lt_mde := abs(effect) < mde]


# --- magnitude-floored trend call + clears-MDE (devils-advocate fix) -------------
# A direction is only a "trend" when |effect| >= MDE_FLOOR * mde; below that it is
# noise, not signal (sub-MDE direction coherence is noise coherence). Rows with no
# MDE get NA, not "below-floor" -- that label asserts a floor was checked, which
# didn't happen here.
MDE_FLOOR <- 0.5
master[, clears_mde := NA]
master[is.finite(mde), clears_mde := abs(effect) >= mde]
master[, trend_call := NA_character_]
master[is.finite(mde), trend_call := fifelse(abs(effect) < MDE_FLOOR * mde, "below-floor",
                                      fifelse(effect > 0, "up", "down"))]
master[is.na(effect), trend_call := "na"]

# --- symmetric panel coverage on every program / gsea row (no set silently protected) ---
master[, `:=`(n_panel = NA_integer_, n_set_total = NA_integer_, panel_cov_frac = NA_real_)]
iset <- master$readout_class %in% c("pathway", "gsea")
mm   <- cov_dt[match(master$feature[iset], cov_dt$set)]
master[iset, `:=`(n_panel = mm$n_panel, n_set_total = mm$n_total,
                  panel_cov_frac = round(mm$n_panel / mm$n_total, 3))]

# --- interpretability + claim-scoping flags (Fix 1 unassigned, Fix 4 MBRT dose) ---
master[, interpretable := !(readout_class == "DE" & unit == "unassigned")]   # unassigned DE rows are grab-bag
master[, dose_confounded := contrast %in% c("MBRT_vs_SBRT", "MBRT_vs_SBRT_4h")]  # MBRT mean dose unrecorded
master[, independent := !(readout_class %in% c("pathway","gsea") & feature == "STING")]  # STING ~50% shares genes with IFN sets

setcolorder(master, c(COLS, "padj_confirmatory",
                      "inference_method", "effect_deseq2_wald", "pvalue_deseq2_wald",
                      "gate_program", "gate_padj", "padj_gated", "padj_ihw",
                      "mde", "mde_scale", "abs_effect_lt_mde", "clears_mde", "trend_call",
                      "n_panel", "n_set_total", "panel_cov_frac",
                      "fit_unit", "n_fits", "n_cells_fit"))
setorder(master, comparison, tier, readout_class, unit, contrast, feature, na.last = TRUE)
fwrite(master, out_p, sep = "\t")

# --- detectability summary (the transparent decision-gate: what is/isn't detectable) ---
det <- master[, .(n = .N, n_clears_mde = sum(clears_mde, na.rm = TRUE),
                  n_trend_up = sum(trend_call == "up", na.rm = TRUE),
                  n_trend_down = sum(trend_call == "down", na.rm = TRUE)),
              by = .(comparison, readout_class, contrast, hypothesis, tier)]
setorder(det, comparison, tier, readout_class, contrast, hypothesis, na.last = TRUE)
fwrite(det, sub("results_master", "detectability_summary", out_p), sep = "\t")

# --- inspect --------------------------------------------------------------------
cat(sprintf("results_master: %d rows | confirmatory %d / exploratory %d\n",
            nrow(master), master[tier=="confirmatory", .N],
            master[tier=="exploratory", .N]))
cat("rows missing effect/pvalue:", master[is.na(effect) | is.na(pvalue), .N], "\n")
cat("\nrows by comparison:\n")
print(master[, .N, by = comparison][order(comparison)])
cat("\nconfirmatory rows by readout_class x hypothesis:\n")
print(dcast(master[tier=="confirmatory"], readout_class ~ hypothesis,
            value.var = "pvalue", fun.aggregate = length))
cat("\nconfirmatory rows by hypothesis x comparison (every column must be a CONFIRMATORY_COHORTS member):\n")
print(dcast(master[tier=="confirmatory"], hypothesis ~ comparison,
            value.var = "pvalue", fun.aggregate = length))
cat(sprintf("\nconfirmatory hits at padj_confirmatory < 0.05: %d\n",
            master[padj_confirmatory < 0.05, .N]))
if (master[padj_confirmatory < 0.05, .N] > 0)
  print(master[padj_confirmatory < 0.05,
               .(comparison, readout_class, unit, feature, contrast, effect = round(effect,3),
                 padj_confirmatory = signif(padj_confirmatory,3))][order(padj_confirmatory)])
cat(sprintf("\nconfirmatory rows with |effect| < MDE (underpowered): %d / %d with an MDE\n",
            master[tier=="confirmatory" & abs_effect_lt_mde==TRUE, .N],
            master[tier=="confirmatory" & is.finite(mde), .N]))
