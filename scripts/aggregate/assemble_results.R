#!/usr/bin/env Rscript
# aggregate.smk T13: tier-tagged master results table with tiered FDR.
# Collates every M02 day-2 differential readout (composition, pseudobulk DE, GSEA,
# pathway program scores, niche composition, spatial mixing, myeloid polarization)
# into one long table, one row per readout x unit x feature x contrast. Each row is
# tagged tier = confirmatory|exploratory against a PRE-REGISTERED family (fixed lookup
# below, not inferred). BH is recomputed across the pooled confirmatory family ->
# padj_confirmatory; exploratory rows keep each analysis's own within-analysis BH ->
# padj_exploratory. Effect size + 95% CI carried through; the matching n=4 MDE from
# power_mde.tsv is joined (composition by cell_type, DE by cell_type, program by
# cell_type x set), and abs_effect_lt_mde flags underpowered confirmatory nulls.
#
# Pre-registered confirmatory family (each across all 3 contrasts):
#  H1 immune activation : composition of immune subtypes; pseudobulk DE within immune
#                         cell types restricted to TypeI/TypeII-IFN + STING genes;
#                         program scores for those 3 sets in immune cells.
#  H2 vascular/oxygen   : composition of endothelial; DE within endothelial restricted
#                         to Angiogenesis genes; Angiogenesis+Hypoxia programs in
#                         endothelial/tumor.
#  H3 stromal sparing   : composition of stromal subtypes; DE within fibroblast/SMC/
#                         adipocyte restricted to Fibrosis + Stromal-stress genes;
#                         those 2 programs in stroma.
# Args: <composition_test> <degs> <gsea> <pathway_test> <niche_test> <mixing_test>
#       <myeloid_test> <power_mde> <pathway_sets.tsv> <out_master>
suppressPackageStartupMessages({ library(data.table) })
a <- commandArgs(trailingOnly = TRUE)
comp_p <- a[1]; de_p <- a[2]; gsea_p <- a[3]; pw_p <- a[4]; ni_p <- a[5]
mx_p <- a[6]; my_p <- a[7]; mde_p <- a[8]; sets_p <- a[9]; cov_p <- a[10]; det_p <- a[11]; out_p <- a[12]

# Panel coverage table for the symmetric per-program coverage columns below.
cov_dt <- fread(cov_p)                                   # set, tier, source, n_total, n_panel, usable, thin

dir.create(dirname(out_p), recursive = TRUE, showWarnings = FALSE)
COLS <- c("readout_class","unit","feature","contrast","effect","effect_type",
          "ci_low","ci_high","se","stat","pvalue","padj_own","hypothesis","tier",
          "n_per_arm","dataset")

# Engine outputs carry sufficient statistics only (estimate/se/df/stat/p) -- correction is
# recomputed HERE, reproducing each source's BH grouping. Engine `unit`/`feature_id`
# semantics differ by engine: lm_engine puts the feature in feature_id (unit="mouse"), while
# count_engine puts the cell-type stratum in `unit` and the gene in feature_id. The adapter
# below reads the correct column per readout. (Schema-harmonization deferred to daylight.)

# --- composition (propeller logit log2FC; lm_engine proportion path) -------------
comp <- fread(comp_p)                                    # engine schema
comp[, padj := p.adjust(p, "BH")]                        # composition.R used global BH
comp_m <- comp[, .(readout_class="composition", unit=feature_id, feature=NA_character_,
  contrast, effect=estimate, effect_type="log2FC_logit",
  ci_low=estimate - qt(0.975, df)*se, ci_high=estimate + qt(0.975, df)*se,
  se, stat, pvalue=p, padj_own=padj,
  hypothesis=NA_character_, n_per_arm=4L, dataset="Mutter_02")]

# --- pseudobulk DE (count_engine; unit=cell_type, feature_id=gene; Wald CI from se) ----
de <- fread(de_p)
de[, padj := p.adjust(p, "BH"), by = .(unit, contrast)]  # DESeq2 grouping: per cell_type x contrast
de_m <- de[, .(readout_class="DE", unit=unit, feature=feature_id, contrast,
  effect=estimate, effect_type="log2FC", ci_low=estimate-1.96*se,
  ci_high=estimate+1.96*se, se, stat, pvalue=p, padj_own=padj,
  hypothesis=NA_character_, n_per_arm=n_samples_used, dataset="Mutter_02")]

# --- pathway program scores (UCell primary; limma estimate, normal CI) -----------
pw <- fread(pw_p)[score_type == "UCell"]
pw_m <- pw[, .(readout_class="pathway", unit=cell_type, feature=pathway_name, contrast,
  effect=estimate, effect_type="limma_score_estimate", ci_low=estimate-1.96*se,
  ci_high=estimate+1.96*se, se, stat=t_stat, pvalue, padj_own=padj_bh,
  hypothesis=NA_character_, n_per_arm=n_samples_per_group, dataset)]

# --- GSEA (Hallmark sweep, always exploratory; NES, no CI) -----------------------
gsea <- fread(gsea_p)
gsea_m <- gsea[, .(readout_class="gsea", unit=cell_type, feature=pathway_name, contrast,
  effect=NES, effect_type="NES", ci_low=NA_real_, ci_high=NA_real_, se=NA_real_,
  stat=NES, pvalue, padj_own=padj_bh, hypothesis=NA_character_,
  n_per_arm=NA_integer_, dataset)]

# --- niche composition (lm_engine proportion path; exploratory) -----------------
ni <- fread(ni_p)
ni[, padj := p.adjust(p, "BH")]                          # niches.R used global BH
ni_m <- ni[, .(readout_class="niche_composition", unit=as.character(feature_id),
  feature=NA_character_, contrast, effect=estimate, effect_type="log2FC_logit",
  ci_low=estimate - qt(0.975, df)*se, ci_high=estimate + qt(0.975, df)*se,
  se, stat, pvalue=p, padj_own=padj, hypothesis=NA_character_, n_per_arm=4L,
  dataset="Mutter_02")]

# --- spatial mixing + myeloid polarization (lm_engine identity path; exploratory) --
mx <- fread(mx_p)
mx[, padj := p.adjust(p, "BH")]
mx_m <- mx[, .(readout_class="spatial_mixing", unit="global", feature=feature_id, contrast,
  effect=estimate, effect_type="limma_estimate",
  ci_low=estimate - qt(0.975, df)*se, ci_high=estimate + qt(0.975, df)*se, se,
  stat, pvalue=p, padj_own=padj, hypothesis=NA_character_,
  n_per_arm=4L, dataset="Mutter_02")]
my <- fread(my_p)
my[, padj := p.adjust(p, "BH")]
my_m <- my[, .(readout_class="myeloid_polarization", unit="Macrophages", feature=feature_id,
  contrast, effect=estimate, effect_type="limma_estimate",
  ci_low=estimate - qt(0.975, df)*se, ci_high=estimate + qt(0.975, df)*se, se,
  stat, pvalue=p, padj_own=padj, hypothesis=NA_character_,
  n_per_arm=4L, dataset="Mutter_02")]

master <- rbindlist(list(comp_m, de_m, pw_m, gsea_m, ni_m, mx_m, my_m),
                    use.names = TRUE)

# --- confirmatory tagging by frozen a-priori claims (not retro pattern-match) -----
# Producers emit hypothesis=NA; the frozen results-tier config declares which
# (readout_class, unit, contrast, feature) rows are each hypothesis's pre-registered
# evidence, and the claims-join fills `hypothesis` for exactly those rows.
source("scripts/aggregate/load_claims.R")
CLAIMS <- load_claims("config/confirmatory_claims.yaml", "config/hypotheses.yaml")
assert_no_multiclaim(CLAIMS)
master <- apply_claims(master, CLAIMS)

master[, tier := fifelse(!is.na(hypothesis), "confirmatory", "exploratory")]

# --- tiered FDR -----------------------------------------------------------------
master[, padj_confirmatory := NA_real_]
master[tier == "confirmatory", padj_confirmatory := p.adjust(pvalue, method = "BH")]
master[, padj_exploratory := NA_real_]
master[tier == "exploratory", padj_exploratory := padj_own]   # keep within-analysis BH

# --- join the matching n=4 MDE from power_mde.tsv -------------------------------
mde <- fread(mde_p)
mde_comp <- mde[readout_class == "composition",   .(unit=cell_type, mde, mde_scale)]
mde_de   <- mde[readout_class == "pseudobulk_DE", .(unit=cell_type, mde, mde_scale)]
mde_prog <- mde[readout_class == "program_score",
                .(key=paste(cell_type, detail), mde, mde_scale)]
master[, `:=`(mde = NA_real_, mde_scale = NA_character_)]
ic <- master$readout_class == "composition"
mm <- mde_comp[match(master$unit[ic], mde_comp$unit)]
master[ic, `:=`(mde = mm$mde, mde_scale = mm$mde_scale)]
id <- master$readout_class == "DE"
mm <- mde_de[match(master$unit[id], mde_de$unit)]
master[id, `:=`(mde = mm$mde, mde_scale = mm$mde_scale)]
ip <- master$readout_class == "pathway"
mm <- mde_prog[match(paste(master$unit, master$feature)[ip], mde_prog$key)]
master[ip, `:=`(mde = mm$mde, mde_scale = mm$mde_scale)]
master[, abs_effect_lt_mde := is.finite(mde) & abs(effect) < mde]

# --- detection-vs-level decomposition (Fix 2): tag every DE row regulation/fraction_shift/ambiguous ---
dtl <- fread(det_p)
master[, `:=`(detection_padj = NA_real_, level_padj = NA_real_,
              mean_among_expr_max = NA_real_, call_class = NA_character_)]
ide <- master$readout_class == "DE"
mm  <- dtl[match(paste(master$feature, master$unit, master$contrast)[ide],
                 paste(dtl$gene, dtl$cell_type, dtl$contrast))]
master[ide, `:=`(detection_padj = mm$detection_padj, level_padj = mm$level_padj,
                 mean_among_expr_max = mm$mean_among_expr_max, call_class = mm$call_class)]

# --- magnitude-floored trend call + clears-MDE (devils-advocate fix) -------------
# A direction is only a "trend" when |effect| >= MDE_FLOOR * mde; below that it is
# noise, not signal (sub-MDE direction coherence is noise coherence).
MDE_FLOOR <- 0.5
master[, clears_mde := is.finite(mde) & abs(effect) >= mde]
master[, trend_call := fifelse(!is.finite(mde) | abs(effect) < MDE_FLOOR * mde, "below-floor",
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
master[, dose_confounded := contrast == "MBRT_vs_SBRT"]                       # MBRT mean dose unrecorded
master[, independent := !(readout_class %in% c("pathway","gsea") & feature == "STING")]  # STING ~50% shares genes with IFN sets

setcolorder(master, c(COLS, "padj_confirmatory", "padj_exploratory",
                      "mde", "mde_scale", "abs_effect_lt_mde", "clears_mde", "trend_call",
                      "n_panel", "n_set_total", "panel_cov_frac"))
setorder(master, tier, readout_class, unit, contrast, feature, na.last = TRUE)
fwrite(master, out_p, sep = "\t")

# --- detectability summary (the transparent decision-gate: what is/isn't detectable) ---
det <- master[, .(n = .N, n_clears_mde = sum(clears_mde, na.rm = TRUE),
                  n_trend_up = sum(trend_call == "up"), n_trend_down = sum(trend_call == "down")),
              by = .(readout_class, contrast, hypothesis, tier)]
setorder(det, tier, readout_class, contrast, hypothesis, na.last = TRUE)
fwrite(det, sub("results_master", "detectability_summary", out_p), sep = "\t")

# --- inspect --------------------------------------------------------------------
cat(sprintf("results_master: %d rows | confirmatory %d / exploratory %d\n",
            nrow(master), master[tier=="confirmatory", .N],
            master[tier=="exploratory", .N]))
cat("rows missing effect/pvalue:", master[is.na(effect) | is.na(pvalue), .N], "\n")
cat("confirmatory rows by readout_class x hypothesis:\n")
print(dcast(master[tier=="confirmatory"], readout_class ~ hypothesis,
            value.var = "pvalue", fun.aggregate = length))
cat(sprintf("\nconfirmatory hits at padj_confirmatory < 0.05: %d\n",
            master[padj_confirmatory < 0.05, .N]))
if (master[padj_confirmatory < 0.05, .N] > 0)
  print(master[padj_confirmatory < 0.05,
               .(readout_class, unit, feature, contrast, effect = round(effect,3),
                 padj_confirmatory = signif(padj_confirmatory,3))][order(padj_confirmatory)])
cat(sprintf("\nconfirmatory rows with |effect| < MDE (underpowered): %d / %d with an MDE\n",
            master[tier=="confirmatory" & abs_effect_lt_mde==TRUE, .N],
            master[tier=="confirmatory" & is.finite(mde), .N]))
