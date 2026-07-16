#!/usr/bin/env Rscript
# aggregate.smk T13: tier-tagged master results table with tiered FDR.
# Collates every M02 day-2 differential readout (composition, pseudobulk DE, GSEA,
# pathway program scores, niche composition, spatial mixing, myeloid polarization)
# into one long table, one row per readout x unit x feature x contrast. Each row is
# tagged tier = confirmatory|exploratory via a frozen a-priori claims-join (the
# results-tier config config/confirmatory_claims.yaml declares each hypothesis's
# pre-registered evidence; producers emit hypothesis=NA and the join fills it, so
# nothing is retro pattern-matched). Primary confirmatory FDR = pooled BH over the
# confirmatory family -> padj_confirmatory. Two auxiliary FDR views ride alongside:
# gatekeeping (gate_padj/padj_gated) and IHW (padj_ihw). Effect size + 95% CI carried
# through; the matching n=4 MDE from power_mde.tsv is joined, and abs_effect_lt_mde
# flags underpowered confirmatory nulls.
#
# The six pre-registered hypotheses and their evidence live in config/hypotheses.yaml
# (biological definition) + config/confirmatory_claims.yaml (readout->confirmatory map).
# Args: <composition_test> <degs> <gsea> <pathway_test> <niche_test> <mixing_test>
#       <myeloid_test> <power_mde> <pathway_sets.tsv> <out_master>
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
comp_p <- a[1]; de_p <- a[2]; gsea_p <- a[3]; pw_p <- a[4]; ni_p <- a[5]
mx_p <- a[6]; my_p <- a[7]; mde_p <- a[8]; sets_p <- a[9]; cov_p <- a[10]; sub_p <- a[11]; ddet_p <- a[12]; out_p <- a[13]

# Panel coverage table for the symmetric per-program coverage columns below.
cov_dt <- fread(cov_p)                                   # set, tier, source, n_total, n_panel, usable, thin

dir.create(dirname(out_p), recursive = TRUE, showWarnings = FALSE)
COLS <- c("readout_class","unit","feature","contrast","effect","effect_type",
          "ci_low","ci_high","se","stat","pvalue","padj_own","hypothesis","tier",
          "n_per_arm","n_samples_used","dataset")

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
  hypothesis=NA_character_, n_per_arm=4L, n_samples_used=NA_integer_,
  dataset="Mutter_02", baseMean=NA_real_)]

# --- substate composition (Fibroblast resting/activated; lm_engine proportion path) --
# Same schema as composition; unit is the sub-state label (e.g. Fibroblast_activated),
# which the myofibroblast_expansion claim targets.
sub <- fread(sub_p)
sub[, padj := p.adjust(p, "BH")]
sub_m <- sub[, .(readout_class="substate_composition", unit=feature_id, feature=NA_character_,
  contrast, effect=estimate, effect_type="log2FC_logit",
  ci_low=estimate - qt(0.975, df)*se, ci_high=estimate + qt(0.975, df)*se,
  se, stat, pvalue=p, padj_own=padj,
  hypothesis=NA_character_, n_per_arm=4L, n_samples_used=NA_integer_,
  dataset="Mutter_02", baseMean=NA_real_)]

# --- pseudobulk DE (count_engine; unit=cell_type, feature_id=gene; Wald CI from se) ----
de <- fread(de_p)
de[, padj := p.adjust(p, "BH"), by = .(unit, contrast)]  # DESeq2 grouping: per cell_type x contrast
de_m <- de[, .(readout_class="DE", unit=unit, feature=feature_id, contrast,
  effect=estimate, effect_type="log2FC", ci_low=estimate-1.96*se,
  ci_high=estimate+1.96*se, se, stat, pvalue=p, padj_own=padj,
  hypothesis=NA_character_, n_per_arm=4L, n_samples_used=n_samples_used,
  dataset="Mutter_02", baseMean=baseMean)]

# --- pathway program scores (UCell primary; limma estimate, normal CI) -----------
pw <- fread(pw_p)[score_type == "UCell"]
pw_m <- pw[, .(readout_class="pathway", unit=cell_type, feature=pathway_name, contrast,
  effect=estimate, effect_type="limma_score_estimate", ci_low=estimate-1.96*se,
  ci_high=estimate+1.96*se, se, stat=t_stat, pvalue, padj_own=padj_bh,
  hypothesis=NA_character_, n_per_arm=4L, n_samples_used=n_samples_per_group,
  dataset, baseMean=NA_real_)]

# --- GSEA (Hallmark sweep, always exploratory; NES, no CI) -----------------------
gsea <- fread(gsea_p)
gsea_m <- gsea[, .(readout_class="gsea", unit=cell_type, feature=pathway_name, contrast,
  effect=NES, effect_type="NES", ci_low=NA_real_, ci_high=NA_real_, se=NA_real_,
  stat=NES, pvalue, padj_own=padj_bh, hypothesis=NA_character_,
  n_per_arm=NA_integer_, n_samples_used=NA_integer_, dataset, baseMean=NA_real_)]

# --- niche composition (lm_engine proportion path; exploratory) -----------------
ni <- fread(ni_p)
ni[, padj := p.adjust(p, "BH")]                          # niches.R used global BH
ni_m <- ni[, .(readout_class="niche_composition", unit=as.character(feature_id),
  feature=NA_character_, contrast, effect=estimate, effect_type="log2FC_logit",
  ci_low=estimate - qt(0.975, df)*se, ci_high=estimate + qt(0.975, df)*se,
  se, stat, pvalue=p, padj_own=padj, hypothesis=NA_character_, n_per_arm=4L,
  n_samples_used=NA_integer_, dataset="Mutter_02", baseMean=NA_real_)]

# --- spatial mixing + myeloid polarization (lm_engine identity path; exploratory) --
mx <- fread(mx_p)
mx[, padj := p.adjust(p, "BH")]
mx_m <- mx[, .(readout_class="spatial_mixing", unit="global", feature=feature_id, contrast,
  effect=estimate, effect_type="limma_estimate",
  ci_low=estimate - qt(0.975, df)*se, ci_high=estimate + qt(0.975, df)*se, se,
  stat, pvalue=p, padj_own=padj, hypothesis=NA_character_,
  n_per_arm=4L, n_samples_used=NA_integer_, dataset="Mutter_02", baseMean=NA_real_)]
my <- fread(my_p)
my[, padj := p.adjust(p, "BH")]
my_m <- my[, .(readout_class="myeloid_polarization", unit="Macrophages", feature=feature_id,
  contrast, effect=estimate, effect_type="limma_estimate",
  ci_low=estimate - qt(0.975, df)*se, ci_high=estimate + qt(0.975, df)*se, se,
  stat, pvalue=p, padj_own=padj, hypothesis=NA_character_,
  n_per_arm=4L, n_samples_used=NA_integer_, dataset="Mutter_02", baseMean=NA_real_)]

# --- differential detection (muscat: fraction of cells expressing gene changed) ------
# Detection-only descriptor (never composition-vs-regulation); always exploratory.
dd <- fread(ddet_p)
dd_m <- dd[, .(readout_class="detection", unit=cell_type, feature=gene, contrast,
  effect=dd_log2fc, effect_type="detection_log2fc", ci_low=NA_real_, ci_high=NA_real_,
  se=NA_real_, stat=NA_real_, pvalue=dd_p, padj_own=dd_padj,
  hypothesis=NA_character_, n_per_arm=4L, n_samples_used=NA_integer_,
  dataset="Mutter_02", baseMean=NA_real_)]

master <- rbindlist(list(comp_m, sub_m, de_m, pw_m, gsea_m, ni_m, mx_m, my_m, dd_m),
                    use.names = TRUE)

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

# --- tiered FDR -----------------------------------------------------------------
# Primary confirmatory FDR = pooled BH over the confirmatory family. Exploratory
# rows keep their within-analysis BH in padj_own (no separate padj_exploratory: it
# was a verbatim copy of padj_own).
master[, padj_confirmatory := NA_real_]
master[tier == "confirmatory", padj_confirmatory := p.adjust(pvalue, method = "BH")]

# --- gatekeeping (secondary FDR view; padj_confirmatory above stays primary) -----
master <- add_gatekeeping(master, CLAIMS, HYP)

# --- IHW auxiliary FDR (keyed on DE baseMean; padj_confirmatory stays primary) ----
master <- add_ihw(master)

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

setcolorder(master, c(COLS, "padj_confirmatory",
                      "gate_program", "gate_padj", "padj_gated", "padj_ihw",
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
