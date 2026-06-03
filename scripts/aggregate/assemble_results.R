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
#       <myeloid_test> <power_mde> <pathway_yaml> <out_master>
suppressPackageStartupMessages({ library(data.table); library(yaml) })
a <- commandArgs(trailingOnly = TRUE)
comp_p <- a[1]; de_p <- a[2]; gsea_p <- a[3]; pw_p <- a[4]; ni_p <- a[5]
mx_p <- a[6]; my_p <- a[7]; mde_p <- a[8]; yaml_p <- a[9]; out_p <- a[10]

IMMUNE   <- c("T cells","NK cells","ILC","Plasma cells","Macrophages","DC",
              "Mast cells","Neutrophils")
STROMA   <- c("Fibroblast","SmoothMuscle","Adipocyte","Endothelial")
STROMA_DE<- c("Fibroblast","SmoothMuscle","Adipocyte")   # endo DE belongs to H2
H1_PROG  <- c("TypeI_interferon","TypeII_interferon","STING")
H2_PROG  <- c("Angiogenesis","Hypoxia")
H3_PROG  <- c("Fibrosis_remodeling","Stromal_stress_senescence")

gs <- lapply(read_yaml(yaml_p), as.character)
H1_GENES <- unique(unlist(gs[H1_PROG]))
H2_GENES <- unique(unlist(gs["Angiogenesis"]))           # H2 DE: Angiogenesis only
H3_GENES <- unique(unlist(gs[H3_PROG]))

dir.create(dirname(out_p), recursive = TRUE, showWarnings = FALSE)
COLS <- c("readout_class","unit","feature","contrast","effect","effect_type",
          "ci_low","ci_high","se","stat","pvalue","padj_own","hypothesis","tier",
          "n_per_arm","dataset")

# --- composition (propeller logit log2FC) ---------------------------------------
comp <- fread(comp_p)
comp_m <- comp[, .(readout_class="composition", unit=cell_type, feature=NA_character_,
  contrast, effect=log2FC_logit, effect_type="log2FC_logit", ci_low=ci_low_log2,
  ci_high=ci_high_log2, se=NA_real_, stat=t_stat, pvalue, padj_own=padj,
  hypothesis=NA_character_, n_per_arm=n_samples_per_group, dataset)]
comp_m[unit %in% IMMUNE,        hypothesis := "H1"]
comp_m[unit == "Endothelial",   hypothesis := "H2;H3"]
comp_m[unit %in% setdiff(STROMA,"Endothelial"), hypothesis := "H3"]

# --- pseudobulk DE (DESeq2 log2FC; Wald CI from lfcSE) ---------------------------
de <- fread(de_p)
de_m <- de[, .(readout_class="DE", unit=cell_type, feature=gene, contrast,
  effect=log2FC, effect_type="log2FC", ci_low=log2FC-1.96*lfcSE,
  ci_high=log2FC+1.96*lfcSE, se=lfcSE, stat, pvalue, padj_own=padj,
  hypothesis=NA_character_, n_per_arm=n_samples_used, dataset)]
de_m[unit %in% IMMUNE    & feature %in% H1_GENES, hypothesis := "H1"]
de_m[unit == "Endothelial" & feature %in% H2_GENES, hypothesis := "H2"]
de_m[unit %in% STROMA_DE & feature %in% H3_GENES, hypothesis := "H3"]

# --- pathway program scores (UCell primary; limma estimate, normal CI) -----------
pw <- fread(pw_p)[score_type == "UCell"]
pw_m <- pw[, .(readout_class="pathway", unit=cell_type, feature=pathway_name, contrast,
  effect=estimate, effect_type="limma_score_estimate", ci_low=estimate-1.96*se,
  ci_high=estimate+1.96*se, se, stat=t_stat, pvalue, padj_own=padj_bh,
  hypothesis=NA_character_, n_per_arm=n_samples_per_group, dataset)]
pw_m[unit %in% IMMUNE              & feature %in% H1_PROG, hypothesis := "H1"]
pw_m[unit %in% c("Endothelial","Tumor") & feature %in% H2_PROG, hypothesis := "H2"]
pw_m[unit %in% STROMA             & feature %in% H3_PROG, hypothesis := "H3"]

# --- GSEA (Hallmark sweep, always exploratory; NES, no CI) -----------------------
gsea <- fread(gsea_p)
gsea_m <- gsea[, .(readout_class="gsea", unit=cell_type, feature=pathway_name, contrast,
  effect=NES, effect_type="NES", ci_low=NA_real_, ci_high=NA_real_, se=NA_real_,
  stat=NES, pvalue, padj_own=padj_bh, hypothesis=NA_character_,
  n_per_arm=NA_integer_, dataset)]

# --- niche composition (propeller; exploratory) ---------------------------------
ni <- fread(ni_p)
ni_m <- ni[, .(readout_class="niche_composition", unit=as.character(niche),
  feature=NA_character_, contrast, effect=log2FC_logit, effect_type="log2FC_logit",
  ci_low=ci_low_log2, ci_high=ci_high_log2, se=NA_real_, stat=t_stat, pvalue,
  padj_own=padj, hypothesis=NA_character_, n_per_arm=n_samples_per_group, dataset)]

# --- spatial mixing + myeloid polarization (limma metric matrices; exploratory) --
mx <- fread(mx_p)
mx_m <- mx[, .(readout_class="spatial_mixing", unit="global", feature=metric, contrast,
  effect=estimate, effect_type="limma_estimate", ci_low, ci_high, se=NA_real_,
  stat=t_stat, pvalue, padj_own=padj, hypothesis=NA_character_,
  n_per_arm=n_samples_per_group, dataset)]
my <- fread(my_p)
my_m <- my[, .(readout_class="myeloid_polarization", unit="Macrophages", feature=metric,
  contrast, effect=estimate, effect_type="limma_estimate", ci_low, ci_high,
  se=NA_real_, stat=t_stat, pvalue, padj_own=padj, hypothesis=NA_character_,
  n_per_arm=n_samples_per_group, dataset)]

master <- rbindlist(list(comp_m, de_m, pw_m, gsea_m, ni_m, mx_m, my_m),
                    use.names = TRUE)
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

setcolorder(master, c(COLS, "padj_confirmatory", "padj_exploratory",
                      "mde", "mde_scale", "abs_effect_lt_mde"))
setorder(master, tier, readout_class, unit, contrast, feature, na.last = TRUE)
fwrite(master, out_p, sep = "\t")

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
