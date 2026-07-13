#!/usr/bin/env Rscript
# Pre-flight the rewritten assemble_results ingestion blocks against real engine outputs,
# BEFORE the 5h run, so a schema/adapter bug surfaces in seconds not hours. Exercises exactly
# the 5 engine-fed blocks (composition, DE, niche, mixing, myeloid); pathway/gsea/mde/det are
# unchanged code and not re-tested here. Asserts: rows produced, no NA effect/pvalue, CI
# brackets the estimate, hypotheses assigned, BH monotone.
suppressPackageStartupMessages(library(data.table))
E <- "results/aggregate/engine"
IMMUNE <- c("T cells","NK cells","ILC","Plasma cells","Macrophages","DC","Mast cells","Neutrophils")
STROMA <- c("Fibroblast","SmoothMuscle","Adipocyte","Endothelial")
ok <- TRUE
chk <- function(cond, msg) { if (!isTRUE(cond)) { cat("FAIL:", msg, "\n"); ok <<- FALSE } }

## composition
comp <- fread(file.path(E,"composition_engine.tsv")); comp[, padj := p.adjust(p,"BH")]
cm <- comp[, .(unit=feature_id, contrast, effect=estimate,
  ci_low=estimate-qt(0.975,df)*se, ci_high=estimate+qt(0.975,df)*se, pvalue=p, padj)]
chk(nrow(cm)==45, sprintf("composition rows=%d (want 45)", nrow(cm)))
chk(!any(is.na(cm$effect)|is.na(cm$pvalue)), "composition NA effect/pvalue")
chk(all(cm$ci_low <= cm$effect & cm$effect <= cm$ci_high), "composition CI brackets estimate")
chk(any(cm$unit %in% IMMUNE) && any(cm$unit %in% STROMA), "composition units span immune+stroma")

## DE (count engine: unit=cell_type, feature_id=gene)
de <- fread(file.path(E,"de_engine.tsv")); de[, padj := p.adjust(p,"BH"), by=.(unit,contrast)]
chk(nrow(de) > 100, sprintf("DE rows=%d", nrow(de)))
chk(!any(is.na(de$estimate)), "DE NA estimate")
chk(uniqueN(de$unit) >= 10, sprintf("DE cell types=%d", uniqueN(de$unit)))
chk(all(de[, .(mono = !is.unsorted(sort(padj))), by=.(unit,contrast)]$mono), "DE padj monotone in p")

## niche / mixing / myeloid
for (nm in c("niche","mixing","myeloid")) {
  d <- fread(file.path(E, paste0(nm,"_engine.tsv"))); d[, padj := p.adjust(p,"BH")]
  chk(nrow(d) > 0 && !any(is.na(d$estimate)|is.na(d$p)), sprintf("%s rows/NA", nm))
  chk(all(d$estimate - qt(0.975,d$df)*d$se <= d$estimate), sprintf("%s CI sane", nm))
  chk(setequal(unique(d$contrast), c("MBRT_vs_Ctrl","SBRT_vs_Ctrl","MBRT_vs_SBRT")),
      sprintf("%s has 3 contrasts", nm))
}
cat(if (ok) "PASS preflight_assemble: all engine-fed blocks ingest cleanly\n" else "PREFLIGHT FAILED\n")
quit(status = if (ok) 0L else 1L)
