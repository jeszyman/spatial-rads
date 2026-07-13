#!/usr/bin/env Rscript
# Task 4 TDD: the lm_engine must reproduce each producer's CURRENT inline limma test
# (estimate + p) before that test is removed. References captured in results/aggregate/_ref/.
# matrix path (mixing, myeloid): estimate = raw logFC. proportion path (niche, substate):
# estimate = log2 logit (reference column log2FC_logit). Live source-compare.
suppressPackageStartupMessages(library(data.table))
S <- "results/aggregate/_lm_smoke"; R <- "results/aggregate/_ref"
cases <- list(
  list(nm="mixing",   ref=file.path(R,"spatial_mixing_test_m02day2.tsv"),
       eng=file.path(S,"mixing_lm.tsv"),   key="metric",    est="estimate"),
  list(nm="myeloid",  ref=file.path(R,"myeloid_m1m2_test_m02day2.tsv"),
       eng=file.path(S,"myeloid_lm.tsv"),  key="metric",    est="estimate"),
  list(nm="niche",    ref=file.path(R,"niche_test_m02day2.tsv"),
       eng=file.path(S,"niche_lm.tsv"),    key="niche",     est="log2FC_logit"),
  list(nm="substate", ref=file.path(R,"composition_substate_test_m02day2.tsv"),
       eng=file.path(S,"substate_lm.tsv"), key="cell_type", est="log2FC_logit"))
fail <- 0L
for (c in cases) {
  ref <- fread(c$ref); eng <- fread(c$eng)
  ref2 <- ref[, .(feature_id = get(c$key), contrast, ref_est = get(c$est), ref_p = pvalue)]
  m <- merge(ref2, eng[, .(feature_id, contrast, estimate, p)], by = c("feature_id","contrast"))
  bad_e <- m[abs(ref_est - estimate)/pmax(abs(ref_est),1e-9) >= 1e-6, .N]
  bad_p <- m[abs(ref_p - p) >= 1e-9, .N]
  status <- if (nrow(m) > 0 && bad_e == 0 && bad_p == 0) "PASS" else {fail <- fail+1L; "FAIL"}
  cat(sprintf("%s %-9s: %d rows | %d est / %d p mismatches\n", status, c$nm, nrow(m), bad_e, bad_p))
}
quit(status = if (fail > 0) 1L else 0L)
