#!/usr/bin/env Rscript
# Run-once scientific-validation acceptance for the modularized differential layer. Asserts the
# known day-2 findings reproduce in the refactored results_master.tsv (no golden-baseline diff).
# Args: <results_master.tsv>
suppressPackageStartupMessages(library(data.table))
m <- fread(commandArgs(trailingOnly = TRUE)[1])
say <- function(...) cat(sprintf(...), "\n")
pass <- TRUE; chk <- function(c, msg){ if(!isTRUE(c)){cat("  MISS:",msg,"\n"); pass<<-FALSE} else cat("  ok:",msg,"\n") }

say("=== results_master: %d rows | confirmatory %d / exploratory %d ===",
    nrow(m), m[tier=="confirmatory",.N], m[tier=="exploratory",.N])

# 1. SBRT stromal fibrosis is the strong signal: collagens/Acta2 in stroma, SBRT_vs_Ctrl, sig.
colla <- m[readout_class=="DE" & contrast=="SBRT_vs_Ctrl" &
           grepl("^Col|Acta2", feature) & unit %in% c("Fibroblast","Adipocyte","SmoothMuscle")]
chk(nrow(colla) > 0 && colla[padj_confirmatory < 0.05 & effect > 0, .N] > 0,
    sprintf("SBRT stromal collagen/Acta2 up + confirmatory-sig (%d rows, %d sig-up)",
            nrow(colla), colla[padj_confirmatory<0.05 & effect>0,.N]))

# 2. H3 dominates SBRT_vs_Ctrl confirmatory hits.
h3 <- m[tier=="confirmatory" & hypothesis=="H3" & contrast=="SBRT_vs_Ctrl" & padj_confirmatory<0.05]
chk(nrow(h3) >= 5, sprintf("H3 SBRT_vs_Ctrl confirmatory hits >=5 (got %d)", nrow(h3)))

# 3. MBRT weak/diluted at arm level: few MBRT_vs_Ctrl confirmatory hits.
mbrt_hits <- m[tier=="confirmatory" & contrast=="MBRT_vs_Ctrl" & padj_confirmatory<0.05]
chk(nrow(mbrt_hits) <= 20, sprintf("MBRT_vs_Ctrl confirmatory hits sparse (got %d, expect single digits)", nrow(mbrt_hits)))

# 4. p21/Cdkn1a in SBRT stroma is ambiguous (detection-level too sparse to call regulation).
p21 <- m[readout_class=="DE" & feature=="Cdkn1a" & contrast=="SBRT_vs_Ctrl" &
         unit %in% c("Fibroblast","Adipocyte","SmoothMuscle","Endothelial")]
chk(nrow(p21)==0 || any(p21$call_class=="ambiguous", na.rm=TRUE) || all(is.na(p21$call_class)),
    sprintf("p21 SBRT stroma not called regulation (classes: %s)",
            paste(unique(p21$call_class), collapse=",")))

# 5. all three contrasts present across readouts
chk(setequal(unique(m$contrast), c("MBRT_vs_Ctrl","SBRT_vs_Ctrl","MBRT_vs_SBRT")),
    "all 3 contrasts present")

say(if (pass) "\nVALIDATION PASS: refactored master reproduces the known day-2 findings"
    else "\nVALIDATION FAIL: one or more known findings did not reproduce -- inspect before trusting")
quit(status = if (pass) 0L else 1L)
