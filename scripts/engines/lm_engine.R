#!/usr/bin/env Rscript
# Linear-model test engine (modularization Task 2). Runs ONE limma moderated-t analysis
# for a family of arm contrasts, reading contrast levels / formula / unit / robust flag
# from the resolved comparison registry -- no embedded literals. Two input paths, one
# transform switch:
#   proportion : per-cell labels -> speckle::getTransformedProps(logit) -> celltype x sample
#                logit-proportion matrix (this IS propeller). Estimate reported on the log2
#                scale (natural-log logit / ln2), matching composition.R.
#   matrix     : a sample x feature matrix (continuous readout: pathway score, mixing metric,
#                myeloid ratio) fit directly (identity). Estimate = raw coefficient.
# A single ~0+condition+slide_id fit over all arms; each pairwise contrast is extracted from
# it (single-model contrast, not a per-pair 2-arm refit). Emits SUFFICIENT STATISTICS ONLY
# (estimate/se/df/stat/p) -- no BH here; correction is the results tier's job.
# Args: <input> <input_kind: proportion|matrix> <cohort> <readout> <comparisons.tsv>
#       <engine_params.yaml> <out.tsv>
suppressPackageStartupMessages({
  library(data.table); library(limma); library(speckle); library(yaml)
})
a <- commandArgs(trailingOnly = TRUE)
input_path <- a[1]; input_kind <- a[2]; coh <- a[3]; readout <- a[4]
comp_path  <- a[5]; params_path <- a[6]; out_path <- a[7]
ln2 <- log(2)

## ---- resolve the contrast family from the registry ----
comp <- fread(comp_path)
fam  <- comp[cohort == coh & kind == "sample" & !is.na(contrast_num_level)]
fam  <- unique(fam, by = c("name", "contrast_num_level", "contrast_den_level"))
stopifnot(nrow(fam) >= 1)
formula_str <- fam$formula[1]
stopifnot(!is.na(formula_str))
unit_val <- fam$unit[1]

## ---- readout-specific eBayes robust flag (reproduces each source script) ----
params <- read_yaml(params_path)
robust <- isTRUE(params$lm_engine$robust[[readout]])

## ---- build the (feature x sample) matrix + the sample design table ----
if (input_kind == "proportion") {
  d <- fread(input_path)                       # per-cell: sample_id, label, condition, slide_id
  props <- getTransformedProps(clusters = d$label, sample = d$sample_id, transform = "logit")
  mat   <- props$TransformedProps              # celltype x sample, natural-log logit
  nonfinite <- apply(mat, 1, function(r) any(!is.finite(r)))
  mat   <- mat[!nonfinite, , drop = FALSE]
  samp  <- unique(d[, .(sample_id, condition, slide_id)])
  est_scale <- ln2                             # report on log2-logit scale (composition.R)
  ftype <- "cell_type_proportion"
} else if (input_kind == "matrix") {
  d <- fread(input_path)                        # long: sample_id, feature_id, value, condition, slide_id
  w <- dcast(d, feature_id ~ sample_id, value.var = "value")
  mat <- as.matrix(w[, -1]); rownames(mat) <- w$feature_id
  mat <- mat[complete.cases(mat), , drop = FALSE]
  samp <- unique(d[, .(sample_id, condition, slide_id)])
  est_scale <- 1                                # raw coefficient (identity readouts)
  ftype <- readout
} else stop("unknown input_kind: ", input_kind)

setkey(samp, sample_id); samp <- samp[colnames(mat)]
stopifnot(identical(samp$sample_id, colnames(mat)))
samp[, condition := factor(condition)]
samp[, slide_id  := factor(slide_id)]

## ---- single fit, all contrasts extracted from it ----
design <- model.matrix(as.formula(formula_str), data = samp)
colnames(design) <- make.names(colnames(design))
fit <- lmFit(mat, design)

con_exprs <- setNames(
  sprintf("condition%s - condition%s", fam$contrast_num_level, fam$contrast_den_level),
  fam$name)
cm <- makeContrasts(contrasts = con_exprs, levels = design)
colnames(cm) <- names(con_exprs)
fit2 <- eBayes(contrasts.fit(fit, cm), robust = robust)
df_total <- fit2$df.residual + fit2$df.prior

## ---- extract sufficient statistics ----
res <- rbindlist(lapply(colnames(cm), function(cn) {
  tt <- topTable(fit2, coef = cn, number = Inf, sort.by = "none")
  se_nat <- fit2$stdev.unscaled[, cn] * sqrt(fit2$s2.post)   # moderated SE, natural-log scale
  data.table(
    comparison   = coh,
    contrast     = cn,
    unit         = unit_val,
    feature_type = ftype,
    feature_id   = rownames(tt),
    estimate     = tt$logFC / est_scale,
    se           = se_nat[rownames(tt)] / est_scale,
    df           = df_total[rownames(tt)],
    stat         = tt$t,
    p            = tt$P.Value)
}))

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
fwrite(res, out_path, sep = "\t")
cat(sprintf("lm_engine[%s/%s]: %d features x %d contrasts | robust=%s | %d rows\n",
            coh, readout, nrow(mat), ncol(cm), robust, nrow(res)))
