if (!exists("PATHWAY_COLS")) source("00_load_data.R")
library(ape)
library(zoo)
library(patchwork)

# --- Load scored object from script 14 ---
# Load scored object from script 14 (local, avoids GCS FUSE)
scored_rds <- file.path(DATA_DIR, "seurat_scored.rds")
if (!file.exists(scored_rds)) stop("Run 14_signature_scoring.R first — need seurat_scored.rds")
cat(sprintf("Loading scored object from: %s\n", scored_rds))
obj <- readRDS(scored_rds)
cat(sprintf("Loaded scored object: %d genes x %d cells\n", nrow(obj), ncol(obj)))

stripe_model <- readRDS(file.path(DATA_DIR, "stripe_model.rds"))
meta <- obj@meta.data %>% as.data.frame()

# Ensure parsed columns
meta$model <- ifelse(grepl("^Tongue", meta$Condition), "tongue", "flank")
meta$treatment <- case_when(
  meta$Condition == "Control" ~ "NT",
  grepl("MBRT", meta$Condition) ~ "MBRT",
  grepl("SBRT", meta$Condition) ~ "SBRT"
)
meta$timepoint_h <- case_when(
  meta$Condition == "Control" ~ 0,
  grepl("1h", meta$Condition) ~ 1,
  grepl("4h", meta$Condition) ~ 4,
  grepl("day2", meta$Condition) ~ 48,
  grepl("day6", meta$Condition) ~ 144,
  grepl("day8", meta$Condition) ~ 192,
  grepl("day10", meta$Condition) ~ 240
)

# ============================================================
# 1. Kinetics: signature score vs timepoint
# ============================================================
cat("Computing signature kinetics...\n")

kinetics <- meta %>%
  group_by(Condition, timepoint_h, treatment, model) %>%
  summarise(
    mean_peak   = mean(peak_up, na.rm = TRUE),
    sd_peak     = sd(peak_up, na.rm = TRUE),
    mean_valley = mean(valley_up, na.rm = TRUE),
    sd_valley   = sd(valley_up, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

kinetics_long <- kinetics %>%
  pivot_longer(cols = c(mean_peak, mean_valley),
               names_to = "signature", values_to = "mean_score",
               names_prefix = "mean_") %>%
  mutate(sd = ifelse(signature == "peak", sd_peak, sd_valley),
         signature = factor(signature, levels = c("peak", "valley"),
                            labels = c("Peak-up", "Valley-up")))

write_tsv(kinetics, file.path(DATA_DIR, "signature_kinetics.tsv"))

# Flank kinetics plot
p_kin_flank <- kinetics_long %>%
  filter(model == "flank", treatment %in% c("MBRT", "SBRT", "NT")) %>%
  ggplot(aes(x = timepoint_h, y = mean_score, color = treatment,
             linetype = signature, group = interaction(treatment, signature))) +
  geom_line(linewidth = 0.8) + geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_score - sd, ymax = mean_score + sd),
                width = 5, alpha = 0.4) +
  scale_color_manual(values = TREATMENT_COLORS) +
  scale_x_continuous(breaks = c(0, 1, 4, 48, 144),
                     labels = c("0", "1h", "4h", "2d", "6d")) +
  labs(title = "Peak/valley signature kinetics (flank)",
       x = "Time post-RT", y = "UCell score (mean ± SD)",
       color = "Treatment", linetype = "Signature",
       caption = "n=1 per condition; signatures exclude DDR genes") +
  theme_bw(base_size = 12)

# Tongue kinetics plot
p_kin_tongue <- kinetics_long %>%
  filter(model == "tongue" | (treatment == "NT")) %>%
  ggplot(aes(x = timepoint_h, y = mean_score, color = treatment,
             linetype = signature, group = interaction(treatment, signature))) +
  geom_line(linewidth = 0.8) + geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_score - sd, ymax = mean_score + sd),
                width = 5, alpha = 0.4) +
  scale_color_manual(values = TREATMENT_COLORS) +
  scale_x_continuous(breaks = c(0, 192, 240),
                     labels = c("0", "8d", "10d")) +
  labs(title = "Peak/valley signature kinetics (tongue)",
       x = "Time post-RT", y = "UCell score (mean ± SD)",
       color = "Treatment", linetype = "Signature",
       caption = "n=1 per condition") +
  theme_bw(base_size = 12)

p_kin_combined <- p_kin_flank / p_kin_tongue + plot_layout(heights = c(2, 1))
ggsave(file.path(PLOT_DIR, "signature_kinetics_bulk.png"),
       plot = p_kin_combined, width = 10, height = 10, dpi = 150)
cat("Saved: signature_kinetics_bulk.png\n")

# ============================================================
# 2. SBRT comparison
# ============================================================
cat("SBRT comparison...\n")

sbrt_compare <- meta %>%
  filter(Condition %in% c("MBRT_4h", "SBRT_4h", "MBRT_day2", "SBRT_day2",
                           "MBRT_day6", "SBRT_day6", "Control")) %>%
  select(Condition, treatment, timepoint_h, peak_up, valley_up) %>%
  pivot_longer(cols = c(peak_up, valley_up),
               names_to = "signature", values_to = "score") %>%
  mutate(signature = gsub("_UCell", "", signature),
         signature = factor(signature, levels = c("peak_up", "valley_up"),
                            labels = c("Peak-up", "Valley-up")),
         Condition = factor(Condition, levels = CONDITION_LEVELS))

p_sbrt <- ggplot(sbrt_compare, aes(x = Condition, y = score, fill = treatment)) +
  geom_boxplot(outlier.size = 0.1, outlier.alpha = 0.1) +
  scale_fill_manual(values = TREATMENT_COLORS) +
  facet_wrap(~signature, scales = "free_y") +
  labs(title = "Peak/valley signatures: MBRT vs SBRT vs Control",
       y = "UCell score", x = NULL,
       caption = "n=1 per condition; signatures exclude DDR genes") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(PLOT_DIR, "signature_sbrt_comparison.png"),
       plot = p_sbrt, width = 14, height = 6, dpi = 150)
cat("Saved: signature_sbrt_comparison.png\n")

# ============================================================
# 3. Spatial maps of signature scores across MBRT timepoints
# ============================================================
cat("Generating spatial maps...\n")

mbrt_conditions <- c("MBRT_1h", "MBRT_4h", "MBRT_day2", "MBRT_day6")
mbrt_meta <- meta %>% filter(Condition %in% mbrt_conditions)

# Shared color scale for comparability across timepoints
peak_range  <- range(mbrt_meta$peak_up, na.rm = TRUE)
valley_range <- range(mbrt_meta$valley_up, na.rm = TRUE)

p_spatial_peak <- mbrt_meta %>%
  mutate(Condition = factor(Condition, levels = mbrt_conditions)) %>%
  ggplot(aes(x = x_slide_mm, y = y_slide_mm, color = peak_up)) +
  geom_point(size = 0.02, alpha = 0.3) +
  scale_color_viridis_c(option = "inferno", limits = peak_range) +
  facet_wrap(~Condition, scales = "free") +
  labs(title = "Peak-up signature score (MBRT flank timepoints)",
       color = "UCell") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom", aspect.ratio = 1)

p_spatial_valley <- mbrt_meta %>%
  mutate(Condition = factor(Condition, levels = mbrt_conditions)) %>%
  ggplot(aes(x = x_slide_mm, y = y_slide_mm, color = valley_up)) +
  geom_point(size = 0.02, alpha = 0.3) +
  scale_color_viridis_c(option = "mako", limits = valley_range) +
  facet_wrap(~Condition, scales = "free") +
  labs(title = "Valley-up signature score (MBRT flank timepoints)",
       color = "UCell") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom", aspect.ratio = 1)

p_spatial <- p_spatial_peak / p_spatial_valley
ggsave(file.path(PLOT_DIR, "signature_spatial_maps.png"),
       plot = p_spatial, width = 16, height = 16, dpi = 150)
cat("Saved: signature_spatial_maps.png\n")

# ============================================================
# 4. Stripe persistence: per-condition stripe fitting + Moran's I
# ============================================================
cat("Stripe persistence analysis...\n")
cat("Each MBRT condition is a different animal/slide — fitting stripe model per condition.\n")

# Need p21 expression for stripe fitting
scored_obj <- obj  # already loaded above
p21_all <- GetAssayData(scored_obj, layer = "data")["Cdkn1a", ]

# Build weight matrix once, reuse for observed + permutations
build_weight_matrix <- function(coords_xy, n_neighbors = 20, max_cells = 10000) {
  n <- nrow(coords_xy)
  if (n > max_cells) {
    idx <- sample(n, max_cells)
    coords_xy <- coords_xy[idx, ]
    n <- max_cells
    attr(coords_xy, "idx") <- idx
  }
  library(RANN)
  nn <- nn2(coords_xy, k = n_neighbors + 1)
  w <- matrix(0, n, n)
  for (i in 1:n) {
    nn_idx <- nn$nn.idx[i, -1]
    nn_dist <- nn$nn.dists[i, -1]
    w[i, nn_idx] <- 1 / nn_dist
  }
  list(w = w, idx = attr(coords_xy, "idx"), n = n)
}

compute_morans_i <- function(values, w_obj) {
  Moran.I(values, w_obj$w)
}

# Per-condition stripe model fitting (replicating script 05 logic)
fit_stripe_model <- function(x_mm, y_mm, p21_vals, fixed_spacing = 1.02) {
  # Scan tilt angles, fix spacing at collimator value
  best_tilt <- NA; best_contrast <- 0; best_offset <- NA
  for (deg in seq(0, 35, by = 1)) {
    rad <- deg * pi / 180
    y_corr <- y_mm + x_mm * tan(rad)
    for (o in seq(min(y_corr) + fixed_spacing/4, min(y_corr) + fixed_spacing * 2, by = 0.05)) {
      centers <- seq(o, max(y_corr) + fixed_spacing, by = fixed_spacing)
      dist <- sapply(y_corr, function(yc) min(abs(yc - centers)))
      pk <- which(dist < fixed_spacing / 4)
      vl <- which(dist > fixed_spacing / 2)
      if (length(pk) < 200 || length(vl) < 200) next
      contrast <- mean(p21_vals[pk], na.rm = TRUE) - mean(p21_vals[vl], na.rm = TRUE)
      if (contrast > best_contrast) {
        best_contrast <- contrast; best_tilt <- deg; best_offset <- o
      }
    }
  }
  if (is.na(best_tilt)) return(NULL)
  rad <- best_tilt * pi / 180
  y_corr <- y_mm + x_mm * tan(rad)
  centers <- seq(best_offset, max(y_corr) + fixed_spacing, by = fixed_spacing)
  # Keep only centers with enough cells
  centers <- centers[sapply(centers, function(c) sum(abs(y_corr - c) < fixed_spacing/4) > 50)]
  list(tilt_deg = best_tilt, spacing_mm = fixed_spacing, beam_centers = centers,
       n_peaks = length(centers), contrast = best_contrast, offset = best_offset)
}

stripe_results <- list()

for (cond in c("MBRT_1h", "MBRT_4h", "MBRT_day2", "MBRT_day6")) {
  cat(sprintf("\n  === %s ===\n", cond))
  sub <- meta %>% filter(Condition == cond)
  sub$p21 <- as.numeric(p21_all[rownames(sub)])

  # Fit stripe model for this condition's tissue section
  cat(sprintf("  Fitting stripe model (%d cells)...\n", nrow(sub)))
  sm <- fit_stripe_model(sub$x_slide_mm, sub$y_slide_mm, sub$p21)

  if (!is.null(sm)) {
    cat(sprintf("  Fit: %d-deg tilt, %.2fmm spacing, %d peaks, contrast=%.4f\n",
                sm$tilt_deg, sm$spacing_mm, sm$n_peaks, sm$contrast))
    rad <- sm$tilt_deg * pi / 180
    sub$y_corr <- sub$y_slide_mm + sub$x_slide_mm * tan(rad)
  } else {
    cat("  WARNING: no stripe model fit for this condition\n")
    sub$y_corr <- sub$y_slide_mm
  }

  # Moran's I on peak-up signature (orientation-agnostic)
  # Build weight matrix ONCE, subsample to 10K for tractability
  cat(sprintf("  Building weight matrix (10K subsample)...\n"))
  coords_xy <- as.matrix(sub[, c("x_slide_mm", "y_slide_mm")])
  w_obj <- tryCatch(
    build_weight_matrix(coords_xy, n_neighbors = 20, max_cells = 10000),
    error = function(e) { cat("  Weight matrix failed:", e$message, "\n"); NULL }
  )

  if (!is.null(w_obj)) {
    vals <- if (!is.null(w_obj$idx)) sub$peak_up[w_obj$idx] else sub$peak_up
    cat(sprintf("  Computing Moran's I on %d cells...\n", w_obj$n))
    mi <- tryCatch(
      compute_morans_i(vals, w_obj),
      error = function(e) list(observed = NA, expected = NA, sd = NA, p.value = NA)
    )
    # Null distribution: 100 permutations (reuse weight matrix)
    cat("  Running 100 permutations...\n")
    null_is <- replicate(100, {
      perm_vals <- sample(vals)
      tryCatch(compute_morans_i(perm_vals, w_obj)$observed, error = function(e) NA)
    })
  } else {
    mi <- list(observed = NA, expected = NA, sd = NA, p.value = NA)
    null_is <- rep(NA, 100)
  }

  stripe_results[[cond]] <- tibble(
    condition = cond,
    morans_i = mi$observed,
    morans_p = mi$p.value,
    null_mean = mean(null_is, na.rm = TRUE),
    null_sd = sd(null_is, na.rm = TRUE),
    fitted_tilt_deg = if (!is.null(sm)) sm$tilt_deg else NA,
    fitted_n_peaks = if (!is.null(sm)) sm$n_peaks else NA,
    p21_contrast = if (!is.null(sm)) sm$contrast else NA,
    n_cells = nrow(sub),
    n_subsampled = min(nrow(sub), 50000)
  )
}

stripe_df <- bind_rows(stripe_results)
write_tsv(stripe_df, file.path(DATA_DIR, "stripe_persistence.tsv"))
cat("\nStripe persistence results:\n")
print(as.data.frame(stripe_df))

# Moran's I plot
p_morans <- ggplot(stripe_df, aes(x = condition, y = morans_i)) +
  geom_point(size = 4, color = "#E74C3C") +
  geom_errorbar(aes(ymin = null_mean - 2*null_sd, ymax = null_mean + 2*null_sd),
                width = 0.2, color = "grey50", linetype = "dashed") +
  geom_point(aes(y = null_mean), shape = 1, size = 3, color = "grey50") +
  labs(title = "Spatial autocorrelation of peak-up signature across timepoints",
       subtitle = "Red = observed Moran's I; grey = null (100 permutations, ±2 SD)\nPer-condition stripe model fit (different animal/slide per timepoint)",
       x = NULL, y = "Moran's I",
       caption = "Subsampled to 50K cells; k=20 NN weights; 1.02mm spacing fixed") +
  theme_bw(base_size = 12)
ggsave(file.path(PLOT_DIR, "stripe_persistence_morans.png"),
       plot = p_morans, width = 8, height = 6, dpi = 150)
cat("Saved: stripe_persistence_morans.png\n")

# ============================================================
# 5. Cell-type stratified kinetics (if applicable)
# ============================================================
ct_sigs_file <- file.path(DATA_DIR, "celltype_signatures.tsv")
if (file.exists(ct_sigs_file)) {
  ct_sigs <- read_tsv(ct_sigs_file, show_col_types = FALSE)
  ct_types <- unique(ct_sigs$cell_type)

  if (length(ct_types) > 0) {
    cat(sprintf("Cell-type stratified kinetics for: %s\n", paste(ct_types, collapse = ", ")))

    ct_kinetics <- meta %>%
      filter(cell_type_validated %in% ct_types, model == "flank") %>%
      group_by(Condition, timepoint_h, treatment, cell_type_validated) %>%
      summarise(
        mean_peak = mean(peak_up, na.rm = TRUE),
        mean_valley = mean(valley_up, na.rm = TRUE),
        n = n(),
        .groups = "drop"
      ) %>%
      pivot_longer(cols = c(mean_peak, mean_valley),
                   names_to = "signature", values_to = "mean_score",
                   names_prefix = "mean_") %>%
      mutate(signature = factor(signature, levels = c("peak", "valley"),
                                labels = c("Peak-up", "Valley-up")))

    p_ct_kin <- ct_kinetics %>%
      filter(treatment %in% c("MBRT", "SBRT")) %>%
      ggplot(aes(x = timepoint_h, y = mean_score, color = treatment,
                 linetype = signature, group = interaction(treatment, signature))) +
      geom_line(linewidth = 0.7) + geom_point(size = 2) +
      scale_color_manual(values = TREATMENT_COLORS) +
      scale_x_continuous(breaks = c(0, 1, 4, 48, 144),
                         labels = c("0", "1h", "4h", "2d", "6d")) +
      facet_wrap(~cell_type_validated, scales = "free_y") +
      labs(title = "Signature kinetics by cell type (flank)",
           x = "Time post-RT", y = "UCell score",
           caption = "n=1 per condition; bulk UCell scores, not cell-type-specific signatures") +
      theme_bw(base_size = 10)
    ggsave(file.path(PLOT_DIR, "signature_kinetics_celltype.png"),
           plot = p_ct_kin, width = 14, height = 10, dpi = 150)
    cat("Saved: signature_kinetics_celltype.png\n")
  }
}

cat("Signature kinetics analysis complete.\n")
