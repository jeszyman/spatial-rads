## 04_celldeath_senescence.R
## Three strands:
##   4a. Apoptosis + senescence signature kinetics (UCell)
##   4b. DAMP / alarmin signature kinetics + boundary analysis
##   4c. Necrosis quantification (Tier 1 priority — central biological output)

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(UCell)
  library(RANN)
})

DATA_DIR <- "data"
PLOT_DIR <- "plots"
OBJECTS_DIR <- "/mnt/data/projects/spatial-rads/analysis/objects/mbrt_vs_sbrt"

obj <- readRDS(file.path(OBJECTS_DIR, "seurat_filtered.rds"))
cat(sprintf("Loaded: %d cells\n", ncol(obj)))

# ============================================================
# Score apoptosis, senescence, DAMP signatures (UCell)
# ============================================================
hallmark_apoptosis_genes <- c(
  "Add1","Aifm3","Ank","Anxa1","App","Atf3","Avpr1a","Bax","Bcap31","Bcl10",
  "Bcl2l1","Bcl2l10","Bcl2l11","Bcl2l2","Bgn","Bid","Bik","Bmf","Bmp2","Bnip3l",
  "Brca1","Btg2","Btg3","Casp1","Casp2","Casp3","Casp4","Casp6","Casp7","Casp8",
  "Casp9","Cd14","Cd2","Cd38","Cd44","Cd69","Cdc25b","Cdk2","Cdkn1a","Cdkn1b",
  "Cflar","Clu","Crebbp","Ctnnb1","Ctsb","Ctsc","Ctsd","Ctso","Ctsv","Ctsz",
  "Cyld","Daxx","Ddit3","Dffa","Diablo","Dnaja1","Dnm1l","Dpyd","Eaf1","Ebag9",
  "Egr3","Emp1","Enpp2","Erbb2","Erbb3","Etf1","Faim2","Fas","Fasl","Fdxr",
  "Fez1","Gadd45a","Gadd45b","Gch1","Gna15","Gpx1","Gpx3","Gpx4","Gsn","Gsr",
  "Gstm1","Gucy2d","H1-0","Hgf","Hmgb2","Hmox1","Hspb1","Ier3","Igf2r","Igfbp6",
  "Il18","Il1a","Il1b","Il6","Irf1","Isg20","Jun","Krt18","Lef1","Lgals3",
  "Lmna","Lum","Madd","Mcl1","Mgmt","Mmp2","Nedd9","Nefh","Pak1","Pdcd4",
  "Pdgfra","Pea15","Plat","Plcb2","Plppr4","Pmaip1","Ppp2r5b","Ppp3r1","Ppt1",
  "Prf1","Psen1","Psen2","Ptk2","Rara","Rela","Rffl","Rhob","Rhot2","Rnasel",
  "Sat1","Satb1","Sc5d","Slc20a1","Smad7","Sod1","Sod2","Sptan1","Tap1",
  "Timp1","Timp2","Timp3","Tnf","Tnfrsf12a","Tnfsf10","Top2a","Tspo","Txnip",
  "Vdac2","Wee1","Xiap")

senmayo_genes <- c(
  "Acvr1b","Ang","Angpt1","Angptl4","Areg","Axl","Bex3","Bmp2","Bmp6","C3",
  "Ccl1","Ccl13","Ccl16","Ccl2","Ccl20","Ccl24","Ccl26","Ccl3","Ccl3l1","Ccl4",
  "Ccl5","Ccl7","Ccl8","Cd55","Cd9","Csf1","Csf2","Csf2rb","Cst10","Ctnnb1",
  "Ctsb","Cxcl1","Cxcl10","Cxcl12","Cxcl16","Cxcl2","Cxcl3","Cxcl8","Cxcr2",
  "Dkk1","Edn1","Egf","Egfr","Ereg","Esm1","Ets2","Fas","Fgf1","Fgf2","Fgf7",
  "Gdf15","Gem","Gmfg","Hgf","Hmgb1","Icam1","Icam3","Igf1","Igfbp1","Igfbp2",
  "Igfbp3","Igfbp4","Igfbp5","Igfbp6","Igfbp7","Il10","Il13","Il15","Il18",
  "Il1a","Il1b","Il2","Il6","Il6st","Il7","Inha","Iqgap2","Itga2","Itpka",
  "Jun","Kitlg","Lcp1","Mif","Mmp1","Mmp10","Mmp12","Mmp13","Mmp14","Mmp2",
  "Mmp3","Mmp9","Nap1l4","Nrg1","Pappa","Pecam1","Pgf","Pigf","Plat","Plau",
  "Plaur","Ppbp","Ptbp1","Ptger2","Ptges","Rps6ka5","Scamp4","Selplg","Sema3f",
  "Serpinb2","Serpine1","Serpine2","Spp1","Spx","Timp2","Tnf","Tnfrsf10c",
  "Tnfrsf11b","Tnfrsf1a","Tnfrsf1b","Tubgcp2","Vegfa","Vegfc","Vgf","Wnt16","Wnt2"
)

damp_genes <- c("Hmgb1","Hmgb2","S100a8","S100a9","Anxa1","Anxa2",
                "Il1b","Il6","Nlrp3","Casp1","Tlr2","Tlr4")

panel_genes <- rownames(obj)
sig_lists <- list(
  APOPTOSIS    = intersect(hallmark_apoptosis_genes, panel_genes),
  SENESCENCE   = intersect(senmayo_genes, panel_genes),
  DAMP         = intersect(damp_genes, panel_genes)
)
cat("Signature panel coverage:\n")
for (sn in names(sig_lists)) {
  cat(sprintf("  %s: %d genes in panel\n", sn, length(sig_lists[[sn]])))
}

if (!"APOPTOSIS_UCell" %in% colnames(obj@meta.data)) {
  cat("\nScoring with UCell (could take 5-10 min)...\n")
  t0 <- Sys.time()
  obj <- AddModuleScore_UCell(obj, features = sig_lists, ncores = 8, name = "_UCell")
  cat(sprintf("UCell done: %.1f min\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

# ============================================================
# 4a. Apoptosis + senescence kinetics
# ============================================================
cat("\n=== 4a. Apoptosis + senescence kinetics ===\n")

# Default exclusion: necrosis_zone at day2/day6
keep_for_main <- !(obj$timepoint_h %in% c(48, 144) & obj$necrosis_zone %in% TRUE)

celldeath_kinetics <- as_tibble(obj@meta.data) %>%
  mutate(use_for_main = keep_for_main) %>%
  filter(use_for_main) %>%
  group_by(cell_type_major, timepoint_h, treatment, Condition) %>%
  summarise(across(c(APOPTOSIS_UCell, SENESCENCE_UCell, DAMP_UCell),
                   list(mean = ~mean(., na.rm = TRUE),
                        sd   = ~sd(., na.rm = TRUE),
                        q75  = ~quantile(., 0.75, na.rm = TRUE))),
            n_cells = n(),
            .groups = "drop")

# Frac high-scorers (above 75th percentile within each timepoint, all cells pooled)
tp_thresholds <- as_tibble(obj@meta.data) %>%
  group_by(timepoint_h) %>%
  summarise(apo_q75 = quantile(APOPTOSIS_UCell, 0.75, na.rm = TRUE),
            sen_q75 = quantile(SENESCENCE_UCell, 0.75, na.rm = TRUE),
            damp_q75 = quantile(DAMP_UCell, 0.75, na.rm = TRUE),
            .groups = "drop")

high_scorer_frac <- as_tibble(obj@meta.data) %>%
  mutate(use_for_main = keep_for_main) %>%
  filter(use_for_main) %>%
  left_join(tp_thresholds, by = "timepoint_h") %>%
  group_by(cell_type_major, timepoint_h, treatment, Condition) %>%
  summarise(frac_apo_high = mean(APOPTOSIS_UCell > apo_q75, na.rm = TRUE),
            frac_sen_high = mean(SENESCENCE_UCell > sen_q75, na.rm = TRUE),
            frac_damp_high = mean(DAMP_UCell > damp_q75, na.rm = TRUE),
            .groups = "drop")

celldeath_kinetics <- celldeath_kinetics %>%
  left_join(high_scorer_frac,
            by = c("cell_type_major", "timepoint_h", "treatment", "Condition"))
write_tsv(celldeath_kinetics, file.path(DATA_DIR, "celldeath_kinetics.tsv"))

# Validation: necrosis_zone vs not — apoptosis/senescence scores
necr_strat <- as_tibble(obj@meta.data) %>%
  filter(!is.na(necrosis_zone)) %>%
  group_by(cell_type_major, treatment, timepoint_h, necrosis_zone) %>%
  summarise(mean_apo = mean(APOPTOSIS_UCell, na.rm = TRUE),
            mean_sen = mean(SENESCENCE_UCell, na.rm = TRUE),
            mean_damp = mean(DAMP_UCell, na.rm = TRUE),
            n_cells = n(),
            .groups = "drop")
write_tsv(necr_strat, file.path(DATA_DIR, "necrosis_zone_vs_signature_scores.tsv"))

# Plot: kinetics line plots
TREATMENT_COLORS <- c("MBRT" = "#E74C3C", "SBRT" = "#3498DB", "NT" = "#7F8C8D")
flank_conditions <- c("Control",
                      "MBRT_1h", "MBRT_4h", "MBRT_day2", "MBRT_day6",
                      "SBRT_4h", "SBRT_day2", "SBRT_day6")

p_apo <- celldeath_kinetics %>%
  filter(Condition %in% flank_conditions, treatment %in% c("MBRT","SBRT","NT")) %>%
  ggplot(aes(x = factor(timepoint_h), y = APOPTOSIS_UCell_mean,
             color = treatment, group = treatment)) +
  geom_line(linewidth = 0.7) + geom_point(size = 2) +
  facet_wrap(~cell_type_major, scales = "free_y", ncol = 3) +
  scale_color_manual(values = TREATMENT_COLORS) +
  labs(title = "Apoptosis (Hallmark) UCell score kinetics",
       x = "Hours post-RT", y = "Mean apoptosis score") + theme_bw()
ggsave(file.path(PLOT_DIR, "apoptosis_kinetics.png"),
       plot = p_apo, width = 12, height = 9, dpi = 150)

p_sen <- celldeath_kinetics %>%
  filter(Condition %in% flank_conditions, treatment %in% c("MBRT","SBRT","NT")) %>%
  ggplot(aes(x = factor(timepoint_h), y = SENESCENCE_UCell_mean,
             color = treatment, group = treatment)) +
  geom_line(linewidth = 0.7) + geom_point(size = 2) +
  facet_wrap(~cell_type_major, scales = "free_y", ncol = 3) +
  scale_color_manual(values = TREATMENT_COLORS) +
  labs(title = "SenMayo senescence UCell score kinetics",
       x = "Hours post-RT", y = "Mean senescence score") + theme_bw()
ggsave(file.path(PLOT_DIR, "senescence_kinetics.png"),
       plot = p_sen, width = 12, height = 9, dpi = 150)

# ============================================================
# 4b. DAMP signature kinetics + boundary
# ============================================================
cat("\n=== 4b. DAMP signature kinetics ===\n")

p_damp <- celldeath_kinetics %>%
  filter(Condition %in% flank_conditions, treatment %in% c("MBRT","SBRT","NT")) %>%
  ggplot(aes(x = factor(timepoint_h), y = DAMP_UCell_mean,
             color = treatment, group = treatment)) +
  geom_line(linewidth = 0.7) + geom_point(size = 2) +
  facet_wrap(~cell_type_major, scales = "free_y", ncol = 3) +
  scale_color_manual(values = TREATMENT_COLORS) +
  labs(title = "DAMP / alarmin signature UCell score kinetics",
       x = "Hours post-RT", y = "Mean DAMP score") + theme_bw()
ggsave(file.path(PLOT_DIR, "damp_kinetics.png"),
       plot = p_damp, width = 12, height = 9, dpi = 150)

# DAMP boundary: distance to nearest necrosis_zone cell (per slide)
cat("Computing distance-to-necrosis-zone (this may take a few minutes)...\n")
slide_col <- if ("Slide" %in% colnames(obj@meta.data)) "Slide" else "slide_ID_numeric"
slides <- unique(obj@meta.data[[slide_col]])

dist_to_necr <- rep(NA_real_, ncol(obj))
for (s in slides) {
  idx <- which(obj@meta.data[[slide_col]] == s)
  if (length(idx) < 50) next
  m_s <- obj@meta.data[idx, c("x_slide_mm", "y_slide_mm")]
  necr_idx <- which(obj@meta.data$necrosis_zone[idx] %in% TRUE)
  if (length(necr_idx) == 0) next
  necr_coords <- as.matrix(m_s[necr_idx, , drop = FALSE])
  query_coords <- as.matrix(m_s)
  nn <- RANN::nn2(necr_coords, query = query_coords, k = 1)
  dist_to_necr[idx] <- nn$nn.dists[, 1]
}
obj$dist_to_necrosis_mm <- dist_to_necr

damp_boundary <- as_tibble(obj@meta.data) %>%
  filter(!is.na(dist_to_necrosis_mm), timepoint_h %in% c(48, 144)) %>%
  mutate(distance_bin = cut(dist_to_necrosis_mm,
                            breaks = c(-Inf, 0.05, 0.1, 0.2, 0.5, 1, Inf),
                            labels = c("<=0.05", "0.05-0.1", "0.1-0.2", "0.2-0.5", "0.5-1.0", ">1.0"))) %>%
  group_by(treatment, timepoint_h, distance_bin) %>%
  summarise(mean_damp = mean(DAMP_UCell, na.rm = TRUE),
            n_cells = n(),
            .groups = "drop")
write_tsv(damp_boundary, file.path(DATA_DIR, "damp_boundary.tsv"))

p_db <- damp_boundary %>%
  filter(!is.na(distance_bin), treatment %in% c("MBRT","SBRT","NT")) %>%
  ggplot(aes(x = distance_bin, y = mean_damp,
             color = treatment, group = treatment)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  facet_wrap(~factor(timepoint_h, levels = c(48,144),
                     labels = c("day2","day6"))) +
  scale_color_manual(values = TREATMENT_COLORS) +
  labs(title = "DAMP score vs distance to necrosis-zone cell",
       x = "Distance bin (mm)", y = "Mean DAMP score") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(PLOT_DIR, "damp_boundary_kinetics.png"),
       plot = p_db, width = 9, height = 5, dpi = 150)

# ============================================================
# 4c. Necrosis quantification (Tier 1 priority)
# ============================================================
cat("\n=== 4c. Necrosis quantification ===\n")

# Per (slide x condition x timepoint) summary
necr_per_slide_cond <- as_tibble(obj@meta.data) %>%
  filter(!is.na(necrosis_zone)) %>%
  rename(slide_id = !!slide_col) %>%
  group_by(slide_id, Condition, treatment, timepoint_h) %>%
  summarise(n_total = n(),
            n_necrosis = sum(necrosis_zone),
            necrosis_pct = 100 * n_necrosis / n_total,
            .groups = "drop")

# Calibrated necrosis flag using Control distribution
ctrl_d20 <- obj@meta.data$local_density_d20[obj@meta.data$Condition == "Control"]
ctrl_thresh <- quantile(ctrl_d20, 0.95, na.rm = TRUE)
cat(sprintf("Control d20 95th percentile threshold: %.4f mm\n", ctrl_thresh))
obj$necrosis_calibrated <- obj$local_density_d20 > ctrl_thresh

necr_calib <- as_tibble(obj@meta.data) %>%
  filter(!is.na(necrosis_calibrated)) %>%
  group_by(Condition, treatment, timepoint_h) %>%
  summarise(n_total = n(),
            n_necr_per_slide = sum(necrosis_zone),
            necr_pct_per_slide = 100 * n_necr_per_slide / n_total,
            n_necr_calib = sum(necrosis_calibrated),
            necr_pct_calib = 100 * n_necr_calib / n_total,
            .groups = "drop")
write_tsv(necr_calib, file.path(DATA_DIR, "necrosis_quantification.tsv"))

cat("\nNecrosis fractions by condition (calibrated to Control 95th pct):\n")
print(necr_calib %>%
        filter(Condition %in% flank_conditions) %>%
        select(Condition, n_total, necr_pct_per_slide, necr_pct_calib) %>%
        arrange(Condition))

# Spatial extent: bin slide into 0.1mm grid; necrotic if >50% of cells in bin are necr
bin_size_mm <- 0.1
extent_results <- list()
for (s in slides) {
  m_s <- as_tibble(obj@meta.data) %>%
    filter(!!sym(slide_col) == s, !is.na(necrosis_zone)) %>%
    mutate(x_bin = round(x_slide_mm / bin_size_mm) * bin_size_mm,
           y_bin = round(y_slide_mm / bin_size_mm) * bin_size_mm)
  if (nrow(m_s) < 50) next

  bins_per_cond <- m_s %>%
    group_by(Condition, treatment, timepoint_h, x_bin, y_bin) %>%
    summarise(n_cells = n(),
              n_necr = sum(necrosis_zone, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(is_necrotic_bin = (n_necr / pmax(n_cells, 1)) >= 0.5 & n_cells >= 3)

  per_cond <- bins_per_cond %>%
    group_by(Condition, treatment, timepoint_h) %>%
    summarise(total_bins = n(),
              necrotic_bins = sum(is_necrotic_bin),
              total_necrotic_area_mm2 = necrotic_bins * bin_size_mm^2,
              .groups = "drop") %>%
    mutate(slide_id = s)

  extent_results[[as.character(s)]] <- per_cond
}
extent_df <- bind_rows(extent_results)
write_tsv(extent_df, file.path(DATA_DIR, "necrosis_spatial_extent.tsv"))

# Connected-component analysis (largest region per condition x timepoint x slide)
# Simple flood-fill on 2D binary grid
flood_fill_max_region <- function(bins_xy) {
  # bins_xy: tibble with x_bin, y_bin (necrotic bins only)
  if (nrow(bins_xy) == 0) return(0L)
  grid_x <- round(bins_xy$x_bin, 4)
  grid_y <- round(bins_xy$y_bin, 4)
  keys <- paste(grid_x, grid_y, sep = "_")
  visited <- rep(FALSE, length(keys))
  max_size <- 0L
  for (i in seq_along(keys)) {
    if (visited[i]) next
    queue <- i; size <- 0L
    while (length(queue) > 0) {
      cur <- queue[1]; queue <- queue[-1]
      if (visited[cur]) next
      visited[cur] <- TRUE; size <- size + 1L
      x <- grid_x[cur]; y <- grid_y[cur]
      # 8-connected neighbors
      for (dx in c(-bin_size_mm, 0, bin_size_mm)) {
        for (dy in c(-bin_size_mm, 0, bin_size_mm)) {
          if (dx == 0 && dy == 0) next
          nbr_key <- paste(round(x+dx, 4), round(y+dy, 4), sep = "_")
          m <- match(nbr_key, keys)
          if (!is.na(m) && !visited[m]) queue <- c(queue, m)
        }
      }
    }
    if (size > max_size) max_size <- size
  }
  max_size
}

cat("\nComputing connected components (largest necrotic region)...\n")
cc_results <- list()
for (s in slides) {
  m_s <- as_tibble(obj@meta.data) %>%
    filter(!!sym(slide_col) == s, !is.na(necrosis_zone)) %>%
    mutate(x_bin = round(x_slide_mm / bin_size_mm, 4) * 1,
           y_bin = round(y_slide_mm / bin_size_mm, 4) * 1)
  if (nrow(m_s) < 50) next
  m_s$x_bin <- round(m_s$x_slide_mm / bin_size_mm) * bin_size_mm
  m_s$y_bin <- round(m_s$y_slide_mm / bin_size_mm) * bin_size_mm

  conds <- unique(m_s$Condition)
  for (cd in conds) {
    bins_cd <- m_s %>% filter(Condition == cd) %>%
      group_by(x_bin, y_bin) %>%
      summarise(n_cells = n(), n_necr = sum(necrosis_zone), .groups = "drop") %>%
      filter(n_cells >= 3, n_necr / n_cells >= 0.5)
    if (nrow(bins_cd) == 0) {
      max_region <- 0
    } else {
      max_region <- flood_fill_max_region(bins_cd %>% select(x_bin, y_bin))
    }
    cc_results[[length(cc_results) + 1]] <- tibble(
      slide_id = s, Condition = cd,
      necrotic_bins = nrow(bins_cd),
      largest_region_bins = max_region,
      largest_region_mm2 = max_region * bin_size_mm^2)
  }
}
cc_df <- bind_rows(cc_results) %>%
  mutate(treatment = case_when(Condition == "Control" ~ "NT",
                               grepl("MBRT", Condition) ~ "MBRT",
                               grepl("SBRT", Condition) ~ "SBRT"),
         timepoint_h = case_when(Condition == "Control" ~ 0,
                                 grepl("1h", Condition) ~ 1,
                                 grepl("4h", Condition) ~ 4,
                                 grepl("day2", Condition) ~ 48,
                                 grepl("day6", Condition) ~ 144))
write_tsv(cc_df, file.path(DATA_DIR, "necrosis_connected_components.tsv"))

# Plot: necrosis kinetics (per-slide percentile + calibrated)
nq_long <- necr_calib %>%
  filter(Condition %in% flank_conditions, treatment %in% c("MBRT","SBRT","NT")) %>%
  select(Condition, treatment, timepoint_h, necr_pct_per_slide, necr_pct_calib) %>%
  pivot_longer(cols = starts_with("necr_pct_"),
               names_to = "method", values_to = "necrosis_pct") %>%
  mutate(method = factor(method,
                         levels = c("necr_pct_per_slide", "necr_pct_calib"),
                         labels = c("Per-slide top 10%", "Control-calibrated (>95th pct of Control)")))

p_nk <- nq_long %>%
  ggplot(aes(x = factor(timepoint_h), y = necrosis_pct,
             color = treatment, group = treatment)) +
  geom_line(linewidth = 0.8) + geom_point(size = 3) +
  facet_wrap(~method) +
  scale_color_manual(values = TREATMENT_COLORS) +
  labs(title = "Necrosis-zone fraction kinetics (4T1 flank)",
       x = "Hours post-RT", y = "% cells in necrosis-zone",
       caption = "Per-slide: top 10% NN-distance per slide. Calibrated: > 95th percentile of Control.") +
  theme_bw()
ggsave(file.path(PLOT_DIR, "necrosis_kinetics.png"),
       plot = p_nk, width = 12, height = 5, dpi = 150)

# Plot: largest contiguous region size by condition
p_cc <- cc_df %>%
  filter(Condition %in% flank_conditions) %>%
  ggplot(aes(x = factor(timepoint_h), y = largest_region_mm2,
             color = treatment, group = treatment)) +
  geom_line(linewidth = 0.8) + geom_point(size = 3) +
  scale_color_manual(values = TREATMENT_COLORS) +
  labs(title = "Largest contiguous necrotic region size",
       x = "Hours post-RT", y = "Largest region (mm²)") +
  theme_bw()
ggsave(file.path(PLOT_DIR, "necrosis_largest_region.png"),
       plot = p_cc, width = 8, height = 5, dpi = 150)

# Cell-type composition of necrosis-zone cells
necr_celltypes <- as_tibble(obj@meta.data) %>%
  filter(!is.na(necrosis_zone), Condition %in% flank_conditions) %>%
  group_by(Condition, treatment, timepoint_h, necrosis_zone) %>%
  count(cell_type_major) %>%
  group_by(Condition, treatment, timepoint_h, necrosis_zone) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()
write_tsv(necr_celltypes, file.path(DATA_DIR, "necrosis_celltype_breakdown.tsv"))

p_nct <- necr_celltypes %>%
  filter(necrosis_zone == TRUE) %>%
  ggplot(aes(x = Condition, y = prop, fill = cell_type_major)) +
  geom_col() +
  scale_fill_brewer(palette = "Set1", name = "Cell type") +
  labs(title = "Cell-type composition of necrosis-zone cells (4T1 flank)",
       x = NULL, y = "Proportion") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(PLOT_DIR, "necrosis_celltype_breakdown.png"),
       plot = p_nct, width = 10, height = 5, dpi = 150)

# Spatial maps per condition × timepoint
flank_meta <- as_tibble(obj@meta.data) %>%
  filter(Condition %in% flank_conditions, !is.na(necrosis_zone)) %>%
  mutate(necrosis_label = ifelse(necrosis_zone, "necrosis_zone", "viable"))

p_spatial <- ggplot(flank_meta,
                    aes(x = x_slide_mm, y = y_slide_mm,
                        color = necrosis_label)) +
  geom_point(size = 0.05, alpha = 0.4) +
  scale_color_manual(values = c("viable" = "grey85", "necrosis_zone" = "red3")) +
  facet_wrap(~Condition, ncol = 4, scales = "free") +
  labs(title = "Necrosis spatial distribution (4T1 flank)",
       x = NULL, y = NULL) +
  theme_void() +
  theme(legend.position = "top",
        strip.text = element_text(size = 8),
        aspect.ratio = 1)
ggsave(file.path(PLOT_DIR, "necrosis_spatial.png"),
       plot = p_spatial, width = 16, height = 10, dpi = 120)

# ============================================================
cat("\n=== Sanity check ===\n")
cat("Necrosis-zone vs viable apoptosis means:\n")
print(necr_strat %>%
        filter(treatment %in% c("MBRT","SBRT","NT"), timepoint_h %in% c(48, 144)) %>%
        select(cell_type_major, treatment, timepoint_h, necrosis_zone, mean_apo, mean_sen) %>%
        arrange(treatment, timepoint_h, cell_type_major, necrosis_zone) %>%
        head(20))

cat("\n=== 04_celldeath_senescence.R complete ===\n")
