## 02_pathway_kinetics.R
## Score curated signatures (Hallmark + SenMayo + antigen-presentation + Yi's pathways)
## per cell, then aggregate per (cell_type_major x timepoint x condition).
## Necrosis_zone cells excluded at day2/day6.

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(UCell)
  library(msigdbr)
  library(BiocParallel)
})

DATA_DIR <- "data"
PLOT_DIR <- "plots"
OBJECTS_DIR <- "/mnt/data/projects/spatial-rads/analysis/objects/mbrt_vs_sbrt"

obj <- readRDS(file.path(OBJECTS_DIR, "seurat_filtered.rds"))
cat(sprintf("Loaded: %d genes x %d cells\n", nrow(obj), ncol(obj)))

# Necrosis-zone exclusion at day2/day6
exclude_necr_at <- c(48, 144)
keep_mask <- !(obj$timepoint_h %in% exclude_necr_at & obj$necrosis_zone %in% TRUE)
obj <- subset(obj, cells = colnames(obj)[keep_mask])
cat(sprintf("After necrosis-zone exclusion at day2/day6: %d cells\n", ncol(obj)))

# ========== Build signature lists ==========
panel_genes <- rownames(obj)
cat(sprintf("Panel size: %d genes\n", length(panel_genes)))

# 1. Hallmark mouse gene sets via msigdbr
hallmark_targets <- c(
  "HALLMARK_APOPTOSIS",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_HYPOXIA",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
)

# msigdbr v10+ uses different API; try both species options
m_df <- tryCatch(
  msigdbr(species = "mouse", collection = "H"),
  error = function(e) {
    cat("msigdbr collection arg failed; trying category arg...\n")
    msigdbr(species = "Mus musculus", category = "H")
  }
)
cat(sprintf("Hallmark sets fetched: %d unique pathways\n",
            length(unique(m_df$gs_name))))
cat(sprintf("Columns in m_df: %s\n", paste(colnames(m_df), collapse = ",")))

# Extract gene symbols column (depends on version)
gene_col <- intersect(c("gene_symbol", "mouse_gene_symbol", "gene_name"), colnames(m_df))[1]
gs_col   <- intersect(c("gs_name", "gs_id", "set_name"), colnames(m_df))[1]
stopifnot(!is.na(gene_col), !is.na(gs_col))

hallmark_lists <- m_df %>%
  dplyr::filter(.data[[gs_col]] %in% hallmark_targets) %>%
  dplyr::select(all_of(c(gs_col, gene_col))) %>%
  dplyr::distinct() %>%
  split(., .[[gs_col]]) %>%
  lapply(function(d) unique(d[[gene_col]]))

# 2. SenMayo (Saul et al. 2022) - hardcoded mouse gene symbols (125 genes)
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

# 3. Curated antigen presentation
antigen_pres_genes <- c("H2-K1","H2-D1","H2-Ab1","H2-Aa","H2-Eb1",
                        "B2m","Tap1","Tap2","Calr","Tapbp","Pdia3","Erap1")

# 4. Curated DAMP / alarmin (preview for script 04 4b — also useful here)
damp_genes <- c("Hmgb1","Hmgb2","S100a8","S100a9","Anxa1","Anxa2",
                "Il1b","Il6","Nlrp3","Casp1","Tlr2","Tlr4")

# Assemble all signatures
signature_lists <- c(
  hallmark_lists,
  list(SENMAYO_SENESCENCE      = senmayo_genes,
       CURATED_ANTIGEN_PRESENT = antigen_pres_genes,
       CURATED_DAMP            = damp_genes)
)

# Coverage report
cov_rows <- lapply(names(signature_lists), function(sn) {
  genes <- signature_lists[[sn]]
  in_panel <- intersect(genes, panel_genes)
  tibble(signature = sn,
         n_genes_total = length(genes),
         n_genes_in_panel = length(in_panel),
         coverage_pct = 100 * length(in_panel) / length(genes),
         genes_in_panel = paste(in_panel, collapse = ","))
})
coverage_df <- bind_rows(cov_rows) %>% arrange(coverage_pct)
write_tsv(coverage_df %>% select(-genes_in_panel),
          file.path(DATA_DIR, "signature_coverage.tsv"))
write_tsv(coverage_df, file.path(DATA_DIR, "signature_gene_lists.tsv"))
cat("\n=== Signature panel coverage ===\n")
print(coverage_df %>% select(signature, n_genes_total, n_genes_in_panel, coverage_pct))

# Filter signatures: keep only those with >=30% coverage AND >= 5 genes in panel
keep_sigs <- coverage_df %>%
  filter(coverage_pct >= 30, n_genes_in_panel >= 5) %>%
  pull(signature)
cat(sprintf("\nSignatures passing >=30%% coverage gate: %d / %d\n",
            length(keep_sigs), length(signature_lists)))
cat("  Kept:", paste(keep_sigs, collapse=", "), "\n")
cat("  Dropped (low coverage):",
    paste(setdiff(names(signature_lists), keep_sigs), collapse=", "), "\n")

# Subset to genes-in-panel for kept signatures
sig_in_panel <- lapply(keep_sigs, function(sn) {
  intersect(signature_lists[[sn]], panel_genes)
})
names(sig_in_panel) <- keep_sigs

# ========== Score signatures with UCell ==========
cat("\nScoring signatures with UCell (this may take 10-20 min on 547k cells)...\n")
t0 <- Sys.time()
obj <- AddModuleScore_UCell(obj,
                            features = sig_in_panel,
                            ncores = 8,
                            name = "_UCell")
cat(sprintf("UCell scoring done: %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# Score columns added: <signature_name>_UCell
sig_cols_ucell <- paste0(keep_sigs, "_UCell")
stopifnot(all(sig_cols_ucell %in% colnames(obj@meta.data)))

# ========== Yi's pre-computed pathway scores (already in metadata) ==========
yi_cols <- c("TypeI_interferon", "TypeII_interferon", "DNA_Damage_Repair", "STING")
yi_present <- yi_cols[yi_cols %in% colnames(obj@meta.data)]
cat(sprintf("\nYi's pathway columns present: %s\n", paste(yi_present, collapse = ",")))

all_score_cols <- c(yi_present, sig_cols_ucell)

# ========== Aggregate per (cell_type_major x timepoint x condition) ==========
agg <- as_tibble(obj@meta.data) %>%
  group_by(cell_type_major, timepoint_h, treatment, Condition) %>%
  summarise(across(all_of(all_score_cols),
                   list(mean = ~mean(., na.rm = TRUE),
                        sd   = ~sd(., na.rm = TRUE))),
            n_cells = n(),
            .groups = "drop")

agg_long <- agg %>%
  pivot_longer(cols = matches("_(mean|sd)$"),
               names_to = c("signature", ".value"),
               names_pattern = "(.+)_(mean|sd)$")
write_tsv(agg_long, file.path(DATA_DIR, "pathway_kinetics.tsv"))
cat(sprintf("Wrote: %s (%d rows)\n",
            file.path(DATA_DIR, "pathway_kinetics.tsv"), nrow(agg_long)))

# ========== Plots: kinetics line plots per cell type ==========
TREATMENT_COLORS <- c("MBRT" = "#E74C3C", "SBRT" = "#3498DB", "NT" = "#7F8C8D")

# Convert wide signature names for prettier panels
nice_label <- function(x) {
  x %>%
    gsub("HALLMARK_", "", .) %>%
    gsub("_UCell$", "", .) %>%
    gsub("_", " ", .) %>%
    str_to_title()
}

agg_long <- agg_long %>%
  mutate(signature_label = nice_label(signature))

MAJOR_TYPES_ORDER <- c("tumor_epithelial", "endothelial", "fibroblast", "smooth_muscle",
                       "T.cell", "B.cell", "myeloid", "NK", "plasma")

# One plot per cell type: signature panels x time, lines per treatment
for (ct in MAJOR_TYPES_ORDER) {
  d <- agg_long %>% filter(cell_type_major == ct)
  if (nrow(d) == 0) next
  p <- d %>%
    ggplot(aes(x = factor(timepoint_h), y = mean, color = treatment, group = treatment)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2, alpha = 0.4) +
    facet_wrap(~signature_label, scales = "free_y", ncol = 4) +
    scale_color_manual(values = TREATMENT_COLORS) +
    labs(title = sprintf("Pathway kinetics: %s", ct),
         x = "Hours post-RT", y = "Mean signature score",
         caption = "n=1 per condition; necrosis-zone cells excluded at day2/day6.") +
    theme_bw() +
    theme(strip.text = element_text(size = 9))
  ggsave(file.path(PLOT_DIR, sprintf("pathway_kinetics_%s.png", ct)),
         plot = p, width = 14, height = 10, dpi = 130)
}

# Top-level synthesis: heatmap of (signature x cell_type x timepoint) — MBRT vs SBRT delta
delta_df <- agg_long %>%
  filter(treatment %in% c("MBRT", "SBRT"), cell_type_major %in% MAJOR_TYPES_ORDER) %>%
  select(cell_type_major, timepoint_h, signature_label, treatment, mean) %>%
  pivot_wider(names_from = treatment, values_from = mean) %>%
  mutate(delta_MBRT_SBRT = MBRT - SBRT) %>%
  filter(!is.na(delta_MBRT_SBRT))

p_delta <- delta_df %>%
  ggplot(aes(x = factor(timepoint_h), y = cell_type_major, fill = delta_MBRT_SBRT)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#3498DB", mid = "white", high = "#E74C3C",
                       midpoint = 0, name = "MBRT - SBRT") +
  facet_wrap(~signature_label, ncol = 4) +
  labs(title = "MBRT - SBRT score delta by cell type and timepoint",
       x = "Hours post-RT", y = NULL,
       caption = "n=1 per condition; necrosis-zone cells excluded at day2/day6.") +
  theme_bw() +
  theme(panel.grid = element_blank(), strip.text = element_text(size = 9))
ggsave(file.path(PLOT_DIR, "pathway_delta_mbrt_vs_sbrt.png"),
       plot = p_delta, width = 14, height = 10, dpi = 130)

write_tsv(delta_df, file.path(DATA_DIR, "pathway_delta_mbrt_vs_sbrt.tsv"))

cat("\n=== Sanity check: DDR + IFN kinetics (tumor_epithelial) ===\n")
sanity <- agg_long %>%
  filter(cell_type_major == "tumor_epithelial",
         signature %in% c("DNA_Damage_Repair", "TypeI_interferon",
                          "HALLMARK_INTERFERON_ALPHA_RESPONSE_UCell",
                          "HALLMARK_APOPTOSIS_UCell"),
         treatment %in% c("MBRT","SBRT","NT")) %>%
  select(signature_label, timepoint_h, treatment, mean) %>%
  pivot_wider(names_from = treatment, values_from = mean) %>%
  arrange(signature_label, timepoint_h)
print(sanity)

cat("\n=== 02_pathway_kinetics.R complete ===\n")
