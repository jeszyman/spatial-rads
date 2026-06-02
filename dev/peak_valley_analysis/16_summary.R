if (!exists("PATHWAY_COLS")) source("00_load_data.R")
library(patchwork)

# --- Load Layer 4 outputs ---
kinetics  <- read_tsv(file.path(DATA_DIR, "signature_kinetics.tsv"), show_col_types = FALSE)
stripe    <- read_tsv(file.path(DATA_DIR, "stripe_persistence.tsv"), show_col_types = FALSE)
scores    <- read_tsv(file.path(DATA_DIR, "signature_scores_all_cells.tsv"), show_col_types = FALSE)
cor_df    <- read_tsv(file.path(DATA_DIR, "scoring_method_correlation.tsv"), show_col_types = FALSE)
peak_sig  <- read_tsv(file.path(DATA_DIR, "peak_signature_bulk.tsv"), show_col_types = FALSE)
valley_sig <- read_tsv(file.path(DATA_DIR, "valley_signature_bulk.tsv"), show_col_types = FALSE)

# --- Summary table ---
summary_tbl <- kinetics %>%
  mutate(Condition = factor(Condition, levels = CONDITION_LEVELS)) %>%
  left_join(stripe %>% select(condition, morans_i, morans_p),
            by = c("Condition" = "condition")) %>%
  select(Condition, timepoint_h, treatment, model,
         mean_peak, sd_peak, mean_valley, sd_valley, n,
         morans_i, morans_p) %>%
  arrange(Condition)

write_tsv(summary_tbl, file.path(DATA_DIR, "layer4_summary.tsv"))
cat("Layer 4 summary table:\n")
print(as.data.frame(summary_tbl))

# --- Composite figure ---
# Panel A: kinetics (flank, MBRT + Control)
kinetics_long <- kinetics %>%
  filter(model == "flank", treatment %in% c("MBRT", "NT")) %>%
  pivot_longer(cols = c(mean_peak, mean_valley),
               names_to = "signature", values_to = "score",
               names_prefix = "mean_") %>%
  mutate(signature = factor(signature, levels = c("peak", "valley"),
                            labels = c("Peak-up", "Valley-up")))

pA <- ggplot(kinetics_long, aes(x = timepoint_h, y = score, color = treatment,
                                linetype = signature,
                                group = interaction(treatment, signature))) +
  geom_line(linewidth = 0.8) + geom_point(size = 2.5) +
  scale_color_manual(values = TREATMENT_COLORS) +
  scale_x_continuous(breaks = c(0, 1, 4, 48, 144),
                     labels = c("0", "1h", "4h", "2d", "6d")) +
  labs(title = "A. Signature kinetics (flank, MBRT vs Control)",
       x = "Time post-RT", y = "UCell score") +
  theme_bw(base_size = 11)

# Panel B: Moran's I
pB <- ggplot(stripe, aes(x = condition, y = morans_i)) +
  geom_col(fill = "#E74C3C", alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "B. Spatial autocorrelation (peak-up)",
       x = NULL, y = "Moran's I") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# Panel C: signature gene counts
sig_info <- tibble(
  category = c("Peak-up genes", "Valley-up genes", "DDR excluded"),
  count = c(nrow(peak_sig), nrow(valley_sig), 18)
)
pC <- ggplot(sig_info, aes(x = category, y = count, fill = category)) +
  geom_col(show.legend = FALSE) +
  scale_fill_manual(values = c("Peak-up genes" = "#E74C3C",
                                "Valley-up genes" = "#3498DB",
                                "DDR excluded" = "grey50")) +
  geom_text(aes(label = count), vjust = -0.3) +
  labs(title = "C. Signature composition", x = NULL, y = "N genes") +
  theme_bw(base_size = 11)

p_composite <- (pA | (pB / pC)) + plot_layout(widths = c(2, 1))
ggsave(file.path(PLOT_DIR, "layer4_composite.png"),
       plot = p_composite, width = 16, height = 8, dpi = 150)
cat("Saved: layer4_composite.png\n")

# --- Print summary stats ---
cat("\n=== Layer 4 Summary ===\n")
cat(sprintf("Signature genes: %d peak-up, %d valley-up (18 DDR excluded)\n",
            nrow(peak_sig), nrow(valley_sig)))
cat(sprintf("UCell vs AMS correlation: peak r=%.3f, valley r=%.3f\n",
            cor_df$spearman_r[1], cor_df$spearman_r[2]))
cat(sprintf("Conditions scored: %d\n", n_distinct(summary_tbl$Condition)))
cat(sprintf("Total cells scored: %d\n", sum(summary_tbl$n)))
cat("Layer 4 summary complete.\n")
