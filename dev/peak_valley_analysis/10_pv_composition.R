if (!exists("PATHWAY_COLS")) source("00_load_data.R")

pv <- read_tsv(file.path(DATA_DIR, "mbrt4h_peak_valley.tsv"), show_col_types = FALSE)

# ---- Cell type proportions by zone ----
comp <- pv %>%
  count(zone, cell_type_validated) %>%
  group_by(zone) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()
write_tsv(comp, file.path(DATA_DIR, "pv_composition.tsv"))

# Stacked bar
p_stack <- ggplot(comp, aes(x = zone, y = prop, fill = cell_type_validated)) +
  geom_col() +
  labs(title = "Cell type composition: Peak vs Valley (MBRT_4h)",
       y = "Proportion", fill = "Cell type") +
  theme_bw()
ggsave(file.path(PLOT_DIR, "pv_composition_stacked.png"), plot = p_stack, width = 8, height = 6, dpi = 150)

# Proportion difference (peak - valley)
diff_df <- comp %>%
  select(zone, cell_type_validated, prop) %>%
  pivot_wider(names_from = zone, values_from = prop, values_fill = 0) %>%
  mutate(diff = peak - valley) %>%
  arrange(desc(abs(diff)))

cat("=== Composition difference (peak - valley) ===\n")
print(diff_df)

p_diff <- ggplot(diff_df, aes(x = reorder(cell_type_validated, diff), y = diff,
                               fill = diff > 0)) +
  geom_col() +
  scale_fill_manual(values = c("TRUE" = ZONE_COLORS["peak"], "FALSE" = ZONE_COLORS["valley"]),
                    labels = c("TRUE" = "Enriched in peaks", "FALSE" = "Enriched in valleys"),
                    guide = guide_legend(title = NULL)) +
  coord_flip() +
  labs(title = "Cell type enrichment: Peak vs Valley", x = NULL,
       y = "Proportion difference (peak - valley)") +
  theme_bw()
ggsave(file.path(PLOT_DIR, "pv_composition_diff.png"), plot = p_diff, width = 10, height = 6, dpi = 150)

cat("Composition analysis complete.\n")
