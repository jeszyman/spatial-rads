## run_all.R
## Master runner for Layer 2 Extended (MBRT vs SBRT, 4T1 flank).
## Sourced from dev/mbrt_vs_sbrt/. Each script self-contained — can also be run individually.
## Tier 1 (must): 00, 01, 02, 03, 04, 05
## Tier 2 (should): 10
## Tier 3 (next): 06, 07, 08
## Tier 4 (stretch): 09

scripts <- c(
  "00_load_and_filter.R",         # ~1 min — apply 3 filters, cache
  "01_degs_kinetics.R",           # ~15 min — DEGs by cell type x timepoint x contrast
  "02_pathway_kinetics.R",        # ~5 min — UCell signatures + Yi pathways
  "03_composition_trajectories.R",# ~30 sec — proportions over time
  "04_celldeath_senescence.R",    # ~10 min — apoptosis/senescence/necrosis quantification
  "05_spatial_nn.R",              # ~3 min — k=20 NN immune fraction kinetics
  ## Tier 3
  "08_set2_validation.R",         # ~10 min — Mutter_02 day2 concordance
  ## Tier 2 (synthesis after Tier 1+3)
  "10_synthesis.R"                # ~1 min — pull summaries
)

run_one <- function(s) {
  cat(sprintf("\n=== %s ===\n", s))
  t0 <- Sys.time()
  source(s, local = new.env())
  cat(sprintf(">>> %s done in %.1f min\n",
              s, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

for (s in scripts) {
  if (file.exists(s)) run_one(s) else cat(sprintf("Skipping (not present): %s\n", s))
}
