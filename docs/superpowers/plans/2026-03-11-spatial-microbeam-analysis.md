# Spatial Microbeam Analysis Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Exploratory analysis of CosMx spatial transcriptomics data comparing MBRT vs SBRT vs Control across timepoints, organized as four progressive layers.

**Architecture:** Development follows the frag-length workup pattern: write standalone R scripts first for fast interactive iteration, then roll working code into org-babel blocks in `spatial-rads.org`. Production pipeline blocks tangle to numbered scripts (`00_`, `01_`, ...); exploratory/diagnostic blocks use `:tangle no`. Persistent R session (`:session spatial`) across all blocks.

**Tech Stack:** R (Seurat, tidyverse, ggplot2, data.table, future, arrow, RANN, patchwork), Python (pandas, pyarrow), org-babel, conda (spatial-rads env)

**Spec:** `docs/superpowers/specs/2026-03-11-spatial-microbeam-analysis-design.md`

---

## Development Workflow

1. **Develop interactively** in standalone `.R` files under `dev/` (fast iteration, no org overhead)
2. **Once working**, roll code into org-babel blocks in `spatial-rads.org`
3. **Production blocks** (`:tangle scripts/NN_name.R`) — numbered sequentially, form the reproducible pipeline
4. **Exploratory blocks** (`:tangle no`) — diagnostics, plots, sensitivity analyses; stay in org but don't tangle
5. **Persistent session** — all blocks share `:session spatial :async yes`
6. **Intermediate outputs** — save as TSV/RDS to `/mnt/gcs/jeszyman/projects/spatial-rads/analysis/` for caching
7. **Explicit reporting** — every block prints key stats (`cat(sprintf(...))`)

## Design Rationale Table

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Cell type labels | Yi's insitutype (Phase 1) | Fast path; validate with markers before trusting |
| Pathway scores | Yi's pre-computed columns | Provenance TBD; recompute in Phase 2 |
| DEG method | Seurat FindMarkers (Wilcoxon) | Standard; rank by effect size, not p-value |
| DEG reporting | log2FC + pct expressed | n=1 pseudo-replication; p-values meaningless |
| Normalization | LogNormalize 1e4, VST 2000 | Standard Seurat workflow for targeted panels |
| PCs | 20 (validate with elbow) | Reasonable for ~1000 gene panel |
| Clustering | Resolution 0.4 | Starting point; adjust if needed |
| Peak/valley classification | Gamma-H2AX IF (protein) | Ground truth; transcript DDR is fallback |
| Batch correction | Assess first, correct if needed | UMAP-by-Slide check in Layer 1 |
| Control model | Treat as flank (TBD) | Confirm with lab |

---

## Conventions

- Before adding/modifying org headings or src blocks, invoke the `org-edit` skill
- Invoke `bioinfo-dev` skill when writing analysis code
- Output paths: `analysis/figures/`, `analysis/objects/`, `analysis/tables/`
- All under `/mnt/gcs/jeszyman/projects/spatial-rads/`
- Existing exploratory code in `* Exploratory analysis` gets reorganized into `* Methods` as it's formalized

---

## Chunk 1: Layer 1 — QC & Landscape

### Task 1: Set up environment and verify data

- [ ] **Step 1: Create output directories**

```bash
mkdir -p /mnt/gcs/jeszyman/projects/spatial-rads/analysis/{figures,objects,tables}
mkdir -p ~/repos/spatial-rads/dev
```

- [ ] **Step 2: Verify input files exist**

```bash
ls -lh /mnt/gcs/jeszyman/projects/spatial-rads/inputs/*.parquet
```

- [ ] **Step 3: Verify conda env has all required packages**

```bash
conda run -n spatial-rads R -e "library(Seurat); library(tidyverse); library(data.table); library(future); library(arrow); library(RANN); library(patchwork); cat('All packages loaded\n')"
```

Install any missing packages before proceeding.

---

### Task 2: Data loading script (dev → org)

**Dev file:** `dev/00_load_data.R`

- [ ] **Step 1: Write standalone R script for data loading**

Write `dev/00_load_data.R` interactively. This loads parquet files, builds the Seurat object, runs QC filtering, and saves a clean object. Include explicit reporting at each step.

```r
library(data.table)
library(Seurat)
library(arrow)

# --- Load counts ---
counts_df <- read_parquet("/mnt/gcs/jeszyman/projects/spatial-rads/inputs/projects_Mutter_01_CosMmR_app_data_raw_tx_counts_matrix.parquet")
counts_dt <- as.data.table(counts_df)

non_gene_cols <- c("Slide", "fov", "cell_id")
gene_cols <- setdiff(colnames(counts_dt), non_gene_cols)
cat(sprintf("Panel size: %d genes\n", length(gene_cols)))
cat(sprintf("Total cells in counts: %d\n", nrow(counts_dt)))
stopifnot(!anyDuplicated(counts_dt$cell_id))

# --- Load metadata ---
meta_df <- read_parquet("/mnt/gcs/jeszyman/projects/spatial-rads/inputs/Analysis_Mutter_01_CosMmR_Mutter_updated_metadata.parquet")
meta_dt <- as.data.table(meta_df)
cat(sprintf("Total cells in metadata: %d\n", nrow(meta_dt)))

# --- Build count matrix (genes x cells) ---
count_mat <- t(as.matrix(counts_dt[, ..gene_cols]))
colnames(count_mat) <- counts_dt$cell_id
rownames(count_mat) <- gene_cols

meta_for_seurat <- as.data.frame(meta_dt)
rownames(meta_for_seurat) <- meta_for_seurat$cell_id
stopifnot(all(colnames(count_mat) %in% rownames(meta_for_seurat)))

# --- Create Seurat object ---
obj <- CreateSeuratObject(
  counts    = count_mat,
  meta.data = meta_for_seurat,
  project   = "CosMx_Mutter"
)
cat(sprintf("Raw Seurat: %d genes x %d cells\n", nrow(obj), ncol(obj)))

# --- QC filtering ---
pre_qc <- table(obj$Condition)
obj <- subset(obj, subset = qcFlagsCell == "Pass")
obj <- subset(obj, subset = nCount_RNA > 20 & nFeature_RNA > 10 & propNegative < 0.5)
post_qc <- table(obj$Condition)

qc_summary <- data.frame(
  condition = names(pre_qc),
  pre_filter = as.integer(pre_qc),
  post_filter = as.integer(post_qc[names(pre_qc)]),
  pct_retained = round(as.integer(post_qc[names(pre_qc)]) / as.integer(pre_qc) * 100, 1)
)
print(qc_summary)
write.csv(qc_summary, "/mnt/gcs/jeszyman/projects/spatial-rads/analysis/tables/qc_summary.csv", row.names = FALSE)
cat(sprintf("Post-QC: %d cells retained\n", ncol(obj)))

# --- Add parsed metadata ---
obj$model <- ifelse(grepl("^Tongue", obj$Condition), "tongue", "flank")
obj$treatment <- dplyr::case_when(
  obj$Condition == "Control" ~ "NT",
  grepl("MBRT", obj$Condition) ~ "MBRT",
  grepl("SBRT", obj$Condition) ~ "SBRT"
)
obj$timepoint_h <- dplyr::case_when(
  obj$Condition == "Control" ~ 0,
  grepl("1h", obj$Condition) ~ 1,
  grepl("4h", obj$Condition) ~ 4,
  grepl("day2", obj$Condition) ~ 48,
  grepl("day6", obj$Condition) ~ 144,
  grepl("day8", obj$Condition) ~ 192,
  grepl("day10", obj$Condition) ~ 240
)

# --- Verify pathway columns exist ---
pathway_cols <- c("TypeI_interferon", "TypeII_interferon", "STING", "DNA_Damage_Repair")
stopifnot(all(pathway_cols %in% colnames(obj@meta.data)))
cat("Pathway score columns confirmed present.\n")

# --- Save ---
saveRDS(obj, "/mnt/gcs/jeszyman/projects/spatial-rads/analysis/objects/seurat_qc.rds")
cat("QC'd Seurat object saved.\n")
```

- [ ] **Step 2: Run interactively, verify outputs**

Expected: ~1000 genes, >80% cell retention per sample, pathway columns present.

- [ ] **Step 3: Roll into org-babel block**

Add to `spatial-rads.org` under `* Methods > ** Layer 1 — QC & Landscape > *** Data loading and QC` as a tangled block:

```
#+begin_src R :session spatial :async yes :tangle scripts/00_load_data.R
... (working code from dev/)
#+end_src
```

---

### Task 3: Normalization, clustering, UMAP (dev → org)

**Dev file:** `dev/01_normalize_cluster.R`

- [ ] **Step 1: Write standalone script**

```r
library(Seurat)
library(future)

obj <- readRDS("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/objects/seurat_qc.rds")

# --- Normalize ---
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 1e4)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj, features = VariableFeatures(obj))
obj <- RunPCA(obj, features = VariableFeatures(obj))

pdf("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/figures/elbow_plot.pdf", width = 6, height = 4)
ElbowPlot(obj, ndims = 30)
dev.off()
cat("Elbow plot saved.\n")

# --- Cluster (parallel) ---
options(future.globals.maxSize = 50 * 1024^3)
plan(multicore, workers = 48)
obj <- FindNeighbors(obj, dims = 1:20)
obj <- FindClusters(obj, resolution = 0.4)
plan(sequential)

# --- UMAP ---
obj <- RunUMAP(obj, dims = 1:20, n.neighbors = 30, min.dist = 0.3)
cat(sprintf("Clusters: %d | Cells: %d\n", length(levels(obj$seurat_clusters)), ncol(obj)))

# --- Save ---
saveRDS(obj, "/mnt/gcs/jeszyman/projects/spatial-rads/analysis/objects/seurat_clustered.rds")
cat("Clustered object saved.\n")
```

- [ ] **Step 2: Run, review elbow plot, verify cluster count**
- [ ] **Step 3: Roll into org-babel block** (`:tangle scripts/01_normalize_cluster.R`)

---

### Task 4: Cell type validation (exploratory, no tangle)

- [ ] **Step 1: Write validation script in `dev/explore_celltype_validation.R`**

```r
library(Seurat)
library(tidyverse)

obj <- readRDS("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/objects/seurat_clustered.rds")

# --- Tabulate Yi's labels ---
ct_table <- sort(table(obj$ImmuneAtlas_ImmGen_Main_cell_Types), decreasing = TRUE)
cat("Cell type frequencies:\n")
print(ct_table)

ct_family <- sort(table(obj$ImmuneAtlas_ImmGen_cellFamily_Main_cell_Types), decreasing = TRUE)
cat("\nCell family frequencies:\n")
print(ct_family)

# --- Check available canonical markers ---
marker_candidates <- c(
  "Cd3e", "Cd3d", "Cd4", "Cd8a", "Cd8b1",
  "Cd68", "Adgre1", "Csf1r", "Itgam",
  "Epcam", "Krt8", "Krt18", "Krt14",
  "Cd19", "Ms4a1",
  "Ncr1", "Klrb1c",
  "Pecam1", "Cdh5"
)
available_markers <- marker_candidates[marker_candidates %in% rownames(obj)]
cat(sprintf("\nAvailable markers: %d / %d\n", length(available_markers), length(marker_candidates)))
print(available_markers)

# --- Dot plot: markers vs cell types ---
DotPlot(obj, features = available_markers,
        group.by = "ImmuneAtlas_ImmGen_Main_cell_Types") +
  RotatedAxis() +
  labs(title = "Marker validation: Yi's cell type labels")
ggsave("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/figures/celltype_marker_validation.pdf",
       width = 14, height = 8)

# --- Characterize "a" population ---
epi_markers <- available_markers[available_markers %in% c("Epcam", "Krt8", "Krt18", "Krt14")]
if (length(epi_markers) > 0) {
  VlnPlot(obj, features = epi_markers,
          group.by = "ImmuneAtlas_ImmGen_Main_cell_Types", pt.size = 0) +
    labs(title = "Is 'a' tumor/epithelial?")
  ggsave("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/figures/a_population_markers.pdf",
         width = 12, height = 6)
}
```

- [ ] **Step 2: Run, review plots, decide on "a" relabeling**

If "a" is confirmed epithelial:
```r
obj$cell_type_validated <- obj$ImmuneAtlas_ImmGen_Main_cell_Types
obj$cell_type_validated[obj$cell_type_validated == "a"] <- "tumor_epithelial"
saveRDS(obj, "/mnt/gcs/jeszyman/projects/spatial-rads/analysis/objects/seurat_clustered.rds")
```

- [ ] **Step 3: Roll into org-babel block** (`:tangle no` — exploratory)

---

### Task 5: Batch assessment and landscape plots (exploratory, no tangle)

- [ ] **Step 1: Write `dev/explore_landscape.R`**

```r
library(Seurat)
library(patchwork)

obj <- readRDS("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/objects/seurat_clustered.rds")

# --- Batch assessment ---
DimPlot(obj, group.by = "Slide", raster = TRUE) +
  labs(title = "UMAP by Slide — batch check")
ggsave("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/figures/umap_by_slide.pdf",
       width = 8, height = 6)

# --- Landscape UMAPs ---
p1 <- DimPlot(obj, group.by = "seurat_clusters", raster = TRUE, label = TRUE) + labs(title = "Clusters")
p2 <- DimPlot(obj, group.by = "Condition", raster = TRUE) + labs(title = "Condition")
p3 <- DimPlot(obj, group.by = "ImmuneAtlas_ImmGen_Main_cell_Types",
              raster = TRUE, label = TRUE, repel = TRUE) + labs(title = "Cell type")
p4 <- DimPlot(obj, group.by = "treatment", raster = TRUE) + labs(title = "Treatment")
(p1 | p2) / (p3 | p4)
ggsave("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/figures/umap_landscape.pdf",
       width = 16, height = 12)

# --- Spatial cell type maps (flank) ---
obj@meta.data %>%
  as.data.frame() %>%
  filter(model == "flank") %>%
  ggplot(aes(x = x_slide_mm, y = y_slide_mm, color = ImmuneAtlas_ImmGen_Main_cell_Types)) +
  geom_point(size = 0.1, alpha = 0.3) +
  facet_wrap(~Condition, scales = "free") +
  labs(title = "Spatial cell type distribution (flank)", color = "Cell type") +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))
ggsave("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/figures/spatial_celltype_flank.pdf",
       width = 16, height = 12)

# --- Cell count summary ---
cell_counts <- obj@meta.data %>%
  as.data.frame() %>%
  count(Condition, ImmuneAtlas_ImmGen_Main_cell_Types) %>%
  pivot_wider(names_from = ImmuneAtlas_ImmGen_Main_cell_Types, values_from = n, values_fill = 0)
write.csv(cell_counts, "/mnt/gcs/jeszyman/projects/spatial-rads/analysis/tables/cell_counts_by_condition_celltype.csv", row.names = FALSE)
cat("Cell count table saved.\n")
```

- [ ] **Step 2: Run, review. Flag if Slide dominates UMAP (needs Harmony).**
- [ ] **Step 3: Roll into org-babel blocks** (`:tangle no`)

- [ ] **Step 4: Commit Layer 1**

```bash
git add spatial-rads.org dev/
git commit -m "feat: Layer 1 QC, clustering, cell type validation, landscape plots"
```

---

## Chunk 2: Layer 2 — Comparative Kinetics

### Task 6: Differential expression (dev → org)

**Dev file:** `dev/02_deg_analysis.R`

- [ ] **Step 1: Write DEG analysis script**

```r
library(Seurat)
library(tidyverse)

obj <- readRDS("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/objects/seurat_clustered.rds")

# --- DEG function (effect-size ranked, NOT p-value based) ---
run_deg <- function(obj, ident.1, ident.2, group.by = "Condition") {
  Idents(obj) <- group.by
  markers <- FindMarkers(obj, ident.1 = ident.1, ident.2 = ident.2,
                         min.pct = 0.1, logfc.threshold = 0.1)
  markers$gene <- rownames(markers)
  markers$comparison <- paste0(ident.1, "_vs_", ident.2)
  markers %>% arrange(desc(abs(avg_log2FC)))
}

# --- MBRT vs Control ---
mbrt_tp <- c("MBRT_1h", "MBRT_4h", "MBRT_day2", "MBRT_day6")
deg_mbrt_ctrl <- map_dfr(mbrt_tp, function(tp) {
  sub <- subset(obj, cells = colnames(obj)[obj$Condition %in% c(tp, "Control")])
  run_deg(sub, tp, "Control")
})

# --- SBRT vs Control ---
sbrt_tp <- c("SBRT_4h", "SBRT_day2", "SBRT_day6")
deg_sbrt_ctrl <- map_dfr(sbrt_tp, function(tp) {
  sub <- subset(obj, cells = colnames(obj)[obj$Condition %in% c(tp, "Control")])
  run_deg(sub, tp, "Control")
})

# --- MBRT vs SBRT at matched timepoints ---
matched <- list(c("MBRT_4h", "SBRT_4h"), c("MBRT_day2", "SBRT_day2"), c("MBRT_day6", "SBRT_day6"))
deg_mbrt_sbrt <- map_dfr(matched, function(pair) {
  sub <- subset(obj, cells = colnames(obj)[obj$Condition %in% pair])
  run_deg(sub, pair[1], pair[2])
})

# --- Combine and save ---
all_degs <- bind_rows(
  mbrt_vs_ctrl = deg_mbrt_ctrl,
  sbrt_vs_ctrl = deg_sbrt_ctrl,
  mbrt_vs_sbrt = deg_mbrt_sbrt,
  .id = "comparison_type"
)
write.csv(all_degs, "/mnt/gcs/jeszyman/projects/spatial-rads/analysis/tables/deg_all_cells.csv", row.names = FALSE)
cat(sprintf("Total DEG rows: %d\n", nrow(all_degs)))

# --- Per-cell-type DEGs (MBRT vs SBRT at 4h) ---
ct_counts <- table(obj$ImmuneAtlas_ImmGen_Main_cell_Types)
major_cts <- names(ct_counts[ct_counts > 1000])
cat(sprintf("Major cell types (>1000 cells): %s\n", paste(major_cts, collapse = ", ")))

deg_by_ct <- map_dfr(major_cts, function(ct) {
  cells <- colnames(obj)[obj$Condition %in% c("MBRT_4h", "SBRT_4h") &
                          obj$ImmuneAtlas_ImmGen_Main_cell_Types == ct]
  if (length(cells) < 50) return(NULL)
  sub <- subset(obj, cells = cells)
  tryCatch(run_deg(sub, "MBRT_4h", "SBRT_4h") %>% mutate(cell_type = ct),
           error = function(e) NULL)
})
write.csv(deg_by_ct, "/mnt/gcs/jeszyman/projects/spatial-rads/analysis/tables/deg_by_celltype_4h.csv", row.names = FALSE)
```

- [ ] **Step 2: Run, review top DEGs for biological plausibility**
- [ ] **Step 3: Roll tangled portion into org** (`:tangle scripts/02_deg_analysis.R`)

---

### Task 7: DEG visualization (exploratory, no tangle)

- [ ] **Step 1: Write `dev/explore_deg_plots.R`**

```r
library(tidyverse)

all_degs <- read.csv("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/tables/deg_all_cells.csv")

# --- Volcano: MBRT vs SBRT ---
all_degs %>%
  filter(comparison_type == "mbrt_vs_sbrt") %>%
  ggplot(aes(x = avg_log2FC, y = -log10(p_val), color = abs(avg_log2FC) > 0.5)) +
  geom_point(alpha = 0.5, size = 0.5) +
  facet_wrap(~comparison) +
  scale_color_manual(values = c("grey60", "red3"), guide = "none") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey40") +
  labs(title = "MBRT vs SBRT DEGs (effect size ranked)",
       x = "log2FC", y = "-log10(p)",
       caption = "n=1 per condition; p-values for visualization only") +
  theme_bw()
ggsave("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/figures/volcano_mbrt_vs_sbrt.pdf",
       width = 12, height = 5)

# --- Top DEGs heatmap across conditions ---
top_genes <- all_degs %>%
  filter(comparison_type == "mbrt_vs_sbrt") %>%
  group_by(comparison) %>%
  slice_max(abs(avg_log2FC), n = 10) %>%
  pull(gene) %>% unique()
cat(sprintf("Top genes for heatmap: %s\n", paste(top_genes, collapse = ", ")))
```

- [ ] **Step 2: Run, review**
- [ ] **Step 3: Roll into org** (`:tangle no`)

---

### Task 8: Pathway score kinetics (exploratory, no tangle)

- [ ] **Step 1: Write `dev/explore_pathway_kinetics.R`**

```r
library(tidyverse)

obj <- readRDS("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/objects/seurat_clustered.rds")
pathway_cols <- c("TypeI_interferon", "TypeII_interferon", "STING", "DNA_Damage_Repair")

# --- By cell family ---
kinetics <- obj@meta.data %>%
  as.data.frame() %>%
  filter(model == "flank") %>%
  group_by(timepoint_h, treatment,
           cell_family = ImmuneAtlas_ImmGen_cellFamily_Main_cell_Types) %>%
  summarise(n_cells = n(),
            across(all_of(pathway_cols), mean, na.rm = TRUE),
            .groups = "drop") %>%
  pivot_longer(cols = all_of(pathway_cols), names_to = "pathway", values_to = "score")

kinetics %>%
  filter(treatment %in% c("MBRT", "SBRT", "NT")) %>%
  ggplot(aes(x = factor(timepoint_h), y = score, color = treatment, group = treatment)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  facet_grid(pathway ~ cell_family, scales = "free_y") +
  labs(x = "Hours post-RT", y = "Mean pathway score", caption = "n=1 per condition") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/figures/pathway_kinetics_by_celltype.pdf",
       width = 20, height = 12)

# --- All cells (simplified) ---
obj@meta.data %>%
  as.data.frame() %>%
  filter(model == "flank", treatment %in% c("MBRT", "SBRT")) %>%
  group_by(timepoint_h, treatment) %>%
  summarise(across(all_of(pathway_cols), mean, na.rm = TRUE), .groups = "drop") %>%
  pivot_longer(cols = all_of(pathway_cols), names_to = "pathway", values_to = "score") %>%
  ggplot(aes(x = factor(timepoint_h), y = score, color = treatment, group = treatment)) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  facet_wrap(~pathway, scales = "free_y") +
  scale_color_manual(values = c("MBRT" = "#E74C3C", "SBRT" = "#3498DB")) +
  labs(x = "Hours post-RT", y = "Mean pathway score", caption = "n=1 per condition") +
  theme_bw()
ggsave("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/figures/pathway_kinetics_all_cells.pdf",
       width = 10, height = 8)

write.csv(kinetics, "/mnt/gcs/jeszyman/projects/spatial-rads/analysis/tables/pathway_kinetics.csv", row.names = FALSE)
```

- [ ] **Step 2: Run, review kinetic patterns**
- [ ] **Step 3: Roll into org** (`:tangle no`)

---

### Task 9: Immune composition (exploratory, no tangle)

- [ ] **Step 1: Write `dev/explore_composition.R`**

```r
library(tidyverse)

obj <- readRDS("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/objects/seurat_clustered.rds")

# --- Proportions ---
composition <- obj@meta.data %>%
  as.data.frame() %>%
  filter(model == "flank") %>%
  count(Condition, timepoint_h, treatment,
        cell_type = ImmuneAtlas_ImmGen_Main_cell_Types) %>%
  group_by(Condition) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# --- Stacked bar ---
cond_levels <- c("Control", "MBRT_1h", "MBRT_4h", "SBRT_4h",
                 "MBRT_day2", "SBRT_day2", "MBRT_day6", "SBRT_day6")
composition %>%
  filter(Condition %in% cond_levels) %>%
  mutate(Condition = factor(Condition, levels = cond_levels)) %>%
  ggplot(aes(x = Condition, y = prop, fill = cell_type)) +
  geom_col() +
  labs(x = NULL, y = "Proportion", fill = "Cell type",
       title = "Cell type composition (flank)") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/figures/composition_stacked.pdf",
       width = 10, height = 6)

# --- Trajectories for key immune types ---
# Adapt cell type names after Task 4 validation
composition %>%
  filter(treatment %in% c("MBRT", "SBRT")) %>%
  ggplot(aes(x = factor(timepoint_h), y = prop, color = treatment, group = treatment)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  facet_wrap(~cell_type, scales = "free_y") +
  labs(x = "Hours post-RT", y = "Proportion", color = "Treatment",
       caption = "n=1 per condition") +
  theme_bw()
ggsave("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/figures/composition_trajectories.pdf",
       width = 14, height = 10)

write.csv(composition, "/mnt/gcs/jeszyman/projects/spatial-rads/analysis/tables/composition.csv", row.names = FALSE)
```

- [ ] **Step 2: Run, review**
- [ ] **Step 3: Roll into org** (`:tangle no`)

---

### Task 10: Spatial neighborhood analysis (exploratory, no tangle)

- [ ] **Step 1: Write `dev/explore_spatial_nn.R`**

```r
library(Seurat)
library(RANN)
library(tidyverse)

obj <- readRDS("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/objects/seurat_clustered.rds")

# --- Nearest-neighbor immune fraction at 4h ---
nn_results <- map_dfr(c("MBRT_4h", "SBRT_4h"), function(cond) {
  meta_sub <- obj@meta.data %>% as.data.frame() %>% filter(Condition == cond)
  coords <- as.matrix(meta_sub[, c("x_slide_mm", "y_slide_mm")])
  ct <- meta_sub$ImmuneAtlas_ImmGen_Main_cell_Types

  nn <- nn2(coords, k = 21)
  nn_idx <- nn$nn.idx[, -1]

  # Fraction of neighbors that are immune (not "a"/tumor)
  immune_frac <- rowMeans(matrix(ct[nn_idx] != "a", nrow = nrow(nn_idx)))

  tibble(condition = cond, cell_type = ct, immune_neighbor_frac = immune_frac)
})

nn_results %>%
  ggplot(aes(x = immune_neighbor_frac, fill = condition)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~cell_type, scales = "free_y") +
  labs(x = "Fraction immune neighbors (k=20)",
       title = "Spatial immune neighborhoods: MBRT vs SBRT at 4h") +
  theme_bw()
ggsave("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/figures/spatial_nn_immune_frac.pdf",
       width = 14, height = 10)
```

- [ ] **Step 2: Run, review**
- [ ] **Step 3: Roll into org** (`:tangle no`)

- [ ] **Step 4: Commit Layer 2**

```bash
git add spatial-rads.org dev/
git commit -m "feat: Layer 2 DEGs, pathway kinetics, composition, spatial neighborhoods"
```

---

## Chunk 3: Layers 3 & 4 — Peak/Valley Analysis (BLOCKED)

> **BLOCKED:** Requires gamma-H2AX IF data from the lab. Tasks below are defined as placeholders — adapt once IF data format is known.

### Task 11: Peak/valley classification from IF

**Dev file:** `dev/03_peak_valley_classify.R`

- [ ] **Step 1: Load IF data, extract per-cell H2AX intensity** (format TBD)
- [ ] **Step 2: Spatial KDE smoothing of H2AX IF signal**
- [ ] **Step 3: Classify cells as peak/valley/boundary via thresholding**
- [ ] **Step 4: Visual validation — do stripes match expected beam geometry?**
- [ ] **Step 5: Roll into org** (`:tangle scripts/03_peak_valley_classify.R`)

---

### Task 12: Peak vs valley differential expression

**Dev file:** `dev/explore_peak_valley_de.R`

- [ ] **Step 1: DEGs between peak and valley cells (global + per cell type)**
- [ ] **Step 2: Characterize peak biology (DDR, apoptosis) and valley biology (bystander, immune priming)**
- [ ] **Step 3: Exclude DDR classification genes to avoid circularity**
- [ ] **Step 4: Roll into org** (`:tangle no`)

---

### Task 13: Peak/valley signature derivation and validation

**Dev file:** `dev/04_signatures.R`

- [ ] **Step 1: Define peak and valley gene signatures from DEGs**
- [ ] **Step 2: Validate — do signatures reconstruct stripe pattern spatially?**
- [ ] **Step 3: Score SBRT 4h with both signatures — does it resemble peaks, valleys, or neither?**
- [ ] **Step 4: Roll into org** (`:tangle scripts/04_signatures.R`)

---

### Task 14: Signature projection across timepoints (Layer 4)

**Dev file:** `dev/explore_signature_projection.R`

- [ ] **Step 1: Score all samples with peak and valley signatures (AddModuleScore)**
- [ ] **Step 2: Temporal trajectory plots — decay rates by cell type**
- [ ] **Step 3: Compare MBRT vs SBRT signature trajectories**
- [ ] **Step 4: Cross-reference with Layer 2 pathway kinetics**
- [ ] **Step 5: Roll into org** (`:tangle no`)

- [ ] **Step 6: Commit Layers 3-4**

```bash
git add spatial-rads.org dev/
git commit -m "feat: Layers 3-4 peak/valley classification, signatures, projection"
```

---

## Summary

| Chunk | Layer | Tasks | Status | Tangled scripts |
|-------|-------|-------|--------|-----------------|
| 1 | Layer 1 — QC & Landscape | 1-5 | Ready | `00_load_data.R`, `01_normalize_cluster.R` |
| 2 | Layer 2 — Comparative Kinetics | 6-10 | Ready (after L1) | `02_deg_analysis.R` |
| 3 | Layers 3-4 — Peak/Valley | 11-14 | **BLOCKED** on IF data | `03_peak_valley_classify.R`, `04_signatures.R` |

**Execution:** Chunk 1 → Chunk 2 → (wait for IF data) → Chunk 3

**After each chunk:** Review outputs with user. Findings inform whether downstream layers are worth pursuing. This is exploratory — let the data guide the next step.

**Dev → org workflow:** Each task follows `dev/*.R` → run interactively → review → roll into org-babel block (tangled or `:tangle no`). Dev files stay in repo as scratch record.
