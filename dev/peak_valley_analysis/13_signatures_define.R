if (!exists("PATHWAY_COLS")) source("00_load_data.R")

# DDR genes used for stripe model fitting — MUST exclude to avoid circularity
DDR_GENES <- c("Cdkn1a", "Gadd45a", "Gadd45b", "Mdm2", "Bax", "Pmaip1", "Bbc3",
               "Ddit3", "Atm", "Atr", "Chek1", "Chek2", "Rad51", "Brca1", "Brca2",
               "Xrcc5", "Xrcc6", "Parp1")
write_tsv(tibble(gene = DDR_GENES, reason = "DDR module — used for stripe model fitting"),
          file.path(DATA_DIR, "excluded_ddr_genes.tsv"))
cat(sprintf("DDR exclusion list: %d genes\n", length(DDR_GENES)))

# --- Load DEGs ---
bulk_degs <- read_tsv(file.path(DATA_DIR, "pv_degs_bulk.tsv"), show_col_types = FALSE)
ct_degs   <- read_tsv(file.path(DATA_DIR, "pv_degs_by_celltype.tsv"), show_col_types = FALSE)

# Verify signature genes exist in Seurat object
# Prefer local copy over GCS mount (FUSE reads are slow for 1.7G object)
local_rds <- file.path(Sys.getenv("HOME"), "data/spatial-rads/seurat_clustered.rds")
gcs_rds <- file.path(ANALYSIS_DIR, "objects", "seurat_clustered.rds")
rds_path <- if (file.exists(local_rds)) local_rds else gcs_rds
cat(sprintf("Loading Seurat object from: %s\n", rds_path))
obj <- readRDS(rds_path)
available_genes <- rownames(obj)

# --- Bulk signatures (secondary) ---
bulk_filtered <- bulk_degs %>%
  filter(!gene %in% DDR_GENES,
         abs(avg_log2FC) > 0.05,
         pct.1 > 0.1 | pct.2 > 0.1,
         gene %in% available_genes)

peak_bulk <- bulk_filtered %>% filter(avg_log2FC > 0) %>% arrange(desc(avg_log2FC))
valley_bulk <- bulk_filtered %>% filter(avg_log2FC < 0) %>% arrange(avg_log2FC)

cat(sprintf("Bulk signatures: %d peak-up, %d valley-up (after DDR exclusion)\n",
            nrow(peak_bulk), nrow(valley_bulk)))

write_tsv(peak_bulk, file.path(DATA_DIR, "peak_signature_bulk.tsv"))
write_tsv(valley_bulk, file.path(DATA_DIR, "valley_signature_bulk.tsv"))

# --- Cell-type-stratified signatures (primary) ---
ct_filtered <- ct_degs %>%
  filter(!gene %in% DDR_GENES,
         abs(avg_log2FC) > 0.10,
         pct.1 > 0.1 | pct.2 > 0.1,
         gene %in% available_genes)

# Find cell types with enough genes for a meaningful signature
ct_gene_counts <- ct_filtered %>%
  group_by(cell_type) %>%
  summarise(n_peak = sum(avg_log2FC > 0),
            n_valley = sum(avg_log2FC < 0),
            n_total = n(), .groups = "drop") %>%
  filter(n_total >= 10) %>%
  arrange(desc(n_total))

cat("\nCell types with >=10 non-DDR genes (|log2FC| > 0.10):\n")
print(as.data.frame(ct_gene_counts))

ct_sigs <- ct_filtered %>%
  filter(cell_type %in% ct_gene_counts$cell_type) %>%
  mutate(direction = ifelse(avg_log2FC > 0, "peak_up", "valley_up")) %>%
  select(gene, cell_type, direction, avg_log2FC, pct.1, pct.2) %>%
  arrange(cell_type, direction, desc(abs(avg_log2FC)))

write_tsv(ct_sigs, file.path(DATA_DIR, "celltype_signatures.tsv"))
cat(sprintf("Cell-type signatures: %d genes across %d cell types\n",
            nrow(ct_sigs), n_distinct(ct_sigs$cell_type)))

# --- Summary ---
cat("\n=== Signature Definition Summary ===\n")
cat(sprintf("DDR genes excluded: %d\n", length(DDR_GENES)))
cat(sprintf("Bulk peak-up genes: %d\n", nrow(peak_bulk)))
cat(sprintf("Bulk valley-up genes: %d\n", nrow(valley_bulk)))
cat(sprintf("Cell-type stratified: %d cell types, %d total gene-celltype pairs\n",
            n_distinct(ct_sigs$cell_type), nrow(ct_sigs)))
cat("Signature definition complete.\n")
