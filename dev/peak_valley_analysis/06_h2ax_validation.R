if (!exists("PATHWAY_COLS")) source("00_load_data.R")
has_officer <- requireNamespace("officer", quietly = TRUE) && requireNamespace("magick", quietly = TRUE)
if (has_officer) { library(officer); library(magick) }

obj <- readRDS(file.path(ANALYSIS_DIR, "objects", "seurat_clustered.rds"))
stripe_model <- readRDS(file.path(DATA_DIR, "stripe_model.rds"))

mbrt4h <- obj@meta.data %>%
  as.data.frame() %>%
  filter(Condition == "MBRT_4h") %>%
  select(cell_id, fov, x_slide_mm, y_slide_mm, cell_type_validated)

rad <- stripe_model$tilt_deg * pi / 180
mbrt4h <- mbrt4h %>% mutate(y_corr = y_slide_mm + x_slide_mm * tan(rad))

# ---- Part A: Extract H2AX image from PPTX ----
pptx_path <- file.path(INPUT_DIR, "Gamma H2AX + FOVs.pptx")
cat(sprintf("Reading PPTX: %s\n", pptx_path))
h2ax_extracted <- FALSE

if (!has_officer) {
  cat("officer/magick not installed — skipping PPTX extraction. Using zone map only.\n")
} else tryCatch({
  pptx <- read_pptx(pptx_path)
  slide_summary <- pptx_summary(pptx)
  n_slides <- if (nrow(slide_summary) > 0 && "slide_index" %in% names(slide_summary)) {
    as.integer(max(slide_summary$slide_index, na.rm = TRUE))
  } else {
    length(pptx)
  }
  cat(sprintf("PPTX has %d slides, %d content elements\n", n_slides, nrow(slide_summary)))

  # Extract images from PPTX media
  img_files <- list.files(
    file.path(dirname(pptx_path), "..", "pptx_media"),
    pattern = "\\.(png|jpg|jpeg|tif|tiff)$", full.names = TRUE, ignore.case = TRUE
  )

  # Alternative: extract directly from pptx zip
  if (length(img_files) == 0) {
    tmpdir <- tempdir()
    unzip(pptx_path, exdir = tmpdir)
    img_files <- list.files(file.path(tmpdir, "ppt", "media"),
                            pattern = "\\.(png|jpg|jpeg|tif|tiff|emf|wmf)$",
                            full.names = TRUE, ignore.case = TRUE)
    cat(sprintf("Extracted %d media files from PPTX\n", length(img_files)))
  }

  if (length(img_files) > 0) {
    # Read the largest image (likely the H2AX slide scan)
    img_sizes <- file.size(img_files)
    h2ax_path <- img_files[which.max(img_sizes)]
    cat(sprintf("Using largest image: %s (%.1f KB)\n", basename(h2ax_path), max(img_sizes) / 1024))

    h2ax_img <- image_read(h2ax_path)
    img_info <- image_info(h2ax_img)
    cat(sprintf("Image dimensions: %d x %d px\n", img_info$width, img_info$height))

    # Save extracted image
    image_write(h2ax_img, file.path(PLOT_DIR, "h2ax_extracted.png"))
    cat("H2AX image saved to plots/h2ax_extracted.png\n")
    h2ax_extracted <- TRUE
  } else {
    cat("WARNING: No images found in PPTX. Side-by-side comparison will use placeholder.\n")
    h2ax_extracted <- FALSE
  }
}, error = function(e) {
  cat(sprintf("PPTX extraction failed: %s\n", e$message))
  cat("Proceeding with zone map only.\n")
  h2ax_extracted <<- FALSE
})

# ---- Part A (cont): Side-by-side zone map ----
fov_centroids <- mbrt4h %>%
  group_by(fov) %>%
  summarise(x = mean(x_slide_mm), y = mean(y_slide_mm), .groups = "drop")

# Classify cells for zone map
mbrt4h <- mbrt4h %>%
  mutate(dist_to_peak = sapply(y_corr, function(yc) min(abs(yc - stripe_model$beam_centers)))) %>%
  mutate(zone = ifelse(dist_to_peak < stripe_model$spacing_mm / 4, "peak", "valley"))

p_zones <- ggplot(mbrt4h, aes(x = x_slide_mm, y = y_slide_mm, color = zone)) +
  geom_point(size = 0.02, alpha = 0.3) +
  scale_color_manual(values = ZONE_COLORS) +
  coord_fixed() +
  labs(title = sprintf("Stripe model: %d-deg tilt, %.2fmm spacing, %d peaks",
                       stripe_model$tilt_deg, stripe_model$spacing_mm, stripe_model$n_peaks),
       color = "Zone") +
  theme_bw(base_size = 14)
for (bc in stripe_model$beam_centers) {
  p_zones <- p_zones + geom_abline(intercept = bc, slope = -tan(rad),
                                    color = "yellow", linewidth = 0.6, alpha = 0.8)
}

# With FOV labels
p_zones_labeled <- p_zones +
  geom_text(data = fov_centroids, aes(x = x, y = y, label = fov),
            color = "black", size = 2, fontface = "bold", inherit.aes = FALSE)
ggsave(file.path(PLOT_DIR, "h2ax_zone_map_labeled.png"), plot = p_zones_labeled, width = 14, height = 10, dpi = 150)

# ---- Part B: Programmatic co-registration attempt ----
if (exists("h2ax_extracted") && h2ax_extracted) {
  tryCatch({
    cat("\nAttempting programmatic co-registration...\n")

    # Convert H2AX image to grayscale for intensity analysis
    h2ax_gray <- image_convert(h2ax_img, colorspace = "gray")
    h2ax_mat <- as.integer(image_data(h2ax_gray, channels = "gray"))
    dim(h2ax_mat) <- c(img_info$height, img_info$width)

    cat(sprintf("Grayscale matrix: %d rows x %d cols\n", nrow(h2ax_mat), ncol(h2ax_mat)))

    # Compute column-averaged intensity profile (horizontal stripes -> vertical profile)
    row_means <- rowMeans(h2ax_mat)

    # Smooth the profile
    row_smooth <- rollmean(row_means, k = max(1, length(row_means) %/% 50), fill = NA, align = "center")

    # Map pixel rows to approximate mm (using tissue extent)
    tissue_rows <- which(row_means > quantile(row_means, 0.1))
    if (length(tissue_rows) > 10) {
      tissue_top <- min(tissue_rows)
      tissue_bot <- max(tissue_rows)
      tissue_span_px <- tissue_bot - tissue_top

      # CosMx tissue y-range
      y_range_mm <- range(mbrt4h$y_slide_mm)
      mm_per_px <- diff(y_range_mm) / tissue_span_px

      intensity_profile <- tibble(
        pixel_row = seq_along(row_smooth),
        intensity = row_smooth,
        y_mm = (pixel_row - tissue_top) * mm_per_px + y_range_mm[1]
      ) %>%
        filter(!is.na(intensity), y_mm >= y_range_mm[1] - 0.5, y_mm <= y_range_mm[2] + 0.5)

      # Plot H2AX intensity profile with beam centers overlaid
      p_intensity <- ggplot(intensity_profile, aes(x = y_mm, y = intensity)) +
        geom_line(color = "brown", linewidth = 0.5) +
        geom_vline(xintercept = stripe_model$beam_centers, color = "red", linetype = "dashed", linewidth = 0.8) +
        annotate("text", x = stripe_model$beam_centers, y = max(intensity_profile$intensity, na.rm = TRUE) * 0.95,
                 label = paste0("P", seq_along(stripe_model$beam_centers)), color = "red", size = 4) +
        labs(x = "y position (mm, approximate)", y = "H2AX intensity (grayscale)",
             title = "H2AX intensity profile vs model beam centers",
             subtitle = "Red dashed = model peak centers; brown = H2AX staining intensity") +
        theme_bw()
      ggsave(file.path(PLOT_DIR, "h2ax_intensity_vs_model.png"), plot = p_intensity, width = 12, height = 5, dpi = 150)

      # Find H2AX intensity peaks
      d <- diff(row_smooth[!is.na(row_smooth)])
      local_maxs <- which(d[-length(d)] > 0 & d[-1] < 0) + 1
      h2ax_peaks_px <- (seq_along(row_smooth)[!is.na(row_smooth)])[local_maxs]
      h2ax_peaks_mm <- (h2ax_peaks_px - tissue_top) * mm_per_px + y_range_mm[1]
      h2ax_peaks_mm <- h2ax_peaks_mm[h2ax_peaks_mm >= y_range_mm[1] & h2ax_peaks_mm <= y_range_mm[2]]

      cat(sprintf("H2AX intensity peaks (mm): %s\n", paste(round(h2ax_peaks_mm, 2), collapse = ", ")))
      cat(sprintf("Model beam centers (mm):   %s\n", paste(round(stripe_model$beam_centers, 2), collapse = ", ")))

      # Compute alignment metric: for each model center, distance to nearest H2AX peak
      if (length(h2ax_peaks_mm) > 0) {
        alignment <- sapply(stripe_model$beam_centers, function(mc) {
          min(abs(mc - h2ax_peaks_mm))
        })
        cat(sprintf("Alignment (mm per center): %s\n", paste(round(alignment, 3), collapse = ", ")))
        cat(sprintf("Mean alignment error: %.3f mm (criterion: <0.2mm)\n", mean(alignment)))
        cat(sprintf("VALIDATION: %s\n", ifelse(mean(alignment) < 0.2, "PASS", "NEEDS REVIEW")))
      }
    }
  }, error = function(e) {
    cat(sprintf("Co-registration failed: %s\n", e$message))
    cat("Side-by-side visual comparison is the fallback.\n")
  })
} else {
  cat("\nH2AX image not available. Using zone map only for validation.\n")
  cat("User must visually compare h2ax_zone_map_labeled.png with PPTX.\n")
}

cat("\nH2AX validation complete. Review plots/h2ax_zone_map_labeled.png\n")
