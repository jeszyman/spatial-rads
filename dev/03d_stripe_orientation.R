library(Seurat)
library(tidyverse)
library(zoo)

obj <- readRDS("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/objects/seurat_clustered.rds")

mbrt4h <- obj@meta.data %>%
  as.data.frame() %>%
  filter(Condition == "MBRT_4h") %>%
  select(cell_id, fov, x_slide_mm, y_slide_mm, cell_type_validated,
         TypeI_interferon, TypeII_interferon, STING, DNA_Damage_Repair)

p21_expr <- GetAssayData(obj, layer = "data")["Cdkn1a", rownames(obj@meta.data[obj@meta.data$Condition == "MBRT_4h", ])]
mbrt4h$p21 <- as.numeric(p21_expr[rownames(mbrt4h)])

# --- X-axis p21 profile (stripes might be vertical in CosMx coords) ---
x_profile <- mbrt4h %>%
  mutate(x_bin = round(x_slide_mm / 0.1) * 0.1) %>%
  group_by(x_bin) %>%
  summarise(mean_p21 = mean(p21, na.rm = TRUE),
            n_cells = n(), .groups = "drop") %>%
  filter(n_cells >= 50) %>%
  arrange(x_bin) %>%
  mutate(p21_smooth = rollmean(mean_p21, k = 5, fill = NA, align = "center"))

cat("=== X-axis p21 profile (0.1mm bins) ===\n")
print(x_profile %>% filter(!is.na(p21_smooth)) %>% select(x_bin, mean_p21, p21_smooth, n_cells), n = 70)

# --- Find peaks/valleys in X profile ---
smooth_x <- x_profile %>% filter(!is.na(p21_smooth))
dx <- diff(smooth_x$p21_smooth)
x_mins <- which(dx[-length(dx)] < 0 & dx[-1] > 0) + 1
x_maxs <- which(dx[-length(dx)] > 0 & dx[-1] < 0) + 1

cat("\nX-axis local minima:\n")
print(smooth_x[x_mins, c("x_bin", "p21_smooth", "n_cells")])
cat("\nX-axis local maxima:\n")
print(smooth_x[x_maxs, c("x_bin", "p21_smooth", "n_cells")])

# --- Comparison: profiles along both axes ---
y_profile <- mbrt4h %>%
  mutate(y_bin = round(y_slide_mm / 0.1) * 0.1) %>%
  group_by(y_bin) %>%
  summarise(mean_p21 = mean(p21, na.rm = TRUE),
            n_cells = n(), .groups = "drop") %>%
  filter(n_cells >= 50) %>%
  arrange(y_bin) %>%
  mutate(p21_smooth = rollmean(mean_p21, k = 5, fill = NA, align = "center"))

# Coefficient of variation of the smoothed profiles (higher = more stripe-like variation)
y_cv <- sd(y_profile$p21_smooth, na.rm = TRUE) / mean(y_profile$p21_smooth, na.rm = TRUE)
x_cv <- sd(x_profile$p21_smooth, na.rm = TRUE) / mean(x_profile$p21_smooth, na.rm = TRUE)
cat(sprintf("\nCoefficient of variation: Y-axis = %.4f, X-axis = %.4f\n", y_cv, x_cv))
cat(sprintf("Axis with more stripe variation: %s\n", ifelse(y_cv > x_cv, "Y", "X")))

# --- Plot both axis profiles side by side ---
library(patchwork)

py <- ggplot(y_profile, aes(x = y_bin)) +
  geom_point(aes(y = mean_p21, size = n_cells), alpha = 0.2) +
  geom_line(aes(y = p21_smooth), color = "red", linewidth = 1.2, na.rm = TRUE) +
  labs(x = "y position (mm)", y = "Mean p21",
       title = sprintf("Y-axis profile (CV=%.3f)", y_cv)) +
  theme_bw()

px <- ggplot(x_profile, aes(x = x_bin)) +
  geom_point(aes(y = mean_p21, size = n_cells), alpha = 0.2) +
  geom_line(aes(y = p21_smooth), color = "blue", linewidth = 1.2, na.rm = TRUE) +
  labs(x = "x position (mm)", y = "Mean p21",
       title = sprintf("X-axis profile (CV=%.3f)", x_cv)) +
  theme_bw()

combo <- py / px
ggsave("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/figures/mbrt4h_p21_xy_profiles.pdf",
       plot = combo, width = 14, height = 10)

# --- Also try rotated axes at 15-degree increments ---
# The stripes in the image are roughly horizontal but might be slightly tilted
angles <- seq(0, 175, by = 15)
angle_cvs <- map_dbl(angles, function(theta) {
  rad <- theta * pi / 180
  # Project onto axis perpendicular to stripes (i.e., along stripe-crossing direction)
  proj <- mbrt4h$x_slide_mm * cos(rad) + mbrt4h$y_slide_mm * sin(rad)
  bin_size <- 0.1
  bins <- round(proj / bin_size) * bin_size

  profile <- tibble(bin = bins, p21 = mbrt4h$p21) %>%
    group_by(bin) %>%
    summarise(mean_p21 = mean(p21, na.rm = TRUE), n = n(), .groups = "drop") %>%
    filter(n >= 50) %>%
    arrange(bin)

  if (nrow(profile) < 7) return(0)
  smoothed <- rollmean(profile$mean_p21, k = 5, fill = NA, align = "center")
  smoothed <- smoothed[!is.na(smoothed)]
  if (length(smoothed) < 3) return(0)
  sd(smoothed) / mean(smoothed)
})

angle_df <- tibble(angle = angles, cv = angle_cvs)
cat("\n=== Stripe direction scan (CV by projection angle) ===\n")
cat("Higher CV = more stripe-like variation along that axis\n")
print(angle_df)
cat(sprintf("\nBest angle: %d degrees (CV=%.4f)\n",
            angle_df$angle[which.max(angle_df$cv)],
            max(angle_df$cv)))

p_angle <- ggplot(angle_df, aes(x = angle, y = cv)) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  geom_vline(xintercept = angle_df$angle[which.max(angle_df$cv)],
             color = "red", linetype = "dashed") +
  labs(x = "Projection angle (degrees)", y = "CV of smoothed p21",
       title = "Stripe direction scan: which axis shows most p21 variation?") +
  theme_bw()
ggsave("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/figures/mbrt4h_stripe_angle_scan.pdf",
       plot = p_angle, width = 10, height = 6)

# --- Best angle profile with peaks/valleys ---
best_angle <- angle_df$angle[which.max(angle_df$cv)]
rad <- best_angle * pi / 180
mbrt4h$proj <- mbrt4h$x_slide_mm * cos(rad) + mbrt4h$y_slide_mm * sin(rad)

best_profile <- mbrt4h %>%
  mutate(proj_bin = round(proj / 0.1) * 0.1) %>%
  group_by(proj_bin) %>%
  summarise(mean_p21 = mean(p21, na.rm = TRUE), n = n(), .groups = "drop") %>%
  filter(n >= 50) %>%
  arrange(proj_bin) %>%
  mutate(p21_smooth = rollmean(mean_p21, k = 5, fill = NA, align = "center"))

smooth_best <- best_profile %>% filter(!is.na(p21_smooth))
db <- diff(smooth_best$p21_smooth)
b_mins <- which(db[-length(db)] < 0 & db[-1] > 0) + 1
b_maxs <- which(db[-length(db)] > 0 & db[-1] < 0) + 1

cat(sprintf("\n=== Best angle (%d deg) profile peaks/valleys ===\n", best_angle))
cat("Valleys:\n")
print(smooth_best[b_mins, c("proj_bin", "p21_smooth", "n")])
cat("Peaks:\n")
print(smooth_best[b_maxs, c("proj_bin", "p21_smooth", "n")])

p_best <- ggplot(best_profile, aes(x = proj_bin)) +
  geom_point(aes(y = mean_p21, size = n), alpha = 0.2) +
  geom_line(aes(y = p21_smooth), color = "red", linewidth = 1.2, na.rm = TRUE) +
  geom_vline(xintercept = smooth_best$proj_bin[b_mins], color = "blue",
             linetype = "dashed") +
  geom_vline(xintercept = smooth_best$proj_bin[b_maxs], color = "darkgreen",
             linetype = "dashed") +
  labs(x = sprintf("Position along %d-deg axis (mm)", best_angle), y = "Mean p21",
       title = sprintf("Best stripe axis (%d deg): p21 profile with peaks(green) / valleys(blue)", best_angle)) +
  theme_bw()
ggsave("/mnt/gcs/jeszyman/projects/spatial-rads/analysis/figures/mbrt4h_p21_best_axis.pdf",
       plot = p_best, width = 14, height = 6)

cat("\nAll stripe orientation plots saved.\n")
