# ===============================
# 1. Package setup
# ===============================
packages <- c("terra", "dplyr", "mapaccuracy", "stringr", "ggplot2", "reshape2")

for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

# ===============================
# 2. Paths
# ===============================
# Input: DIR to all water change Geotiffs, desired input and output DIR - change this path to your own input and output path
raster_dir <- "C:/Users/"
val_dir    <- "C:/Users/"
out_dir    <- "C:/Users/"

dir.create(out_dir, showWarnings = FALSE)

# ===============================
# 3. Periods
# ===============================
wanted_periods <- c(
  "2017Q2_to_2017Q3",
  "2018Q2_to_2018Q3",
  "2019Q3_to_2019Q4",
  "2020Q1_to_2020Q2",
  "2021Q3_to_2021Q4",
  "2022Q2_to_2022Q3",
  "2023Q3_to_2023Q4",
  "2024Q1_to_2024Q2",
  "2025Q2_to_2025Q3"
)

raster_files <- list.files(
  raster_dir,
  pattern = paste(wanted_periods, collapse = "|"),
  full.names = TRUE
)

if (length(raster_files) == 0) stop("No rasters found — check names")

# ===============================
# 4. Load validation data
# ===============================
validation_all <- read.csv(
  file.path(val_dir, "All_combined_change_layers_validation_points.csv")
)
validation_all$ref_class <- as.character(validation_all$ref_class)
validation_all$map_class <- as.character(validation_all$map_class)
validation_all$layer     <- as.character(validation_all$layer)

# ===============================
# 5. Class label mapping
# ===============================
# Classes: 0 = Non-water, 1 = Stable water, 3 = Water gain, 4 = Water loss
# Note: class 2 does not exist in this dataset.
class_labels <- c("0" = "Non-water", "1" = "Stable water", "3" = "Water gain", "4" = "Water loss")
known_labels <- c("Non-water", "Stable water", "Water gain", "Water loss")

# ===============================
# 6. Colour palette
# ===============================
COL_PRIMARY   <- "#1f77b4"   # matplotlib tab:blue
COL_SECONDARY <- "#ff7f0e"   # matplotlib tab:orange
COL_GREY      <- "#7f7f7f"   # matplotlib tab:grey
COL_DARK      <- "#2c2c2c"   # near-black
COL_WHITE     <- "#ffffff"

# UA / PA fill colours
ua_pa_colours <- c("UA" = COL_PRIMARY, "PA" = COL_SECONDARY)

# Class colours — thematic water mapping palette
class_colours <- c(
  "Non-water"    = COL_GREY,
  "Stable water" = COL_PRIMARY,
  "Water gain"   = "#17becf",
  "Water loss"   = COL_SECONDARY
)

# ===============================
# 7. ggplot theme
# ===============================
base_theme <- theme_classic(base_size = 11) +
  theme(
    text              = element_text(color = COL_DARK),
    axis.text         = element_text(color = COL_DARK),
    axis.text.x       = element_text(angle = 45, hjust = 1),
    axis.line         = element_line(color = "#cccccc"),
    panel.background  = element_rect(fill = COL_WHITE, color = NA),
    plot.background   = element_rect(fill = COL_WHITE, color = NA),
    panel.grid.major  = element_line(color = "#e5e5e5"),
    panel.grid.minor  = element_blank(),
    legend.background = element_rect(fill = COL_WHITE, color = NA),
    legend.position   = "bottom",
    strip.background  = element_rect(fill = "#f0f4f8", color = NA),
    strip.text        = element_text(color = COL_DARK, face = "bold")
  )

# ===============================
# 8. Olofsson function
# ===============================

run_olofsson <- function(raster_path, val_df) {
  
  change_map    <- rast(raster_path)
  change_map    <- as.factor(change_map)
  
  strata        <- as.data.frame(freq(change_map))
  strata        <- strata[!is.na(strata$value), ]
  
  nh            <- strata$count
  names(nh)     <- as.character(strata$value)
  nh            <- nh[sort(names(nh))]
  
  # Derive pixel area from raster resolution
  res_m         <- res(change_map)
  pixel_area_m2 <- prod(res_m)
  area_km2      <- (nh * pixel_area_m2) / 1e6
  
  # Subset validation to classes present in nh only
  valid_classes <- names(nh)
  val_df        <- val_df[val_df$map_class %in% valid_classes &
                            val_df$ref_class %in% valid_classes, ]
  
  acc <- olofsson(
    val_df$ref_class,
    val_df$map_class,
    nh
  )
  
  # Raw sample-point count confusion matrix
  all_cls      <- sort(unique(c(val_df$map_class, val_df$ref_class)))
  count_matrix <- table(
    Reference = factor(val_df$ref_class, levels = all_cls),
    Map       = factor(val_df$map_class, levels = all_cls)
  )
  
  list(
    OA            = acc$OA,
    UA            = acc$UA,
    PA            = acc$PA,
    SEua          = acc$SEua,
    SEpa          = acc$SEpa,
    SEoa          = acc$SEoa,
    nh            = nh,
    area_km2      = area_km2,
    pixel_area_m2 = pixel_area_m2,
    ni            = table(val_df$map_class),
    matrix        = acc$matrix,
    count_matrix  = count_matrix
  )
}

# ===============================
# 9. Main loop
# ===============================

summary_results  <- list()
period_raw_store <- list()
area_results     <- list()

for (r in raster_files) {
  
  period_to     <- str_extract(basename(r), "\\d{4}Q\\d_to_\\d{4}Q\\d")
  period_csv    <- str_replace(period_to, "_to_", "_")
  layer_pattern <- paste0("change_", period_csv)
  
  message("Processing ", period_to)
  
  val_period <- validation_all %>%
    dplyr::filter(stringr::str_detect(layer, layer_pattern))
  
  if (nrow(val_period) == 0) {
    warning("No validation points for ", period_to); next
  }
  
  res <- run_olofsson(r, val_period)
  period_raw_store[[period_to]] <- res
  
  # Rename class codes to labels
  ua_named   <- res$UA;   names(ua_named)   <- class_labels[names(ua_named)]
  pa_named   <- res$PA;   names(pa_named)   <- class_labels[names(pa_named)]
  seua_named <- res$SEua; names(seua_named) <- class_labels[names(seua_named)]
  sepa_named <- res$SEpa; names(sepa_named) <- class_labels[names(sepa_named)]
  
  metrics_df <- data.frame(
    Class  = names(ua_named),
    UA     = round(ua_named   * 100, 2),
    PA     = round(pa_named   * 100, 2),
    UA_CI  = round(qnorm(0.975) * seua_named * 100, 2),
    PA_CI  = round(qnorm(0.975) * sepa_named * 100, 2)
  )
  
  write.csv(metrics_df,
            file.path(out_dir, paste0("Accuracy_metrics_", period_to, ".csv")),
            row.names = FALSE)
  
  # Olofsson proportional confusion matrix
  write.csv(round(res$matrix * 100, 2),
            file.path(out_dir, paste0("Confusion_matrix_", period_to, ".csv")))
  
  # Raw sample-count confusion matrix
  write.csv(as.data.frame.matrix(res$count_matrix),
            file.path(out_dir, paste0("Confusion_matrix_counts_", period_to, ".csv")))
  
  summary_results[[period_to]] <- data.frame(
    Period = period_to,
    OA     = round(res$OA * 100, 2),
    OA_CI  = round(qnorm(0.975) * res$SEoa * 100, 2)
  )
  
  # Surface area per class
  area_km2_named        <- res$area_km2
  names(area_km2_named) <- class_labels[names(area_km2_named)]
  area_row              <- as.data.frame(t(round(area_km2_named, 4)))
  area_row$Period       <- period_to
  area_results[[period_to]] <- area_row
}

# ===============================
# 10. Export per-period summary
# ===============================
final_summary <- dplyr::bind_rows(summary_results)
write.csv(final_summary,
          file.path(out_dir, "Olofsson_overall_accuracy_summary.csv"),
          row.names = FALSE)

# ===============================
# 11. Surface area table
# ===============================
area_table <- dplyr::bind_rows(area_results)
area_table <- area_table[, c("Period", intersect(known_labels, colnames(area_table)))]

write.csv(area_table,
          file.path(out_dir, "Surface_area_km2_per_class_per_period.csv"),
          row.names = FALSE)

message("Surface area table written.")
print(area_table)

# ===============================
# 12. Pooled / Combined accuracy across all periods
# ===============================
oa_vec   <- sapply(period_raw_store, function(x) x$OA)
seoa_vec <- sapply(period_raw_store, function(x) x$SEoa)
area_vec <- sapply(period_raw_store, function(x) sum(x$nh))

# A) Macro-average
macro_OA <- mean(oa_vec)
macro_SE <- sqrt(sum(seoa_vec^2)) / length(seoa_vec)
macro_CI <- qnorm(0.975) * macro_SE

# B) Area-weighted
w         <- area_vec / sum(area_vec)
pooled_OA <- sum(w * oa_vec)
pooled_SE <- sqrt(sum((w * seoa_vec)^2))
pooled_CI <- qnorm(0.975) * pooled_SE

combined_summary <- data.frame(
  Method   = c("Macro-average", "Area-weighted pooled"),
  OA_pct   = round(c(macro_OA, pooled_OA) * 100, 2),
  OA_CI_95 = round(c(macro_CI, pooled_CI) * 100, 2)
)

write.csv(combined_summary,
          file.path(out_dir, "Olofsson_combined_accuracy.csv"),
          row.names = FALSE)

message("Combined accuracy estimates written.")
print(combined_summary)

# ===============================
# 13. Plots
# ===============================
plots_dir <- file.path(out_dir, "Plots")
dir.create(plots_dir, showWarnings = FALSE)

# Save with white background
save_plot <- function(plot, filename, width, height) {
  ggsave(file.path(plots_dir, filename),
         plot = plot, width = width, height = height, dpi = 300, bg = COL_WHITE)
}

# Build shared class_metrics_long
class_files <- list.files(out_dir, pattern = "Accuracy_metrics_.*\\.csv$", full.names = TRUE)

class_metrics_all <- do.call(rbind, lapply(class_files, function(f) {
  df <- read.csv(f)
  df$Period <- str_extract(basename(f), "\\d{4}Q\\d_to_\\d{4}Q\\d")
  df
}))

# CSVs already store label strings — filter to the 4 known labels
class_metrics_all <- class_metrics_all[class_metrics_all$Class %in% known_labels, ]
class_metrics_all$Class_label <- factor(class_metrics_all$Class, levels = known_labels)

class_metrics_long <- reshape2::melt(
  class_metrics_all,
  id.vars      = c("Class", "Class_label", "Period"),
  measure.vars = c("UA", "PA")
)

ci_long <- reshape2::melt(
  class_metrics_all,
  id.vars      = c("Class", "Class_label", "Period"),
  measure.vars = c("UA_CI", "PA_CI")
)
ci_long$variable <- gsub("_CI", "", ci_long$variable)
names(ci_long)[names(ci_long) == "value"] <- "CI"

class_metrics_long <- merge(
  class_metrics_long, ci_long,
  by = c("Class", "Class_label", "Period", "variable")
)

# Overall Accuracy over time
oa_plot <- ggplot(final_summary, aes(x = Period, y = OA, group = 1)) +
  geom_ribbon(aes(ymin = OA - OA_CI, ymax = OA + OA_CI),
              fill = COL_PRIMARY, alpha = 0.18) +
  geom_line(color = COL_PRIMARY, linewidth = 1.2) +
  geom_point(color = COL_SECONDARY, size = 2.5) +
  ylim(0, 100) +
  base_theme +
  labs(
    title   = "Overall Accuracy per Period (Olofsson)",
    y       = "Overall Accuracy (%)",
    x       = "Period",
    caption = "Ribbon = 95% confidence interval"
  )

save_plot(oa_plot, "Overall_Accuracy_over_time.png", 9, 5)

# Class-wise UA & PA bar chart
class_plot <- ggplot(class_metrics_long,
                     aes(x = Class_label, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = value - CI, ymax = value + CI),
    position = position_dodge(0.8), width = 0.25, color = COL_DARK, linewidth = 0.4
  ) +
  scale_fill_manual(values = ua_pa_colours, name = "Metric") +
  facet_wrap(~Period, nrow = 3) +
  ylim(0, 110) +
  base_theme +
  labs(
    y       = "Accuracy (%)",
    x       = "Class",
    title   = "Class-wise User & Producer Accuracy per Period",
    caption = "Error bars = 95% confidence interval"
  )

save_plot(class_plot, "Classwise_Accuracy_per_Period.png", 14, 10)

# Box-and-whisker — UA & PA across all 9 periods, per class
box_plot <- ggplot(class_metrics_long,
                   aes(x = Class_label, y = value, fill = variable)) +
  geom_boxplot(
    position       = position_dodge(0.75),
    width          = 0.60,
    outlier.shape  = 16,            # solid circle
    outlier.size   = 2.2,
    outlier.colour = COL_DARK,
    color          = COL_DARK,
    linewidth      = 0.5
  ) +
  scale_fill_manual(values = ua_pa_colours, name = "Metric") +
  ylim(0, 110) +
  base_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(
    title   = "UA & PA Distribution across all Periods (n = 9 per class)",
    y       = "Accuracy (%)",
    x       = "Class",
    caption = paste0(
      "Centre line = median  |  Box = IQR  |  Whiskers = 1.5\u00d7IQR\n",
      "Filled dots outside whiskers = outliers"
    )
  )

save_plot(box_plot, "UA_PA_Boxplot_all_periods.png", 9, 6)

# Confusion matrices
count_conf_files <- list.files(out_dir,
                               pattern = "Confusion_matrix_counts_.*\\.csv$",
                               full.names = TRUE)

for (f in count_conf_files) {
  cm <- read.csv(f, row.names = 1, check.names = FALSE)
  
  # Strip leading "X" that read.csv adds to numeric column names
  colnames(cm) <- gsub("^X", "", colnames(cm))
  rownames(cm) <- gsub("^X", "", rownames(cm))
  
  rownames(cm) <- ifelse(rownames(cm) %in% names(class_labels),
                         class_labels[rownames(cm)], rownames(cm))
  colnames(cm) <- ifelse(colnames(cm) %in% names(class_labels),
                         class_labels[colnames(cm)], colnames(cm))
  
  # Keep only the 4 known classes
  lv <- known_labels[known_labels %in% rownames(cm)]
  cm <- cm[lv, lv, drop = FALSE]
  
  cm_long        <- reshape2::melt(as.matrix(cm))
  names(cm_long) <- c("Reference", "Map", "Count")
  cm_long$Reference <- factor(cm_long$Reference, levels = rev(lv))
  cm_long$Map       <- factor(cm_long$Map,       levels = lv)
  
  max_count <- max(cm_long$Count, na.rm = TRUE)
  period_label <- gsub("Confusion_matrix_counts_|\\.csv", "", basename(f))
  
  p <- ggplot(cm_long, aes(x = Map, y = Reference, fill = Count)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(
      aes(label = Count),
      color    = ifelse(cm_long$Count > max_count * 0.5, "white", COL_DARK),
      size     = 4,
      fontface = "bold"
    ) +
    scale_fill_gradient(low = "#d4e9f7", high = COL_PRIMARY,
                        name = "Sample\npoints") +
    base_theme +
    theme(
      axis.text.x     = element_text(angle = 30, hjust = 1),
      legend.position = "right"
    ) +
    labs(
      title    = paste0("Confusion Matrix (sample counts): ", period_label),
      subtitle = "Cell values = number of validation sample points",
      x        = "Map Class",
      y        = "Reference Class"
    )
  
  save_plot(p,
            paste0("Confusion_Matrix_counts_", period_label, ".png"),
            7, 5.5)
}

# Surface area
area_long <- reshape2::melt(
  area_table,
  id.vars       = "Period",
  variable.name = "Class",
  value.name    = "Area_km2"
)
area_long$Class <- factor(area_long$Class, levels = known_labels)

area_plot <- ggplot(area_long, aes(x = Period, y = Area_km2, fill = Class)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(values = class_colours, name = "Class") +
  base_theme +
  labs(
    title   = "Surface Area per Class per Period",
    y       = expression("Area (km"^2*")"),
    x       = "Period",
    caption = "Pixel area derived from raster resolution"
  )

save_plot(area_plot, "Surface_area_per_class_stacked.png", 10, 5.5)

# Grouped version for easier per-class comparison
area_plot_dodge <- ggplot(area_long, aes(x = Period, y = Area_km2, fill = Class)) +
  geom_bar(stat = "identity", position = position_dodge(0.75), width = 0.65) +
  scale_fill_manual(values = class_colours, name = "Class") +
  base_theme +
  labs(
    title   = "Surface Area per Class per Period (grouped)",
    y       = expression("Area (km"^2*")"),
    x       = "Period"
  )

save_plot(area_plot_dodge, "Surface_area_per_class_grouped.png", 11, 5.5)

# Combined accuracy summary
combined_plot <- ggplot(combined_summary,
                        aes(x = Method, y = OA_pct, fill = Method)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_errorbar(aes(ymin = OA_pct - OA_CI_95, ymax = OA_pct + OA_CI_95),
                width = 0.15, color = COL_DARK, linewidth = 0.7) +
  scale_fill_manual(values = c("Macro-average"        = COL_GREY,
                               "Area-weighted pooled" = COL_PRIMARY)) +
  ylim(0, 100) +
  base_theme +
  theme(legend.position = "none",
        axis.text.x     = element_text(angle = 0, hjust = 0.5)) +
  labs(
    title   = "Combined Overall Accuracy (all periods)",
    y       = "Overall Accuracy (%)",
    x       = "",
    caption = "Error bars = 95% confidence interval"
  )

save_plot(combined_plot, "Combined_Accuracy_Summary.png", 6, 5)

message("All plots saved to: ", plots_dir)
message("Olofsson validation pipeline complete.")