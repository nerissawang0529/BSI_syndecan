# === Clear environment and load packages ===
rm(list = ls())

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, ComplexHeatmap, circlize, rms, scales, grid)

# === Define marker rename map (to standardized abbreviations) ===
marker_label_map <- c(
  # Endothelial
  "Ang1" = "ANG1",
  "Ang2" = "ANG2",
  "CX3CL1" = "Fractalkine",
  "Thrombomodulin" = "sThrombomodulin",
  "ESM1" = "Endocan",
  
  # Inflammation
  "CD163" = "CD163",
  "MMP8" = "MMP-8",
  "NGAL" = "NGAL",
  "Procalcitonin" = "Procalcitonin",
  "RAGE" = "RAGE",
  "TenascinC" = "Tenascin C",
  "TREM1" = "sTREM-1",
  "IL10" = "IL-10",
  "IL18" = "IL-18",
  "IL1ra" = "IL-1RA",
  "IL6" = "IL-6",
  "IL8" = "IL-8",
  
  # Coagulation
  "Tissue factor" = "Tissue Factor",
  "Ddimer" = "D-dimer",
  "Platelets" = "Platelet count",
  "PT" = "PT",
  "Antitrombin" = "Antithrombin",
  "PROC" = "Protein-C"
)

# === Load and rename raw data ===
merged_data <- read.csv("original_data/final_data.csv")
colnames(merged_data)[colnames(merged_data) == "CoagulationFactorIII"] <- "Tissue factor"
colnames(merged_data)[colnames(merged_data) == "Platelets_value_1"] <- "Platelets"
colnames(merged_data)[colnames(merged_data) == "PT_Max_24h"] <- "PT"

# === Define marker domains ===
markers_end <- c("Syndecan1", "Ang1", "Ang2", "CX3CL1", "Thrombomodulin", "ESM1")
markers_inf <- c("Syndecan1", "CD163", "MMP8", "NGAL", "Procalcitonin", "RAGE", "TenascinC", "TREM1",
                 "IL10", "IL18", "IL1ra", "IL6", "IL8")
markers_coa <- c("Syndecan1", "Tissue factor", "Ddimer", "Platelets", "PT", "Antitrombin", "PROC")

# === Log10 transform and Z-score ===
all_markers <- unique(c(markers_end, markers_inf, markers_coa))
merged_data <- merged_data %>%
  select(all_of(all_markers)) %>%
  mutate(across(everything(), ~ log10(.))) %>%
  mutate(across(everything(), ~ scale(.)[, 1]))

# === Create syndecan-1 value sequence for prediction ===
syndecan_seq <- seq(min(merged_data$Syndecan1, na.rm = TRUE),
                    max(merged_data$Syndecan1, na.rm = TRUE),
                    length.out = 100)

# === Function to calculate slope vector ===
get_slope_vector <- function(marker) {
  fit <- lm(as.formula(paste0("`", marker, "` ~ rcs(Syndecan1, 3)")), data = merged_data)
  pred <- predict(fit, newdata = data.frame(Syndecan1 = syndecan_seq))
  slope <- diff(pred) / diff(syndecan_seq)
  return(slope)
}

# === Heatmap generator function ===
make_heatmap <- function(marker_list, heatmap_name, col_fun = NULL) {
  slope_matrix <- do.call(rbind, lapply(marker_list[-1], get_slope_vector))
  rownames(slope_matrix) <- marker_label_map[marker_list[-1]]  # rename markers
  
  x_ticks <- pretty(syndecan_seq, n = 5)
  x_ticks <- x_ticks[x_ticks >= min(syndecan_seq) & x_ticks <= max(syndecan_seq)]
  x_tick_positions <- round(scales::rescale(x_ticks, from = range(syndecan_seq), to = c(1, 99)))
  col_labels <- rep("", 99)
  col_labels[x_tick_positions] <- format(x_ticks, digits = 2)
  colnames(slope_matrix) <- col_labels
  
  Heatmap(
    slope_matrix,
    name = heatmap_name,
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = FALSE,
    column_names_rot = 0,
    border = FALSE,
    row_names_gp = gpar(fontsize = 10, fontfamily = "sans"),
    column_names_gp = gpar(fontsize = 10, fontfamily = "sans"),
    heatmap_legend_param = list(
      title = "slope",
      title_gp = gpar(fontsize = 12, fontfamily = "sans", fontface = "plain"),
      labels_gp = gpar(fontsize = 10, fontfamily = "sans")
    ),
    height = unit(nrow(slope_matrix) * 0.6, "cm"),
    width = unit(30, "mm"),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.rect(x = x, y = y, width = width, height = height,
                gp = gpar(col = NA, fill = fill))
      if (i < nrow(slope_matrix)) {
        grid.lines(x = unit(c(0, 1), "npc"), y = unit(c(y - height/2, y - height/2), "npc"), gp = gpar(col = "black", lwd = 1.2))
      }
      if (j == 1) {
        grid.lines(x = unit(c(x - width/2, x - width/2), "npc"), y = unit(c(y - height/2, y + height/2), "npc"), gp = gpar(col = "black", lwd = 1.2))
      }
      if (j == ncol(slope_matrix)) {
        grid.lines(x = unit(c(x + width/2, x + width/2), "npc"), y = unit(c(y - height/2, y + height/2), "npc"), gp = gpar(col = "black", lwd = 1.2))
      }
      if (i == 1) {
        grid.lines(x = unit(c(x - width/2, x + width/2), "npc"), y = unit(c(y + height/2, y + height/2), "npc"), gp = gpar(col = "black", lwd = 1.2))
      }
      if (i == nrow(slope_matrix)) {
        grid.lines(x = unit(c(x - width/2, x + width/2), "npc"), y = unit(c(y - height/2, y - height/2), "npc"), gp = gpar(col = "black", lwd = 1.2))
      }
    }
  )
}

# === Generate and draw heatmaps ===
slope_limit <- quantile(abs(do.call(rbind, lapply(unique(c(markers_end[-1], markers_inf[-1], markers_coa[-1])), get_slope_vector))), probs = 0.95, na.rm = TRUE)
shared_col_fun <- colorRamp2(c(-slope_limit, 0, slope_limit), c("#006d2c", "#fdd9b5", "#b35806"))

ht1 <- make_heatmap(markers_end, "Endothelial", col_fun = shared_col_fun)
ht2 <- make_heatmap(markers_inf, "Inflammation", col_fun = shared_col_fun)
ht3 <- make_heatmap(markers_coa, "Coagulation", col_fun = shared_col_fun)

pdf("figures/Figure_continous_hedges_slope_zscore.pdf", width = 4, height = 12)
draw(ht1 %v% ht2 %v% ht3,
     column_title_gp = gpar(fontsize = 12, fontfamily = "sans", fontface = "plain", hjust = 0),
     padding = unit(c(0.1, 0.1, 0.3, 0.1), "cm"),
     merge_legend = TRUE)

decorate_heatmap_body("Coagulation", {
  x_ticks <- pretty(syndecan_seq, n = 5)
  x_ticks <- x_ticks[x_ticks >= min(syndecan_seq) & x_ticks <= max(syndecan_seq)]
  x_pos <- scales::rescale(x_ticks, from = range(syndecan_seq), to = c(0, 1))
  
  grid.xaxis(at = x_pos, label = format(x_ticks, digits = 2),
             gp = gpar(fontsize = 10, fontfamily = "sans"))
  grid.text("Syndecan-1 (Z-score)", y = unit(-1, "cm"),
            gp = gpar(fontsize = 12, fontfamily = "sans"))
})
dev.off()

# === Load patchwork and build spline plots with standard labels ===
pacman::p_load(patchwork)

color_map <- c(
  set_names(rep("#006d2c", length(markers_end[-1])), marker_label_map[markers_end[-1]]),
  set_names(rep("#fdd9b5", length(markers_inf[-1])), marker_label_map[markers_inf[-1]]),
  set_names(rep("#b35806", length(markers_coa[-1])), marker_label_map[markers_coa[-1]])
)

all_markers <- c(markers_end[-1], markers_inf[-1], markers_coa[-1])

full_data <- merged_data %>%
  select(Syndecan1, all_of(all_markers))

get_y_range <- function(marker, data) {
  rng <- range(data[[marker]], na.rm = TRUE)
  if (diff(rng) == 0) rng <- c(-1, 1)
  return(rng)
}

plot_single_marker_z <- function(marker, color = "black", y_lim = c(-2, 2)) {
  df <- full_data %>% select(Syndecan1, all_of(marker)) %>% drop_na()
  colnames(df) <- c("Syndecan1_z", "Marker_z")
  
  fit <- lm(Marker_z ~ rcs(Syndecan1_z, 3), data = df)
  new_df <- tibble(Syndecan1_z = seq(min(df$Syndecan1_z), max(df$Syndecan1_z), length.out = 200))
  preds <- predict(fit, newdata = new_df, se.fit = TRUE)
  new_df <- new_df %>%
    mutate(pred = preds$fit, lower = preds$fit - 1.96 * preds$se.fit, upper = preds$fit + 1.96 * preds$se.fit)
  
  ggplot(df, aes(x = Syndecan1_z, y = Marker_z)) +
    geom_point(size = 1.2, alpha = 0.5, color = "black") +
    geom_ribbon(data = new_df, aes(y = pred, ymin = lower, ymax = upper), fill = color, alpha = 0.3) +
    geom_line(data = new_df, aes(y = pred), color = color, size = 1) +
    labs(title = marker_label_map[[marker]], x = NULL, y = NULL) +
    coord_cartesian(ylim = y_lim) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(size = 12, hjust = 0, family = "sans", face = "plain"),
      axis.text = element_text(size = 10, family = "sans"),
      plot.margin = margin(5, 5, 5, 5),
      aspect.ratio = 1
    )
}

plot_list_z <- lapply(all_markers, function(marker) {
  y_range <- get_y_range(marker, full_data)
  plot_single_marker_z(marker, color = color_map[[marker_label_map[[marker]]]], y_lim = y_range)
})

ncol <- 5
n_missing <- ncol - (length(plot_list_z) %% ncol)
if (n_missing < ncol) {
  plot_list_z <- c(plot_list_z, rep(list(patchwork::plot_spacer()), n_missing))
}

final_plot_z <- wrap_plots(plot_list_z, ncol = ncol) +
  plot_annotation(
    title = "Spline fits: Syndecan-1 (Z-score) vs host response markers (Z-score)",
    subtitle = "Color-coded by domain: Endothelial (green), Inflammation (beige), Coagulation (brown)",
    theme = theme(
      plot.title = element_text(size = 12, family = "sans", hjust = 0),
      plot.subtitle = element_text(size = 11, family = "sans", hjust = 0)
    )
  )

ggsave("figures/spline_markers_zscore_lines.svg", final_plot_z, width = 15, height = 16)
