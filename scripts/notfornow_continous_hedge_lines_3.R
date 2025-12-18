pacman::p_load(tidyverse, patchwork, rms)

# === Marker groups and color map ===
markers_end <- c("Ang1", "Ang2", "CX3CL1", "Thrombomodulin", "ESM1")
markers_inf <- c("CD163", "MMP8", "NGAL", "Procalcitonin", "RAGE", "TenascinC", "TREM1",
                 "IL10", "IL18", "IL1ra", "IL6", "IL8")
markers_coa <- c("Tissue factor", "Ddimer", "Platelets", "PT", "Antitrombin", "PROC")

color_map <- c(
  set_names(rep("#006d2c", length(markers_end)), markers_end),   # green: endothelial markers
  set_names(rep("#fdd9b5", length(markers_inf)), markers_inf),   # beige: inflammation markers
  set_names(rep("#b35806", length(markers_coa)), markers_coa)    # brown: coagulation markers
)
all_markers <- c(markers_end, markers_inf, markers_coa)

# === Data preprocessing ===
merged_data <- read.csv("original_data/final_data.csv")
colnames(merged_data)[colnames(merged_data) == "CoagulationFactorIII"] <- "Tissue factor"
colnames(merged_data)[colnames(merged_data) == "Platelets_value_1"] <- "Platelets"
colnames(merged_data)[colnames(merged_data) == "PT_Max_24h"] <- "PT"

# === Winsorize extreme values at 99th percentile (upper limit only) ===
winsorize_upper <- function(x) {
  cap <- quantile(x, 0.99, na.rm = TRUE)
  pmin(x, cap)
}
full_data <- merged_data %>%
  select(Syndecan1, all_of(all_markers)) %>%
  mutate(across(all_of(all_markers), winsorize_upper))

# === Function to get y-axis range dynamically ===
get_y_range <- function(marker, data) {
  rng <- range(data[[marker]], na.rm = TRUE)
  if (diff(rng) == 0) rng <- c(0, 1)
  return(rng)
}

# === Single marker plot function (log-log spline fit + back-transform) ===
plot_single_marker <- function(marker, color = "black", y_lim = c(0, 7)) {
  df <- full_data %>% select(Syndecan1, all_of(marker)) %>% drop_na()
  df <- df %>% filter(Syndecan1 > 1e-3) %>%  # avoid log(0) error
    mutate(log_syn = log10(Syndecan1), log_marker = log10(.data[[marker]]))
  
  fit <- lm(log_marker ~ rcs(log_syn, 3), data = df)
  
  new_df <- tibble(Syndecan1 = seq(min(df$Syndecan1), max(df$Syndecan1), length.out = 200)) %>%
    mutate(log_syn = log10(Syndecan1))
  preds <- predict(fit, newdata = new_df, se.fit = TRUE)
  
  new_df <- new_df %>%
    mutate(pred = 10^preds$fit - 1,
           lower = 10^(preds$fit - 1.96 * preds$se.fit) - 1,
           upper = 10^(preds$fit + 1.96 * preds$se.fit) - 1)
  
  ggplot(df, aes(x = Syndecan1, y = .data[[marker]])) +
    geom_point(size = 1.2, alpha = 0.5, color = "black") +
    geom_ribbon(data = new_df, aes(y = pred, ymin = lower, ymax = upper), fill = color, alpha = 0.3) +
    geom_line(data = new_df, aes(y = pred), color = color, size = 1) +
    labs(title = marker, x = NULL, y = NULL) +
    coord_cartesian(ylim = y_lim) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(size = 12, hjust = 0, family = "sans", face = "plain"),
      axis.text = element_text(size = 10, family = "sans"),
      plot.margin = margin(5, 5, 5, 5),
      aspect.ratio = 1
    )
}

# === Batch plotting ===
plot_list <- lapply(all_markers, function(marker) {
  y_range <- get_y_range(marker, full_data)
  plot_single_marker(marker, color = color_map[[marker]], y_lim = y_range)
})

# === Add blank plots to make the layout evenly divisible by 4 columns ===
ncol <- 4
n_missing <- ncol - (length(plot_list) %% ncol)
if (n_missing < ncol) {
  plot_list <- c(plot_list, rep(list(patchwork::plot_spacer()), n_missing))
}

# === Combine all plots and export as SVG ===
final_plot <- wrap_plots(plot_list, ncol = ncol) +
  plot_annotation(
    title = "Spline fits: Syndecan-1 vs host response markers (original scale, log-log model)",
    subtitle = "Color-coded by domain: Endothelial (green), Inflammation (beige), Coagulation (brown)",
    theme = theme(
      plot.title = element_text(size = 12, family = "sans", hjust = 0),
      plot.subtitle = element_text(size = 11, family = "sans", hjust = 0)
    )
  )

ggsave("figures/FigureS_spline_markers_loglog_bestfix——1.svg", final_plot, width = 15, height = 16)
