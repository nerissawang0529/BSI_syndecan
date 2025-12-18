pacman::p_load(tidyverse, patchwork, rms, scales)

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

# === 统一坐标轴格式函数 ===
format_axis_exp <- function(x) {
  # 只显示为10^x格式
  parse(text = paste0("10^", round(log10(x))))
}

# === Single marker plot function ===
plot_single_marker <- function(marker, color = "black") {
  df <- full_data %>% 
    select(Syndecan1, all_of(marker)) %>% 
    drop_na() %>%
    filter(Syndecan1 > 1e-3) %>%
    mutate(log_syn = log10(Syndecan1), 
           log_marker = log10(.data[[marker]] + 1e-3))  # Avoid log(0)
  
  fit <- lm(log_marker ~ rcs(log_syn, 3), data = df)
  
  new_df <- tibble(
    Syndecan1 = seq(min(df$Syndecan1), max(df$Syndecan1), length.out = 200)
  ) %>%
    mutate(log_syn = log10(Syndecan1))
  
  preds <- predict(fit, newdata = new_df, se.fit = TRUE)
  
  new_df <- new_df %>%
    mutate(pred = 10^preds$fit - 1e-3,
           lower = 10^(preds$fit - 1.96 * preds$se.fit) - 1e-3,
           upper = 10^(preds$fit + 1.96 * preds$se.fit) - 1e-3)
  
  y_range <- range(df[[marker]], na.rm = TRUE)
  x_range <- range(df$Syndecan1, na.rm = TRUE)
  
  # 取合适的断点用于指数显示
  y_breaks <- 10^seq(floor(log10(y_range[1]+1)), ceiling(log10(y_range[2]+1)), by = 1)
  x_breaks <- 10^seq(floor(log10(x_range[1]+1)), ceiling(log10(x_range[2]+1)), by = 1)
  
  ggplot(df, aes(x = Syndecan1, y = .data[[marker]])) +
    geom_point(size = 1.2, alpha = 0.5, color = "black") +
    geom_ribbon(data = new_df, aes(y = pred, ymin = lower, ymax = upper), 
                fill = color, alpha = 0.3) +
    geom_line(data = new_df, aes(y = pred), color = color, size = 1) +
    labs(title = marker, x = NULL, y = NULL) +
    scale_y_continuous(
      trans = "log10",
      breaks = y_breaks,
      labels = format_axis_exp,
      limits = c(min(y_breaks), max(y_breaks)),
      n.breaks = 5
    ) +
    scale_x_continuous(
      trans = "log10",
      breaks = x_breaks,
      labels = format_axis_exp,
      limits = c(min(x_breaks), max(x_breaks)),
      n.breaks = 5
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5, face = "plain"),
      axis.text = element_text(size = 10),
      axis.text.y = element_text(angle = 0, hjust = 1),
      plot.margin = margin(5, 5, 5, 5),
      aspect.ratio = 1
    )
}

# === Batch plotting ===
plot_list <- lapply(all_markers, function(marker) {
  plot_single_marker(marker, color = color_map[[marker]])
})

# === Arrange and save plots ===
ncol <- 4
n_plots <- length(plot_list)
n_rows <- ceiling(n_plots / ncol)

final_plot <- wrap_plots(plot_list, ncol = ncol) +
  plot_annotation(
    title = "Spline fits: Syndecan-1 vs host response markers",
    subtitle = "Color-coded by domain: Endothelial (green), Inflammation (beige), Coagulation (brown)",
    theme = theme(
      plot.title = element_text(size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5)
    )
  )

ggsave("figures/FigureS_spline_markers_loglog.svg", 
       final_plot, 
       width = 15, 
       height = 4 * n_rows,
       limitsize = FALSE)
