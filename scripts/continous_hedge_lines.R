pacman::p_load(tidyverse, patchwork, rms, scales)

# === Marker groups and color map ===
markers_end <- c("Ang1", "Ang2", "CX3CL1", "Thrombomodulin", "ESM1")
markers_inf <- c("CD163", "MMP8", "NGAL", "Procalcitonin", "RAGE", "TenascinC", "TREM1",
                 "IL10", "IL18", "IL1ra", "IL6", "IL8")
markers_coa <- c("Tissue factor", "Ddimer", "Platelets", "PT", "Antitrombin", "PROC")

# Create color mapping
color_map <- c(
  set_names(rep("#006d2c", length(markers_end)), markers_end),   # green: endothelial markers
  set_names(rep("#fdd9b5", length(markers_inf)), markers_inf),   # beige: inflammation markers
  set_names(rep("#b35806", length(markers_coa)), markers_coa)    # brown: coagulation markers
)

# Combine all markers
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

# === Improved y-axis label formatting ===
format_y_labels <- function(x) {
  sapply(x, function(val) {
    if (is.na(val)) return("")
    if (val >= 1e3) {
      exp <- round(log10(val), 1)
      return(parse(text = paste0("10^", exp)))
    }
    as.character(round(val, 1))
  })
}


format_x_labels <- function(x) {
  sapply(x, function(val) {
    if (is.na(val)) return("")
    if (val >= 1e3) {
      exp <- round(log10(val))
      return(parse(text = paste0("10^", exp)))
    }
    as.character(round(val, 1))
  })
}

marker_label_map <- c(
  # Endothelial
  "Ang1" = "ANG1",
  "Ang2" = "ANG2",
  "CX3CL1" = "Fractalkine",
  "Thrombomodulin" = "sThrombomodulin",
  "ESM1" = "Endocan",
  "Syndecan1" = "Syndecan-1",
  
  # Inflammation
  "CD163" = "CD163",
  "IL10" = "IL-10",
  "IL18" = "IL-18",
  "IL1ra" = "IL-1RA",
  "IL6" = "IL-6",
  "IL8" = "IL-8",
  "MMP8" = "MMP-8",
  "NGAL" = "NGAL",
  "Procalcitonin" = "Procalcitonin",
  "RAGE" = "RAGE",
  "TenascinC" = "Tenascin C",
  "TREM1" = "sTREM-1",
  
  # Coagulation
  "Tissue factor" = "Tissue Factor",
  "Ddimer" = "D-dimer",
  "Platelets" = "Platelet count",
  "PT" = "PT",
  "Antitrombin" = "Antithrombin",
  "PROC" = "Protein-C"
)






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
  
  ggplot(df, aes(x = Syndecan1, y = .data[[marker]])) +
    geom_point(size = 1.2, alpha = 0.5, color = "black") +
    geom_ribbon(data = new_df, aes(y = pred, ymin = lower, ymax = upper), 
                fill = color, alpha = 0.3) +
    geom_line(data = new_df, aes(y = pred), color = color, size = 1) +
    labs(title = marker_label_map[[marker]], x = NULL, y = NULL) +
    scale_x_continuous(
      labels = format_x_labels, 
      n.breaks = 5
    ) +
    scale_y_continuous(
      labels = format_y_labels,
      limits = c(0, max(y_range) * 1.05),
      n.breaks = 5
    ) +
    theme_classic(base_size = 9) +
    theme(
      plot.title = element_text(size = 9, hjust = 0.5, face = "plain"),
      axis.text = element_text(size = 8),
      axis.text.y = element_text(angle = 0, hjust = 1),
      plot.margin = margin(2, 2, 2, 2),
      aspect.ratio = 0.8  
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

final_plot <- wrap_plots(plot_list, ncol = ncol, byrow = TRUE) +
  plot_annotation(
    #title = "Spline fits: Syndecan-1 vs host response markers",
    #subtitle = "Color-coded by domain: Endothelial (green), Inflammation (beige), Coagulation (brown)",
    theme = theme(
      plot.title = element_text(size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5)
    )
  )


ggsave("figures/FigureS_spline_markers_loglog.svg",
       final_plot,
       width = 8.27,  
       height = 11.69, 
       units = "in",   
       limitsize = FALSE)
