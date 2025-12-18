rm(list = ls())

# === 加载依赖包 ===
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, ComplexHeatmap, circlize, rms, scales, grid)

# === 数据准备 ===
merged_data <- read.csv("original_data/final_data.csv")
colnames(merged_data)[colnames(merged_data) == "CoagulationFactorIII"] <- "Tissue factor"
colnames(merged_data)[colnames(merged_data) == "Platelets_value_1"] <- "Platelets"
colnames(merged_data)[colnames(merged_data) == "PT_Max_24h"] <- "PT"

# === Marker 组定义 ===
markers_end <- c("Syndecan1", "Ang1", "Ang2", "CX3CL1", "Thrombomodulin", "ESM1")
markers_inf <- c("Syndecan1", "CD163", "MMP8", "NGAL", "Procalcitonin", "RAGE", "TenascinC", "TREM1",
                 "IL10", "IL18", "IL1ra", "IL6", "IL8")
markers_coa <- c("Syndecan1", "Tissue factor", "Ddimer", "Platelets", "PT", "Antitrombin", "PROC")

# === 合并所有 marker 并做 log(x + 1) ===
all_markers <- unique(c(markers_end, markers_inf, markers_coa))
merged_data <- merged_data %>%
  select(all_of(all_markers)) %>%
  mutate(across(all_of(all_markers[-1]), ~ log10(. + 1)))

# === 构建 Syndecan-1 拟合序列 ===
syndecan_seq <- seq(min(merged_data$Syndecan1, na.rm = TRUE),
                    max(merged_data$Syndecan1, na.rm = TRUE),
                    length.out = 100)

# === 斜率计算函数 ===
get_slope_vector <- function(marker) {
  fit <- lm(reformulate("rcs(Syndecan1, 3)", paste0("log10(`", marker, "` + 1)")), data = merged_data)
  pred_log <- predict(fit, newdata = data.frame(Syndecan1 = syndecan_seq))
  pred_raw <- 10^(pred_log) - 1
  slope <- diff(pred_raw) / diff(syndecan_seq)
  return(slope)
}

# === 计算 slope 并提取 95% 分位值范围 ===
all_marker_list <- unique(c(markers_end[-1], markers_inf[-1], markers_coa[-1]))
all_slope_matrix <- do.call(rbind, lapply(all_marker_list, get_slope_vector))
slope_limit <- quantile(abs(all_slope_matrix), probs = 0.95, na.rm = TRUE)

# === 统一配色方案 ===
shared_col_fun <- colorRamp2(
  c(-slope_limit, 0, slope_limit),
  c("#006d2c", "#fdd9b5", "#b35806")
)

# === 热图构建函数（字体样式调整）===
make_heatmap <- function(marker_list, heatmap_name, col_fun = NULL) {
  slope_matrix <- do.call(rbind, lapply(marker_list[-1], get_slope_vector))
  rownames(slope_matrix) <- marker_list[-1]
  
  # 构造列标签
  x_ticks <- pretty(syndecan_seq, n = 5)
  x_ticks <- x_ticks[x_ticks >= min(syndecan_seq) & x_ticks <= max(syndecan_seq)]
  x_tick_positions <- round(scales::rescale(x_ticks,
                                            from = range(syndecan_seq),
                                            to = c(1, 99)))
  col_labels <- rep("", 99)
  col_labels[x_tick_positions] <- format(x_ticks, big.mark = ",", scientific = FALSE)
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
        grid.lines(x = unit(c(0, 1), "npc"),
                   y = unit(c(y - height/2, y - height/2), "npc"),
                   gp = gpar(col = "black", lwd = 1.2))
      }
      if (j == 1) {
        grid.lines(x = unit(c(x - width/2, x - width/2), "npc"),
                   y = unit(c(y - height/2, y + height/2), "npc"),
                   gp = gpar(col = "black", lwd = 1.2))
      }
      if (j == ncol(slope_matrix)) {
        grid.lines(x = unit(c(x + width/2, x + width/2), "npc"),
                   y = unit(c(y - height/2, y + height/2), "npc"),
                   gp = gpar(col = "black", lwd = 1.2))
      }
      if (i == 1) {
        grid.lines(x = unit(c(x - width/2, x + width/2), "npc"),
                   y = unit(c(y + height/2, y + height/2), "npc"),
                   gp = gpar(col = "black", lwd = 1.2))
      }
      if (i == nrow(slope_matrix)) {
        grid.lines(x = unit(c(x - width/2, x + width/2), "npc"),
                   y = unit(c(y - height/2, y - height/2), "npc"),
                   gp = gpar(col = "black", lwd = 1.2))
      }
    }
  )
}

# === 构建热图 ===
ht1 <- make_heatmap(markers_end, "Endothelial", col_fun = shared_col_fun)
ht2 <- make_heatmap(markers_inf, "Inflammation", col_fun = shared_col_fun)
ht3 <- make_heatmap(markers_coa, "Coagulation", col_fun = shared_col_fun)

# === 导出 PDF（含 X轴标注与字体统一）===
pdf("figures/Figure_continous_hedgesg.pdf", width = 4, height = 12)

draw(ht1 %v% ht2 %v% ht3,
     column_title_gp = gpar(fontsize = 12, fontfamily = "sans", fontface = "plain", hjust = 0),
     padding = unit(c(0.1, 0.1, 0.3, 0.1), "cm"),
     merge_legend = TRUE)

decorate_heatmap_body("Coagulation", {
  x_ticks <- pretty(syndecan_seq, n = 5)
  x_ticks <- x_ticks[x_ticks >= min(syndecan_seq) & x_ticks <= max(syndecan_seq)]
  x_pos <- scales::rescale(x_ticks, from = range(syndecan_seq), to = c(0, 1))
  
  grid.xaxis(at = x_pos,
             label = format(x_ticks, big.mark = ",", scientific = FALSE),
             gp = gpar(fontsize = 10, fontfamily = "sans"))
  
  grid.text("Syndecan-1", y = unit(-1, "cm"),
            gp = gpar(fontsize = 12, fontfamily = "sans"))
})

dev.off()




# === 加载依赖包 ===
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, patchwork, rms)

# === Marker 列表 ===
markers_end <- c("Ang1", "Ang2", "CX3CL1", "Thrombomodulin", "ESM1")
markers_inf <- c("CD163", "MMP8", "NGAL", "Procalcitonin", "RAGE", "TenascinC", "TREM1",
                 "IL10", "IL18", "IL1ra", "IL6", "IL8")
markers_coa <- c("Tissue factor", "Ddimer", "Platelets", "PT", "Antitrombin", "PROC")

# === 配色方案 ===
color_map <- c(
  set_names(rep("#006d2c", length(markers_end)), markers_end),
  set_names(rep("#fdd9b5", length(markers_inf)), markers_inf),
  set_names(rep("#b35806", length(markers_coa)), markers_coa)
)

# === 合并 marker 列表（不含 Syndecan1）===
all_markers <- c(markers_end, markers_inf, markers_coa)

# === 数据读取与预处理 ===
merged_data <- read.csv("original_data/final_data.csv")
colnames(merged_data)[colnames(merged_data) == "CoagulationFactorIII"] <- "Tissue factor"
colnames(merged_data)[colnames(merged_data) == "Platelets_value_1"] <- "Platelets"
colnames(merged_data)[colnames(merged_data) == "PT_Max_24h"] <- "PT"

full_data <- merged_data %>%
  select(Syndecan1, all_of(all_markers)) %>%
  mutate(across(-Syndecan1, ~ log10(. + 1)))

# === 统一 y 轴范围 ===
get_y_range <- function(marker_list, data) {
  data %>%
    select(all_of(marker_list)) %>%
    pivot_longer(everything()) %>%
    drop_na() %>%
    summarise(min = min(value), max = max(value)) %>%
    unlist()
}
y_range_all <- get_y_range(all_markers, full_data)

# === 单图绘图函数 ===
plot_single_marker <- function(marker, color = "black", y_lim = c(0, 7)) {
  df <- full_data %>%
    select(Syndecan1, all_of(marker)) %>%
    drop_na()
  
  fit <- lm(reformulate("rcs(Syndecan1, 3)", response = marker), data = df)
  new_df <- tibble(Syndecan1 = seq(min(df$Syndecan1), max(df$Syndecan1), length.out = 200))
  preds <- predict(fit, newdata = new_df, se.fit = TRUE)
  new_df <- new_df %>%
    mutate(pred = preds$fit,
           lower = preds$fit - 1.96 * preds$se.fit,
           upper = preds$fit + 1.96 * preds$se.fit)
  
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

# === 批量生成图 ===
plot_list <- lapply(all_markers, function(marker) {
  plot_single_marker(marker, color = color_map[[marker]], y_lim = y_range_all)
})

# === 填充空白图，确保总数为 4 的倍数 ===
ncol <- 4
n_missing <- ncol - (length(plot_list) %% ncol)
if (n_missing < ncol) {
  plot_list <- c(plot_list, rep(list(patchwork::plot_spacer()), n_missing))
}

# === 拼接成网格 ===
final_plot <- wrap_plots(plot_list, ncol = ncol) +
  plot_annotation(
    title = "Spline fits: Syndecan-1 vs host response markers",
    subtitle = "Color-coded by domain: Endothelial (green), Inflammation (beige), Coagulation (brown)",
    theme = theme(
      plot.title = element_text(size = 12, family = "sans", hjust = 0, face = "plain"),
      plot.subtitle = element_text(size = 11, family = "sans", hjust = 0)
    )
  )

# === 显示图像于 R 中 ===
print(final_plot)

# === 导出 SVG 图像 ===
ggsave("figures/FigureS_spline_markers.svg", final_plot, width = 15, height = 16)
