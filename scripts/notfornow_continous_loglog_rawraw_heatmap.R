#continous_loglog_rawraw

pacman::p_load(tidyverse, ComplexHeatmap, circlize)

# === 计算每个 marker 的 spline 斜率 ===
spline_outputs <- lapply(all_markers, function(marker) {
  df <- full_data %>% select(Syndecan1, all_of(marker)) %>% drop_na()
  df <- df %>% filter(Syndecan1 > 1e-3) %>%
    mutate(log_syn = log10(Syndecan1), log_marker = log10(.data[[marker]] + 1))
  
  fit <- lm(log_marker ~ rms::rcs(log_syn, 3), data = df)
  
  new_df <- tibble(Syndecan1 = seq(min(df$Syndecan1), max(df$Syndecan1), length.out = 200)) %>%
    mutate(log_syn = log10(Syndecan1))
  preds <- predict(fit, newdata = new_df, se.fit = TRUE)
  
  new_df <- new_df %>%
    mutate(pred = 10^preds$fit - 1) %>%
    select(Syndecan1, pred)
  
  delta <- data.frame(
    x = diff(new_df$Syndecan1),
    y = diff(new_df$pred)
  )
  delta$slope <- delta$y / delta$x
  return(delta$slope)
})
names(spline_outputs) <- all_markers

# === 构建斜率矩阵 ===
slope_matrix <- do.call(rbind, spline_outputs)
colnames(slope_matrix) <- paste0("Pt", seq_len(ncol(slope_matrix)))
rownames(slope_matrix) <- all_markers

# === 限制斜率范围在 [-20, 30] 之间 ===
# 将斜率限制在更窄范围（例如 -2 到 5）
slope_clipped <- pmin(pmax(slope_matrix, -10), 10)

# 颜色映射范围相应调整
col_fun <- colorRamp2(c(-2, 0, 5), c("#006d2c", "#fdd9b5", "#b35806"))
# 更小的颜色映射区间，提升色彩敏感度
col_fun <- colorRamp2(c(-30, 0, 30), c("#006d2c", "#fdd9b5", "#b35806"))


# === 构建 marker 类型注释 ===
marker_class <- case_when(
  all_markers %in% markers_end ~ "Endothelial",
  all_markers %in% markers_inf ~ "Inflammation",
  all_markers %in% markers_coa ~ "Coagulation"
)
row_anno <- rowAnnotation(
  Class = marker_class,
  col = list(Class = c(
    "Endothelial" = "#006d2c",
    "Inflammation" = "#fdd9b5",
    "Coagulation" = "#b35806"
  )),
  show_annotation_name = FALSE
)

# === 创建热图对象（颜色为绿-米-棕）===
ht <- Heatmap(
  slope_clipped,
  name = "Slope",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 10),
  left_annotation = row_anno,
  column_title = "Slope heatmap (Syndecan-1 gradient, clipped -20 to 30)",
  column_title_gp = gpar(fontsize = 12, fontface = "bold")
)

# === 显示热图（RStudio 绘图窗口）===
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
grid::grid.text("Slope heatmap clipped to [-20, 30]",
                x = 0.01, y = 0.97, just = "left",
                gp = gpar(fontsize = 14, fontface = "bold"))

# === 保存 PDF ===
pdf("figures/FigureS_spline_slope_heatmap_clipped_20_30_colored.pdf", width = 10, height = 8)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
grid::grid.text("Slope heatmap clipped to [-20, 30]",
                x = 0.01, y = 0.97, just = "left",
                gp = gpar(fontsize = 14, fontface = "bold"))
dev.off()



# 生成一个包含 max, min, mean slope 的数据框
slope_summary <- data.frame(
  Marker = rownames(slope_matrix),
  Max_Slope = apply(slope_matrix, 1, max, na.rm = TRUE),
  Min_Slope = apply(slope_matrix, 1, min, na.rm = TRUE),
  Mean_Slope = apply(slope_matrix, 1, mean, na.rm = TRUE)
)

# 查看前几行
head(slope_summary)
