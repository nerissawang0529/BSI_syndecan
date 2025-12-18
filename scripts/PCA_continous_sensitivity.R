rm(list = ls())
pacman::p_load(tidyverse, ggfortify, ggrepel, readr, rio, naniar, grDevices, scales, grid, patchwork)

# 1) 读最终数据
merged_data <- read.csv("original_data/final_data.csv")

# 2) 读没有 Antithrombin / PROC 的 ID（刚才导出的）
missing_ids <- read.csv("original_data/AT_PROC_missing_ids.csv")

# 确保 ID 类型一致
merged_data$ICU_ID_from_datasource <- as.character(merged_data$ICU_ID_from_datasource)
missing_ids$ICU_ID_from_datasource <- as.character(missing_ids$ICU_ID_from_datasource)

# 3) 删掉这些 missing IDs（只保留原始有 AT/PROC 的 80% 患者）
merged_data_cc <- merged_data %>%
  filter(!ICU_ID_from_datasource %in% missing_ids$ICU_ID_from_datasource)

cat("Complete-case sample size:", nrow(merged_data_cc), "\n")

# 4) 下面直接用你原来的 PCA 代码，只把 merged_data 换成 merged_data_cc

## ---------- Global color scale for Syndecan-1 (用全体样本算范围，和主图一致) ----------
sdc1_all  <- log10(merged_data$Syndecan1)
sdc1_all  <- sdc1_all[is.finite(sdc1_all)]
q_sdc1    <- quantile(sdc1_all, probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE)
vals_sdc1 <- scales::rescale(q_sdc1, from = range(q_sdc1))

pca_div_cols <- c(
  "#2D5F9A",
  "#B4C5E4",
  "#FFFFFF",
  "#F4B7A9",
  "#B22222"
)

# ============================================================
# --- Coagulation Activation PCA (complete-case) -------------
# ============================================================
colnames(merged_data_cc)[colnames(merged_data_cc) == "CoagulationFactorIII"] <- "Tissue factor"
colnames(merged_data_cc)[colnames(merged_data_cc) == "Platelets_value_1"]    <- "Platelets"
colnames(merged_data_cc)[colnames(merged_data_cc) == "PT_Max_24h"]           <- "PT"

markers <- c("Syndecan1", "Antitrombin", "PROC", "Tissue factor", "Ddimer", "Platelets", "PT")
markers_data <- merged_data_cc[, markers] %>%
  mutate(across(all_of(markers), ~ log10(.))) %>%
  na.omit()

endo_pca <- prcomp(markers_data %>% select(-Syndecan1), center = TRUE, scale. = TRUE)

arrow_df <- as.data.frame(endo_pca$rotation[, 1:2])
arrow_df$label <- rownames(arrow_df)
arrow_df <- arrow_df %>%
  mutate(PC1 = PC1 * 4, PC2 = PC2 * 4) %>%
  mutate(label = recode(label,
                        "Tissue factor" = "Tissue Factor",
                        "Ddimer"        = "D-dimer",
                        "Platelets"     = "Platelet count",
                        "PT"            = "PT",
                        "Antitrombin"   = "Antithrombin",
                        "PROC"          = "Protein C"))

SDC1 <- markers_data$Syndecan1
PC1  <- endo_pca$x[, 1]
PC2  <- endo_pca$x[, 2]
cor_PC1 <- cor.test(SDC1, PC1)
cor_PC2 <- cor.test(SDC1, PC2)

label_text <- paste0(
  "PC1 vs Syndecan-1: r = ", round(cor_PC1$estimate, 2),
  ", p = ", format(cor_PC1$p.value, scientific = TRUE, digits = 2), "\n",
  "PC2 vs Syndecan-1: r = ", round(cor_PC2$estimate, 2),
  ", p = ", format(cor_PC2$p.value, scientific = TRUE, digits = 2)
)

p_coa_cc <- ggplot(as.data.frame(endo_pca$x), aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = markers_data$Syndecan1), shape = 21, size = 4, colour = "black") +
  scale_fill_gradientn(
    colours = pca_div_cols,
    values  = vals_sdc1,
    limits  = range(q_sdc1),
    oob     = scales::squish,
    name    = "Syndecan-1\n(low → high)"
  ) +
  geom_segment(data = arrow_df,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.25, "cm")),
               color = "black", linewidth = 1) +
  geom_text(data = arrow_df,
            aes(x = PC1 * 1.1, y = PC2 * 1.1, label = label),
            color = "black", size = 4, family = "sans") +
  annotate("text", x = min(PC1), y = min(PC2), label = label_text,
           hjust = 0, vjust = 0, size = 4, family = "sans") +
  theme_classic(base_size = 14, base_family = "sans") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  labs(title = "Coagulation Activation (complete-case)",
       x = "PC1", y = "PC2")

p_coa_cc
ggsave(
  filename = "figures/PCA_Coagulation_complete_case.svg",
  plot = p_coa_cc,
  device = "svg",
  width = 8,
  height = 6
)



















library(tidyverse)

data <- read.csv("original_data/final_data.csv")

# 分组：和 Hedges 一样，用三个 tertile
data$Syndecan_group <- cut(
  data$Syndecan1,
  breaks = quantile(data$Syndecan1, probs = seq(0, 1, length = 4), na.rm = TRUE),
  include.lowest = TRUE,
  labels = c("1", "2", "3")
)

# log10 这两个 marker
data <- data %>%
  mutate(
    log_Antithrombin = log10(Antitrombin),
    log_PROC         = log10(PROC)
  )

# ---- 统计学检验：Kruskal–Wallis ----
kw_AT  <- kruskal.test(log_Antithrombin ~ Syndecan_group, data = data)
kw_PROC <- kruskal.test(log_PROC ~ Syndecan_group, data = data)

p_AT   <- signif(kw_AT$p.value, 3)
p_PROC <- signif(kw_PROC$p.value, 3)

cat("Kruskal–Wallis p (Antithrombin) =", p_AT, "\n")
cat("Kruskal–Wallis p (Protein C)   =", p_PROC, "\n")

# ---- 画图：把 p 值放在 subtitle 里 ----
p_AT_plot <- ggplot(
  data %>% filter(!is.na(log_Antithrombin), !is.na(Syndecan_group)),
  aes(x = Syndecan_group, y = log_Antithrombin)
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
  labs(
    title    = "Antithrombin across Syndecan-1 tertiles",
    subtitle = paste0("Kruskal–Wallis p = ", p_AT),
    x = "Syndecan-1 tertiles",
    y = "log10(Antithrombin)"
  ) +
  theme_classic(base_size = 12)

p_PROC_plot <- ggplot(
  data %>% filter(!is.na(log_PROC), !is.na(Syndecan_group)),
  aes(x = Syndecan_group, y = log_PROC)
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
  labs(
    title    = "Protein C across Syndecan-1 tertiles",
    subtitle = paste0("Kruskal–Wallis p = ", p_PROC),
    x = "Syndecan-1 tertiles",
    y = "log10(Protein C)"
  ) +
  theme_classic(base_size = 12)

# 如果你想拼在一起（需要 patchwork）
# library(patchwork)
p_AT_plot + p_PROC_plot
