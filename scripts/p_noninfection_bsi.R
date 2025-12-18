load("original_data/df_biomarker_BSI_second_time_from_Hessel.RData")
data <- df_biomarker_wide_erik3

#box plot
library(dplyr)
library(ggplot2)

# 1. Prepare and relabel groups
data_bsi_filtered <- data_bsi %>%
  filter(Syndecan1 > 0) %>%
  mutate(
    log10_Syndecan1 = log10(Syndecan1),
    FinalGroup_pathegon = case_when(
      FinalGroup_pathegon == "Mixed_pathegon" ~ "Mix",
      FinalGroup_pathegon == "Other_pathogens" ~ "Other",
      TRUE ~ FinalGroup_pathegon
    )
  )

# 2. Custom factor level order: alphabetical except Mix and Other at end
ordered_levels <- c(
  sort(setdiff(unique(data_bsi_filtered$FinalGroup_pathegon), c("Mix", "Other"))),
  "Mix",
  "Other"
)
data_bsi_filtered$FinalGroup_pathegon <- factor(
  data_bsi_filtered$FinalGroup_pathegon,
  levels = ordered_levels
)

# 3. Non-infectious median (no subtitle)
noninfectious_median <- data %>%
  filter(Microbe_groups == "Non-infectious", `Syndecan-1/CD138 (29)` > 0) %>%
  summarise(median_log10 = median(log10(`Syndecan-1/CD138 (29)`), na.rm = TRUE)) %>%
  pull(median_log10)
#    median       Q1       Q3
#    5750.708 3719.035 9366.429

# 使用完整的 data，确保 Non-infectious 组不为空
noninfectious_group <- data %>%
  filter(Microbe_groups == "Non-infectious", `Syndecan-1/CD138 (29)` > 0) %>%
  mutate(log10_Syndecan1 = log10(`Syndecan-1/CD138 (29)`)) %>%
  pull(log10_Syndecan1)



# 📂 Load your data
data <- read.csv("original_data/final_data.csv")

# 对 data$Syndecan1 取 log10
log10_data_syndecan1 <- log10(data$Syndecan1[data$Syndecan1 > 0])  # 避免 log(0) 错误

# 进行两组比较，使用 Wilcoxon 秩和检验
test_result <- wilcox.test(noninfectious_group, log10_data_syndecan1, exact = FALSE)

# 查看 p 值
test_result$p.value
