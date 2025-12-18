rm(list = ls())
#ICUA_ARDS Shock_on_admis ICUA_AKI

## install.packages("pacman")
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)
install.packages("svglite")

# Read data
data_bsi <- read.csv("original_data/final_data.csv")

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
      FinalGroup_pathegon == "Mixed_pathegon" ~ "Mixed",
      FinalGroup_pathegon == "Other_pathogens" ~ "Other",
      TRUE ~ FinalGroup_pathegon
    )
  )

# 2. Custom factor level order: alphabetical except Mix and Other at end
ordered_levels <- c(
  sort(setdiff(unique(data_bsi_filtered$FinalGroup_pathegon), c("Mixed", "Other"))),
  "Mixed",
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

# 
custom_colors <- c(
  "E.coli"        = "#d9ead3",  # 柔和鼠尾草绿 sage green
  "Enterobacter"  = "#c4dfc3",  # 稍深一点的绿色
  "Enterococcus"  = "#b7d3b0",  # 绿中带灰
  "Klebsiella"    = "#fbe9d1",  # 非常柔的米黄沙橙
  "Pseudomonas"   = "#f4d8b0",  # 淡橙色
  "S.aureus"      = "#eec4a5",  # 淡珊瑚橙
  "Streptococcus" = "#e3a587",  # 灰橙调（muted coral）
  "Mixed"           = "#d9d9d9",  # light grey
  "Other"         = "#bdbdbd"   # darker grey
)
custom_colors <- c(
  "E.coli"        = "#3cb371",  # medium sea green
  "Enterobacter"  = "#66c2a5",  # muted teal
  "Enterococcus"  = "#a6dba0",  # light green
  "Klebsiella"    = "#fdd9b5",  # light orange-beige (同 Figure 1 boxplot)
  "Pseudomonas"   = "#fdae6b",  # orange
  "S.aureus"      = "#fd8d3c",  # dark orange
  "Streptococcus" = "#e6550d",  # strong orange
  "Mix"           = "#999999",  # medium gray
  "Other"         = "#cccccc"   # light gray
)
custom_colors <- c(
  "E.coli"        = "#d9ead3",  # 柔和浅绿色
  "Enterobacter"  = "#cfe2d4",  # 更淡
  "Enterococcus"  = "#c2d9c0",  # 绿灰调
  "Klebsiella"    = "#fce5cd",  # 沙橙淡调
  "Pseudomonas"   = "#e49e62",  # 明亮橙
  "S.aureus"      = "#b35806",  # 柔橙
  "Streptococcus" = "#a63603",  # muted coral
  "Mixed"           = "#e0e0e0",  # 淡灰
  "Other"         = "#c9c9c9"   # 中灰
)
label_bacteria <- function(x) dplyr::recode(x,
                                            "E.coli"        = "<i>E. coli</i>",
                                            "Enterobacter"  = "<i>Enterobacter</i>",
                                            "Enterococcus"  = "<i>Enterococcus</i>",
                                            "Klebsiella"    = "<i>Klebsiella</i>",
                                            "Pseudomonas"   = "<i>Pseudomonas</i>",
                                            "S.aureus"      = "<i>S. aureus</i>",
                                            "Streptococcus" = "<i>Streptococcus</i>",
                                            "Mixed"         = "Mixed",
                                            "Other"         = "Other"
)

# 4. Plot
# once per machine
install.packages("ggtext")

# in your script
library(ggtext)
p_pathogenbox <- ggplot(data_bsi_filtered, aes(x = FinalGroup_pathegon, y = log10_Syndecan1, fill = FinalGroup_pathegon)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.2, size = 1.2, alpha = 0.5, color = "black") +
  geom_hline(yintercept = noninfectious_median, linetype = "dashed", color = "gray30", size = 1) +
  scale_fill_manual(values = custom_colors) +
  scale_x_discrete(labels = label_bacteria) +
  labs(x = NULL, y = expression(log[10]*"(Syndecan-1 concentration (pg/ml))")) +
  theme_classic(base_size = 14) +
  theme(
    text = element_text(family = "sans"),
    axis.text.x = element_markdown(angle = 45, hjust = 1, size = 10, family = "sans"),
    axis.text.y = element_text(size = 10, family = "sans"),
    axis.title.y = element_text(size = 14, margin = margin(r = 10), family = "sans"),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.tag = element_text(face = "bold", size = 16)
  ) +
  labs(tag = "B")
p_pathogenbox




# 📦 Load required libraries
library(ggplot2)
library(dplyr)

# 📂 Load your data
data <- read.csv("original_data/final_data.csv")

# ✅ Define infection source
data$FinalSource <- "unknown"
data$FinalSource[data$CNS_source == 1] <- "CNS"
data$FinalSource[data$abdominal == 1] <- "Abdominal"
data$FinalSource[data$respiratory == 1] <- "Respiratory"
data$FinalSource[data$urinary == 1] <- "Urinary"
data$FinalSource[data$cardiovascular == 1] <- "Cardiovascular"
data$FinalSource[data$skin == 1] <- "Skin"
data$FinalSource[data$other_source == 1] <- "Other"

# 🧮 Fit linear model with all specified covariates
fit <- lm(data = data, log(Syndecan1) ~ 
            FinalGroup_pathegon + 
            FinalSource +
            SOFAtot_new +
            age_yrs + gender + 
            ckd + Cerebrovascular_disease + 
            Past_myocardial_infarction + 
            diabetes + Immune_deficiency)

# 📊 Extract variance explained
var_explained <- anova(fit)
var_explained$PctExp <- (var_explained$"Sum Sq" / sum(var_explained$"Sum Sq")) * 100

# 🧱 Create dataframe for plotting
var_df <- data.frame(Variable = rownames(var_explained), Variance = var_explained$PctExp)

# ✍️ Rename variables for plotting
var_df$Variable <- recode(var_df$Variable,
                          "SOFAtot_new" = "SOFA Score",
                          "FinalGroup_pathegon" = "Pathogen Type",
                          "FinalSource" = "Infection Site",
                          "age_yrs" = "Age",
                          "gender" = "Sex",
                          "ckd" = "CKD",
                          "Cerebrovascular_disease" = "Cerebrovascular Disease",
                          "Past_myocardial_infarction" = "Past MI",
                          "diabetes" = "Diabetes",
                          "Immune_deficiency" = "Immune Deficiency")

# 🟡 Categorize variables for plotting
var_df$Category <- ifelse(var_df$Variable == "Residuals", "Residuals",
                          ifelse(var_df$Variable == "SOFA Score", "Disease Severity",
                                 ifelse(var_df$Variable %in% c("Pathogen Type", "Infection Site"), "Pathogen Info", 
                                        "Comorbidities")))

var_df$Category <- factor(var_df$Category,
                          levels = c("Residuals", "Disease Severity", "Comorbidities", "Pathogen Info"))


complication_vars <- c("CKD", "Cerebrovascular Disease", "Past MI", "Diabetes", "Immune Deficiency")
complication_sum <- var_df %>%
  filter(Variable %in% complication_vars) %>%
  summarise(Variance = sum(Variance)) %>%
  mutate(Variable = "Comorbidities")


var_df_filtered <- var_df %>%
  filter(Variable %in% c("SOFA Score", "Pathogen Type", "Infection Site", "Age", "Sex", "Residuals")) %>%
  bind_rows(complication_sum)

library(ggplot2)
library(dplyr)

# 
var_df_filtered$Variable <- factor(var_df_filtered$Variable, levels = rev(c(
  "SOFA Score", "Pathogen Type", "Infection Site", "Comorbidities", "Age", "Sex", "Residuals"
)))

# 
bar_colors <- c(
  "SOFA Score"     = "#b3c2a2",  # desaturated sage green / 鼠尾草绿+灰
  "Pathogen Type"  = "#d4a07a",  # desaturated sand orange / 沙橙+灰
  "Infection Site" = "#e2cfc3",  # 粉橘灰 / 低对比过渡色
  "Comorbidities"   = "#bcbcbc",  # 中性灰（冷调灰）
  "Age"            = "#a9bcd0",  # muted steel blue / 钢蓝灰
  "Gender"         = "#d6dce4",  # desaturated mist blue / 薄雾蓝
  "Residuals"      = "#f5f5f5"   # background very light grey
)

bar_colors <- c(
  "SOFA Score"     = "#228b22",  # strong green
  "Pathogen Type"  = "#ff7f0e",  # orange
  "Infection Site" = "#fdd9b5",  # light orange
  "Complication"   = "#969696",  # medium gray
  "Age"            = "#a9a9a9",  # gray
  "Gender"         = "#d9d9d9",  # light gray
  "Residuals"      = "#f5f5f5"   # very light gray
)
bar_colors <- c(
  "SOFA Score"     = "#228b22",   # 绿（同 Figure 1 曲线）
  "Pathogen Type"  = "#ff7f0e",   # 橙（同 Figure 1 曲线）
  "Infection Site" = "#fdd9b5",   # 浅橙
  "Complication"   = "#bdbdbd",   # 灰
  "Age"            = "#cccccc",   # 灰白
  "Gender"         = "#dddddd",   # 更浅灰
  "Residuals"      = "#f5f5f5"    # 背景灰
)
bar_colors <- c(
  "SOFA Score"     = "#a6dba0",  # 淡绿色，对应 Figure 2A 中 Enterococcus 等绿色菌
  "Pathogen Type"  = "#f4b183",  # 柔和橙色，对应 Klebsiella/Pseudomonas/S.aureus
  "Infection Site" = "#fdebd0",  # 非常淡的米色（淡沙橙）
  "Complication"   = "#c8c8c8",  # 中灰（统一暖调）
  "Age"            = "#dcdcdc",  # 灰白
  "Gender"         = "#e6e6e6",  # 更浅灰白
  "Residuals"      = "#f5f5f5"   # 背景灰
)
bar_colors <- c(
  "SOFA Score"     = "#b3c2a2",  # greenish tan
  "Pathogen Type"  = "#e49e62",  # orange-tan
  "Infection Site" = "#fdd9b5",  # pale peach
  "Complication"   = "#d9d9d9",  # gray
  "Age"            = "#e5e5e5",  # light gray
  "Gender"         = "#f0f0f0",  # very light gray
  "Residuals"      = "#f7f7f7"   # background
)
bar_colors <- c(
  "SOFA Score"     = "#b3c2a2",  # 绿（与Fig1一致的绿色系）
  "Pathogen Type"  = "#e49e62",  # 橙（与Fig1一致）
  "Infection Site" = "#fdd9b5",  # 浅橙
  "Comorbidities"  = "#f5f5f5",  # 中灰
  "Age"            = "#cfcfcf",  # 略深一点的灰
  "Gender"         = "#d8d8d8",  # 浅灰
  "Residuals"      = "#d0d0d0"   # **中灰**（比之前明显深）
)
bar_colors <- c(
  "SOFA Score"     = "#9dad8f",  # 深一些的柔和绿（比#b3c2a2更暗）
  "Pathogen Type"  = "#cc844f",  # 深一些的暖橙（比#e49e62更暗）
  "Infection Site" = "#f2c89a",  # 浅橙下调亮度
  "Comorbidities"  = "#e0e0e0",  # 中灰（与Mix一致）
  "Age"            = "#bcbcbc",  # 稍深灰
  "Gender"         = "#c7c7c7",  # 中浅灰
  "Residuals"      = "#b0b0b0"   # 更深的灰，稳重
)
bar_colors <- c(
  "SOFA Score"     = "#9dad8f",  # 深柔和绿
  "Pathogen Type"  = "#cc844f",  # 深暖橙
  "Infection Site" = "#f2c89a",  # 柔浅橙
  "Comorbidities"  = "#a6c0b8",  # 柔和蓝绿（与整体协调，增加冷色平衡）
  "Age"            = "#bcbcbc",  # 稍深灰
  "Sex"         = "#c7c7c7",  # 中浅灰
  "Residuals"      = "#b0b0b0"   # 稳重深灰
)

#"#f7f7f7", "#b3c2a2", "#fdd9b5", "#e49e62", "#b35806"
#"Unadjusted" = "#006d2c", "Adjusted" = "#a63603"


var_df_filtered$Color <- bar_colors[as.character(var_df_filtered$Variable)]


p_varplot <- ggplot(var_df_filtered, aes(x = Variance, y = Variable, fill = Color)) +
  geom_col(width = 0.85) +  
  scale_fill_identity() +
  scale_y_discrete(expand = c(0.01, 0.01)) +
  labs(x = "Variance Explained (%)", y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    text = element_text(family = "sans"),
    axis.text.y = element_text(size = 10, family = "sans"),
    axis.text.x = element_text(size = 10, family = "sans"),
    axis.title.x = element_text(size = 10, family = "sans", margin = margin(t = 8)),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey80", linewidth = 0.4),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 10),
    plot.tag = element_text(face = "bold", size = 16)
  ) +
  geom_text(data = var_df_filtered %>% filter(Variable == "SOFA Score"),
            aes(label = "***", x = Variance + 1.5),
            size = 5, color = "black") +
  labs(tag = "B")


p_varplot
# ---- unify text style ----
txt_family <- "sans"
txt_size   <- 10          # 与 axis.text.* 一致的 pt 大小
txt_color  <- "grey20"    # 与 axis.text.y 颜色一致

p_varplot <- p_varplot +
  geom_text(aes(label = round(Variance, 2)),
            hjust  = -0.1,
            size   = txt_size / .pt,   # pt -> ggplot size
            family = txt_family,
            color  = txt_color) +
  theme(
    axis.text.y = element_text(size = txt_size, family = txt_family, colour = txt_color),
    axis.text.x = element_text(size = txt_size, family = txt_family, colour = txt_color)
  )

# ---- combine ----
library(patchwork)
p_combined <- p_varplot / p_pathogenbox +
  plot_layout(ncol = 1, heights = c(1.2, 1)) +
  plot_annotation(tag_levels = 'A', tag_suffix = "")

p_combined

# 
ggsave("figures/Figure2_combined.svg", p_combined, width = 12, height = 8)


## =========================
## Figure 2 (A–B) horizontal combined (simplified)
## =========================
library(patchwork)
library(ggplot2)
library(ggtext)

# --- unified text settings ---
txt_family <- "sans"
txt_size   <- 10
txt_color  <- "grey20"

# --- Panel A (variance explained) ---
pA <- p_varplot +
  labs(tag = "A") +
  theme(
    axis.text.y  = element_text(size = txt_size, family = txt_family, colour = txt_color),
    axis.text.x  = element_text(size = txt_size, family = txt_family, colour = txt_color),
    axis.title.x = element_text(size = 10, margin = margin(t = 8), family = txt_family),
    plot.margin  = margin(5, 10, 5, 12)
  )

# --- Panel B (pathogen boxplot, narrower) ---
pB <- p_pathogenbox +
  labs(tag = "B") +
  theme(
    plot.margin = margin(5, 4, 5, 0),
    axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1, size = 10, family = txt_family),
    axis.text.y = element_text(size = txt_size, family = txt_family),
    axis.title.y = element_text(size = 12, family = txt_family, margin = margin(r = 10))
  ) +
  coord_cartesian(clip = "off")

# --- Combine horizontally (make B narrower) ---
p_fig2 <- pA | pB +
  plot_layout(widths = c(2.3, 1.0)) +
  plot_annotation(tag_levels = 'A', tag_suffix = "")

# --- Preview & Save ---
p_fig2
ggsave("figures/Figure2_horizontal.svg", p_fig2, width = 12, height = 4, dpi = 300)
ggsave("figures/Figure2_horizontal.pdf", p_fig2, width = 12, height = 4, dpi = 300)








#the statistic####
#Syndecan-1 concentrations were log₁₀-transformed prior to statistical analysis to account for right-skewed distribution.
library(dplyr)
library(FSA)  # for dunnTest()

# 1. Prepare the data: filter out Syndecan1 <= 0 and log-transform
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

# Order groups alphabetically, with "Mix" and "Other" last
ordered_levels <- c(
  sort(setdiff(unique(data_bsi_filtered$FinalGroup_pathegon), c("Mix", "Other"))),
  "Mix",
  "Other"
)
data_bsi_filtered$FinalGroup_pathegon <- factor(
  data_bsi_filtered$FinalGroup_pathegon,
  levels = ordered_levels
)

# 2. Kruskal-Wallis test (overall significance)
kruskal_result <- kruskal.test(log10_Syndecan1 ~ FinalGroup_pathegon, data = data_bsi_filtered)
print(kruskal_result)

# 3. Post-hoc Dunn's test with Bonferroni correction
dunn_result <- dunnTest(log10_Syndecan1 ~ FinalGroup_pathegon, data = data_bsi_filtered, method = "bonferroni")
print(dunn_result)





#the following code is not for the paper, is for check the sofa for non-infection
# Load required libraries
library(dplyr)
library(ggplot2)
library(FSA)        # for Dunn's test
library(rstatix)    # optional for tidy output

# 1. Summary statistics per group
data %>%
  group_by(Microbe_groups) %>%
  summarise(
    n = length(SOFAtot_new),
    mean = mean(SOFAtot_new, na.rm = TRUE),
    median = median(SOFAtot_new, na.rm = TRUE),
    sd = sd(SOFAtot_new, na.rm = TRUE),
    IQR = IQR(SOFAtot_new, na.rm = TRUE)
  )

# 2. Boxplot of SOFAtot by Microbe_groups
ggplot(data, aes(x = Microbe_groups, y = SOFAtot_new)) +
  geom_boxplot(fill = "skyblue") +
  theme_minimal() +
  labs(title = "SOFAtot by Microbe_groups", x = "Microbe Group", y = "SOFA Total") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 3. Kruskal-Wallis test
kruskal_result <- kruskal.test(SOFAtot_new ~ Microbe_groups, data = data)
print(kruskal_result)

# 4. Post-hoc Dunn test (if Kruskal-Wallis is significant)
dunn_result <- dunnTest(SOFAtot_new ~ Microbe_groups, data = data, method = "bonferroni")
print(dunn_result)
