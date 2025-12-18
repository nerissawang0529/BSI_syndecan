### -----------------------------------------------------------
### 0. 数据准备（与你原来一致）
### -----------------------------------------------------------
data <- read.csv("original_data/final_data.csv")

data_2 <- data[,c("Syndecan1","mortality_d30","gender","ckd",
                  "Cerebrovascular_disease","Past_myocardial_infarction",
                  "diabetes","age_yrs","SOFAtot_new")]

mort <- data_2[!is.na(data$mortality_d30), ]
mort$mortality_d30 <- as.factor(mort$mortality_d30)

library(rms)
d <- datadist(mort)
options(datadist="d")

### -----------------------------------------------------------
### 1. 三个模型：unadj / no-SOFA / SOFA
### -----------------------------------------------------------

# model 1: unadjusted
model_unadj <- lrm(mortality_d30 ~ rcs(Syndecan1, 3), data = mort)

# model 2: adjusted for age+sex+comorbidities (no SOFA)
model_noSOFA <- lrm(mortality_d30 ~ rcs(Syndecan1, 3) +
                      gender + ckd + Cerebrovascular_disease +
                      Past_myocardial_infarction + diabetes + age_yrs,
                    data = mort)

# model 3: adjusted for age+sex+comorbidities + SOFA
model_SOFA <- lrm(mortality_d30 ~ rcs(Syndecan1, 3) +
                    gender + ckd + Cerebrovascular_disease +
                    Past_myocardial_infarction + diabetes + age_yrs +
                    SOFAtot_new,
                  data = mort)


### -----------------------------------------------------------
### 2. Predict 三条曲线
### -----------------------------------------------------------

synd_median <- median(mort$Syndecan1, na.rm = TRUE)

p1 <- Predict(model_unadj, Syndecan1, ref.zero = FALSE)
p2 <- Predict(model_noSOFA, Syndecan1, ref.zero = FALSE)
p3 <- Predict(model_SOFA, Syndecan1, ref.zero = FALSE)

# Function to convert to OR relative to median
process_pred <- function(x, model_name){
  ref_y <- x$yhat[which.min(abs(x$Syndecan1 - synd_median))]
  data.frame(
    Syndecan1 = x$Syndecan1,
    yhat  = exp(x$yhat  - ref_y),
    lower = exp(x$lower - ref_y),
    upper = exp(x$upper - ref_y),
    Model = model_name
  )
}

df1 <- process_pred(p1, "Unadjusted")
df2 <- process_pred(p2, "No-SOFA adjusted")
df3 <- process_pred(p3, "Add SOFA adjusted")

df_all <- rbind(df1, df2, df3)

df_all$Model <- factor(df_all$Model,
                       levels = c("Unadjusted",
                                  "No-SOFA adjusted",
                                  "Add SOFA adjusted"))


### -----------------------------------------------------------
### 3. 三条曲线的颜色
### -----------------------------------------------------------

curve_colors <- c(
  "Unadjusted"        = "#006d2c",  # dark green
  "No-SOFA adjusted"  = "#1f78b4",  # blue
  "Add SOFA adjusted"     = "#a63603"   # brown
)


### -----------------------------------------------------------
### 4. 绘图
### -----------------------------------------------------------

library(ggplot2)

p_three <- ggplot(df_all, aes(x = Syndecan1, y = yhat,
                              color = Model, fill = Model)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.12, color = NA) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = synd_median, linetype = "dashed") +
  scale_color_manual(values = curve_colors) +
  scale_fill_manual(values = curve_colors) +
  scale_y_continuous("30-day Mortality OR (95% CI)",
                     limits = c(0, 20),
                     breaks = c(1, 2, 4, 8, 16, 20)) +
  scale_x_continuous("Syndecan-1 concentration (pg/ml)") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    axis.title = element_text(size = 12, family = "Arial"),
    axis.text = element_text(size = 10, family = "Arial"),
    legend.title = element_blank(),
    legend.text = element_text(size = 11, family = "Arial")
  )

p_three








### -----------------------------------------------------------
### 0. 数据准备
### -----------------------------------------------------------
data <- read.csv("original_data/final_data.csv")

data_2 <- data[,c("Syndecan1","mortality_d30","gender","ckd",
                  "Cerebrovascular_disease","Past_myocardial_infarction",
                  "diabetes","age_yrs","SOFAtot_new",
                  "APACHE_IV_Score")]   # <<< 新增 APACHE 总分

mort <- data_2[!is.na(data$mortality_d30), ]
mort$mortality_d30 <- as.factor(mort$mortality_d30)

library(rms)
d <- datadist(mort)
options(datadist="d")

### -----------------------------------------------------------
### 1. 五个模型
### -----------------------------------------------------------

# 1) Unadjusted
model_unadj <- lrm(mortality_d30 ~ rcs(Syndecan1, 3), data = mort)

# 2) Adjust clinical only
model_clinical <- lrm(mortality_d30 ~ rcs(Syndecan1, 3) +
                        gender + ckd + Cerebrovascular_disease +
                        Past_myocardial_infarction + diabetes + age_yrs,
                      data = mort)

# 3) Clinical + SOFA
model_SOFA <- lrm(mortality_d30 ~ rcs(Syndecan1, 3) +
                    gender + ckd + Cerebrovascular_disease +
                    Past_myocardial_infarction + diabetes + age_yrs +
                    SOFAtot_new,
                  data = mort)

# 4) Clinical + APACHE
model_APACHE <- lrm(mortality_d30 ~ rcs(Syndecan1, 3) +
                      gender + ckd + Cerebrovascular_disease +
                      Past_myocardial_infarction + diabetes + age_yrs +
                      APACHE_IV_Score,
                    data = mort)

# 5) Clinical + SOFA + APACHE
model_SOFA_APACHE <- lrm(mortality_d30 ~ rcs(Syndecan1, 3) +
                           gender + ckd + Cerebrovascular_disease +
                           Past_myocardial_infarction + diabetes + age_yrs +
                           SOFAtot_new + APACHE_IV_Score,
                         data = mort)


### -----------------------------------------------------------
### 2. Predict + OR 转换
### -----------------------------------------------------------

synd_median <- median(mort$Syndecan1, na.rm = TRUE)

p1 <- Predict(model_unadj, Syndecan1, ref.zero = FALSE)
p2 <- Predict(model_clinical, Syndecan1, ref.zero = FALSE)
p3 <- Predict(model_SOFA, Syndecan1, ref.zero = FALSE)
p4 <- Predict(model_APACHE, Syndecan1, ref.zero = FALSE)
p5 <- Predict(model_SOFA_APACHE, Syndecan1, ref.zero = FALSE)

process_pred <- function(x, model_name){
  ref_y <- x$yhat[which.min(abs(x$Syndecan1 - synd_median))]
  data.frame(
    Syndecan1 = x$Syndecan1,
    yhat  = exp(x$yhat  - ref_y),
    lower = exp(x$lower - ref_y),
    upper = exp(x$upper - ref_y),
    Model = model_name
  )
}

df1 <- process_pred(p1, "Unadjusted")
df2 <- process_pred(p2, "Clinical adjusted")
df3 <- process_pred(p3, "Clinical + SOFA")
df4 <- process_pred(p4, "Clinical + APACHE")
df5 <- process_pred(p5, "Clinical + SOFA + APACHE")

df_all <- rbind(df1, df2, df3, df4, df5)

df_all$Model <- factor(df_all$Model,
                       levels = c("Unadjusted",
                                  "Clinical adjusted",
                                  "Clinical + SOFA",
                                  "Clinical + APACHE",
                                  "Clinical + SOFA + APACHE"))


### -----------------------------------------------------------
### 3. 颜色
### -----------------------------------------------------------

### -----------------------------------------------------------
### 3. 颜色（重新定义一下，更拉开差异）
### -----------------------------------------------------------

curve_colors <- c(
  "Unadjusted"              = "#006d2c",  # dark green
  "Clinical adjusted"       = "#1f78b4",  # blue
  "Clinical + SOFA"         = "#e6550d",  # orange
  "Clinical + APACHE"       = "#756bb1",  # purple
  "Clinical + SOFA + APACHE"= "#636363"   # dark grey
)

### -----------------------------------------------------------
### 4. 绘图：只画线，不画 ribbon
### -----------------------------------------------------------

library(ggplot2)

p_five <- ggplot(df_all, aes(x = Syndecan1, y = yhat, color = Model)) +
  geom_line(size = 1.3) +                         # 只画线
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = synd_median, linetype = "dashed") +
  scale_color_manual(values = curve_colors) +
  scale_y_continuous("30-day Mortality OR (95% CI)",
                     limits = c(0, 8),
                     breaks = c(1, 2, 4, 8)) +
  scale_x_continuous("Syndecan-1 concentration (pg/ml)") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.title = element_text(size = 12, family = "Arial"),
    axis.text = element_text(size = 10, family = "Arial"),
    legend.text = element_text(size = 11, family = "Arial")
  )

p_five





## =============================
## Figure 2A: Variance explained
## Charlson_score version
## =============================

library(ggplot2)
library(dplyr)

# 📂 Load your data
data <- read.csv("original_data/final_data.csv")

# ================================
# 1. Define infection source
# ================================
data$FinalSource <- "unknown"
data$FinalSource[data$CNS_source == 1]          <- "CNS"
data$FinalSource[data$abdominal == 1]           <- "Abdominal"
data$FinalSource[data$respiratory == 1]         <- "Respiratory"
data$FinalSource[data$urinary == 1]             <- "Urinary"
data$FinalSource[data$cardiovascular == 1]      <- "Cardiovascular"
data$FinalSource[data$skin == 1]                <- "Skin"
data$FinalSource[data$other_source == 1]        <- "Other"

# ===============================================
# 2. Fit linear model with Charlson_score
# ===============================================
fit <- lm(
  data = data,
  log(Syndecan1) ~ 
    FinalGroup_pathegon +      # Pathogen type
    FinalSource +              # Infection site
    SOFAtot_new +              # Disease severity
    Charlson_score +           # ⭐ consolidated comorbidity index
    age_yrs + gender           # Demographics
)

# ===============================================
# 3. Variance decomposition
# ===============================================
var_explained <- anova(fit)
var_explained$PctExp <- (var_explained$"Sum Sq" / sum(var_explained$"Sum Sq")) * 100

var_df <- data.frame(
  Variable = rownames(var_explained),
  Variance = var_explained$PctExp
)

# ================================
# 4. Clean variable names
# ================================
var_df$Variable <- recode(var_df$Variable,
                          "SOFAtot_new"         = "SOFA Score",
                          "FinalGroup_pathegon" = "Pathogen Type",
                          "FinalSource"         = "Infection Site",
                          "Charlson_score"      = "Charlson score",
                          "age_yrs"             = "Age",
                          "gender"              = "Sex")

# ================================
# 5. Filter only relevant rows
# ================================
var_df_filtered <- var_df %>%
  filter(Variable %in% c(
    "SOFA Score",
    "Pathogen Type",
    "Infection Site",
    "Charlson score",
    "Age",
    "Sex",
    "Residuals"
  ))

# ================================
# 6. Ordering
# ================================
var_df_filtered$Variable <- factor(
  var_df_filtered$Variable,
  levels = rev(c("SOFA Score", "Pathogen Type", "Infection Site", 
                 "Charlson score", "Age", "Sex", "Residuals"))
)

# ================================
# 7. Colors (Charlson included)
# ================================
bar_colors <- c(
  "SOFA Score"     = "#9dad8f",  # 深柔和绿
  "Pathogen Type"  = "#cc844f",  # 深暖橙
  "Infection Site" = "#f2c89a",  # 浅橙
  "Charlson score" = "#a6c0b8",  # 蓝绿色（合并症）
  "Age"            = "#bcbcbc",  # 灰
  "Sex"            = "#c7c7c7",  # 浅灰
  "Residuals"      = "#b0b0b0"   # 深灰
)

var_df_filtered$Color <- bar_colors[as.character(var_df_filtered$Variable)]

# ================================
# 8. Plot
# ================================
txt_family <- "sans"
txt_size   <- 10
txt_color  <- "grey20"

p_varplot <- ggplot(var_df_filtered, aes(x = Variance, y = Variable, fill = Color)) +
  geom_col(width = 0.85) +
  scale_fill_identity() +
  scale_y_discrete(expand = c(0.01, 0.01)) +
  labs(x = "Variance Explained (%)", y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    text = element_text(family = txt_family),
    axis.text.y  = element_text(size = txt_size, family = txt_family, colour = txt_color),
    axis.text.x  = element_text(size = txt_size, family = txt_family, colour = txt_color),
    axis.title.x = element_text(size = 10, margin = margin(t = 8)),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey80", linewidth = 0.4),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 10),
    plot.tag = element_text(face = "bold", size = 16)
  ) +
  geom_text(
    aes(label = round(Variance, 2)),
    hjust = -0.1,
    size = txt_size / .pt,
    family = txt_family,
    colour = txt_color
  ) +
  labs(tag = "A")

p_varplot

# Save if needed
ggsave("figures/Figure2A_Charlson.svg", p_varplot, width = 6, height = 4)
ggsave("figures/Figure2A_Charlson.pdf", p_varplot, width = 6, height = 4, dpi = 300)







## ===============================================
## Figure 2B: APACHE instead of SOFA
## ===============================================

fit_apache <- lm(
  data = data,
  log(Syndecan1) ~ 
    FinalGroup_pathegon +
    FinalSource +
    APACHE_IV_Score +     # ← 替换 SOFA
    Charlson_score +
    age_yrs + gender
)

var_explained_apache <- anova(fit_apache)
var_explained_apache$PctExp <- (var_explained_apache$"Sum Sq" / 
                                  sum(var_explained_apache$"Sum Sq")) * 100

var_df_apache <- data.frame(
  Variable = rownames(var_explained_apache),
  Variance = var_explained_apache$PctExp
)

var_df_apache$Variable <- recode(var_df_apache$Variable,
                                 "APACHE_IV_Score"    = "APACHE Score",
                                 "FinalGroup_pathegon" = "Pathogen Type",
                                 "FinalSource"         = "Infection Site",
                                 "Charlson_score"      = "Charlson score",
                                 "age_yrs"             = "Age",
                                 "gender"              = "Sex")

var_df_apache_filtered <- var_df_apache %>%
  filter(Variable %in% c(
    "APACHE Score", "Pathogen Type", "Infection Site",
    "Charlson score", "Age", "Sex", "Residuals"
  ))

var_df_apache_filtered$Variable <- factor(
  var_df_apache_filtered$Variable,
  levels = rev(c("APACHE Score","Pathogen Type","Infection Site",
                 "Charlson score","Age","Sex","Residuals"))
)

# same color palette except replacing SOFA color by APACHE color
bar_colors_apache <- c(
  "APACHE Score"   = "#cc844f",  # 橙色（与你原来 pathogen 类似，但保持明显）
  "Pathogen Type"  = "#cc844f",
  "Infection Site" = "#f2c89a",
  "Charlson score" = "#a6c0b8",
  "Age"            = "#bcbcbc",
  "Sex"            = "#c7c7c7",
  "Residuals"      = "#b0b0b0"
)

var_df_apache_filtered$Color <- bar_colors_apache[as.character(var_df_apache_filtered$Variable)]

p_varplot_apache <- ggplot(var_df_apache_filtered, aes(x = Variance, y = Variable, fill = Color)) +
  geom_col(width = 0.85) +
  scale_fill_identity() +
  scale_y_discrete(expand = c(0.01, 0.01)) +
  labs(x = "Variance Explained (%)", y = NULL, tag = "B") +
  theme_minimal(base_size = 14) +
  theme(
    text = element_text(family = txt_family),
    axis.text.y = element_text(size = txt_size, colour = txt_color),
    axis.text.x = element_text(size = txt_size, colour = txt_color),
    axis.title.x = element_text(size = txt_size),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey80", linewidth = 0.4)
  ) +
  geom_text(aes(label = round(Variance, 2)),
            hjust = -0.1, size = txt_size / .pt,
            family = txt_family, colour = txt_color)

p_varplot_apache

ggsave("figures/Figure2B_APACHE_Charlson.pdf", p_varplot_apache, width = 6, height = 4)
ggsave("figures/Figure2B_APACHE_Charlson.svg", p_varplot_apache, width = 6, height = 4)




## ===============================================
## Figure 2C: SOFA + APACHE
## ===============================================

fit_both <- lm(
  data = data,
  log(Syndecan1) ~ 
    FinalGroup_pathegon +
    FinalSource +
    SOFAtot_new + 
    APACHE_IV_Score +
    Charlson_score +
    age_yrs + gender
)

var_explained_both <- anova(fit_both)
var_explained_both$PctExp <- (var_explained_both$"Sum Sq" / 
                                sum(var_explained_both$"Sum Sq")) * 100

var_df_both <- data.frame(
  Variable = rownames(var_explained_both),
  Variance = var_explained_both$PctExp
)

var_df_both$Variable <- recode(var_df_both$Variable,
                               "SOFAtot_new"        = "SOFA Score",
                               "APACHE_IV_Score"    = "APACHE Score",
                               "FinalGroup_pathegon" = "Pathogen Type",
                               "FinalSource"         = "Infection Site",
                               "Charlson_score"      = "Charlson score",
                               "age_yrs"             = "Age",
                               "gender"              = "Sex")

var_df_both_filtered <- var_df_both %>%
  filter(Variable %in% c(
    "SOFA Score", "APACHE Score",
    "Pathogen Type", "Infection Site",
    "Charlson score", "Age", "Sex", "Residuals"
  ))

var_df_both_filtered$Variable <- factor(
  var_df_both_filtered$Variable,
  levels = rev(c("SOFA Score","APACHE Score",
                 "Pathogen Type","Infection Site",
                 "Charlson score","Age","Sex","Residuals"))
)

# colors with both SOFA + APACHE visible
bar_colors_both <- c(
  "SOFA Score"     = "#9dad8f",   # green
  "APACHE Score"   = "#cc844f",   # orange
  "Pathogen Type"  = "#cc844f",
  "Infection Site" = "#f2c89a",
  "Charlson score" = "#a6c0b8",
  "Age"            = "#bcbcbc",
  "Sex"            = "#c7c7c7",
  "Residuals"      = "#b0b0b0"
)

var_df_both_filtered$Color <- bar_colors_both[as.character(var_df_both_filtered$Variable)]

p_varplot_both <- ggplot(var_df_both_filtered, aes(x = Variance, y = Variable, fill = Color)) +
  geom_col(width = 0.85) +
  scale_fill_identity() +
  scale_y_discrete(expand = c(0.01, 0.01)) +
  labs(x = "Variance Explained (%)", y = NULL, tag = "C") +
  theme_minimal(base_size = 14) +
  theme(
    text = element_text(family = txt_family),
    axis.text.y = element_text(size = txt_size, colour = txt_color),
    axis.text.x = element_text(size = txt_size, colour = txt_color),
    axis.title.x = element_text(size = txt_size),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey80", linewidth = 0.4)
  ) +
  geom_text(aes(label = round(Variance, 2)),
            hjust = -0.1, size = txt_size / .pt,
            family = txt_family, colour = txt_color)

p_varplot_both

ggsave("figures/Figure2C_SOFA_APACHE_Charlson.pdf", p_varplot_both, width = 6, height = 4)
ggsave("figures/Figure2C_SOFA_APACHE_Charlson.svg", p_varplot_both, width = 6, height = 4)



