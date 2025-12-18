rm(list = ls())

##install.packages("pacman")
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)
library(ggrepel)
install.packages("rstatix")
library(rstatix)
library(dplyr)
install.packages("survminer")
library(survminer)


#ready data
data <- read.csv("original_data/final_data.csv")
data$Syndecan_group <- cut(data$Syndecan1,breaks = quantile(data$Syndecan1, probs = seq(0, 1, length = 4), na.rm = TRUE),
                           include.lowest = TRUE,labels = c(1, 2, 3))

#For KM line
KM <- data[,c("time", "event","Syndecan_group")]
KM$time <- as.numeric(KM$time)
KM$event <- as.numeric(KM$event)
KM$Syndecan_group <- factor(KM$Syndecan_group, levels = c(1, 2, 3), labels = c("Group 1", "Group 2", "Group 3"))
KM$event <- ifelse(is.na(KM$event), 0, KM$event)
surv_object <- survival::survfit(Surv(time, event) ~ Syndecan_group, data = KM)

p <- ggsurvplot(
  surv_object,
  data = KM,
  pval = TRUE,
  pval.method = TRUE,
  xlab = "Time (days)",
  ylab = "Survival Probability",
  risk.table = TRUE,
  palette = c("#228b22","#377eb8", "#ff7f0e"), 
  lwd = 1, # Adjust the line width here
  ylim = c(0.7, 1),
  legend.labs = c("Lowest tertile", "Middle tertile", "Highest tertile") # Update legend labels
)

# Customize the theme to use Calibri font and set font sizes explicitly
font_size <- 14 # Set the font size for general text
axis_font_size <- 14 # Set a larger size for axis labels
axis_text_size <- 14 # Set a size for axis ticks

# Customize main plot
p$plot <- p$plot + theme(
  text = element_text(family = "Arial", size = font_size),
  axis.title = element_text(family = "Arial", size = axis_font_size),
  axis.text = element_text(family = "Arial", size = axis_text_size),
  legend.text = element_text(family = "Arial", size = font_size),
  legend.title = element_text(family = "Arial", size = font_size),
  plot.title = element_text(family = "Arial", size = font_size, hjust = 0.5)
)

# Customize risk table
p$table <- p$table + theme(
  text = element_text(family = "Arial", size = font_size),
  axis.text = element_text(family = "Arial", size = axis_text_size),
  axis.title = element_text(family = "Arial", size = axis_font_size)
)

# Combine main plot and risk table using patchwork
library(patchwork)
combined_plot <- p$plot / p$table + plot_layout(heights = c(3, 1))

# Export to SVG
ggsave("figures/mortality_groups_combined.svg", combined_plot, width = 8, height = 8)
ggsave("figures/mortality_groups_combined.pdf", combined_plot, width = 8, height = 8, dpi = 300)

print(combined_plot)
dev.off()




# Syndecan-1 concentration as the x-axis, mortality as the y-axis
data_2 <- data[,c("Syndecan1","mortality_d30","gender","ckd","Cerebrovascular_disease",
                  "Past_myocardial_infarction","diabetes","age_yrs","SOFAtot_new")]
mortality <- data_2[!is.na(data$mortality_d30), ]
class(mortality$Syndecan1)
mort <- mortality

class(mort$mortality_d30)
mort$mortality_d30 <- as.factor(mort$mortality_d30)
install.packages("rms")  # Only if not already installed
library(rms)

d <- datadist(mort) 
options(datadist="d")

unadj_syndecan <- lrm(mortality_d30 ~ rcs(Syndecan1,3), data= mort)
summary(unadj_syndecan)
anova(unadj_syndecan)

unadj_syndecan <- update(unadj_syndecan)

x1 <- Predict(unadj_syndecan, Syndecan1, ref.zero = T, fun=exp)
x1
plot(x1)


plot(x1, ylab="30-day Mortality (%)")
#ggplot(x1) werkt ook en wellicht wat vriendelijker om assen en dergelijke aan te passen.
x1$yhat
library("scales")

y_breaks = c(1,4)
x3 <- ggplot(x1) + 
  xlab("Syndecan concentration (pg/ml)") +
  ylab("30-day Mortality OR (95% Cl)") + 
  ggtitle("Unadjusted") + 
  theme_bw() +
  theme(
    aspect.ratio = 0.5/1,
    legend.position = "none", 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 16),  # Increase x-axis label size
    axis.title.y = element_text(size = 16),  # Increase y-axis label size
    axis.text.x = element_text(size = 7),   # Increase x-axis tick mark size
    axis.text.y = element_text(size = 10)    # Increase y-axis tick mark size
  ) +
  geom_hline(yintercept= 1, linetype='solid', col = 'black') +
  scale_y_continuous(breaks = c(1.0, 2.0, 4.0, 8.0, 16.0, 20.0))+
  coord_cartesian(ylim = c(0, 20))
x3


#adjusted model with SOFA
# #
adj_syndecan <- lrm(mortality_d30 ~ rcs(Syndecan1,3) + gender + ckd + Cerebrovascular_disease + Past_myocardial_infarction + diabetes + 
                      age_yrs + SOFAtot_new , data= mort)
summary(adj_syndecan) ##why after this code, the result shows warnings? https://groups.google.com/g/unmarked/c/SEQxGBVT9-k?pli=1 can ignore it.
anova(adj_syndecan)

x1 <- Predict(adj_syndecan, Syndecan1, ref.zero = T, fun=exp)
summary(x1)
# 

plot(x1)
# 
plot(x1, ylab="30-day Mortality (%)")
#ggplot(x1) werkt ook en wellicht wat vriendelijker om assen en dergelijke aan te passen.
x1$yhat
library("scales")
# 
y_breaks = c(1,4)
x2 <- ggplot(x1) + 
  xlab("Syndecan concentration (pg/ml)") +
  ylab("30-day Mortality OR (95% Cl)") + 
  ggtitle("Adjusted Model (SOFA and Selected Demographics, Comorbidities)") + 
  theme_bw() +
  theme(
    aspect.ratio = 0.5/1,
    legend.position = "none", 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 16),  # Increase x-axis label size
    axis.title.y = element_text(size = 16),  # Increase y-axis label size
    axis.text.x = element_text(size = 7),   # Increase x-axis tick mark size
    axis.text.y = element_text(size = 10)    # Increase y-axis tick mark size
  ) +
  geom_hline(yintercept= 1, linetype='solid', col = 'black') +
  geom_vline(xintercept = 2143.6, linetype = 'dashed', col = 'black')+
  scale_y_continuous(breaks = c(1.0, 2.0, 4.0, 8.0, 16.0,20.0))+
  coord_cartesian(ylim = c(0, 20))

x2

# 中位数
synd_median <- median(mort$Syndecan1, na.rm = TRUE)

# Unadjusted 和 Adjusted 模型
model_unadj <- lrm(mortality_d30 ~ rcs(Syndecan1, 3), data = mort)
model_adj   <- lrm(mortality_d30 ~ rcs(Syndecan1, 3) + gender + ckd + 
                     Cerebrovascular_disease + Past_myocardial_infarction + 
                     diabetes + age_yrs + SOFAtot_new, data = mort)

# 预测 (不使用 exp)，手动归一化 OR（ref = median）
x1 <- Predict(model_unadj, Syndecan1, ref.zero = FALSE)
x2 <- Predict(model_adj,   Syndecan1, ref.zero = FALSE)

ref_y1 <- x1$yhat[which.min(abs(x1$Syndecan1 - synd_median))]
ref_y2 <- x2$yhat[which.min(abs(x2$Syndecan1 - synd_median))]

x1$yhat_adj <- exp(x1$yhat - ref_y1)
x1$lower_adj <- exp(x1$lower - ref_y1)
x1$upper_adj <- exp(x1$upper - ref_y1)
x1$Model <- "Unadjusted"

x2$yhat_adj <- exp(x2$yhat - ref_y2)
x2$lower_adj <- exp(x2$lower - ref_y2)
x2$upper_adj <- exp(x2$upper - ref_y2)
x2$Model <- "Adjusted"

# 合并数据
df_all <- rbind(
  data.frame(Syndecan1 = x1$Syndecan1, yhat = x1$yhat_adj,
             lower = x1$lower_adj, upper = x1$upper_adj, Model = x1$Model),
  data.frame(Syndecan1 = x2$Syndecan1, yhat = x2$yhat_adj,
             lower = x2$lower_adj, upper = x2$upper_adj, Model = x2$Model)
)

# 获取 p 值
p_unadj <- anova(model_unadj)["Syndecan1", "P"]
p_adj   <- anova(model_adj)["Syndecan1", "P"]

# 自定义颜色
mortality_colors <- c("Unadjusted" = "#006d2c", "Adjusted" = "#a63603")


# Ensure Model factor level order: Unadjusted first
df_all$Model <- factor(df_all$Model, levels = c("Unadjusted", "Adjusted"))

# Plot
p_mortality <- ggplot(df_all, aes(x = Syndecan1, y = yhat, color = Model, fill = Model)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  geom_vline(xintercept = synd_median, linetype = "dashed", color = "black") +
  annotate("text", x = min(df_all$Syndecan1) + 1000, y = 17,
           label = paste0("P (Unadj) = ", signif(p_unadj, 3)),
           color = mortality_colors["Unadjusted"], size = 4, hjust = 0) +
  annotate("text", x = min(df_all$Syndecan1) + 1000, y = 15,
           label = paste0("P (Adjusted) = ", signif(p_adj, 3)),
           color = mortality_colors["Adjusted"], size = 4, hjust = 0) +
  scale_color_manual(values = mortality_colors) +
  scale_fill_manual(values = mortality_colors) +
  scale_y_continuous("30-day Mortality OR (95% CI)", 
                     limits = c(0, 20), breaks = c(1, 2, 4, 8, 16, 20)) +
  scale_x_continuous("Syndecan-1 concentration (pg/ml)") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title = element_text(size = 12, family = "Arial"),
    axis.text = element_text(size = 10, family = "Arial"),
    plot.title = element_text(size = 12, hjust = 0, family = "Arial"),
    plot.margin = margin(5, 5, 5, 5)
  ) #+
#ggtitle("Relationship between syndecan-1 levels and 30-day mortality")

p_mortality

#box plot for ARDS, shock, AKI
#ICUA_ARDS,Shock_24h,AKI_24h box plot
# 📦 
library(ggplot2)
library(dplyr)
library(ggpubr)
library(gridExtra)

# 📂 
data <- read.csv("original_data/final_data.csv")

data$Syndecan1 <- as.numeric(as.character(data$Syndecan1))

# 
data$log10_sdc1 <- log10(data$Syndecan1)

# 📦 Load required libraries
library(ggplot2)
library(ggpubr)
library(gridExtra)

# 🎨 Define custom colors for each panel

#boxplot_colors <- c("gray80", "#e5f5e0", "#fee8c8")  # neutral / green / orange
boxplot_colors <- c("#f7f7f7", "#b3c2a2", "#fdd9b5") 
#"#f7f7f7", "#b3c2a2", "#fdd9b5", "#e49e62", "#b35806"
# 📏 Set unified y-axis limits and annotation height
y_limit <- c(3, 5.2)
y_pos <- 5.1

## ====== Common title size for all panels ======
title_size <- 14   # adjust if you want (e.g., 15–16)

## If not already defined above:
mortality_colors <- c("Unadjusted" = "#006d2c", "Adjusted" = "#a63603")

## ====== PANEL A: mortality spline (consistent title size) ======
p_mortality <- ggplot(df_all, aes(x = Syndecan1, y = yhat, color = Model, fill = Model)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  geom_vline(xintercept = synd_median, linetype = "dashed", color = "black") +
  annotate("text", x = min(df_all$Syndecan1) + 1000, y = 17,
           label = paste0("P (Unadj) = ", signif(p_unadj, 3)),
           color = mortality_colors["Unadjusted"], size = 4, hjust = 0) +
  annotate("text", x = min(df_all$Syndecan1) + 1000, y = 15,
           label = paste0("P (Adjusted) = ", signif(p_adj, 3)),
           color = mortality_colors["Adjusted"], size = 4, hjust = 0) +
  scale_color_manual(values = mortality_colors) +
  scale_fill_manual(values = mortality_colors) +
  scale_y_continuous("30-day Mortality OR (95% CI)",
                     limits = c(0, 20), breaks = c(1, 2, 4, 8, 16, 20)) +
  scale_x_continuous("Syndecan-1 concentration (pg/ml)") +
  #ggtitle("Relationship between syndecan-1 levels and 30-day mortality") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title = element_text(size = 12, family = "Arial"),
    axis.text  = element_text(size = 10, family = "Arial"),
    plot.title = element_text(size = title_size, hjust = 0, family = "Arial"),
    plot.margin = margin(5, 5, 5, 5)
  )

## ====== (Optional) Panels B–D use the same title_size ======
make_boxplot <- function(var, title, signif_label = NULL) {
  ggplot(data, aes(x = factor(.data[[var]]), y = log10_sdc1)) +
    geom_boxplot(fill = NA, color = "black", width = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.4, color = "black", size = 1.2) +
    labs(x = NULL,
         y = expression(log[10]~"(Syndecan-1 concentration (pg/ml))"),
         title = title) +
    scale_x_discrete(labels = c("0" = "No", "1" = "Yes")) +
    coord_cartesian(ylim = c(3, 5.2)) +
    theme_classic(base_size = 12) +
    theme(
      axis.text    = element_text(size = 10, family = "Arial"),
      axis.title.y = element_text(size = 12, family = "Arial"),
      plot.title   = element_text(size = title_size, hjust = 0, family = "Arial"),
      plot.margin  = margin(t = 5, r = 5, b = 5, l = 20)
    ) +
    if (!is.null(signif_label)) {
      stat_pvalue_manual(
        data = data.frame(group1="0", group2="1", y.position=5.1, p.adj.signif = signif_label),
        label = "p.adj.signif", size = 5
      )
    } else {
      stat_compare_means(method = "wilcox.test", label = "p.signif",
                         comparisons = list(c("0","1")), label.y = 5.1,
                         size = 5, tip.length = 0.02)
    }
}

# Short titles per your tutor
p1 <- make_boxplot("ICUA_ARDS", "ARDS")
p2 <- make_boxplot("Shock_24h", "Shock")
p3 <- make_boxplot("AKI_24h",   "AKI", "***")

# 4. --- Combine all ---
library(ggpubr)
final_plot <- ggarrange(
  p_mortality, 
  ggarrange(p2, p3, p1, ncol = 3, labels = c("B", "C", "D")),
  nrow = 2,
  heights = c(1.1, 1),
  labels = "AUTO"  # 自动 A, B, C, D
)
final_plot
# 保存
ggsave("figures/Figure1_combined.svg", final_plot, width = 12, height = 8)


## =========================
## Figure A–D (uniform axes)
## =========================
rm(list = ls()); options(stringsAsFactors = FALSE)

# ---- packages ----
need <- c("tidyverse","rms","ggpubr","patchwork")
to_install <- need[!sapply(need, requireNamespace, quietly = TRUE)]
if (length(to_install)) install.packages(to_install)
invisible(lapply(need, library, character.only = TRUE))

# ---- data ----
dat <- read.csv("original_data/final_data.csv")

## ========= A) Mortality spline (rms::lrm + Predict) =========
mort <- dat %>%
  transmute(
    Syndecan1,
    mortality_d30 = as.factor(mortality_d30),
    gender, ckd, Cerebrovascular_disease, Past_myocardial_infarction,
    diabetes, age_yrs, SOFAtot_new
  ) %>%
  filter(!is.na(Syndecan1), !is.na(mortality_d30))

d <- datadist(mort); options(datadist = "d")

# models
model_unadj <- lrm(mortality_d30 ~ rcs(Syndecan1, 3), data = mort)
model_adj   <- lrm(mortality_d30 ~ rcs(Syndecan1, 3) + gender + ckd +
                     Cerebrovascular_disease + Past_myocardial_infarction +
                     diabetes + age_yrs + SOFAtot_new, data = mort)

# predictions → OR relative to the sample median
synd_median <- median(mort$Syndecan1, na.rm = TRUE)
x1 <- Predict(model_unadj, Syndecan1, ref.zero = FALSE)
x2 <- Predict(model_adj,   Syndecan1, ref.zero = FALSE)

ref_y1 <- x1$yhat[which.min(abs(x1$Syndecan1 - synd_median))]
ref_y2 <- x2$yhat[which.min(abs(x2$Syndecan1 - synd_median))]

x1 <- transform(x1, yhat = exp(yhat - ref_y1), lower = exp(lower - ref_y1),
                upper = exp(upper - ref_y1), Model = "Unadjusted")
x2 <- transform(x2, yhat = exp(yhat - ref_y2), lower = exp(lower - ref_y2),
                upper = exp(upper - ref_y2), Model = "Adjusted")

df_all <- rbind(
  data.frame(Syndecan1 = x1$Syndecan1, yhat = x1$yhat,
             lower = x1$lower, upper = x1$upper, Model = x1$Model),
  data.frame(Syndecan1 = x2$Syndecan1, yhat = x2$yhat,
             lower = x2$lower, upper = x2$upper, Model = x2$Model)
)
df_all$Model <- factor(df_all$Model, levels = c("Unadjusted","Adjusted"))

# p-values
p_unadj <- anova(model_unadj)["Syndecan1","P"]
p_adj   <- anova(model_adj)["Syndecan1","P"]

# ---- shared style numbers ----
title_size   <- 14
y_title_size <- 10
y_tick_size  <- 9
x_title_size <- 12
x_tick_size  <- 9
mortality_colors <- c("Unadjusted" = "#006d2c", "Adjusted" = "#a63603")

# ---- Panel A (axes styled exactly like B/C/D) ----
pA <- ggplot(df_all, aes(Syndecan1, yhat, color = Model, fill = Model)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  geom_vline(xintercept = synd_median, linetype = "dashed", color = "black") +
  annotate("text", x = min(df_all$Syndecan1) + 1000, y = 17,
           label = paste0("P (Unadj) = ", signif(p_unadj, 3)),
           color = mortality_colors["Unadjusted"], size = 4, hjust = 0) +
  annotate("text", x = min(df_all$Syndecan1) + 1000, y = 15,
           label = paste0("P (Adjusted) = ", signif(p_adj, 3)),
           color = mortality_colors["Adjusted"], size = 4, hjust = 0) +
  scale_color_manual(values = mortality_colors) +
  scale_fill_manual(values = mortality_colors) +
  scale_y_continuous("30-day Mortality OR (95% CI)",
                     limits = c(0, 20), breaks = c(1, 2, 4, 8, 16, 20)) +
  scale_x_continuous("Syndecan-1 concentration (pg/ml)") +
  labs(title = "30-day mortality") +
  theme_classic(base_size = 12) +
  theme(
    panel.grid   = element_blank(),
    legend.position = "none",
    plot.title   = element_text(size = title_size, hjust = 0, family = "Arial"),
    axis.title.y = element_text(size = y_title_size, family = "Arial"),
    axis.text.y  = element_text(size = y_tick_size,  family = "Arial"),
    axis.title.x = element_text(size = 12, family = "Arial"),   # match Panels B–D
    axis.text.x  = element_text(size = 9,  family = "Arial"),
    plot.margin  = margin(5,12,5,12)
  )

## ========= B–D) Boxplots =========
dat$Syndecan1  <- as.numeric(dat$Syndecan1)
dat$log10_sdc1 <- log10(dat$Syndecan1)

make_boxplot <- function(var, title, signif_label = NULL) {
  p <- ggplot(dat, aes(x = factor(.data[[var]]), y = log10_sdc1)) +
    geom_boxplot(fill = NA, color = "black", width = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.4, color = "black", size = 1.2) +
    labs(x = NULL,
         y = expression(log[10]~"(Syndecan-1 concentration (pg/ml))"),
         title = title) +
    scale_x_discrete(labels = c("0" = "No", "1" = "Yes")) +
    coord_cartesian(ylim = c(3, 5.2)) +
    theme_classic(base_size = 12) +
    theme(
      axis.title.y = element_text(size = y_title_size, family = "Arial"),
      axis.text.y  = element_text(size = y_tick_size,  family = "Arial"),
      axis.text.x  = element_text(size = x_tick_size,  family = "Arial"),
      plot.title   = element_text(size = title_size, hjust = 0, family = "Arial"),
      plot.margin  = margin(5,6,5,2)
    )
  if (is.null(signif_label)) {
    p + ggpubr::stat_compare_means(method = "wilcox.test", label = "p.signif",
                                   comparisons = list(c("0","1")), label.y = 5.1,
                                   size = 5, tip.length = 0.02)
  } else {
    p + ggpubr::stat_pvalue_manual(
      data = data.frame(group1="0", group2="1", y.position=5.1, p.adj.signif = signif_label),
      label = "p.adj.signif", size = 5
    )
  }
}

pB <- make_boxplot("Shock_24h", "Shock")
pC <- make_boxplot("AKI_24h",   "AKI", "***")
pD <- make_boxplot("ICUA_ARDS", "ARDS")

# balanced outer margins
pB <- pB + theme(plot.margin = margin(5,6,5,2))
pC <- pC + theme(plot.margin = margin(5,6,5,2))
pD <- pD + theme(plot.margin = margin(5,12,5,2))

## ========= Arrange A–D on one row =========
final_row <- (pA | pB | pC | pD) +
  plot_layout(widths = c(1.4, 0.9, 0.9, 0.9)) +
  plot_annotation(tag_levels = "A")   # adds A, B, C, D

final_row




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


##combine the above figures
## =========================
## Final composite: Figure 1 (A–D) + Figure 2 (E–F)
## =========================
library(patchwork)
library(ggplot2)

# --- Adjust row1 (A–D) ---
row1 <- final_row +
  plot_annotation(title = NULL) +
  theme(plot.margin = margin(5, 5, 0, 5))  # tighter bottom gap

# --- Adjust row2 (E–F) ---
row2 <- p_fig2 +
  plot_annotation(title = NULL) +
  theme(plot.margin = margin(0, 5, 5, 5))

# --- Combine vertically ---
fig_all <- row1 / row2 +
  plot_layout(heights = c(1.2, 1)) +      # Row1 slightly taller
  plot_annotation(tag_levels = 'A', tag_suffix = "") &
  theme(
    plot.tag = element_text(face = "bold", size = 14, family = "sans"),
    text = element_text(family = "sans")
  )

# --- Preview ---
fig_all

# --- Save ---
ggsave("figures/Figure1_total_combined.svg", fig_all, width = 12, height = 9, dpi = 300)

