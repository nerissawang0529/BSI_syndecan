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














#adjusted model without SOFA#############
# #
adj_syndecan <- lrm(mortality_d30 ~ rcs(Syndecan1,3) + gender + ckd + Cerebrovascular_disease + Past_myocardial_infarction + diabetes + 
                      age_yrs, data= mort)
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
  ggtitle("Adjusted Model (Selected Demographics and Comorbidities)") + 
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


#cox model####
# For Cox proportional hazard model 
## For unadjusted_cox
###calculate
unadjusted_cox <- coxph(Surv(time, event) ~ Syndecan_group, data = KM)
summary(unadjusted_cox)
anova(unadjusted_cox)

#continuous
unadjusted_cox_continuous <- coxph(Surv(time, event) ~ log10(Syndecan1), data = data)
summary(unadjusted_cox_continuous)
anova(unadjusted_cox_continuous)

###get the result
unadjusted_results <- summary(unadjusted_cox)$coefficients
print(unadjusted_results)

#continuous
unadjusted_results_continuous <- summary(unadjusted_cox_continuous)$coefficients
print(unadjusted_results_continuous)


## For adjusted_cox
###calculate
data <- data %>%
  mutate(
    Group = factor(Syndecan_group, levels = c('1', '2', '3')),
    gender = as.character(gender),
    ckd = as.character(ckd),
    Cerebrovascular_disease = as.character(Cerebrovascular_disease),
    Past_myocardial_infarction = as.character(Past_myocardial_infarction),
    diabetes = as.character(diabetes)
  )
str(data)

adjusted_cox <- coxph(Surv(time, event) ~ gender + ckd + Cerebrovascular_disease + Past_myocardial_infarction + diabetes + age_yrs + SOFAtot + Syndecan_group, data = data)
summary(adjusted_cox)
anova(adjusted_cox)

#continuous
adjusted_cox_continuous <- coxph(Surv(time, event) ~ gender + ckd + Cerebrovascular_disease + Past_myocardial_infarction + diabetes + age_yrs + SOFAtot 
                                 + Syndecan1, data = data)
summary(adjusted_cox_continuous)
anova(adjusted_cox_continuous)


###get the result
adjusted_results <- summary(adjusted_cox)$coefficients
print(adjusted_results)

#continuous
adjusted_results_continuous <- summary(adjusted_cox_continuous)$coefficients
print(adjusted_results_continuous)



#not adjust the severity
adjusted_cox_no_severity <- coxph(Surv(time, event) ~ gender + ckd + Cerebrovascular_disease + Past_myocardial_infarction + diabetes + age_yrs + Syndecan_group, data = data)
summary(adjusted_cox_no_severity)
anova(adjusted_cox_no_severity)

#continuous
adjusted_cox_no_severity_continuous <- coxph(Surv(time, event) ~ gender + ckd + Cerebrovascular_disease + Past_myocardial_infarction + diabetes + 
                                               age_yrs + Syndecan1, data = data)
summary(adjusted_cox_no_severity_continuous)
anova(adjusted_cox_no_severity_continuous)


###get the result
adjusted_results_no_severity <- summary(adjusted_cox_no_severity)$coefficients
print(adjusted_results_no_severity)
anova(adjusted_results_no_severity)

#continuous
adjusted_cox_no_severity_continuous <- summary(adjusted_cox_no_severity_continuous)$coefficients
print(adjusted_cox_no_severity_continuous)



##test assumption
#####Proportional Hazards Assumption
# 1. Test Proportional Hazards Assumption
ph_test <- cox.zph(adjusted_cox)
print(ph_test)
#this part is the most important. to check if the values are significant.

##test for multilinear
vif_model <- lm(rep(1, nrow(data)) ~ Syndecan_group + gender + ckd +
                  Cerebrovascular_disease + Past_myocardial_infarction +
                  diabetes + age_yrs + SOFAtot, data = data)
vif(vif_model)
