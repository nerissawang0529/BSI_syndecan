rm(list = ls())

## install.packages("pacman")
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)
install.packages("svglite")

# Read data
data <- read.csv("Documents/BSI/R_code/original_data/clinical_marker_scource_pathogen.csv")

#box plot
#box plots for syndecan from healthy and patients

my_colors <- c("Non-infectious" = "#ADD8E6",    # Light blue
               "Bacteremia" = "#E76F51",       # Blue
               "Fungal" = "#BC8F8F",           # Purple
               "Viral" = "#FFDAB9")            # Gray

# box plots
p <- ggplot(data, aes(x = factor(Microbe_groups, levels = c("Non-infectious", "Bacteremia", "Fungal", "Viral"), labels = c("Non-infectious", "Bacteremia", "Fungal", "Viral")), 
                                        y = Syndecan.1.CD138..29., 
                                        fill = factor(Microbe_groups, levels = c("Non-infectious", "Bacteremia", "Fungal", "Viral"), labels = c("Non-infectious", "Bacteremia", "Fungal", "Viral")))) +
  geom_boxplot(outlier.shape = NA, size = 0.1) +  # Adjust boxplot line thickness
  scale_fill_manual(values = my_colors, name = NULL) +  
  labs(x = NULL, y = "Syndecan-1 (pg/ml)") +  
  theme_minimal() +
  theme(
    axis.title.y = element_text(size = 30),   
    axis.text.x = element_text(size = 12),     
    axis.text.y = element_text(size = 12),     
    plot.title = element_text(size = 16),      
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5),  # Adjust panel border thickness
    axis.line = element_line(color = "black", size = 0.5),  # Adjust axis line thickness
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0, 40000))
p

##export figure
file_path <- "Documents/BSI/R_code/figures/syndecan_4groups.svg" 
ggsave(file_path, plot = p, width = 8, height = 8)

#make p valuble
#calculate the p value between groups
library(dunn.test)
kruskal_test_result <- kruskal.test(Syndecan.1.CD138..29. ~ Microbe_groups, data = data)
dunn_result <- dunn.test(data$Syndecan.1.CD138..29., 
                         g = data$Microbe_groups, 
                         method = "BH")


#correct for severity
# Linear model with log(Syndecan) as the outcome and group, SOFAtot, and APACHE_IV_Score.x as predictors
linear_model_correct <- lm(log(Syndecan.1.CD138..29.) ~ SOFAtot + APACHE_IV_Score.x + Microbe_groups, data = data)
anova(linear_model_correct)
coef(linear_model_correct)
coef(linear_model_correct)[2]
coef(linear_model_correct)[3]
syndecan_adjusted <- log(data$Syndecan.1.CD138..29.)- coef(linear_model_correct)[2] * data$SOFAtot - coef(linear_model_correct)[3] * data$APACHE_IV_Score.x
data$syndecan_adjusted <-syndecan_adjusted

# box plots
p <- ggplot(data, aes(x = factor(Microbe_groups, levels = c("Non-infectious", "Bacteremia", "Fungal", "Viral"), labels = c("Non-infectious", "Bacteremia", "Fungal", "Viral")), 
                      y = syndecan_adjusted, 
                      fill = factor(Microbe_groups, levels = c("Non-infectious", "Bacteremia", "Fungal", "Viral"), labels = c("Non-infectious", "Bacteremia", "Fungal", "Viral")))) +
  geom_boxplot(outlier.shape = NA, size = 0.1) +  # Adjust boxplot line thickness
  scale_fill_manual(values = my_colors, name = NULL) +  
  labs(x = NULL, y = "Syndecan-1 adjusted severity") +  
  theme_minimal() +
  theme(
    axis.title.y = element_text(size = 30),   
    axis.text.x = element_text(size = 12),     
    axis.text.y = element_text(size = 12),     
    plot.title = element_text(size = 16),      
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5),  # Adjust panel border thickness
    axis.line = element_line(color = "black", size = 0.5),  # Adjust axis line thickness
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(6, 10))
p

##export figure
file_path <- "Documents/BSI/R_code/figures/syndecan_4groups_adjust.svg" 
ggsave(file_path, plot = p, width = 8, height = 8)
