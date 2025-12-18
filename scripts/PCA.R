rm(list = ls())

##install.packages("pacman")
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)
library(ggrepel)
install.packages("rstatix")
library(rstatix)
library(dplyr)
install.packages("writexl")
library(writexl)
# === Rename PCA marker names to formatted abbreviations ===
rename_markers <- function(pca_obj) {
  new_names <- rownames(pca_obj$rotation)
  new_names <- gsub("^Ang1$", "ANG1", new_names)
  new_names <- gsub("^Ang2$", "ANG2", new_names)
  new_names <- gsub("^CX3CL1$", "Fractalkine", new_names)
  new_names <- gsub("^Thrombomodulin$", "sThrombomodulin", new_names)
  new_names <- gsub("^ESM1$", "Endocan", new_names)
  new_names <- gsub("^Syndecan1$", "Syndecan-1", new_names)
  new_names <- gsub("^IL10$", "IL-10", new_names)
  new_names <- gsub("^IL18$", "IL-18", new_names)
  new_names <- gsub("^IL1ra$", "IL-1RA", new_names)
  new_names <- gsub("^IL6$", "IL-6", new_names)
  new_names <- gsub("^IL8$", "IL-8", new_names)
  new_names <- gsub("^MMP8$", "MMP-8", new_names)
  new_names <- gsub("^Ddimer$", "D-dimer", new_names)
  new_names <- gsub("^Platelets_value_1$", "Platelet count", new_names)
  new_names <- gsub("^PT_Max_24h$", "PT", new_names)
  new_names <- gsub("^CoagulationFactorIII$", "Tissue Factor", new_names)
  new_names <- gsub("^PROC$", "Protein-C", new_names)
  new_names <- gsub("^TREM1$", "sTREM-1", new_names)
  rownames(pca_obj$rotation) <- new_names
  return(pca_obj)
}


#ready data
data <- read.csv("Original_data/merged_data.csv")
# List of markers (y-axis variables)
markers <- c("Antitrombin", "PROC", "Ang1", "Ang2", "CD163", "CoagulationFactorIII",
             "CX3CL1", "Ddimer", "ESM1", "IL10", "IL18", "IL1ra", "IL6", 
             "IL8", "MMP8", "NGAL", "Procalcitonin", "RAGE", "TenascinC", "Thrombomodulin",
             "TREM1", "Platelets_value_1", "PT_Max_24h")
markers <- c("Ang1", "Ang2", "CD163", "CoagulationFactorIII",
             "CX3CL1", "Ddimer", "ESM1", "IL10", "IL18", "IL1ra", "IL6", 
             "IL8", "MMP8", "NGAL", "Procalcitonin", "RAGE", "TenascinC", "Thrombomodulin",
             "TREM1", "Platelets_value_1", "PT_Max_24h")
#log the markers data
data[markers] <- log10(data[markers])

#Grouped the patients based on syndecan
data$Syndecan_group <- cut(data$Syndecan1,breaks = quantile(data$Syndecan1, probs = seq(0, 1, length = 4), na.rm = TRUE),
                           include.lowest = TRUE,labels = c(1, 2, 3))


#Systemic Inflammation and Organ Damage
endo <- data[c("CD163", "MMP8", "NGAL", "Procalcitonin", "RAGE", "TenascinC", "TREM1","IL10", "IL18", "IL1ra", "IL6", "IL8","Syndecan1","Syndecan_group")]

#PCA
endo <- na.omit(endo) 
endo_pca <- endo %>% select(c("CD163", "MMP8", "NGAL", "Procalcitonin", "RAGE", "TenascinC", "TREM1","IL10", "IL18", "IL1ra", "IL6", "IL8",
                  "Syndecan_group"))
endo_pca$Group <- as.factor(endo_pca$Syndecan_group)
e.pca <- prcomp(endo_pca[ ,1:12], center = TRUE,scale. = T) 
e.pca <- rename_markers(e.pca)
summary(e.pca)

#*-1 for PC1 and PC2, just for understanding
#e.pca$x[,1] <- -e.pca$x[,1]
e.pca$x[,2] <- -e.pca$x[,2]

#e.pca$rotation[,1] <- -e.pca$rotation[,1]
e.pca$rotation[,2] <- -e.pca$rotation[,2]


#PCA figure
colors <- c("#228b22","#377eb8", "#ff7f0e")  
#run ggbiplot.script first

e.plot1 <- ggbiplot(e.pca, ellipse = FALSE, obs.scale = 1.5, var.scale = 1, var.axes = TRUE, group = endo_pca$Group,
                    circle = FALSE, varname.size = 4, alpha = 0, varname.adjust = c(2)) +
  scale_color_manual(name = "Group", values = colors) +
  geom_point(aes(colour = endo_pca$Group), size = 1, alpha = 0) +
  stat_ellipse(aes(colour = endo_pca$Group), size = 1, type = "norm", level = 0.25) +
  coord_cartesian(xlim = c(-6, 6), ylim = c(-2, 2)) +
  theme_classic(base_size = 14, base_family = "sans") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title = element_text(size = 12, family = "sans"),
    axis.text = element_text(size = 10, family = "sans", color = "black"),
    legend.title = element_text(size = 12, family = "sans"),
    legend.text = element_text(size = 10, family = "sans"),
    plot.margin = margin(5, 5, 5, 5),
    legend.position = "none",
    plot.title = element_blank()  # 
  ) +
  scale_y_continuous(breaks = seq(-3, 3, by = 1))


e.plot1
file_path <- "figures/Inflammation_PCA_groups.svg" 
ggsave(file_path, plot = e.plot1, width = 8, height = 8)


##PCA box p value
## test plots => first get data
pcdat <- as.data.frame(e.pca[["x"]])
pcdat$group <- as.factor(endo$Syndecan_group)
##ANOVA test
anova_results <- summary(aov(PC1 ~ group, data = pcdat))
anova_results
p_value_anova_PC1 <- anova_results[[1]][["Pr(>F)"]][1]
#PC1 6.08e-12 ***
#PC2 0.0368 *

#need to install in a fresh R sesssion
install.packages("rstatix")
library("rstatix")
#store the tukey_hsd resluts
if (p_value_anova_PC1 < 0.05) {
  tukey_results <- tukey_hsd(aov(PC1 ~ group, data = pcdat))
  # Convert to data frame and store
  tukey_results_df <- as.data.frame(tukey_results)
} else {
  tukey_results_df <- data.frame()  # Empty data frame if ANOVA is not significant
}
tukey_results_df
# PC2 group1_3 *
# PC1 group1_2 ***; group1_3 ***; group2_3 *

#spearman
# Spearman to make the relationship between syndecan and PC2
cor_syndecan_PC1 <- cor.test(endo$Syndecan1, pcdat$PC1, method = "spearman")
cor_syndecan_PC1 <- cor.test(endo$Syndecan1, pcdat$PC1, method = "pearson")

spearman_rho <- cor_syndecan_PC1$estimate
spearman_rho
#PC2 0.2124463 
#PC1 0.5827811   
spearman_p <- cor_syndecan_PC1$p.value
spearman_p
#PC2 0.003481517
#PC1 5.80122e-24
# Plot boxplot for PC1

##for the min and max
min_value <- min(pcdat$PC1, na.rm = TRUE) 
max_value <- max(pcdat$PC1, na.rm = TRUE)
# Plot boxplot for PC1
PC1_boxplot <- ggplot(pcdat, aes(x = group, y = PC1, fill = group)) +
  geom_boxplot(outlier.shape = NA) +  
  scale_fill_manual(values = colors) +  
  theme_bw() +  
  theme(legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),  
        axis.title.y = element_text(size = 15), 
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),  # Removes major grid lines
        panel.grid.minor = element_blank()   # Removes minor grid lines
  ) +  
  coord_cartesian(ylim = c(min_value, max_value))


print(PC1_boxplot)
file_path <- "figures/inflammation_PC1_groups.svg" 
ggsave(file_path, plot = PC1_boxplot, width = 4, height = 4)

# Extract the contributions of each biomarker to the first two principal components
biomarker_contributions <- as.data.frame(e.pca$rotation[, 1:2])
colnames(biomarker_contributions) <- c("Contribution_to_PC1", "Contribution_to_PC2")
biomarker_contributions$Percentage_Contribution_to_PC1 <- (biomarker_contributions$Contribution_to_PC1^2) * 100
biomarker_contributions$Percentage_Contribution_to_PC2 <- (biomarker_contributions$Contribution_to_PC2^2) * 100
biomarker_contributions <- biomarker_contributions %>%
  select(-Contribution_to_PC1, -Contribution_to_PC2)
biomarker_contributions

biomarker_contributions$marker <- rownames(biomarker_contributions)
write_xlsx(biomarker_contributions, path = "figures/inflammation_contributions_3groups.xlsx")



#Endothelial Cell Activation and Dysfunction
endo <- data[c("Ang1", "Ang2", "CX3CL1", "Thrombomodulin", "ESM1",
               "Syndecan1","Syndecan_group")]

#PCA
endo <- na.omit(endo) 
endo_pca <- endo %>% select(c("Ang1", "Ang2", "CX3CL1", "Thrombomodulin", "ESM1",
                              "Syndecan_group"))
endo_pca$Group <- as.factor(endo_pca$Syndecan_group)
e.pca <- prcomp(endo_pca[ ,1:5], center = TRUE,scale. = T) 
e.pca <- rename_markers(e.pca)

summary(e.pca)

#*-1 for PC1 and PC2, just for understanding
e.pca$x[,1] <- -e.pca$x[,1]
#e.pca$x[,2] <- -e.pca$x[,2]

e.pca$rotation[,1] <- -e.pca$rotation[,1]
#e.pca$rotation[,2] <- -e.pca$rotation[,2]


#PCA figure
colors <- c("#228b22","#377eb8", "#ff7f0e")  
#run ggbiplot.script first


pc1_var <- round(summary(e.pca)$importance[2, 1] * 100, 1)
pc2_var <- round(summary(e.pca)$importance[2, 2] * 100, 1)

e.plot1 <- ggbiplot(e.pca, ellipse = FALSE, obs.scale = 1.5, var.scale = 1, var.axes = TRUE,
                    group = endo_pca$Group, circle = FALSE, varname.size = 4, alpha = 0,
                    varname.adjust = c(2)) +
  scale_color_manual(name = "Group", values = colors) +
  geom_point(aes(colour = endo_pca$Group), size = 1, alpha = 0) +
  stat_ellipse(aes(colour = endo_pca$Group), size = 1, type = "norm", level = 0.25) +
  coord_cartesian(xlim = c(-4, 4), ylim = c(-2, 2)) +
  scale_y_continuous(breaks = seq(-3, 3, by = 1)) +
  ggtitle("Endothelial Cell Activation and Dysfunction") +
  labs(
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)")
  ) +
  theme_bw(base_family = "sans") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title = element_text(size = 12, family = "sans"),
    axis.text.x = element_text(size = 10, family = "sans", color = "black"),
    axis.text.y = element_text(size = 10, family = "sans", color = "black"),
    plot.title = element_text(size = 12, hjust = 0, family = "sans"),
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5)
  )
e.plot1

file_path <- "figures/Endothelial_PCA_groups.svg" 
ggsave(file_path, plot = e.plot1, width = 8, height = 8)




##PCA box p value
## test plots => first get data
pcdat <- as.data.frame(e.pca[["x"]])
pcdat$group <- as.factor(endo$Syndecan_group)
##ANOVA test
anova_results <- summary(aov(PC2 ~ group, data = pcdat))
anova_results
p_value_anova_PC2 <- anova_results[[1]][["Pr(>F)"]][1]
#PC2 0.151
#PC1 ***

#need to install in a fresh R sesssion
install.packages("rstatix")
library("rstatix")
#store the tukey_hsd resluts
if (p_value_anova_PC2 < 0.05) {
  tukey_results <- tukey_hsd(aov(PC2 ~ group, data = pcdat))
  # Convert to data frame and store
  tukey_results_df <- as.data.frame(tukey_results)
} else {
  tukey_results_df <- data.frame()  # Empty data frame if ANOVA is not significant
}
tukey_results_df
# PC2 
# PC1 group1_2****group1_3****group2_3 ***

#spearman
# Spearman to make the relationship between syndecan and PC2
cor_syndecan_PC1 <- cor.test(endo$Syndecan1, pcdat$PC1, method = "pearson")

spearman_rho <- cor_syndecan_PC1$estimate
spearman_rho
#PC2 -0.1143659
#PC1 0.6954044  
spearman_p <- cor_syndecan_PC1$p.value
spearman_p
#PC2 0.1171177
#PC1 1.218985e-28
# Plot boxplot for PC2

##for the min and max
min_value <- min(pcdat$PC2, na.rm = TRUE) 
max_value <- max(pcdat$PC2, na.rm = TRUE)
# Plot boxplot for PC2
PC2_boxplot <- ggplot(pcdat, aes(x = group, y = PC2, fill = group)) +
  geom_boxplot(outlier.shape = NA) +  
  scale_fill_manual(values = colors) +  
  theme_bw() +  
  theme(legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),  
        axis.title.y = element_text(size = 15), 
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),  # Removes major grid lines
        panel.grid.minor = element_blank()   # Removes minor grid lines
  ) +  
  coord_cartesian(ylim = c(min_value, max_value))


print(PC2_boxplot)
file_path <- "figures/endothelial_PC2_groups.svg" 
ggsave(file_path, plot = PC1_boxplot, width = 4, height = 4)



# Extract the contributions of each biomarker to the first two principal components
biomarker_contributions <- as.data.frame(e.pca$rotation[, 1:2])
colnames(biomarker_contributions) <- c("Contribution_to_PC1", "Contribution_to_PC2")
biomarker_contributions$Percentage_Contribution_to_PC1 <- (biomarker_contributions$Contribution_to_PC1^2) * 100
biomarker_contributions$Percentage_Contribution_to_PC2 <- (biomarker_contributions$Contribution_to_PC2^2) * 100
biomarker_contributions <- biomarker_contributions %>%
  select(-Contribution_to_PC1, -Contribution_to_PC2)
biomarker_contributions

biomarker_contributions$marker <- rownames(biomarker_contributions)
write_xlsx(biomarker_contributions, path = "figures/Endothelial Cell Activation and Dysfunction_contributions_3groups.xlsx")



#Coagulation Activation
endo <- data[c("CoagulationFactorIII", "Ddimer", "Platelets_value_1",  "PT_Max_24h", "Antitrombin", "PROC",
               "Syndecan1","Syndecan_group")]
endo <- data[c("CoagulationFactorIII", "Ddimer", "Platelets_value_1",  "PT_Max_24h", 
               "Syndecan1","Syndecan_group")]
#PCA
endo <- na.omit(endo) 
endo_pca <- endo %>% select(c("CoagulationFactorIII", "Ddimer", "Platelets_value_1",  "PT_Max_24h", 
                              "Syndecan_group"))
endo_pca$Group <- as.factor(endo_pca$Syndecan_group)
e.pca <- prcomp(endo_pca[ ,1:4], center = TRUE,scale. = T) 
e.pca <- rename_markers(e.pca)
summary(e.pca)

#*-1 for PC1 and PC2, just for understanding
#e.pca$x[,1] <- -e.pca$x[,1]
e.pca$x[,2] <- -e.pca$x[,2]

#e.pca$rotation[,1] <- -e.pca$rotation[,1]
e.pca$rotation[,2] <- -e.pca$rotation[,2]


#PCA figure
colors <- c("#228b22","#377eb8", "#ff7f0e")  
#run ggbiplot.script first

# 
pc1_var <- round(summary(e.pca)$importance[2, 1] * 100, 1)
pc2_var <- round(summary(e.pca)$importance[2, 2] * 100, 1)

# 
e.plot1 <- ggbiplot(e.pca, ellipse = FALSE, obs.scale = 1.5, var.scale = 1,
                    var.axes = TRUE, group = endo_pca$Group, circle = FALSE,
                    varname.size = 4, alpha = 0, varname.adjust = c(2)) +
  scale_color_manual(name = "Group", values = colors) +
  geom_point(aes(colour = endo_pca$Group), size = 1, alpha = 0) +
  stat_ellipse(aes(colour = endo_pca$Group), size = 1, type = "norm", level = 0.25) +
  coord_cartesian(xlim = c(-4, 4), ylim = c(-2, 2)) +
  scale_y_continuous(breaks = seq(-3, 3, by = 1)) +
  labs(
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)")
  ) +
  theme_bw(base_family = "sans") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title = element_text(size = 12, family = "sans"),
    axis.text.x = element_text(size = 10, family = "sans", color = "black"),
    axis.text.y = element_text(size = 10, family = "sans", color = "black"),
    plot.title = element_text(size = 12, hjust = 0, family = "sans"),
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5)
  )
e.plot1

file_path <- "figures/coagulation_PCA_groups.svg" 
ggsave(file_path, plot = e.plot1, width = 8, height = 8)


##PCA box p value
## test plots => first get data
pcdat <- as.data.frame(e.pca[["x"]])
pcdat$group <- as.factor(endo$Syndecan_group)
##ANOVA test
anova_results <- summary(aov(PC2 ~ group, data = pcdat))
anova_results
p_value_anova_PC2 <- anova_results[[1]][["Pr(>F)"]][1]
#PC2 ***
#PC1 ***

#need to install in a fresh R sesssion
install.packages("rstatix")
library("rstatix")
#store the tukey_hsd resluts
if (p_value_anova_PC2 < 0.05) {
  tukey_results <- tukey_hsd(aov(PC2 ~ group, data = pcdat))
  # Convert to data frame and store
  tukey_results_df <- as.data.frame(tukey_results)
} else {
  tukey_results_df <- data.frame()  # Empty data frame if ANOVA is not significant
}
tukey_results_df
# PC2 group1_2 **; group1_3 ****; group2_3 ns
# PC1 group1_2 ns; group1_3 ***; group2_3 *

#spearman
# Spearman to make the relationship between syndecan and PC2
cor_syndecan_PC2 <- cor.test(endo$Syndecan1, pcdat$PC2, method = "pearson")

spearman_rho <- cor_syndecan_PC2$estimate
spearman_rho
#PC1 -0.3344602  
#PC1 0.3284
spearman_p <- cor_syndecan_PC2$p.value
spearman_p
#PC1 2.554406e-06
#PC1 3.961332e-06
# Plot boxplot for PC1

##for the min and max
min_value <- min(pcdat$PC2, na.rm = TRUE) 
max_value <- max(pcdat$PC2, na.rm = TRUE)
# Plot boxplot for PC1
PC2_boxplot <- ggplot(pcdat, aes(x = group, y = PC2, fill = group)) +
  geom_boxplot(outlier.shape = NA) +  
  scale_fill_manual(values = colors) +  
  theme_bw() +  
  theme(legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),  
        axis.title.y = element_text(size = 15), 
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),  # Removes major grid lines
        panel.grid.minor = element_blank()   # Removes minor grid lines
  ) +  
  coord_cartesian(ylim = c(min_value, max_value))


print(PC2_boxplot)
file_path <- "figures/coagulation_PC2_groups.svg" 
ggsave(file_path, plot = PC2_boxplot, width = 4, height = 4)


# Extract the contributions of each biomarker to the first two principal components
biomarker_contributions <- as.data.frame(e.pca$rotation[, 1:2])
colnames(biomarker_contributions) <- c("Contribution_to_PC1", "Contribution_to_PC2")
biomarker_contributions$Percentage_Contribution_to_PC1 <- (biomarker_contributions$Contribution_to_PC1^2) * 100
biomarker_contributions$Percentage_Contribution_to_PC2 <- (biomarker_contributions$Contribution_to_PC2^2) * 100
biomarker_contributions <- biomarker_contributions %>%
  select(-Contribution_to_PC1, -Contribution_to_PC2)
biomarker_contributions

biomarker_contributions$marker <- rownames(biomarker_contributions)
write_xlsx(biomarker_contributions, path = "figures/Coagulation Activation_contributions_3groups.xlsx")





#1.moved some NA each figure. need to fix it