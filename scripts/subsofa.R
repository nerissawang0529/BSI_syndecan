rm(list = ls())
install.packages("rms")
library(rms)
install.packages("FSA") 
library(FSA)
install.packages("VGAM")

## install.packages("pacman")
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)
install.packages("dplyr")
install.packages("FSA")
library(dplyr)
library(mice)
library(tidyr)
library(FSA)


#sensitivy analysis 1. ALL patients (N=190) 2. only those with known single pathogen (N=150)

#ready data
data <- read.csv("original_data/final_data.csv")

#I compared the new subsofa from Hessel and the old subsofa, they almost the same. the NAs had been solved by Hessel

sofa_data <- data[,c("MARSID","Syndecan1","SOFAtot_new","SOFA_Circulation_new","SOFA_Coagulation_new","SOFA_Liver_new","SOFA_Renal_new","SOFA_Respiration_new")]

dd <- datadist(sofa_data)
options(datadist = "dd")

#proportional odds model

##1.SOFA_Circulation_new
###1.1the original score are 0 1 3 4 
# Convert SOFA_Circulation_new to an ordered factor (without omitting '2')
sofa_data$SOFA_Circulation_new <- factor(sofa_data$SOFA_Circulation_new, 
                                         levels = c(0, 1, 3, 4), 
                                         ordered = TRUE)
# Fit the proportional odds model with rms
fit_Circulation_original <- orm(SOFA_Circulation_new ~ Syndecan1, data = sofa_data)
# View the model summary
summary(fit_Circulation_original)
anova(fit_Circulation_original)

###1.2sensitive make score into 0 1 2 3
# Recode SOFA_Circulation_new (omit '2' and map 4 → 3)
sofa_data <- sofa_data %>%
  mutate(SOFA_Circulation_sens = case_when(
    SOFA_Circulation_new == 0 ~ 0,
    SOFA_Circulation_new == 1 ~ 1,
    SOFA_Circulation_new == 3 ~ 2,  # 3 becomes 2
    SOFA_Circulation_new == 4 ~ 3   # 4 becomes 3
  ))

# Convert to an ordered factor
sofa_data$SOFA_Circulation_sens <- factor(sofa_data$SOFA_Circulation_sens,
                                          levels = c(0, 1, 2, 3),
                                          ordered = TRUE)
# Fit the proportional odds model with rms
fit_Circulation_sens <- orm(SOFA_Circulation_sens ~ Syndecan1, data = sofa_data)
# View the model summary
summary(fit_Circulation_sens)
anova(fit_Circulation_sens)

#2.SOFA_Coagulation_new
sofa_data$SOFA_Coagulation_new <- factor(sofa_data$SOFA_Coagulation_new, 
                                         levels = c(0, 1, 2, 3, 4), 
                                         ordered = TRUE)
fit_Coagulation <- orm(SOFA_Coagulation_new ~ Syndecan1, data = sofa_data)
summary(fit_Coagulation)
anova(fit_Coagulation)

#3.SOFA_Liver_new
sofa_data$SOFA_Liver_new <- factor(sofa_data$SOFA_Liver_new, 
                                         levels = c(0, 1, 2, 3, 4), 
                                         ordered = TRUE)
fit_Liver <- orm(SOFA_Liver_new ~ Syndecan1, data = sofa_data)
summary(fit_Liver)
anova(fit_Liver)

#4.SOFA_Renal_new
sofa_data$SOFA_Renal_new <- factor(sofa_data$SOFA_Renal_new, 
                                   levels = c(0, 1, 2, 3, 4), 
                                   ordered = TRUE)
fit_Renal <- orm(SOFA_Renal_new ~ Syndecan1, data = sofa_data)
summary(fit_Renal)
anova(fit_Renal)


#5.SOFA_Respiration_new
sofa_data$SOFA_Respiration_new <- factor(sofa_data$SOFA_Respiration_new, 
                                   levels = c(0, 1, 2, 3, 4), 
                                   ordered = TRUE)
fit_Respiration <- orm(SOFA_Respiration_new ~ Syndecan1, data = sofa_data)
summary(fit_Respiration)
anova(fit_Respiration)

# Function to extract OR, CI, and p-value from an orm model
extract_orm_results <- function(model, model_name) {
  sum_model <- summary(model)  # Get model summary
  anova_model <- anova(model)  # Get p-value from anova
  
  # Extract effect estimate, SE, and confidence intervals
  effect <- sum_model[1, "Effect"]
  lower_ci <- sum_model[1, "Lower 0.95"]
  upper_ci <- sum_model[1, "Upper 0.95"]
  
  # Extract p-value for syndecan_1 from anova table
  p_value <- anova_model[which(rownames(anova_model) == "Syndecan1"), "P"]
  
  return(data.frame(
    Predictor = model_name,
    OR = exp(effect),  # Convert log-odds to OR
    Lower = exp(lower_ci),
    Upper = exp(upper_ci),
    p_value = p_value
  ))
}

# Extract results for all SubSOFA models
forest_data <- rbind(
  extract_orm_results(fit_Circulation_original, "SOFA Circulation"),
  extract_orm_results(fit_Coagulation, "SOFA Coagulation"),
  extract_orm_results(fit_Liver, "SOFA Liver"),
  extract_orm_results(fit_Renal, "SOFA Renal"),
  extract_orm_results(fit_Respiration, "SOFA Respiration")
)

# Format p-values for display
forest_data$p_text <- ifelse(forest_data$p_value < 0.001, "< 0.001", 
                             sprintf("p = %.3f", forest_data$p_value))

# Print extracted data
print(forest_data)

# Format p-values consistently (remove "p = ")
forest_data$p_text <- ifelse(forest_data$p_value < 0.001, "< 0.001", 
                             sprintf("%.3f", forest_data$p_value))

# Set base theme with font size 12 and sans-serif family
theme_set(theme_minimal(base_size = 12, base_family = "sans"))

# Create forest plot and store it in object `p`
p <- ggplot(forest_data, aes(y = Predictor, x = OR)) +
  geom_point(size = 3) +  # Point for odds ratio
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2) +  # Horizontal error bars
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +  # Reference line at OR = 1
  geom_text(aes(x = Upper + 0.2, label = p_text), hjust = 0, size = 4, family = "sans") +  # Add p-value text
  xlab("Odds Ratio (95% CI)") +
  ylab("SubSOFA Components") +
  ggtitle("Forest Plot of SubSOFA (Proportional Odds Model)") +
  theme(
    axis.text.y = element_text(size = 12, family = "sans"),
    axis.text.x = element_text(size = 12, family = "sans"),
    axis.title.x = element_text(size = 12, family = "sans"),
    axis.title.y = element_text(size = 12, family = "sans"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12, family = "sans")
  )

# Display the plot
print(p)

# Create the output directory if it doesn't exist
dir.create("figures", showWarnings = FALSE)

# Save the plot as an SVG file
ggsave(filename = "figures/sub_sofa.svg", plot = p,
       width = 8, height = 4, units = "in", dpi = 300)





###Correlation of syndecan-1 and sofa
# Convert all relevant columns to numeric
numeric_columns <- c("SOFA_Respiration_new", "SOFA_Coagulation_new", "SOFA_Liver_new",
                     "SOFA_Circulation_new", "SOFA_Renal_new", "SOFAtot_new", "syndecan_1")

sofa_data[numeric_columns] <- lapply(sofa_data[numeric_columns], function(x) as.numeric(as.character(x)))

# Remove rows with missing values before plotting
sofa_data_clean <- sofa_data %>% filter(!is.na(SOFAtot_new) & !is.na(syndecan_1))

# Scatter plot: Syndecan-1 vs. SOFA total score
p1 <- ggplot(sofa_data_clean, aes(x = syndecan_1, y = SOFAtot_new)) +
  geom_point(color = "black", alpha = 0.6) +  # Scatter points
  geom_smooth(method = "lm", color = "purple", fill = "purple", alpha = 0.2) +  # Regression line
  theme_minimal() +
  labs(x = "Syndecan-1", y = "SOFA score") +
  theme(axis.title = element_text(face = "bold")) 

print(p1)  # Display the scatter plot
# Calculate Spearman correlation between Syndecan-1 and SOFA total score
cor_test <- cor.test(sofa_data_clean$syndecan_1, sofa_data_clean$SOFAtot_new, 
                     method = "spearman", use = "complete.obs")

# Extract values
spearman_rho <- cor_test$estimate  # Correlation coefficient
p_value <- cor_test$p.value  # p-value

# Print the results
cat("Spearman’s rho:", round(spearman_rho, 2), "\n")
cat("p-value:", format.pval(p_value, digits = 3, eps = 0.001), "\n")






# Convert all relevant columns to numeric
numeric_columns <- c("SOFA_Respiration_new", "SOFA_Coagulation_new", "SOFA_Liver_new",
                     "SOFA_Circulation_new", "SOFA_Renal_new", "SOFAtot_new", "syndecan_1")

sofa_data[numeric_columns] <- lapply(sofa_data[numeric_columns], function(x) as.numeric(as.character(x)))

# Remove rows with missing values before correlation analysis
sofa_data_clean <- sofa_data %>% filter(!is.na(SOFAtot_new) & !is.na(syndecan_1))

# Extract SOFA subscores
sofa_subscores <- sofa_data_clean[, c("SOFA_Respiration_new", "SOFA_Coagulation_new", 
                                      "SOFA_Liver_new", "SOFA_Circulation_new", "SOFA_Renal_new")]

# Calculate Spearman correlations and p-values
cor_results <- sapply(sofa_subscores, function(x) {
  if (is.numeric(x)) {
    test <- cor.test(x, sofa_data_clean$syndecan_1, method = "spearman", use = "complete.obs")
    return(c(rho = test$estimate, p_value = test$p.value))
  } else {
    return(c(rho = NA, p_value = NA))
  }
})

# Convert results into a tidy dataframe
cor_results <- as.data.frame(t(cor_results))  # Transpose to make it readable
colnames(cor_results) <- c("rho", "p_value")
cor_results$Component <- rownames(cor_results)

# Add significance stars based on p-values
cor_results$significance <- ifelse(cor_results$p_value < 0.001, "***",
                                   ifelse(cor_results$p_value < 0.01, "**",
                                          ifelse(cor_results$p_value < 0.05, "*", "")))

# Ensure rho is numeric
cor_results$rho <- as.numeric(cor_results$rho)

# Sort by rho in descending order
cor_results <- cor_results %>% arrange(-rho)  # Use "-" instead of desc()

# Rename components
rename_map <- c("SOFA_Respiration_new" = "Respiration",
                "SOFA_Renal_new" = "Renal",
                "SOFA_Liver_new" = "Liver",
                "SOFA_Coagulation_new" = "Coagulation",
                "SOFA_Circulation_new" = "Circulation")

cor_results$Component <- rename_map[cor_results$Component]

# Create heatmap with updated ordering and labels
p2 <- ggplot(cor_results, aes(x = "Syndecan-1", y = reorder(Component, rho), fill = rho)) +
  geom_tile() +
  geom_text(aes(label = significance), color = "black", size = 6) +  # Add significance stars
  scale_fill_gradient(low = "white", high = "darkred", limits = c(0, 0.4)) +  
  labs(fill = "rho") +
  theme_minimal() +
  theme(axis.text.x = element_text(face = "bold"), 
        axis.text.y = element_text(face = "bold"),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())

# Combine plots
final_plot <- ggarrange(p1, p2, labels = c("A", "B"), ncol = 2)

# Display plot
print(final_plot)

print(cor_values)



#Correlation with continus varioubles 
data$Bilirubin          #liver
data$P_F_ratio_Min_24h  #respiration


#Coagulation: Platelets_value_1/Coagulation
#check if linear or non-linear of the valuebles and syndecan-1
# Load required libraries
data_platelet <- data[,c("syndecan_original","Platelets_value_1")]
# Ensure there are no missing values (rms does not handle NAs well)
data_platelet <- na.omit(data_platelet)
# Set up datadist
d <- datadist(data_platelet)  # Define datadist based on your dataset
options(datadist = "d")  # Tell rms to use this datadist

# Fit a linear regression model
linear_model <- ols(Platelets_value_1 ~ syndecan_original, data= data_platelet)
# Fit a restricted cubic spline (RCS) model with 3 knots
spline_model <- ols(Platelets_value_1 ~ rcs(syndecan_original, 3), data = data_platelet)
# Compare models using ANOVA to test for non-linearity
anova(spline_model)
# Compare models using AIC (lower AIC = better model)
AIC(linear_model) #2412.841
AIC(spline_model) #2407.978
#for platelet, the linear model is better

# Fit a spline model
spline_model_plate <- ols(Platelets_value_1 ~ rcs(syndecan_original, 3), data = data_platelet)# Adjust reference level for Syndecan-1 at 7.670242
median(data$syndecan_original) #13398.26
d$limits["Adjust to", "syndecan_1"] <- 13398.26   

# Update the model to reflect new reference
linear_model <- update(spline_model_plate)

# Predict Platelet Count while setting reference at 13398.26
x1 <- Predict(spline_model_plate, syndecan_original, ref.zero = TRUE)

# Summary of predictions
summary(x1)

# Create the plot
x3 <- ggplot(x1) + 
  xlab("Syndecan-1 (pg/ml)") +
  ylab("Platelet Count") + 
  ggtitle("Association Between Syndecan-1 and Platelet Count") + 
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
  geom_hline(yintercept= 0, linetype='solid', col = 'black') +  # Reference at zero difference
  geom_vline(xintercept = 13398.26, linetype = 'dashed', col = 'black') +  # Reference level
  scale_y_continuous() +
  coord_cartesian(ylim = c(min(x1$yhat) - 30, max(x1$yhat) + 30))  # Adjust Y-axis limits

# Display the plot
x3


#Creatinine, renal
data_Creatinine <- data[,c("syndecan_original","Creatinine_value_1")]
# Ensure there are no missing values (rms does not handle NAs well)
data_Creatinine <- na.omit(data_Creatinine)
# Set up datadist
d <- datadist(data_Creatinine)  # Define datadist based on your dataset
options(datadist = "d")  # Tell rms to use this datadist

# Fit a linear regression model
linear_model_creatinine <- ols(Creatinine_value_1 ~ syndecan_original, data= data_Creatinine)
# Fit a restricted cubic spline (RCS) model with 3 knots
spline_model_creatinine <- ols(Creatinine_value_1 ~ rcs(syndecan_original, 3), data = data_Creatinine)
# Compare models using ANOVA to test for non-linearity
anova(spline_model_creatinine)
# Compare models using AIC (lower AIC = better model)
AIC(linear_model_creatinine) #2472.857
AIC(spline_model_creatinine) #2465.097
#for platelet, the linear model is better

# Fit a linear model
spline_model_creatinine <- ols(Creatinine_value_1 ~ rcs(syndecan_original, 3), data = data_Creatinine)
median(data$syndecan_original) #13398.26
d$limits["Adjust to", "syndecan_1"] <- 13398.26 

# Update the model to reflect new reference
linear_model <- update(spline_model_creatinine)

# Predict Platelet Count while setting reference at 9.514
x1 <- Predict(spline_model_creatinine, syndecan_original, ref.zero = TRUE)

# Summary of predictions
summary(x1)

# Create the plot
x3 <- ggplot(x1) + 
  xlab("Syndecan-1 (pg/ml)") +
  ylab("Creatinine (µmol/L)") + 
  ggtitle("Association Between Syndecan-1 and Creatinine (Sofa Renal)") + 
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
  geom_hline(yintercept= 0, linetype='solid', col = 'black') +  # Reference at zero difference
  geom_vline(xintercept = 13398.26, linetype = 'dashed', col = 'black') +  # Reference level
  scale_y_continuous() +
  coord_cartesian(ylim = c(min(x1$yhat) - 30, max(x1$yhat) + 30))  # Adjust Y-axis limits

# Display the plot
x3


#Bilirubin for renal
#imputate Bilirubin
