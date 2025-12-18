rm(list = ls())
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)
p_load(rms)

merged_data <- read.csv("original_data/merged_data.csv")
# List of markers (y-axis variables)
markers <- c("Syndecan1","Antitrombin", "PROC", "Ang1", "Ang2", "CD163", "CoagulationFactorIII",
             "CX3CL1", "Ddimer", "ESM1", "IL10", "IL18", "IL1ra", "IL6", 
             "IL8", "MMP8", "NGAL", "Procalcitonin", "RAGE", "TenascinC", "Thrombomodulin",
             "TREM1", "Platelets_value_1", "PT_Max_24h")

# Standardize biomarkers

merged_data <- merged_data %>%
  mutate(across(all_of(markers), ~ log(. + 1))) # Add 1 to avoid log(0) issues
merged_data[markers] <- scale(merged_data[markers])

# Create an empty dataframe to store AIC and Nonlinear P-values
results_df <- data.frame(Biomarker = markers, Linear_AIC = NA, Nonlinear_AIC = NA, NonlinearP = NA)

# Loop through each biomarker
for (marker in markers) {
  
  # Filter out missing values for the current biomarker
  data <- merged_data %>% filter(!is.na(!!sym(marker)))
  
  # Define models
  formula_linear <- as.formula(paste(marker, " ~ Syndecan1"))
  formula_nonlinear <- as.formula(paste(marker, " ~ rcs(Syndecan1, 3)")) # Restricted cubic spline with 3 knots
  
  # Fit linear model
  model_linear <- ols(formula_linear, data = data)
  linear_aic <- AIC(model_linear)
  
  # Fit nonlinear model
  model_nonlinear <- ols(formula_nonlinear, data = data)
  nonlinear_aic <- AIC(model_nonlinear)
  
  # Get Nonlinear P-value from ANOVA test
  anova_result <- anova(model_nonlinear)
  nonlinear_p <- anova_result[3, "P"]  # The p-value for nonlinearity (3rd row in ANOVA)
  
  # Store results
  results_df[results_df$Biomarker == marker, "Linear_AIC"] <- linear_aic
  results_df[results_df$Biomarker == marker, "Nonlinear_AIC"] <- nonlinear_aic
  results_df[results_df$Biomarker == marker, "NonlinearP"] <- nonlinear_p
}

# Adjust p-values using Benjamini-Hochberg (FDR correction)
results_df$NonlinearP <- p.adjust(results_df$NonlinearP, method = "BH")
