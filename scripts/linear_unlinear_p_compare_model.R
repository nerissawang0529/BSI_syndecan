rm(list = ls())
# install.packages("pacman")
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)

# Read the dataset
merged_data <- read.csv("original_data/merged_data.csv")
colnames(merged_data)[colnames(merged_data) == "CoagulationFactorIII"] <- "Tissue factor"
colnames(merged_data)[colnames(merged_data) == "Platelets_value_1"] <- "Platelets"
colnames(merged_data)[colnames(merged_data) == "PT_Max_24h"] <- "PT"

# List of biomarkers (excluding Syndecan1, which is the predictor)
markers <- c("Syndecan1","Antitrombin", "PROC", "Ang1", "Ang2", "CD163", "Tissue factor",
             "CX3CL1", "Ddimer", "ESM1", "IL10", "IL18", "IL1ra", "IL6", 
             "IL8", "MMP8", "NGAL", "Procalcitonin", "RAGE", "TenascinC", "Thrombomodulin",
             "TREM1", "Platelets", "PT")

markers_data <- merged_data[,c("Syndecan1","Antitrombin", "PROC", "Ang1", "Ang2", "CD163", "Tissue factor",
                               "CX3CL1", "Ddimer", "ESM1", "IL10", "IL18", "IL1ra", "IL6", 
                               "IL8", "MMP8", "NGAL", "Procalcitonin", "RAGE", "TenascinC", "Thrombomodulin",
                               "TREM1", "Platelets", "PT","SOFAtot_new.x","APACHE_IV_Score.x",
                               "gender","ckd","Cerebrovascular_disease","diabetes","age_yrs","Immune_deficiency","Past_myocardial_infarction")]

# Log-transform  the biomarkers
markers_data <- markers_data %>%
  mutate(across(all_of(markers), ~ log(. + 1))) # Adding 1 to avoid log(0) issues

merged_markers_log_studydata <- markers_data


# Load necessary library
library(rms)  # for the rcs function

# Initialize a data frame to store the results
results_2 <- data.frame(Marker = character(), 
                        Model = character(), 
                        P_Value = numeric(), 
                        P_Value_Significance = character(), 
                        P_Value_Adjusted = numeric(), 
                        P_Value_Adjusted_Significance = character(), 
                        P_Value_adjusted_without_Mews = numeric(), 
                        P_Value_adjusted_without_Mews_Significance = character(), 
                        stringsAsFactors = FALSE)
# Function to determine significance stars
get_significance_stars <- function(p_value) {
  if (p_value < 0.001) {
    return("***")
  } else if (p_value < 0.01) {
    return("**")
  } else if (p_value < 0.05) {
    return("*")
  } else {
    return(" ")
  }
}

# Loop through each marker
for (marker in markers) {
  
  # Fit linear model
  linear_model <- lm(merged_markers_log_studydata[[marker]] ~ Syndecan1, data = merged_markers_log_studydata)
  
  # Fit spline (non-linear) model
  spline_model <- lm(merged_markers_log_studydata[[marker]]  ~ rcs(Syndecan1, 3), data = merged_markers_log_studydata)
  
  # Compare models using ANOVA
  model_comparison <- anova(linear_model, spline_model)
  
  # Determine the better model based on ANOVA
  if (model_comparison$`Pr(>F)`[2] < 0.05) {
    # Nonlinear model is better
    chosen_model <- "Nonlinear"
    # Extract the p-value from the ANOVA comparison
    p_value <- model_comparison$`Pr(>F)`[2]
    
    # Adjust the nonlinear model for confounders (including MEWS_score)
    adjusted_model <- lm(merged_markers_log_studydata[[marker]] ~ 
                           gender + age_yrs + ckd + Cerebrovascular_disease + 
                           Past_myocardial_infarction + diabetes + Immune_deficiency + 
                           SOFAtot_new.x + APACHE_IV_Score.x + rcs(Syndecan1, 3),
                         data = merged_markers_log_studydata)
    
    # Adjust the nonlinear model for confounders (excluding MEWS_score)
    adjusted_model_without_mews <- lm(merged_markers_log_studydata[[marker]] ~ 
                                        gender + age_yrs + ckd + Cerebrovascular_disease + 
                                        Past_myocardial_infarction + diabetes + Immune_deficiency + 
                                        rcs(Syndecan1, 3),
                                      data = merged_markers_log_studydata)
    
  } else {
    # Linear model is better
    chosen_model <- "Linear"
    # Extract the p-value for the linear model's first term
    p_value <- anova(linear_model)["Syndecan1","Pr(>F)"]
    # Adjust the linear model for confounders (including MEWS_score)
    adjusted_model <- lm(merged_markers_log_studydata[[marker]] ~ 
                           gender + age_yrs + ckd + Cerebrovascular_disease + 
                           Past_myocardial_infarction + diabetes + Immune_deficiency + 
                           SOFAtot_new.x + APACHE_IV_Score.x + Syndecan1,
                         data = merged_markers_log_studydata)
    
    # Adjust the linear model for confounders (excluding MEWS_score)
    adjusted_model_without_mews <- lm(merged_markers_log_studydata[[marker]] ~ 
                                        gender + age_yrs + ckd + Cerebrovascular_disease + 
                                        Past_myocardial_infarction + diabetes + Immune_deficiency + 
                                        Syndecan1,
                                      data = merged_markers_log_studydata)
  }
  
  # Extract the adjusted p-value for the marker of interest
  #p_value_adjusted <- anova(adjusted_model)["rcs(Syndecan1, 3)","Pr(>F)"]
  
  # Extract the adjusted p-value for the model without MEWS_score
  #p_value_adjusted_without_mews <- anova(adjusted_model_without_mews)["rcs(Syndecan1, 3)","Pr(>F)"]
  
  # Get significance stars for each p-value
  p_value_significance <- get_significance_stars(p_value)
  #p_value_adjusted_significance <- get_significance_stars(p_value_adjusted)
  #p_value_adjusted_without_mews_significance <- get_significance_stars(p_value_adjusted_without_mews)
  
  # Store the results
  results_2 <- rbind(results_2, data.frame(Marker = marker, 
                                           Model = chosen_model, 
                                           P_Value = p_value, 
                                           P_Value_Significance = p_value_significance, 
                                           #P_Value_Adjusted = p_value_adjusted, 
                                           #P_Value_Adjusted_Significance = p_value_adjusted_significance, 
                                           #P_Value_adjusted_without_Mews = p_value_adjusted_without_mews, 
                                           #P_Value_adjusted_without_Mews_Significance = p_value_adjusted_without_mews_significance, 
                                           stringsAsFactors = FALSE))
}


# Write the data frame to a CSV file
destination_folder <- "original_data" 
export_file_name <- "linear_nonlinear_p_value.csv" 
write.csv(results_2, file.path(destination_folder, export_file_name), row.names = FALSE)
cat("File exported to:", file.path(destination_folder, export_file_name), "\n")
