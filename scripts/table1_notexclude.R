#use the sofa score from Hessel

rm(list = ls())

## install.packages("pacman")
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)
install.packages("dplyr")
install.packages("FSA")
install.packages("rms")
library(dplyr)
library(mice)
library(tidyr)
library(FSA)
library(rms)

# Read data
data <- read.csv("original_data/final_data.csv")

# Cutting the sorted data into three equal parts (1, 2, 3)
data$Syndecan_group <- cut(data$Syndecan1, breaks = quantile(data$Syndecan1, probs = seq(0, 1, length = 4), na.rm = TRUE),
                                             include.lowest = TRUE,
                                             labels = c(1, 2, 3))

#table
#pick column what i want
# Load required packages
pacman::p_load(tableone, dplyr)

# Create analysis dataset
data_table1 <- data[, c(
  "MARSID", "Syndecan1", "age_yrs", "gender", "BMI",
  "COPD", "Chronic_heart_failure", "Past_myocardial_infarction", "Cerebrovascular_disease",
  "mneoplasm", "ckd", "diabetes", "ICUA_no_complication",
  "Platelets_value_1", "WBC_2_1", "Creatinine_value_1", "Blood_Urea_Nitrogen_value_1", "Bilirubin",
  "SOFAtot_new", "MEWS_score", "APACHE_IV_Score",
  "ICUA_ARDS", "AKI_24h", "Shock_24h",
  "length_of_stay", "Death_in_hospital.x", "death_in_ICU_correct", "mort30", "mort90",
  "FinalGroup_pathegon", "FinalGroup_gram","Syndecan_group","Immune_deficiency"
)]

# Define variable list for Table 1
allvars <- c(
  # 0.syndecan-1
  "Syndecan1", 
  # 1. Demographic characteristics
  "age_yrs", "gender", "BMI",
  # 2. Comorbidities
  "COPD", "Chronic_heart_failure", "Past_myocardial_infarction", "Cerebrovascular_disease",
  "mneoplasm", "Immune_deficiency", "ckd", "diabetes", "ICUA_no_complication",
  # 3. Laboratory parameters
 "Platelets_value_1", "WBC_2_1",
  "Creatinine_value_1", "Blood_Urea_Nitrogen_value_1", "Bilirubin",
  # 4. Clinical severity scores
  "SOFAtot_new", "MEWS_score", "APACHE_IV_Score",
  # 5. Organ dysfunction indicators
  "ICUA_ARDS", "AKI_24h", "Shock_24h",
  # 6. Clinical outcomes
  "length_of_stay", "Death_in_hospital.x", "death_in_ICU_correct", "mort30", "mort90",
  # 7. Pathogen-related characteristics
  "FinalGroup_pathegon", "FinalGroup_gram",
  "Syndecan_group"
)

# Define categorical variables
catvars <- c(
  "gender", "COPD", "Chronic_heart_failure", "Past_myocardial_infarction", 
  "Cerebrovascular_disease", "mneoplasm","Immune_deficiency","ckd", "diabetes", "ICUA_no_complication",
  "ICUA_ARDS", "AKI_24h", "Shock_24h",
  "Death_in_hospital.x", "death_in_ICU_correct", "mort30", "mort90",
  "FinalGroup_pathegon", "FinalGroup_gram","Syndecan_group"
)

# Define non-normally distributed continuous variables
nonnormal <- c(
  "Syndecan1", "age_yrs", "BMI", "Platelets_value_1", "WBC_2_1", 
  "Creatinine_value_1", "Blood_Urea_Nitrogen_value_1", "Bilirubin",
  "SOFAtot_new", "MEWS_score", "APACHE_IV_Score", "length_of_stay"
)

# Generate Table 1
tab1 <- CreateTableOne(
  vars       = allvars,
  data       = data_table1,
  #strata = "Syndecan_group",
  factorVars = catvars,
  test       = TRUE
)

# Print Table 1
print(tab1, 
      nonnormal = nonnormal, 
      quote = FALSE, 
      noSpaces = TRUE, 
      smd = TRUE, 
      missing = TRUE)



# Convert the CreateTableOne object to a data frame

tab1_df <- as.data.frame(print(tab1, nonnormal = nonnormal, quote = FALSE, noSpaces = TRUE, smd = TRUE, missing = TRUE))
tab1_df$tablename <- rownames(tab1_df) 
tab1_df <- tab1_df[, c("tablename", setdiff(names(tab1_df), "tablename"))]


# Write the data frame to a CSV file
destination_folder <- "original_data" 
export_file_name <- "table1_syndecan_continous.csv" 
write.csv(tab1_df, file.path(destination_folder, export_file_name), row.names = FALSE)
cat("File exported to:", file.path(destination_folder, export_file_name), "\n")

