rm(list = ls())

#the data is from Hessel

## install.packages("pacman")
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)


#data for BSI 
data <- read.csv("original_data/clinical_marker_scource.csv")
subset_data <- data[data$Microbe_groups == "Non-infectious", ]
subset_data <- subset_data %>%
  mutate(
    adm_diag_group = case_when(
      grepl("respiratory|pneumonia|asthma|emphysema|bronchitis", APACHE_IV_adm_diag_primary, ignore.case = TRUE) ~ "Respiratory failure",
      grepl("neuro|cerebrovascular|stroke|subarachnoid|intracranial|coma", APACHE_IV_adm_diag_primary, ignore.case = TRUE) ~ "Acute neurological disease (CVA, haemorrhage)",
      grepl("cardiac arrest", APACHE_IV_adm_diag_primary, ignore.case = TRUE) ~ "Cardiac arrest",
      grepl("CHF|cardio|heart failure", APACHE_IV_adm_diag_primary, ignore.case = TRUE) ~ "Cardiac failure",
      grepl("surgery", APACHE_IV_adm_diag_primary, ignore.case = TRUE) & grepl("gastro|GI", APACHE_IV_adm_diag_primary, ignore.case = TRUE) ~ "Gastrointestinal surgery",
      grepl("trauma|injury", APACHE_IV_adm_diag_primary, ignore.case = TRUE) ~ "Trauma",
      grepl("aneurysm", APACHE_IV_adm_diag_primary, ignore.case = TRUE) ~ "Dissecting aortic aneurysm",
      grepl("diabetic ketoacidosis", APACHE_IV_adm_diag_primary, ignore.case = TRUE) ~ "Diabetic ketoacidosis",
      grepl("acid", APACHE_IV_adm_diag_primary, ignore.case = TRUE) ~ "Acid-base electrolyte disturbance",
      grepl("transplant", APACHE_IV_adm_diag_primary, ignore.case = TRUE) ~ "Complications after kidney transplantation",
      grepl("overdose|toxin|poison|drug", APACHE_IV_adm_diag_primary, ignore.case = TRUE) ~ "Drug overdose",
      grepl("endarterectomy", APACHE_IV_adm_diag_primary, ignore.case = TRUE) ~ "Endarterectomy",
      grepl("hemorrhage|haemorrhage", APACHE_IV_adm_diag_primary, ignore.case = TRUE) & grepl("GI|gastro", APACHE_IV_adm_diag_primary, ignore.case = TRUE) ~ "Gastrointestinal haemorrhage",
      grepl("hepatic", APACHE_IV_adm_diag_primary, ignore.case = TRUE) ~ "Hepatic failure",
      grepl("drowning", APACHE_IV_adm_diag_primary, ignore.case = TRUE) ~ "Near drowning",
      grepl("embolus|embolism", APACHE_IV_adm_diag_primary, ignore.case = TRUE) ~ "Pulmonary embolism",
      grepl("lymph node|retroperitoneal", APACHE_IV_adm_diag_primary, ignore.case = TRUE) ~ "Retroperitoneal lymph node dissection",
      TRUE ~ "Unspecified"
    )
  )

#table
#pick column what i want
# Load required packages
pacman::p_load(tableone, dplyr)


# Define variable list for Table 1
allvars <- c(
  # 1. Demographic characteristics
  "age_yrs", "gender", "BMI",
  # 2. Comorbidities
  "COPD", "Chronic_heart_failure", "Past_myocardial_infarction", "Cerebrovascular_disease",
  "mneoplasm", "Immune_deficiency", "ckd", "diabetes", "ICUA_no_complication",
  # 3. Laboratory parameters
  "Platelets_value_1", "WBC_2_1",
  "Creatinine_value_1", "Blood_Urea_Nitrogen_value_1", "Bilirubin",
  # 4. Clinical severity scores
  "MEWS_score", "APACHE_IV_Score", "SOFAtot",
  # 5. Organ dysfunction indicators
  "ICUA_ARDS", "AKI_24h", "Shock_24h",
  # 6. Clinical outcomes
  "length_of_stay", "Death_in_hospital.x", "death_in_ICU_correct", "mort30", "mort90",
  "adm_diag_group"
)

# Define categorical variables
catvars <- c(
  "gender", "COPD", "Chronic_heart_failure", "Past_myocardial_infarction", 
  "Cerebrovascular_disease", "mneoplasm","Immune_deficiency","ckd", "diabetes", "ICUA_no_complication",
  "ICUA_ARDS", "AKI_24h", "Shock_24h",
  "Death_in_hospital.x", "death_in_ICU_correct", "mort30", "mort90","adm_diag_group"
)

# Define non-normally distributed continuous variables
nonnormal <- c(
  "Syndecan1", "age_yrs", "BMI", "Platelets_value_1", "WBC_2_1", 
  "Creatinine_value_1", "Blood_Urea_Nitrogen_value_1", "Bilirubin",
  "SOFAtot_new", "MEWS_score", "APACHE_IV_Score", "SOFAtot","length_of_stay"
)

# Generate Table 1
tab1 <- CreateTableOne(
  vars       = allvars,
  data       = subset_data,
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
export_file_name <- "table1_noninfection.csv" 
write.csv(tab1_df, file.path(destination_folder, export_file_name), row.names = FALSE)
cat("File exported to:", file.path(destination_folder, export_file_name), "\n")
