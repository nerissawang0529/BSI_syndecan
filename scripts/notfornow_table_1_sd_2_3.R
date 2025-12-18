rm(list = ls())

# install.packages("pacman")
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)


# Read data
data <- read.csv("original_data/merged_data_new_group.csv")

#table
#pick column what i want
data_table1 <- data[,c("ICU_ID_from_datasource","syndecan_original","age_yrs","gender","BMI","Microbe_groups","COPD","Chronic_heart_failure","Past_myocardial_infarction","Cerebrovascular_disease","mneoplasm","ckd","diabetes","Platelets_value_1","WBC_2_1","Creatinine_value_1","Blood_Urea_Nitrogen_value_1","Bilirubin","SOFAtot","ICUA_no_complication","ICUA_infection","MEWS_score","Charlson_with_age","APACHE_IV_Score.x","length_of_stay",
                       "mort30","mort90","Shock_on_admis","ICUA_no_complication","ICUA_infection","Death_in_hospital.x","death_in_ICU_correct","FinalGroup_pathegon","sd_group_based_healthy","sd_group_based_healthy_3",
                       "FinalGroup_gram")]
##  
allvars <- c("age_yrs","gender","BMI",
             "COPD","Chronic_heart_failure","Past_myocardial_infarction","Cerebrovascular_disease","mneoplasm","ckd","diabetes","ICUA_no_complication",
             "syndecan_original","Platelets_value_1","WBC_2_1","Creatinine_value_1","Blood_Urea_Nitrogen_value_1","Bilirubin",
             "ICUA_infection","FinalGroup_pathegon",
             "SOFAtot","MEWS_score","APACHE_IV_Score.x","Shock_on_admis",
             "length_of_stay","Death_in_hospital.x","death_in_ICU_correct","mort30","mort90","sd_group_based_healthy","sd_group_based_healthy_3","FinalGroup_gram")

catvars <- c("gender","COPD",	"Chronic_heart_failure","Past_myocardial_infarction","Cerebrovascular_disease",	"mneoplasm","ckd","diabetes","ICUA_no_complication","ICUA_infection","Shock_on_admis","Death_in_hospital.x","death_in_ICU_correct","mort30","mort90")

nonnormal <- c("syndecan_original","age_yrs","BMI","Platelets_value_1", "WBC_2_1","Creatinine_value_1", "Blood_Urea_Nitrogen_value_1","Bilirubin","SOFAtot","MEWS_score","APACHE_IV_Score.x","length_of_stay")




tab1     <- CreateTableOne(
  vars        = allvars, 
  data        = data_table1, 
  strata = "sd_group_based_healthy_3",
  factorVars  = catvars,
  test        = TRUE)
print(tab1, nonnormal = nonnormal, quote = FALSE, noSpaces = TRUE, smd = T, missing = T)

# Convert the CreateTableOne object to a data frame

tab1_df <- as.data.frame(print(tab1, nonnormal = nonnormal, quote = FALSE, noSpaces = TRUE, smd = TRUE, missing = TRUE))
tab1_df$tablename <- rownames(tab1_df) 
tab1_df <- tab1_df[, c("tablename", setdiff(names(tab1_df), "tablename"))]


# Write the data frame to a CSV file
destination_folder <- "original_data" 
export_file_name <- "table1_group_sd3.csv" 
write.csv(tab1_df, file.path(destination_folder, export_file_name), row.names = FALSE)
cat("File exported to:", file.path(destination_folder, export_file_name), "\n")
