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
data <- read.csv("original_data/clinical_marker_scource_pathogen.csv")


#delete the outliers for endothelial 
#4298 14507  1535
data <- data %>% 
  filter(!ICU_ID_from_datasource %in% c(273, 4570, 1038))

#delete the patients Gram_negative and Gram_positive at the same time
data_2 <- data
data_2 <- data_2[!(data_2$Gram_negative == 1 & data_2$Gram_positive == 1), ]

categories <- c("Gram_negative", "Gram_positive")
category_var <- factor(NA, levels = categories)
category_var[data_2$Gram_negative == 1] <- "Gram_negative"
category_var[data_2$Gram_positive == 1] <- "Gram_positive"

data_2$Gram_group <- category_var
data_2<- data_2 %>% filter(!is.na(Gram_group))


####7 bacteria group#####
data_bacteria <- data_2 %>%
  filter(!(microbenum %in% c(NA, 2, 3, 4)))

# Create category variable
categories_2 <- c("Klebsiella", "S.aureus", "CoNS", "Staph_unknown", "Enterococcus", "E.coli", "Streptococcus")
category_var_2 <- factor(NA, levels = categories_2)
category_var_2[data_bacteria$Klebsiella == 1] <- "Klebsiella"
category_var_2[data_bacteria$S.aureus == 1] <- "S.aureus"
category_var_2[data_bacteria$CoNS == 1] <- "CoNS"
category_var_2[data_bacteria$Staph_unknown == 1] <- "Staph_unknown"
category_var_2[data_bacteria$Enterococcus == 1] <- "Enterococcus"
category_var_2[data_bacteria$E.coli == 1] <- "E.coli"
category_var_2[data_bacteria$Streptococcus == 1] <- "Streptococcus"

data_bacteria$seven_group <- category_var_2
data_bacteria<- data_bacteria %>% filter(!is.na(seven_group))


# Cutting the sorted data into three equal parts (1, 2, 3)
data_bacteria$Syndecan_group_bacteria <- cut(data_bacteria$`ANG_ratio`,
                                             breaks = quantile(data_bacteria$`Syndecan.1.CD138..29.`, probs = seq(0, 1, length = 4), na.rm = TRUE),
                                             include.lowest = TRUE,
                                             labels = c(1, 2, 3))



#table
#pick column what i want
data_table1_seven <- data_bacteria[,c("ICU_ID_from_datasource","Syndecan.1.CD138..29.","age_yrs","gender","BMI","Microbe_groups","COPD","Chronic_heart_failure","Past_myocardial_infarction","Cerebrovascular_disease","mneoplasm","ckd","diabetes","Platelets_value_1","WBC_2_1","Creatinine_value_1","Blood_Urea_Nitrogen_value_1","Bilirubin","SOFAtot","ICUA_no_complication","ICUA_infection","MEWS_score","Charlson_with_age","APACHE_IV_Score.x","length_of_stay",
                                      "mort30","mort90","Shock_on_admis","ICUA_no_complication","ICUA_infection","Death_in_hospital.x","death_in_ICU_correct","Klebsiella","S.aureus","CoNS","Staph_unknown","Enterococcus","E.coli","Streptococcus","microbenum","Syndecan_group","Syndecan_group_bacteria",
                                      "CX3CL1.Fractalkine..46.","Cystatin.C..75.","ESM.1..22.","Ang.1..64.","Ang.2..26.","ANG_ratio","Gram_group")]
##  
allvars <- c("age_yrs","gender","BMI",
             "COPD","Chronic_heart_failure","Past_myocardial_infarction","Cerebrovascular_disease","mneoplasm","ckd","diabetes","ICUA_no_complication",
             "Syndecan.1.CD138..29.","CX3CL1.Fractalkine..46.","Cystatin.C..75.","ESM.1..22.","Ang.1..64.","Ang.2..26.","ANG_ratio",
             "Platelets_value_1","WBC_2_1","Creatinine_value_1","Blood_Urea_Nitrogen_value_1","Bilirubin",
             "ICUA_infection","Microbe_groups","Klebsiella","S.aureus","CoNS","Staph_unknown","Enterococcus","E.coli","Streptococcus",
             "SOFAtot","MEWS_score","APACHE_IV_Score.x","Shock_on_admis",
             "length_of_stay","Death_in_hospital.x","death_in_ICU_correct","mort30","mort90","Syndecan_group_bacteria","Gram_group")

catvars <- c("gender","COPD",	"Chronic_heart_failure","Past_myocardial_infarction","Cerebrovascular_disease",	"mneoplasm","ckd","diabetes","ICUA_no_complication","ICUA_infection","Syndecan_group_bacteria","Shock_on_admis","Death_in_hospital.x","death_in_ICU_correct","mort30","mort90","Klebsiella","S.aureus","CoNS","Staph_unknown","Enterococcus","E.coli","Streptococcus")

nonnormal <- c("Syndecan.1.CD138..29.","CX3CL1.Fractalkine..46.","Cystatin.C..75.","ESM.1..22.","Ang.1..64.","Ang.2..26.","ANG_ratio","age_yrs","BMI","Platelets_value_1", "WBC_2_1","Creatinine_value_1", "Blood_Urea_Nitrogen_value_1","Bilirubin","SOFAtot","MEWS_score","APACHE_IV_Score.x","length_of_stay")




tab1     <- CreateTableOne(
  vars        = allvars, 
  data        = data_table1_seven, 
  strata = "Syndecan_group_bacteria",
  factorVars  = catvars,
  test        = TRUE)
print(tab1, nonnormal = nonnormal, quote = FALSE, noSpaces = TRUE, smd = T, missing = T)

# Convert the CreateTableOne object to a data frame

tab1_df <- as.data.frame(print(tab1, nonnormal = nonnormal, quote = FALSE, noSpaces = TRUE, smd = TRUE, missing = TRUE))
tab1_df$tablename <- rownames(tab1_df) 
tab1_df <- tab1_df[, c("tablename", setdiff(names(tab1_df), "tablename"))]

# Write the data frame to a CSV file
#export clinical_marker_unique data
destination_folder <- "Documents/BSI/R_code/original_data" 
export_file_name <- "table1_Stratified_by_syndecan_group.csv" 
write.csv(tab1_df, file.path(destination_folder, export_file_name), row.names = FALSE)
cat("File exported to:", file.path(destination_folder, export_file_name), "\n")



#correct for severity
#####corrected for severity####
#checked if linearity with sofa/APACHE_IV_Score.x
linear_model_seven <- lm(log(CX3CL1.Fractalkine..46.) ~ SOFAtot, data = data_bacteria)
spline_model_seven <- lm(log(CX3CL1.Fractalkine..46.) ~ rcs(SOFAtot, 3), data = data_bacteria)

# Compare models using ANOVA
model_comparison_mews <- anova(linear_model_seven, spline_model_seven)
print("Model comparison for MEWS:")
print(model_comparison_mews)
#Syndecan.1.CD138..29. ,CX3CL1.Fractalkine..46.,Cystatin.C..75.,ESM.1..22., Ang.1..64.,Ang.2..26.,ANG_ratio,and APACHE_IV_Score.x/SOFAtot, linearity is better.

#the relationship is linearity. 
# Linear model with log(Syndecan) as the outcome and group, SOFAtot, and APACHE_IV_Score.x as predictors
linear_model_correct_seven <- lm(log(ANG_ratio) ~ SOFAtot + APACHE_IV_Score.x + Syndecan_group_bacteria, data = data_bacteria)
# Print the summary of the model
anova(linear_model_correct_seven)