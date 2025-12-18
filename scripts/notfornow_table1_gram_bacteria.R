rm(list = ls())

## install.packages("pacman")
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)
library(dplyr)
library(FSA)
library(rms)

# Read data
data <- read.csv("Documents/BSI/R_code/original_data/clinical_marker_scource_pathogen.csv")
data_bacteremia <- data %>% 
  filter(Microbe_groups == "Bacteremia")

# Write the data frame to a CSV file
destination_folder <- "Documents/BSI/R_code/original_data" 
export_file_name <- "data_bacteremia.csv" 
write.csv(data_bacteremia, file.path(destination_folder, export_file_name), row.names = FALSE)
cat("File exported to:", file.path(destination_folder, export_file_name), "\n")



#make syndecan groups
data_bacteremia$Syndecan_group <- cut(data_bacteremia$`Syndecan.1.CD138..29.`,breaks = quantile(data_bacteremia$`Syndecan.1.CD138..29.`, probs = seq(0, 1, length = 4), na.rm = TRUE),
                                      include.lowest = TRUE,labels = c(1, 2, 3))
#ratio
#data_bacteremia$Syndecan_group <- cut(data_bacteremia$`ANG2_1`,breaks = quantile(data_bacteremia$`ANG2_1`, probs = seq(0, 1, length = 4), na.rm = TRUE),
#                                     include.lowest = TRUE,labels = c(1, 2, 3))



#Group for Gram####
data_bacteremia$Gram_group <- with(data_bacteremia, ifelse(GN == 0 & GP == 1, "Gram_positive",
                                     ifelse(GN == 1 & GP == 0, "Gram_negative",
                                            ifelse(GN == 1 & GP == 1, "Both", "Not"))))

#table
#pick column what i want
data_table1_gram_syndecangroup <- data_bacteremia[,c("MARSID","Syndecan.1.CD138..29.","age_yrs","gender","BMI","Microbe_groups","COPD","Chronic_heart_failure","Past_myocardial_infarction","Cerebrovascular_disease","mneoplasm","ckd","diabetes","Platelets_value_1","WBC_2_1","Creatinine_value_1","Blood_Urea_Nitrogen_value_1","Bilirubin","SOFAtot","ICUA_no_complication","ICUA_infection","MEWS_score","Charlson_with_age","APACHE_IV_Score.x","length_of_stay",
                              "mort30","mort90","Shock_on_admis","ICUA_no_complication","ICUA_infection","Death_in_hospital.x","death_in_ICU_correct","FinalGroup_pathegon","Klebsiella", "S.aureus", "Enterococcus", "E.coli", "Streptococcus", "Pseudomonas", "Enterobacter", "Gram_group",
                              "CX3CL1.Fractalkine..46.","ESM.1..22.","Ang.1..64.","Ang.2..26.","ANG2_1","Syndecan_group","respiratory", "abdominal", "cardiovascular", "urinary", "CNS_source", "skin", "other_source", "unknown_source")]

##  
allvars <- c("age_yrs","gender","BMI",
             "COPD","Chronic_heart_failure","Past_myocardial_infarction","Cerebrovascular_disease","mneoplasm","ckd","diabetes","ICUA_no_complication",
             "Syndecan.1.CD138..29.","CX3CL1.Fractalkine..46.","ESM.1..22.","Ang.1..64.","Ang.2..26.","ANG2_1",
             "Platelets_value_1","WBC_2_1","Creatinine_value_1","Blood_Urea_Nitrogen_value_1","Bilirubin",
             "ICUA_infection","Microbe_groups","FinalGroup_pathegon","Klebsiella", "S.aureus", "Enterococcus", "E.coli", "Streptococcus", "Pseudomonas", "Enterobacter",
             "SOFAtot","MEWS_score","APACHE_IV_Score.x","Shock_on_admis",
             "length_of_stay","Death_in_hospital.x","death_in_ICU_correct","mort30","mort90","Gram_group","Syndecan_group","respiratory", "abdominal", "cardiovascular", "urinary", "CNS_source", "skin", "other_source", "unknown_source")

catvars <- c("gender","COPD",	"Chronic_heart_failure","Past_myocardial_infarction","Cerebrovascular_disease",	"mneoplasm","ckd","diabetes","ICUA_no_complication","ICUA_infection","Gram_group","Syndecan_group","Shock_on_admis","Death_in_hospital.x","death_in_ICU_correct","mort30","mort90","FinalGroup_pathegon","Klebsiella", "S.aureus", "Enterococcus", "E.coli", "Streptococcus", "Pseudomonas", "Enterobacter","respiratory", "abdominal", "cardiovascular", "urinary", "CNS_source", "skin", "other_source", "unknown_source")

nonnormal <- c("Syndecan.1.CD138..29.","CX3CL1.Fractalkine..46.","ESM.1..22.","Ang.1..64.","Ang.2..26.","ANG2_1","age_yrs","BMI","Platelets_value_1", "WBC_2_1","Creatinine_value_1", "Blood_Urea_Nitrogen_value_1","Bilirubin","SOFAtot","MEWS_score","APACHE_IV_Score.x","length_of_stay")

##table for gram
tab1     <- CreateTableOne(
  vars        = allvars, 
  data        = data_table1_gram_syndecangroup, 
  strata = "Gram_group",
  factorVars  = catvars,
  test        = TRUE)
print(tab1, nonnormal = nonnormal, quote = FALSE, noSpaces = TRUE, smd = T, missing = T)

# Convert the CreateTableOne object to a data frame

tab1_df <- as.data.frame(print(tab1, nonnormal = nonnormal, quote = FALSE, noSpaces = TRUE, smd = TRUE, missing = TRUE))
tab1_df$tablename <- rownames(tab1_df) 
tab1_df <- tab1_df[, c("tablename", setdiff(names(tab1_df), "tablename"))]

# Write the data frame to a CSV file
destination_folder <- "Documents/BSI/R_code/original_data" 
export_file_name <- "table1_Gram_positive_negative_both.csv" 
write.csv(tab1_df, file.path(destination_folder, export_file_name), row.names = FALSE)
cat("File exported to:", file.path(destination_folder, export_file_name), "\n")


###table for pathegon
tab2     <- CreateTableOne(
  vars        = allvars, 
  data        = data_table1_gram_syndecangroup, 
  strata = "Syndecan_group",
  factorVars  = catvars,
  test        = TRUE)
print(tab2, nonnormal = nonnormal, quote = FALSE, noSpaces = TRUE, smd = T, missing = T)

# Convert the CreateTableOne object to a data frame

tab2_df <- as.data.frame(print(tab2, nonnormal = nonnormal, quote = FALSE, noSpaces = TRUE, smd = TRUE, missing = TRUE))
tab2_df$tablename <- rownames(tab2_df) 
tab2_df <- tab2_df[, c("tablename", setdiff(names(tab2_df), "tablename"))]

# Write the data frame to a CSV file
destination_folder <- "Documents/BSI/R_code/original_data" 
export_file_name <- "table1_syndecan_bacteria.csv" 
write.csv(tab2_df, file.path(destination_folder, export_file_name), row.names = FALSE)
cat("File exported to:", file.path(destination_folder, export_file_name), "\n")


#exclude other_bacteria == 0
#data_table1_gram_syndecangroup_non_other_bacteria <- filter(data_table1_gram_syndecangroup,other_bacteria == 0 )
#tab3     <- CreateTableOne(vars        = allvars,  data        = data_table1_gram_syndecangroup_non_other_bacteria, strata = "Syndecan_group",factorVars  = catvars,test        = TRUE)
#print(tab3, nonnormal = nonnormal, quote = FALSE, noSpaces = TRUE, smd = T, missing = T)


#correct for severity
#checked if linearity with sofa/APACHE_IV_Score.x
linear_model_gram <- lm(log(Syndecan.1.CD138..29.) ~ SOFAtot, data = data_table1_gram)
spline_model_gram <- lm(log(Syndecan.1.CD138..29.) ~ rcs(SOFAtot, 3), data = data_table1_gram)

# Compare models using ANOVA
model_comparison_mews <- anova(linear_model_gram, spline_model_gram)
print("Model comparison for MEWS:")
print(model_comparison_mews)
#Syndecan.1.CD138..29. ,CX3CL1.Fractalkine..46.,Cystatin.C..75.,ESM.1..22., Ang.1..64.,Ang.2..26.,ANG_ratio,and APACHE_IV_Score.x/SOFAtot, linearity is better.

#the relationship is linearity. 
# Linear model with log(Syndecan) as the outcome and group, SOFAtot, and APACHE_IV_Score.x as predictors
linear_model_correct_gram <- lm(log(Syndecan.1.CD138..29.) ~ SOFAtot + APACHE_IV_Score.x + Gram_group, data = data_table1_gram_syndecangroup)
# Print the summary of the model
anova(linear_model_correct_gram)
#####

#correct for bacteria
correct_bacteria <- glm(Enterobacter ~ SOFAtot + APACHE_IV_Score.x + Syndecan.1.CD138..29., data = data_table1_gram_syndecangroup, family = binomial)
anova(correct_bacteria)

#Syndecan_group
#install.packages("nnet")
#library(nnet)
#bacteria_correct <- multinom(Pseudomonas~ SOFAtot + APACHE_IV_Score.x + Syndecan.1.CD138..29., data = data_table1_gram_syndecangroup )

#table for syndecan tertile
#make the syndecan-1 group
