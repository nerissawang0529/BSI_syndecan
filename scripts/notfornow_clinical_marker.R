rm(list = ls())

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


#this is the BSI ICU_admition_data
load("Documents/BSI/R_code/original_data/df_biomarker_BSI_second_time_from_Hessel.RData")
df_biomarker_wide_erik3$ANG2_1 <- df_biomarker_wide_erik3$`Ang-2 (26)`/df_biomarker_wide_erik3$`Ang-1 (64)`

#Ready studydata(MARS database)
studydata <- import("Documents/BSI/R_code/original_data/MARS database.xlsx")

studydata$Index_ICU_admittance_datetime
double <- studydata[duplicated(studydata$ICU_ID_from_datasource), ]
check <- filter(studydata, studydata$ICU_ID_from_datasource == "13080")
good_one <- filter(check, check$hours_of_sepsis_before_admission == -18.5)
studydata <- studydata[!studydata$ICU_ID_from_datasource %in% check$ICU_ID_from_datasource, ]
studydata <- rbind(studydata, good_one)
names(studydata)[names(studydata) == "ICU_ID_from_datasource"] <- "ID"

#combine the df_biomarker_wide_erik3 and studydata
##There is no difference of APACHE_IV_Score.x/y,APACHE_IV_Score.x/y,Index_Hospital_Name.x/y,death_in_ICU.x/y

clinical_marker <- merge(df_biomarker_wide_erik3, studydata, by.x = "MARSID", by.y = "ID", all = FALSE)
clinical_marker_unique <- clinical_marker %>%
  distinct(MARSID, .keep_all = TRUE)

# 3286  1063  3613 14497, there IDs have no source information, excluded them
exclude_values <- c(3286, 1063, 3613, 14497)
clinical_marker_unique <- clinical_marker_unique[!clinical_marker_unique$MARSID %in% exclude_values, ]

# 6611 1967 2791 were outliers for markers data for both all marker and also for marker except Cystatin.C..75.
exclude_values_2 <- c(6611, 1967, 2791)
clinical_marker_unique <- clinical_marker_unique[!clinical_marker_unique$MARSID %in% exclude_values_2, ]

# 5433, the syndecan-1 value is 'NA', delete it for now
clinical_marker_unique <- clinical_marker_unique[!clinical_marker_unique$MARSID == 5433, ]

# 1906 2129 2405 2632 2914 3652  4596 10006 10369 10737 11688 13805 are delete from 'ss01 <- ss[ss$culture_day == 0 | ss$culture_day == 1 | ss$culture_day == -1, ]'
exclude_values_3 <- c(1906,2129,2405,2632,2914,3652,4596,10006,10369,10737,11688,13805)
clinical_marker_unique <- clinical_marker_unique[!clinical_marker_unique$MARSID %in% exclude_values_3, ]

# 12328, 13323, 13700, 8493, 8690 are ambiguous (vague) bacteria decriptions (n=5) in the code pathegone_exclude_CNS
exclude_values_4 <- c(12328, 13323, 13700, 8493, 8690)
clinical_marker_unique <- clinical_marker_unique[!clinical_marker_unique$MARSID %in% exclude_values_4, ]

#Exclude CoNS in the code pathegone_exclude_CNS
exclude_values_5 <- c(12034, 13104, 13501, 13744, 13906, 1535, 189, 2877, 3471, 364, 5305, 5408, 5764, 5841, 6312, 6437, 8563, 8758, 8811)
clinical_marker_unique <- clinical_marker_unique[!clinical_marker_unique$MARSID %in% exclude_values_5, ]

#Exclude likelihood none
exclude_values_6 <- c(1056,  2724 , 4908 , 6386,  6675,  8495,  8991, 11399)
clinical_marker_unique <- clinical_marker_unique[!clinical_marker_unique$MARSID %in% exclude_values_6, ]


#MEWS_score
clinical_marker_unique$MEWS_score <- NULL  #there is a original MEWS_score in the database, i deleted this one and make a new one below.
##For imputation, because when we imputate values for MEWS score, the value about outcome should not be included
##in the original dataframe there are 251 variables, it is not easy to figure out which is the outcome, so I only choose the values which will be used for the table1 (exclud the values for outcome)
clinical_marker_unique_imputation <- clinical_marker_unique[,c("Use_of_vasoactive_medic_24h", "ABPs_First_24h", "ABPs_Max_24h", 
                                     "ABPs_Min_24h", "HR_First_24h", "HR_Min_24h", 
                                     "Mechanical_ventilation_24h", "Mechanical_vent_admission_24h", 
                                     "Respiratory_rate_Max_24h", "Respiratory_rate_Min_24h", 
                                     "Core_temp_Max_24h", "Core_temp_Min_24h", "Escore_min_24h",
                                     "MARSID", "Site_of_infection_primary", 
                                     "Patient_age_at_ICU_Admission", "Gender_is_male", "Patient_BMI", "Diabetes", "Malignancy", 
                                     "Myocardial_infarction", "Cerebrovascular_disease", "Chronic_renal_insufficiency", 
                                     "Congestive_heart_failure", "Platelet_count_Min_24h", "White_cell_count_Max_24h", 
                                     "Creatinin_Max_24h", "Urea_Max_24h", "Bilirubin_Max_24h", "Mechanical_ventilation_24h", 
                                     "Mechanical_vent_admission_24h", "Respiratory_rate_Max_24h", "HR_Max_24h", "ABPs_Min_24h", 
                                     "Core_temp_Max_24h", "Escore_min_24h", "Vscore_min_24h")]    # "Index_ICU_admittance_datetime", "hours_of_sepsis_before_admission" this value not be included here

##turn the format of different values
str(clinical_marker_unique_imputation)
clinical_marker_unique_imputation <- clinical_marker_unique_imputation %>%
  mutate(
    # into num
    Escore_min_24h = as.numeric(Escore_min_24h), HR_Max_24h = as.numeric(HR_Max_24h),ABPs_Min_24h = as.numeric(ABPs_Min_24h),Core_temp_Max_24h = as.numeric(Core_temp_Max_24h),
    Escore_min_24h = as.numeric(Escore_min_24h),Vscore_min_24h = as.numeric(Vscore_min_24h),Escore_min_24h.1= as.numeric(Escore_min_24h.1),
    
    # into factor
    MARSID = as.character(MARSID),Use_of_vasoactive_medic_24h = as.factor(Use_of_vasoactive_medic_24h),
    Mechanical_ventilation_24h = as.factor(Mechanical_ventilation_24h),Mechanical_vent_admission_24h = as.factor(Mechanical_vent_admission_24h),
    Site_of_infection_primary = as.factor(Site_of_infection_primary),Gender_is_male = as.factor(Gender_is_male),
    Diabetes = as.factor(Diabetes), Malignancy = as.factor(Malignancy),
    Myocardial_infarction = as.factor(Myocardial_infarction),Cerebrovascular_disease = as.factor(Cerebrovascular_disease),
    Chronic_renal_insufficiency = as.factor(Chronic_renal_insufficiency),Congestive_heart_failure = as.factor(Congestive_heart_failure),
    Mechanical_vent_admission_24h.1 = as.factor(Mechanical_vent_admission_24h.1), Mechanical_ventilation_24h.1 = as.factor(Mechanical_ventilation_24h.1)
  )
str(clinical_marker_unique_imputation)

##check if there is NA in the values for MEWS score
na_counts <- sapply(clinical_marker_unique_imputation[, c("Use_of_vasoactive_medic_24h", "ABPs_First_24h", "ABPs_Max_24h", 
                                             "ABPs_Min_24h", "HR_First_24h", "HR_Min_24h", 
                                             "Mechanical_ventilation_24h", "Mechanical_vent_admission_24h", 
                                             "Respiratory_rate_Max_24h", "Respiratory_rate_Min_24h", 
                                             "Core_temp_Max_24h", "Core_temp_Min_24h", "Escore_min_24h")], 
                    function(x) sum(is.na(x)))

print(na_counts) #HR_Min_24h;  Respiratory_rate_Min_24h; Core_temp_Max_24h;  Core_temp_Min_24h; ABPs_First_24h,ABPs_Max_24h,ABPs_Min_24h

##imputation
set.seed(0529)
### Create a new dataframe and set the ID column to 0
clinical_marker_unique_imputation_check <- clinical_marker_unique_imputation #for check if the final imputation is right
clinical_marker_unique_imputation_MARSID_0 <- clinical_marker_unique_imputation
clinical_marker_unique_imputation_MARSID_0$MARSID <- 0 #ID can't be used in the imputation modle

### Initialize the mice object
imp <- mice(clinical_marker_unique_imputation_MARSID_0, pri=FALSE, maxit=0)
##Warning message:Number of logged events: 7 
predM <- imp$predictorMatrix
meth <- imp$method

### Specify the variables to impute
vars_to_impute <- c('HR_Min_24h', 'Respiratory_rate_Min_24h', 'Core_temp_Max_24h', 'Core_temp_Min_24h','ABPs_First_24h','ABPs_Max_24h','ABPs_Min_24h')

### Update the predictor matrix and imputation methods
predM[ , vars_to_impute] <- 1
meth[vars_to_impute] <- "cart"

### Perform the imputation process
clinical_marker_unique_imputation_imp <- mice(data=clinical_marker_unique_imputation_MARSID_0, m=20, maxit=20, meth=meth, pred=predM, printFlag=FALSE)
imputation_long <- complete(clinical_marker_unique_imputation_imp, action="long", include = FALSE)

### Add the original ID column to the long format data
imputation_long$OriginalID <- rep(clinical_marker_unique_imputation$MARSID, each = 20)

### Calculate the median imputed values for each original ID
median_imputed <- aggregate(imputation_long[, vars_to_impute], 
                            by = list(OriginalID = imputation_long$OriginalID), 
                            FUN = function(x) median(x, na.rm = TRUE))

### Sort by the original ID order
median_imputed <- median_imputed[order(match(median_imputed$OriginalID, clinical_marker_unique_imputation$MARSID)), ]
### Ensure the ID order is consistent
stopifnot(all(median_imputed$OriginalID == clinical_marker_unique_imputation$MARSID))

### Update the original dataset with the imputed median values
for (var in vars_to_impute) {
  is_na <- is.na(clinical_marker_unique_imputation[[var]])
  clinical_marker_unique_imputation[is_na, var] <- median_imputed[is_na, var]
}

### Restore the original ID column
clinical_marker_unique_imputation$MARSID <- median_imputed$OriginalID

## correct MEWS
clinical_marker_unique_imputation$MEWS_sys_bp <- with(clinical_marker_unique_imputation, ifelse(!is.na(Use_of_vasoactive_medic_24h) & Use_of_vasoactive_medic_24h == 1, 3, 
                                                                      ifelse(is.na(ABPs_First_24h), NA,
                                                                             ifelse(ABPs_First_24h <= 70, 3,
                                                                                    ifelse(ABPs_First_24h %in% 71:80, 2, 
                                                                                           ifelse(ABPs_First_24h %in% 81:100, 1,
                                                                                                  ifelse(ABPs_First_24h %in% 101:199, 0,
                                                                                                         ifelse(ABPs_First_24h >= 200, 2, NA))))))))

clinical_marker_unique_imputation$MEWS_sys_bp <-  with(clinical_marker_unique_imputation, ifelse(!is.na(ABPs_First_24h), MEWS_sys_bp,
                                                                       ifelse(!is.na(ABPs_Max_24h) & ABPs_Max_24h <= 70, 3,
                                                                              ifelse(!is.na(ABPs_Min_24h) & ABPs_Min_24h <= 70, 3,
                                                                                     ifelse(!is.na(ABPs_Min_24h) & ABPs_Min_24h %in% 71:80, 2,
                                                                                            ifelse(!is.na(ABPs_Max_24h) & ABPs_Max_24h %in% 71:80, 2,
                                                                                                   ifelse(!is.na(ABPs_Max_24h) & ABPs_Max_24h >=200, 2,
                                                                                                          ifelse(!is.na(ABPs_Min_24h) & ABPs_Min_24h >= 200, 2,
                                                                                                                 ifelse(!is.na(ABPs_Min_24h) & ABPs_Min_24h %in% 81:100, 1,
                                                                                                                        ifelse(!is.na(ABPs_Max_24h) & ABPs_Max_24h %in% 81:100, 1,
                                                                                                                               ifelse(!is.na(ABPs_Min_24h) & ABPs_Min_24h %in% 101:199 | !is.na(ABPs_Max_24h) & ABPs_Max_24h %in% 101:199, 0, MEWS_sys_bp)))))))))))

## correct
clinical_marker_unique_imputation$MEWS_hrt_rate <- with(clinical_marker_unique_imputation, ifelse(is.na(HR_First_24h), NA,
                                                                        ifelse(HR_First_24h < 40, 2,
                                                                               ifelse(HR_First_24h %in% 40:50, 1, 
                                                                                      ifelse(HR_First_24h %in% 51:100, 0,
                                                                                             ifelse(HR_First_24h %in% 101:110, 1,
                                                                                                    ifelse(HR_First_24h %in% c(111:129,116.1), 2,
                                                                                                           ifelse(HR_First_24h >=130, 3, NA))))))))
clinical_marker_unique_imputation$MEWS_hrt_rate <-  with(clinical_marker_unique_imputation, ifelse(is.na(HR_First_24h) & !is.na(HR_Min_24h) & HR_Min_24h <50, 1, 
                                                                         ifelse(is.na(HR_First_24h) & !is.na(HR_Min_24h) & HR_Min_24h <100, 0, MEWS_hrt_rate)))

## correct
clinical_marker_unique_imputation$MEWS_resp_rate <- with(clinical_marker_unique_imputation, ifelse(Mechanical_ventilation_24h == 1 | Mechanical_vent_admission_24h == 1, 3, 
                                                                         ifelse(is.na(Respiratory_rate_Max_24h) & is.na(Respiratory_rate_Min_24h), NA,
                                                                                ifelse(!is.na(Respiratory_rate_Min_24h) & Respiratory_rate_Min_24h >= 30, 3,
                                                                                       ifelse(!is.na(Respiratory_rate_Max_24h) & Respiratory_rate_Max_24h >= 30, 3,
                                                                                              ifelse(!is.na(Respiratory_rate_Min_24h) & Respiratory_rate_Min_24h < 9, 2,
                                                                                                     ifelse(!is.na(Respiratory_rate_Max_24h) & Respiratory_rate_Max_24h %in% 21:29, 2,
                                                                                                            ifelse(!is.na(Respiratory_rate_Max_24h) & Respiratory_rate_Max_24h %in% 15:20, 1,
                                                                                                                   ifelse(!is.na(Respiratory_rate_Min_24h) & Respiratory_rate_Min_24h %in% 15:20, 1,
                                                                                                                          ifelse(!is.na(Respiratory_rate_Min_24h) & Respiratory_rate_Min_24h %in% 9:14 | !is.na(Respiratory_rate_Max_24h) & Respiratory_rate_Max_24h %in% 9:14, 0, NA))))))))))

## correct all 0 are within 35 and 38.4
clinical_marker_unique_imputation$MEWS_temp <- with(clinical_marker_unique_imputation,ifelse(is.na(Core_temp_Max_24h) & is.na(Core_temp_Min_24h), NA,
                                                                   ifelse(!is.na(Core_temp_Min_24h) & Core_temp_Min_24h < 35, 2, 
                                                                          ifelse(!is.na(Core_temp_Max_24h) & Core_temp_Max_24h < 35, 2, 
                                                                                 ifelse(!is.na(Core_temp_Max_24h) & Core_temp_Max_24h >= 38.5, 2, 
                                                                                        ifelse(!is.na(Core_temp_Min_24h) & Core_temp_Min_24h >= 38.5, 2, 0))))))


## correct
clinical_marker_unique_imputation$MEWS_conciousness <- with(clinical_marker_unique_imputation, ifelse(Mechanical_ventilation_24h == 1 | Mechanical_vent_admission_24h == 1, 3, 
                                                                            ifelse(is.na(Escore_min_24h), NA,
                                                                                   ifelse(Escore_min_24h == 4, 0,
                                                                                          ifelse(Escore_min_24h == 3, 1,
                                                                                                 ifelse(Escore_min_24h == 2, 2,
                                                                                                        ifelse(Escore_min_24h == 1, 3, NA)))))))


clinical_marker_unique_imputation$MEWS_score <- with(clinical_marker_unique_imputation, MEWS_resp_rate + 
                                          MEWS_hrt_rate + 
                                          MEWS_sys_bp + MEWS_conciousness + MEWS_temp)

selected_data <- clinical_marker_unique_imputation[, c("MARSID", "MEWS_score")]

clinical_marker_unique <- merge(clinical_marker_unique, selected_data, by = "MARSID", all.x = TRUE)

#Demographics
##age
names(clinical_marker_unique)[names(clinical_marker_unique) == "Patient_age_at_ICU_Admission"] <- "age_yrs"

##gender
colnames(clinical_marker_unique)[colnames(clinical_marker_unique) == "Gender_is_male"] <- "gender"
clinical_marker_unique$gender <- ifelse(clinical_marker_unique$gender == 1, "Male", "Female")

##Patient_BMI
names(clinical_marker_unique)[names(clinical_marker_unique) == "Patient_BMI"] <- "BMI" 


#Comorbidities
##COPD
clinical_marker_unique$COPD

##diabetes
names(clinical_marker_unique)[names(clinical_marker_unique) == "Diabetes"] <- "diabetes"

##malignancy
names(clinical_marker_unique)[names(clinical_marker_unique) == "Malignancy"] <- "mneoplasm" 

##Past_myocardial_infarction
clinical_marker_unique$Myocardial_infarction
names(clinical_marker_unique)[names(clinical_marker_unique) == "Myocardial_infarction"] <- "Past_myocardial_infarction"

#Cerebrovascular_disease
clinical_marker_unique$Cerebrovascular_disease

##ckd
clinical_marker_unique$Chronic_renal_insufficiency
names(clinical_marker_unique)[names(clinical_marker_unique) == "Chronic_renal_insufficiency"] <- "ckd"

##Chronic_heart_failure
clinical_marker_unique$Congestive_heart_failure
names(clinical_marker_unique)[names(clinical_marker_unique) == "Congestive_heart_failure"] <- "Chronic_heart_failure"

#lab values
names(clinical_marker_unique)[names(clinical_marker_unique) == "Platelet_count_Min_24h"] <- "Platelets_value_1"
names(clinical_marker_unique)[names(clinical_marker_unique) == "White_cell_count_Max_24h"] <- "WBC_2_1"
names(clinical_marker_unique)[names(clinical_marker_unique) == "Creatinin_Max_24h"] <- "Creatinine_value_1"
names(clinical_marker_unique)[names(clinical_marker_unique) == "Urea_Max_24h"] <- "Blood_Urea_Nitrogen_value_1"
names(clinical_marker_unique)[names(clinical_marker_unique) == "Bilirubin_Max_24h"] <- "Bilirubin"

#Clinical course
##length_of_stay
names(clinical_marker_unique)[names(clinical_marker_unique) == "Length_of_Hospital_stay"] <- "length_of_stay" 


##mortality_d30
names(clinical_marker_unique)[names(clinical_marker_unique) == "mortality_30_day_correct"] <- "mortality_d30"



#for KM line####
##death_date, admission_date
names(clinical_marker_unique)[names(clinical_marker_unique) == "Date_of_death"] <- "death_date"
clinical_marker_unique$death_date <- as.Date(clinical_marker_unique$death_date, format = "%Y-%m-%d")
names(clinical_marker_unique)[names(clinical_marker_unique) == "Index_ICU_admittance_datetime"] <- "admission_date"
clinical_marker_unique$admission_date <- as.Date(clinical_marker_unique$admission_date, format = "%Y-%m-%d")

clinical_marker_unique$admission_date
clinical_marker_unique$Date_last_confirmed_alive <- as.Date(clinical_marker_unique$Date_last_confirmed_alive,  format = "%Y-%m-%d")

clinical_marker_unique$death_date
clinical_marker_unique$admission_date
clinical_marker_unique$Days_survived_confirmed
clinical_marker_unique$days_survived <- difftime(clinical_marker_unique$death_date, clinical_marker_unique$admission_date, units = "days")
clinical_marker_unique$Days_survived_known <- difftime(clinical_marker_unique$Date_last_confirmed_alive, clinical_marker_unique$admission_date, units = "days")


## mortality with calculated using date last confirmed alive
## note that MARS was linked with CBS so if date last known alive is missing, patient was alive.
## Indeed these 6 patients were discharge with reason was ready for discharge + discharged to ward same hospital + discharged home afterwards
x <- filter(clinical_marker_unique, is.na(clinical_marker_unique$mortality_d30))
x$discharge_reason
x$Discharge_location
x$length_of_stay

## KM
clinical_marker_unique$time <- ifelse((!is.na(clinical_marker_unique$days_survived) & clinical_marker_unique$days_survived <= 30), clinical_marker_unique$days_survived,
                         ifelse((!is.na(clinical_marker_unique$days_survived) & clinical_marker_unique$days_survived >= 30), 31, NA))
clinical_marker_unique$time <- ifelse(is.na(clinical_marker_unique$time) & clinical_marker_unique$Days_survived_confirmed >= 30, 31, clinical_marker_unique$time)
clinical_marker_unique$time <- ifelse(clinical_marker_unique$mortality_d30 == 0, 31, clinical_marker_unique$time)
clinical_marker_unique$event <- clinical_marker_unique$mortality_d30

## censoring
clinical_marker_unique$time <- ifelse(is.na(clinical_marker_unique$time), clinical_marker_unique$Days_survived_confirmed, clinical_marker_unique$time)
clinical_marker_unique$event <- ifelse(is.na(clinical_marker_unique$mortality_d30), 0, clinical_marker_unique$event)
#####

#combine the pathegon
#pathegon_result <- read.csv("Documents/BSI/R_code/original_data/pathegon_result.csv")
#clinical_marker_unique <- merge(pathegon_result, clinical_marker_unique,by.x  = "ICU_ID_from_datasource",by.y = "MARSID",all = TRUE) 

#order by syndecan-1
#clinical_marker_unique <- clinical_marker_unique[order(clinical_marker_unique$`Syndecan-1/CD138 (29)`), ]

# Cutting the sorted data into three equal parts (1, 2, 3)
#clinical_marker_unique$Syndecan_group <- cut(clinical_marker_unique$`Syndecan-1/CD138 (29)`,breaks = quantile(clinical_marker_unique$`Syndecan-1/CD138 (29)`, probs = seq(0, 1, length = 4), na.rm = TRUE),
#                                                    include.lowest = TRUE,labels = c(1, 2, 3))


#export clinical_marker_unique data
destination_folder <- "Documents/BSI/R_code/original_data" 
export_file_name <- "clinical_marker_unique.csv" 
write.csv(clinical_marker_unique, file.path(destination_folder, export_file_name), row.names = FALSE)
cat("File exported to:", file.path(destination_folder, export_file_name), "\n")


#log marker data
clinical_marker_unique_log <- clinical_marker_unique %>%
  mutate_at(vars(24:44,50), ~ log(.))

#export clinical_marker_unique data
destination_folder <- "Documents/BSI/R_code/original_data" 
export_file_name <- "clinical_marker_unique_log.csv" 
write.csv(clinical_marker_unique_log, file.path(destination_folder, export_file_name), row.names = FALSE)
cat("File exported to:", file.path(destination_folder, export_file_name), "\n")


#pick column what i want
#clinical_marker_unique_table1 <- clinical_marker_unique[,c("ICU_ID_from_datasource","Syndecan-1/CD138 (29)","age_yrs","gender","BMI","Microbe_groups","COPD","Chronic_heart_failure","Past_myocardial_infarction","Cerebrovascular_disease","mneoplasm","ckd","diabetes","Platelets_value_1","WBC_2_1","Creatinine_value_1","Blood_Urea_Nitrogen_value_1","Bilirubin","SOFAtot","ICUA_no_complication","ICUA_infection","MEWS_score","Charlson_with_age","APACHE_IV_Score.x","length_of_stay",
#                              "mort30","mort90","Shock_on_admis","ICUA_no_complication","ICUA_infection","Death_in_hospital.x","death_in_ICU_correct","Klebsiella","S.aureus","CoNS","Staph_unknown","Enterococcus","E.coli","Streptococcus","microbenum","Syndecan_group")]


##  
#allvars <- c("Syndecan-1/CD138 (29)","age_yrs","gender","BMI",
#             "COPD","Chronic_heart_failure","Past_myocardial_infarction","Cerebrovascular_disease","mneoplasm","ckd","diabetes","ICUA_no_complication",
#             "Platelets_value_1","WBC_2_1","Creatinine_value_1","Blood_Urea_Nitrogen_value_1","Bilirubin",
#             "ICUA_infection","Microbe_groups","Klebsiella","S.aureus","CoNS","Staph_unknown","Enterococcus","E.coli","Streptococcus",
#             "SOFAtot","MEWS_score","APACHE_IV_Score.x","Shock_on_admis",
#             "length_of_stay","Death_in_hospital.x","death_in_ICU_correct","mort30","mort90")

#catvars <- c("gender","COPD",	"Chronic_heart_failure","Past_myocardial_infarction","Cerebrovascular_disease",	"mneoplasm","ckd","diabetes","ICUA_no_complication","ICUA_infection","Microbe_groups","Shock_on_admis","Death_in_hospital.x","death_in_ICU_correct","mort30","mort90","Klebsiella","S.aureus","CoNS","Staph_unknown","Enterococcus","E.coli","Streptococcus")

#nonnormal <- c("Syndecan-1/CD138 (29)","age_yrs","BMI","Platelets_value_1", "WBC_2_1","Creatinine_value_1", "Blood_Urea_Nitrogen_value_1","Bilirubin","SOFAtot","MEWS_score","APACHE_IV_Score.x","length_of_stay")

#tab1     <- CreateTableOne(vars = allvars, data = clinical_marker_unique_table1, strata = "Syndecan_group",factorVars  = catvars,test        = TRUE)
#print(tab1, nonnormal = nonnormal, quote = FALSE, noSpaces = TRUE, smd = T, missing = T)


#Grouped by pathegon
#

#Klebsiella_syndecan <- clinical_marker_unique_table1 %>% filter(Klebsiella == 1) %>% .$`Syndecan-1/CD138 (29)` 
#S.aureus_syndecan <- clinical_marker_unique_table1 %>% filter(S.aureus == 1) %>% .$`Syndecan-1/CD138 (29)` 
#CoNS_syndecan <- clinical_marker_unique_table1 %>% filter(CoNS == 1) %>% .$`Syndecan-1/CD138 (29)` 
#Staph_unknown_syndecan <- clinical_marker_unique_table1 %>% filter(Staph_unknown == 1) %>% .$`Syndecan-1/CD138 (29)` 
#Enterococcus_syndecan <- clinical_marker_unique_table1 %>% filter(Enterococcus == 1) %>% .$`Syndecan-1/CD138 (29)` 
#E.coli_syndecan <- clinical_marker_unique_table1 %>% filter(E.coli == 1) %>% .$`Syndecan-1/CD138 (29)` 
#Streptococcus_syndecan<- clinical_marker_unique_table1 %>% filter(Streptococcus == 1) %>% .$`Syndecan-1/CD138 (29)` 


#gathered_data <-data.frame(pathegon_type = c(rep("KL", length(Klebsiella_syndecan)), 
#                                              rep("SA", length(S.aureus_syndecan)),
#                                              rep("CO", length(CoNS_syndecan)),
#                                              rep("ST", length(Staph_unknown_syndecan)),
#                                              rep("EN", length(Enterococcus_syndecan)),
#                                              rep("EC", length(E.coli_syndecan)),
#                                              rep("STR", length(Streptococcus_syndecan))),
#                           syndecan_value = c(Klebsiella_syndecan, S.aureus_syndecan, 
#                                              CoNS_syndecan, Staph_unknown_syndecan, 
#                                              Enterococcus_syndecan, E.coli_syndecan, Streptococcus_syndecan))

#aov(syndecan_value ~ pathegon_type, data = gathered_data) %>% summary()
#library(rstatix)
#gathered_data %>% pairwise_t_test(syndecan_value ~ pathegon_type, p.adjust.method = "bonferroni")


# Perform the Kruskal-Wallis test
#kruskal_test_result <- kruskal.test(syndecan_value ~ pathegon_type, data = gathered_data)
#print(kruskal_test_result)
# If the Kruskal-Wallis test is significant, perform Dunn's post-hoc test
#dunn_test_result <- dunnTest(syndecan_value ~ pathegon_type, data = gathered_data, method = "bonferroni")
#print(dunn_test_result)



#Questions about table1
##1.need to correct for immune_sup 
##6.imputation for MEWS, has warning