#MARS script to impute SOFA

#packages
library(stringi)
library(ggplot2)
library(ggpubr)
library(haven)
library(dplyr)
library(readxl)
library(reshape)
library(stringr)
library(dplyr)
library(tidyr)
library(broom)
library(ipw)
library(survival)
library(splines)
library(zoo)
library(mice)
library(pammtools)
options(max.print=10000000)

# Longitudinal part, dailys
sofas <- read_sas("D:/MARS brondata/MARS r dataframes/daily_sepsis_rel_organ_fail.sas7bdat")

admission <- read_sas("D:/MARS brondata/MARS r dataframes/icu_admission.sas7bdat")

patient <- read_sas("D:/MARS brondata/MARS r dataframes/patient.sas7bdat")

episode_of_sepsis <- read_sas("D:/MARS brondata/MARS r dataframes/episode_of_sepsis.sas7bdat")

assessments <- read_sas("D:/MARS brondata/MARS r dataframes/assessment_dates.sas7bdat")

sofas2 <- merge(sofas, assessments[,c("ICU_Admission_TK", "Assessment_Dates_TK","Assessment_date")], by="Assessment_Dates_TK")

sofas2 <- merge(sofas2, admission[,c("ICU_ID_from_datasource", "ICU_Admission_TK",
                                     "Admission_source","Admission_type","Patient_length","Patient_weight","death_in_ICU",               
                                     "Patient_age_at_ICU_Admission","Planned_admission","Primary_specialty", "Index_ICU_discharge_datetime")], by="ICU_Admission_TK")

hospital <- read_sas("D:/MARS brondata/MARS r dataframes/hospital_admission.sas7bdat")
hospital$Batch_id <- NULL
hospital$Original_TK<- NULL
#hospital$Index_ICU_admittance_datetime <- NULL
followup <- read_sas("D:/MARS brondata/MARS r dataframes/follow_up.sas7bdat")
followup$Batch_id <- NULL
followup$batch_id_follow_up <- NULL
followup$Original_TK<- NULL
#followup$Index_ICU_admittance_datetime <- NULL
patient <- merge(patient, hospital[,c("Hospital_Admission_TK", "Patient_TK", "Index_ICU_admittance_datetime")], by=c("Patient_TK", "Index_ICU_admittance_datetime"), all.x = T)
patient <- merge(patient, admission[,c("Hospital_Admission_TK", "ICU_ID_from_datasource","Index_ICU_admittance_datetime")], by=c("Hospital_Admission_TK", "Index_ICU_admittance_datetime"), all.x = T)

hospital_followup <- unique(merge(hospital, followup, by=c("Hospital_Admission_TK","Index_ICU_admittance_datetime"), all.x = T))
hospital_followup$batch_id <- hospital_followup$QA_Status <- hospital_followup$QA_Status_datetime <- NULL
hospital_followup <- merge(hospital_followup, admission[, c("ICU_Admission_TK", "Hospital_Admission_TK", "Index_ICU_admittance_datetime")], 
                           by=c("Hospital_Admission_TK", "Index_ICU_admittance_datetime"), all.x = T)
hospital_followup$Index_ICU_admittance_datetime<-NULL

#SOFA scores
sofas2$SOFA_Coagulation <- substr(sofas2$SOFA_Coagulation, 0, 1)
sofas2$SOFA_Liver <- substr(sofas2$SOFA_Liver, 0, 1)
sofas2$SOFA_Renal <- substr(sofas2$SOFA_Renal, 0, 1)
sofas2$SOFA_Respiration <- substr(sofas2$SOFA_Respiration, 0, 1)
sofas2$SOFA_Circulation <- substr(sofas2$SOFA_Circulation, 0, 1)
sofas2$SOFA_CNS <- substr(sofas2$SOFA_CNS, 0, 1)

sofas2$SOFA_CNS <- as.numeric(sofas2$SOFA_CNS)
sofas2$SOFA_Coagulation <- as.numeric(sofas2$SOFA_Coagulation)
sofas2$SOFA_Liver <- as.numeric(sofas2$SOFA_Liver)
sofas2$SOFA_Renal <- as.numeric(sofas2$SOFA_Renal)
sofas2$SOFA_Respiration <- as.numeric(sofas2$SOFA_Respiration)
sofas2$SOFA_Circulation <- as.numeric(sofas2$SOFA_Circulation)

sofas2$SOFA_CNS[sofas2$SOFA_CNS=='9']<-NA

sofas2$SOFAtot<-sofas2$SOFA_Coagulation+sofas2$SOFA_Liver+sofas2$SOFA_Renal+sofas2$SOFA_Respiration+
  sofas2$SOFA_Circulation

sofas2$Acute_kidney_injury_score <- substr(sofas2$Acute_kidney_injury_score, 0, 1)
sofas2$Acute_kidney_injury_score[sofas2$Acute_kidney_injury_score=='9']<-NA
sofas2$Acute_kidney_injury_score <- as.numeric(sofas2$Acute_kidney_injury_score)
summary(sofas2$Acute_kidney_injury_score)

#set time
Sys.getlocale("LC_TIME") -> loc # First save your current locale
Sys.setlocale("LC_TIME", "English") # Set to an english locale (so that "Thu" and "Nov" are recognized) 
sofas2$timediff_admin_rifle<-round(difftime(strptime(sofas2$Assessment_date,format="%Y-%m-%d"),
                                            strptime(sofas2$Index_ICU_admittance_datetime,format="%Y-%m-%d"),units="days"))

#include std daily lab
daily_lab <- read_sas("D:/MARS brondata/MARS r dataframes/daily_standard_lab_results.sas7bdat")
daily_lab <- merge(daily_lab, assessments[,c("ICU_Admission_TK", "Assessment_Dates_TK","Assessment_date")], by="Assessment_Dates_TK")
daily_lab <- merge(daily_lab, admission[,c("ICU_ID_from_datasource", "ICU_Admission_TK")], by="ICU_Admission_TK")
daily_lab$Original_TK <- daily_lab$Batch_id <- daily_lab$Assessment_datetime<- NULL
daily_lab$Index_ICU_admittance_datetime<- NULL

sofas3 <- merge(sofas2, daily_lab, by=c("ICU_ID_from_datasource", "ICU_Admission_TK", "Assessment_Dates_TK","Assessment_date"))

baseline_lab <- read_sas("D:/MARS brondata/MARS r dataframes/baseline_laboratory_results.sas7bdat")
baseline_lab$Original_TK <- baseline_lab$Index_ICU_admittance_datetime <- baseline_lab$Batch_id <- NULL

baseline_period <- read_sas("D:/MARS brondata/MARS r dataframes/baseline_period.sas7bdat")
baseline_period$Batch_id <- baseline_period$Original_TK <- baseline_period$Index_ICU_admittance_datetime <- NULL

baseline_lab <- merge(baseline_lab, baseline_period, by="Baseline_Period_TK")

sofas4 <- merge(sofas3, baseline_lab, by="ICU_Admission_TK")

sofas4$eerste_daily <- ave(sofas4$Assessment_datetime, factor(sofas4$ICU_ID_from_datasource), #hier een kleine exercitie om het verschil tussen "baseline blauwe variablen" en de "eerste daily gele" variabelen
                           FUN= function(x) min(x,na.rm = T))

sofas4$eerste_daily2 <- as.character(strptime(sofas4$eerste_daily,format="%Y-%m-%d"))
sofas4$eerste_daily3 <- "06:00:00"
sofas4$eerste_daily4 <- paste(sofas4$eerste_daily2, sofas4$eerste_daily3, " ")
sofas4$eerste_daily2<- sofas4$eerste_daily3 <- NULL
sofas4$eerste_daily4 <- strptime(sofas4$eerste_daily4,format="%Y-%m-%d %H:%M:%OS")
sofas4$diff_baseline_eerste_daily <- round(difftime(strptime(sofas4$eerste_daily4,format="%Y-%m-%d %H:%M:%OS"),
                                                    strptime(sofas4$Baseline_start_datetime,format="%Y-%m-%d %H:%M:%OS"),units="hours"))

ggplot(sofas4 %>% distinct(ICU_ID_from_datasource, .keep_all = TRUE), aes(diff_baseline_eerste_daily))+
  geom_density()+
  geom_vline(aes(xintercept=mean(diff_baseline_eerste_daily)),
             color="blue", linetype="dashed", size=1)

#chronic comorbidity
baseline_chronic_co <- read_sas("D:/MARS brondata/MARS r dataframes/baseline_chronic_comorbidity.sas7bdat")
baseline_chronic_co <- merge(baseline_chronic_co, baseline_period, by="Baseline_Period_TK")
baseline_chronic_co$Batch_id <- NULL
baseline_chronic_co$Original_TK<- NULL
baseline_chronic_co$Index_ICU_admittance_datetime <- NULL
baseline_chronic_co$Baseline_start_datetime <- NULL
baseline_chronic_co$Baseline_Period_TK <- NULL

sofas5 <- merge(sofas4, baseline_chronic_co, by="ICU_Admission_TK")
sofas5$ICU_ID_from_datasource <- as.numeric(sofas5$ICU_ID_from_datasource)

sofas5$diff_admin_daily <- as.numeric(round(difftime(strptime(sofas5$Assessment_datetime,format="%Y-%m-%d"),
                                                     strptime(sofas5$Index_ICU_admittance_datetime,format="%Y-%m-%d"),units="days")))

#include daily treatments
daily_treatment <- read_sas("D:/MARS brondata/MARS r dataframes/daily_treatments_and_events.sas7bdat")
daily_treatment <- merge(daily_treatment, assessments[,c("ICU_Admission_TK", "Assessment_Dates_TK","Assessment_date")], by="Assessment_Dates_TK")
daily_treatment <- merge(daily_treatment, admission[,c("ICU_ID_from_datasource", "ICU_Admission_TK")], by="ICU_Admission_TK")
daily_treatment$Original_TK <- daily_treatment$Batch_id <- daily_treatment$Assessment_datetime<- NULL
daily_treatment$Index_ICU_admittance_datetime<- NULL
daily_treatment$Treatments_And_Events_TK <- NULL

sofas6 <- merge(sofas5, daily_treatment, by=c("ICU_ID_from_datasource", "ICU_Admission_TK", "Assessment_Dates_TK","Assessment_date"),all.x = T)

#include daily monitor resp
daily_monitor_r <- read_sas("D:/MARS brondata/MARS r dataframes/daily_monitor_data_respiratory.sas7bdat")
daily_monitor_r <- merge(daily_monitor_r, assessments[,c("ICU_Admission_TK", "Assessment_Dates_TK","Assessment_date")], by="Assessment_Dates_TK")
daily_monitor_r <- merge(daily_monitor_r, admission[,c("ICU_ID_from_datasource", "ICU_Admission_TK")], by="ICU_Admission_TK")
daily_monitor_r$Original_TK <- daily_monitor_r$Batch_id <- daily_monitor_r$Assessment_datetime<- NULL
daily_monitor_r$Index_ICU_admittance_datetime<- NULL
daily_monitor_r$Treatments_And_Events_TK <- NULL

sofas7 <- merge(sofas6, daily_monitor_r, by=c("ICU_ID_from_datasource", "ICU_Admission_TK", "Assessment_Dates_TK","Assessment_date"),all.x = T)

#include daily monitor hemo
daily_monitor_h <- read_sas("D:/MARS brondata/MARS r dataframes/daily_monitor_data_hemodynamic.sas7bdat")
daily_monitor_h <- merge(daily_monitor_h, assessments[,c("ICU_Admission_TK", "Assessment_Dates_TK","Assessment_date")], by="Assessment_Dates_TK")
daily_monitor_h <- merge(daily_monitor_h, admission[,c("ICU_ID_from_datasource", "ICU_Admission_TK")], by="ICU_Admission_TK")
daily_monitor_h$Original_TK <- daily_monitor_h$Batch_id <- daily_monitor_h$Assessment_datetime<- NULL
daily_monitor_h$Index_ICU_admittance_datetime<- NULL
daily_monitor_h$Treatments_And_Events_TK <- NULL

sofas8 <- merge(sofas7, daily_monitor_h, by=c("ICU_ID_from_datasource", "ICU_Admission_TK", "Assessment_Dates_TK","Assessment_date"),all.x = T)

#include daily sepsis status & #include daily additional lab
daily_sepsis <- read_sas("D:/MARS brondata/MARS r dataframes/daily_sepsis_status.sas7bdat")
daily_sepsis <- merge(daily_sepsis, assessments[,c("ICU_Admission_TK", "Assessment_Dates_TK","Assessment_date")], by="Assessment_Dates_TK")
daily_sepsis <- merge(daily_sepsis, admission[,c("ICU_ID_from_datasource", "ICU_Admission_TK")], by="ICU_Admission_TK")
daily_sepsis$Original_TK <- daily_sepsis$Batch_id <- daily_sepsis$Assessment_datetime<- NULL
daily_sepsis$Index_ICU_admittance_datetime<- NULL
daily_sepsis$Treatments_And_Events_TK <- NULL

daily_add <- read_sas("D:/MARS brondata/MARS r dataframes/daily_additional_lab_results.sas7bdat")
daily_add <- merge(daily_add, assessments[,c("ICU_Admission_TK", "Assessment_Dates_TK","Assessment_date")], by="Assessment_Dates_TK")
daily_add <- merge(daily_add, admission[,c("ICU_ID_from_datasource", "ICU_Admission_TK")], by="ICU_Admission_TK")
daily_add$Original_TK <- daily_add$Batch_id <- daily_add$Assessment_datetime<- NULL
daily_add$Index_ICU_admittance_datetime<- NULL
daily_add$Additional_Lab_Results_TK <- NULL

sofas8 <- merge(sofas8, daily_sepsis, by=c("ICU_ID_from_datasource", "ICU_Admission_TK", "Assessment_Dates_TK","Assessment_date"),all.x = T)
sofas8 <- merge(sofas8, daily_add, by=c("ICU_ID_from_datasource", "ICU_Admission_TK", "Assessment_Dates_TK","Assessment_date"),all.x = T)

#vasopressor, zelf gemaakt in MARS
setwd("D:/MARS brondata/MARS/")
vasopressor<-read.csv("MARS_Daily_vasopressor_FU.csv",sep=";", dec = ".", header = T, stringsAsFactors = FALSE ) #eigen gemaakt variable in MARS door Fabrice
vasopressor$Assessment_date<-strptime(vasopressor$Assessment_date,"%d/%m/%Y")
head(vasopressor$Assessment_date)
head(sofas8$Assessment_date)
vasopressor$Assess_day_ICU<- NULL
sofas8$Assessment_date<-strptime(sofas8$Assessment_date,format="%Y-%m-%d")
head(sofas8$Assessment_date)
sofas9 <- merge(sofas8, vasopressor, by=c("ICU_ID_from_datasource", "Assessment_date"),all.x = T)
sofas9 <- unique(sofas9)

#include baseline acute diagnoses
baseline_acute_diag <- read_sas("D:/MARS brondata/MARS r dataframes/baseline_acute_diagnoses.sas7bdat")
baseline_acute_diag <- merge(baseline_acute_diag, baseline_period, by="Baseline_Period_TK")
baseline_acute_diag$Batch_id <- NULL
baseline_acute_diag$Original_TK<- NULL
baseline_acute_diag$Index_ICU_admittance_datetime <- NULL
baseline_acute_diag$Baseline_start_datetime <- NULL
baseline_acute_diag$Baseline_Period_TK <- NULL

sofas10 <- merge(sofas9, baseline_acute_diag, by="ICU_Admission_TK",all.x = T)

#include baseline acute diagnoses
baseline_chronic_meds <- read_sas("D:/MARS brondata/MARS r dataframes/baseline_chronic_medication_use.sas7bdat")
baseline_chronic_meds <- merge(baseline_chronic_meds, baseline_period, by="Baseline_Period_TK")
baseline_chronic_meds$Batch_id <- NULL
baseline_chronic_meds$Original_TK<- NULL
baseline_chronic_meds$Index_ICU_admittance_datetime <- NULL
baseline_chronic_meds$Baseline_Period_TK <- NULL
baseline_chronic_meds$Medications_free_text1<-NULL
baseline_chronic_meds$Medications_free_text2<-NULL

sofas10 <- merge(sofas10, baseline_chronic_meds, by="ICU_Admission_TK",all.x = T)

#include baseline observations treatments
baseline_obs_treat <- read_sas("D:/MARS brondata/MARS r dataframes/baseline_observations_treatments.sas7bdat")
baseline_obs_treat <- merge(baseline_obs_treat, baseline_period, by="Baseline_Period_TK")
baseline_obs_treat$Batch_id <- NULL
baseline_obs_treat$Original_TK<- NULL
baseline_obs_treat$Index_ICU_admittance_datetime <- NULL
baseline_obs_treat$Baseline_start_datetime <- NULL
baseline_obs_treat$Baseline_Period_TK <- NULL

sofas11 <- merge(sofas10, baseline_obs_treat, by="ICU_Admission_TK",all.x = T)

#include baseline monitoring
baseline_monitor <- read_sas("D:/MARS brondata/MARS r dataframes/baseline_monitor_data.sas7bdat")
baseline_monitor <- merge(baseline_monitor, baseline_period, by="Baseline_Period_TK")
baseline_monitor$Batch_id <- NULL
baseline_monitor$Original_TK<- NULL
baseline_monitor$Index_ICU_admittance_datetime <- NULL
baseline_monitor$Baseline_start_datetime <- NULL
baseline_monitor$Baseline_Period_TK <- NULL

sofas12 <- merge(sofas11, baseline_monitor, by="ICU_Admission_TK",all.x = T)

#create SOFA score >2 (sepsis3), and impute missings
sofas13 <- sofas12[order( sofas12$diff_admin_daily, sofas12$ICU_ID_from_datasource),] 
#impute SOFA score
sofa13_mis_sofa<-subset(sofas13, is.na(sofas13$SOFAtot), select = c("ICU_ID_from_datasource", 
                                                                    "ICU_Admission_TK",
                                                                    "Assessment_datetime", 
                                                                    "Acute_kidney_injury_score",
                                                                    "Creatinin_Max" ,
                                                                    "Urea_Max"     ,
                                                                    "Use_of_mech_ventilation",
                                                                    "Use_of_renal_rep_therapy",
                                                                    "Platelet_count_Min",
                                                                    "SOFA_Circulation",
                                                                    "SOFA_CNS",
                                                                    "SOFA_Coagulation",
                                                                    "SOFA_Liver",
                                                                    "SOFA_Renal",
                                                                    "SOFA_Respiration",
                                                                    "Vasopressor",
                                                                    "Platelet_count_Min_24h",
                                                                    "SOFAtot",
                                                                    "ABPm_Min",
                                                                    "FiO2",
                                                                    "SpO2",
                                                                    "Bilirubin_Max",
                                                                    "Creatinin_Max_24h",
                                                                    "Urea_Max_24h",
                                                                    "P_F_ratio_Min_24h", "Bilirubin_Max_24h", 
                                                                    "Mechanical_ventilation_24h", "ABPm_Min_24h",
                                                                    "FiO2_at_highest_A_aDO2", "PaO2",
                                                                    "PaO2_at_highest_A_aDO2",
                                                                    "Total_daily_urine_output",
                                                                    "Urine_output_Sum_24h",
                                                                    "diff_admin_daily"))

ggplot(subset(sofa13_mis_sofa,sofa13_mis_sofa$diff_admin_daily<31), aes(diff_admin_daily))+
  geom_bar()+
  scale_y_continuous(name = "Missing SOFA score", limits=c(0, 700) )

colnames(sofa13_mis_sofa)
num <- c(1:3, 5:34) 
sofa13_mis_sofa[num] <- lapply(sofa13_mis_sofa[num], as.numeric)
#sofa coagulation
sofa13_mis_sofa$SOFA_Coagulation<-ifelse((sofa13_mis_sofa$Platelet_count_Min>=150), 
                                         "0", sofa13_mis_sofa$SOFA_Coagulation)
sofa13_mis_sofa$SOFA_Coagulation<-ifelse((sofa13_mis_sofa$Platelet_count_Min<150), 
                                         "1", sofa13_mis_sofa$SOFA_Coagulation)
sofa13_mis_sofa$SOFA_Coagulation<-ifelse((sofa13_mis_sofa$Platelet_count_Min<100), 
                                         "2", sofa13_mis_sofa$SOFA_Coagulation)
sofa13_mis_sofa$SOFA_Coagulation<-ifelse((sofa13_mis_sofa$Platelet_count_Min<50), 
                                         "3", sofa13_mis_sofa$SOFA_Coagulation)
sofa13_mis_sofa$SOFA_Coagulation<-ifelse((sofa13_mis_sofa$Platelet_count_Min<20), 
                                         "4", sofa13_mis_sofa$SOFA_Coagulation)
table(sofa13_mis_sofa$SOFA_Coagulation, exclude = NULL) 
#sofa13_mis_sofa$SOFA_Coagulation<-ifelse(is.na(sofa13_mis_sofa$SOFA_Coagulation), 0, sofa13_mis_sofa$SOFA_Coagulation)
table(sofa13_mis_sofa$SOFA_Coagulation, exclude = NULL) 
sofa13_mis_sofa$SOFA_Coagulation<-as.numeric(sofa13_mis_sofa$SOFA_Coagulation)
#sofa respiration
table(sofa13_mis_sofa$P_F_ratio_Min_24h, exclude = NULL)
sofa13_mis_sofa$PaO2_FIO2_ratio<-round(sofa13_mis_sofa$PaO2/(sofa13_mis_sofa$FiO2/100))
table(sofa13_mis_sofa$PaO2_FIO2_ratio, exclude = NULL)
sofa13_mis_sofa$PaO2_FIO2_ratio_highest<-round(sofa13_mis_sofa$PaO2_at_highest_A_aDO2/(sofa13_mis_sofa$FiO2_at_highest_A_aDO2/100))
table(sofa13_mis_sofa$PaO2_FIO2_ratio_highest, exclude = NULL)
sofa13_mis_sofa$min_P_F_ratio <- pmin(sofa13_mis_sofa$P_F_ratio_Min_24h, sofa13_mis_sofa$PaO2_FIO2_ratio, sofa13_mis_sofa$PaO2_FIO2_ratio_highest, na.rm=T  )
sofa13_mis_sofa$SOFA_Respiration<-ifelse((sofa13_mis_sofa$min_P_F_ratio>=400), 
                                         "0", sofa13_mis_sofa$SOFA_Respiration)
sofa13_mis_sofa$SOFA_Respiration<-ifelse((sofa13_mis_sofa$min_P_F_ratio<400), 
                                         "1", sofa13_mis_sofa$SOFA_Respiration)
sofa13_mis_sofa$SOFA_Respiration<-ifelse((sofa13_mis_sofa$min_P_F_ratio<300), 
                                         "2", sofa13_mis_sofa$SOFA_Respiration)
sofa13_mis_sofa$SOFA_Respiration<-ifelse((sofa13_mis_sofa$min_P_F_ratio<200 & sofa13_mis_sofa$Use_of_mech_ventilation=="1"), 
                                         "3", sofa13_mis_sofa$SOFA_Respiration)
sofa13_mis_sofa$SOFA_Respiration<-ifelse((sofa13_mis_sofa$min_P_F_ratio<100 & sofa13_mis_sofa$Use_of_mech_ventilation=="1"), 
                                         "4", sofa13_mis_sofa$SOFA_Respiration)
table(sofa13_mis_sofa$SOFA_Respiration, exclude = NULL)
#sofa13_mis_sofa$SOFA_Respiration<-ifelse(is.na(sofa13_mis_sofa$SOFA_Respiration), 0, sofa13_mis_sofa$SOFA_Respiration)
table(sofa13_mis_sofa$SOFA_Respiration, exclude = NULL)
#sofa liver
sofa13_mis_sofa$max_bili <- pmax(sofa13_mis_sofa$Bilirubin_Max, sofa13_mis_sofa$Bilirubin_Max_24h, na.rm=T  )
table(sofa13_mis_sofa$max_bili, exclude = NULL)
sofa13_mis_sofa$SOFA_Liver<-ifelse((sofa13_mis_sofa$max_bili>204), 
                                   "4", sofa13_mis_sofa$SOFA_Liver)
sofa13_mis_sofa$SOFA_Liver<-ifelse((sofa13_mis_sofa$max_bili<=204), 
                                   "3", sofa13_mis_sofa$SOFA_Liver)
sofa13_mis_sofa$SOFA_Liver<-ifelse((sofa13_mis_sofa$max_bili<=101), 
                                   "2", sofa13_mis_sofa$SOFA_Liver)
sofa13_mis_sofa$SOFA_Liver<-ifelse((sofa13_mis_sofa$max_bili<=32), 
                                   "1", sofa13_mis_sofa$SOFA_Liver)
sofa13_mis_sofa$SOFA_Liver<-ifelse((sofa13_mis_sofa$max_bili<=20), 
                                   "0", sofa13_mis_sofa$SOFA_Liver)
table(sofa13_mis_sofa$SOFA_Liver, exclude = NULL) 
#sofa13_mis_sofa$SOFA_Liver<-ifelse(is.na(sofa13_mis_sofa$SOFA_Liver), 0, sofa13_mis_sofa$SOFA_Liver)
table(sofa13_mis_sofa$SOFA_Liver, exclude = NULL) 
#sofa renal
sofa13_mis_sofa$max_creat <- pmax(sofa13_mis_sofa$Creatinin_Max, sofa13_mis_sofa$Creatinin_Max_24h, na.rm=T  )
sofa13_mis_sofa$min_urine <- pmin(sofa13_mis_sofa$Total_daily_urine_output, sofa13_mis_sofa$Urine_output_Sum_24h, na.rm=T  )
table(sofa13_mis_sofa$max_creat, exclude = NULL)
table(sofa13_mis_sofa$min_urine, exclude = NULL)

sofa13_mis_sofa$SOFA_Renal<-ifelse((sofa13_mis_sofa$max_creat>440), 
                                   4, NA)
sofa13_mis_sofa$SOFA_Renal<-ifelse((sofa13_mis_sofa$max_creat<=440), 
                                   3, sofa13_mis_sofa$SOFA_Renal)
sofa13_mis_sofa$SOFA_Renal<-ifelse((sofa13_mis_sofa$max_creat<=299), 
                                   2, sofa13_mis_sofa$SOFA_Renal)
sofa13_mis_sofa$SOFA_Renal<-ifelse((sofa13_mis_sofa$max_creat<=170), 
                                   1, sofa13_mis_sofa$SOFA_Renal)
sofa13_mis_sofa$SOFA_Renal<-ifelse((sofa13_mis_sofa$max_creat<110), 
                                   0, sofa13_mis_sofa$SOFA_Renal)

sofa13_mis_sofa$min_urine_cat<-ifelse(sofa13_mis_sofa$min_urine<500, 3, NA)
sofa13_mis_sofa$min_urine_cat<-ifelse(sofa13_mis_sofa$min_urine<200, 4, sofa13_mis_sofa$min_urine_cat)
table(sofa13_mis_sofa$min_urine_cat)

sofa13_mis_sofa$SOFA_Renal[which(sofa13_mis_sofa$min_urine_cat==3)] <- 3 
sofa13_mis_sofa$SOFA_Renal[which(sofa13_mis_sofa$min_urine_cat==4)] <- 4 

table(sofa13_mis_sofa$SOFA_Renal, exclude = NULL)
#sofa13_mis_sofa$SOFA_Renal<-ifelse(is.na(sofa13_mis_sofa$SOFA_Renal), 0, sofa13_mis_sofa$SOFA_Renal)
table(sofa13_mis_sofa$SOFA_Renal, exclude = NULL)


#sofa cardiovascular
meds_icu <- read_sas("D:/MARS brondata/MARS r dataframes/event_medications.sas7bdat")
meds_icu <- merge(meds_icu, admission[,c("ICU_Admission_TK", "ICU_ID_from_datasource")], by="ICU_Admission_TK", all.x=T) #need the index
meds_icu$ICU_ID_from_datasource <- as.numeric(meds_icu$ICU_ID_from_datasource)
meds_icu_MARS_ALL <- read_excel("D:/AKI grant/meds_icu_MARS_ALL.xlsx", sheet = "meds_icu_new")

meds_icu$diff_admin_meds <- round(difftime(as.POSIXct(meds_icu$Medication_startdatetime,format="%Y-%m-%d"),
                                           as.POSIXct(meds_icu$Index_ICU_admittance_datetime,format="%Y-%m-%d"),units="days"))

meds_icu$time_gift_min <- round(difftime(as.POSIXct(meds_icu$medication_stopdatetime,format="%Y-%m-%d %M:%H:%OS"),
                                         as.POSIXct(meds_icu$Medication_startdatetime,format="%Y-%m-%d %M:%H:%OS"),units="min"))
table(meds_icu$diff_admin_meds)
head(meds_icu$time_gift_min)

meds_icu <- merge(meds_icu, admission[,c("ICU_ID_from_datasource","Patient_weight" )], by="ICU_ID_from_datasource", all.y = T) #37 have missing on weight, but they all are less then 3 hours or so in the ICU, dont bother

drugs_vasopressor <- subset(meds_icu_MARS_ALL, meds_icu_MARS_ALL$Vaso_alpha=="y" | 
                              meds_icu_MARS_ALL$Vaso_alpha_beta=="y")$Var1
drugs_vasopressor <- drugs_vasopressor[!(drugs_vasopressor %in% c("Methyldopa", "Fenylefrine OK/HC dosis","Epinefrine e.d","Fenylefrine e.d.", 
                                                                  "fenylefrine","Fenylefrine d.d." ))] #eig vreemd dat fenylefrine niet bij de sofa cardio hoort
drugs_vasopressor

meds_icu$vasopressor <- ifelse(meds_icu$medication_generic_drug_name %in% drugs_vasopressor, 1, 0)
meds_icu_longitudinal_vasopressor <- meds_icu[,c("ICU_ID_from_datasource","diff_admin_meds","medication_generic_drug_name", "vasopressor","time_gift_min", "medication_dose", "medication_unit", "Patient_weight")]
meds_icu_longitudinal_vasopressor <- meds_icu_longitudinal_vasopressor[meds_icu_longitudinal_vasopressor$vasopressor==1 & meds_icu_longitudinal_vasopressor$diff_admin_meds>=0,] #only this treatment
meds_icu_longitudinal_vasopressor <- meds_icu_longitudinal_vasopressor[!(meds_icu_longitudinal_vasopressor$medication_unit=="drop"),] #exclude drops

meds_icu_longitudinal_vasopressor$ug_kg_min <- (meds_icu_longitudinal_vasopressor$medication_dose*1000) / meds_icu_longitudinal_vasopressor$Patient_weight / as.numeric(meds_icu_longitudinal_vasopressor$time_gift_min)
meds_icu_longitudinal_vasopressor$sofa_cardio_2 <- ifelse( grepl("dobu",meds_icu_longitudinal_vasopressor$medication_generic_drug_name,ignore.case = T) | 
                                                             (grepl("dopa",meds_icu_longitudinal_vasopressor$medication_generic_drug_name,ignore.case = T) & 
                                                                meds_icu_longitudinal_vasopressor$ug_kg_min<=5) , 1, 0)
meds_icu_longitudinal_vasopressor$sofa_cardio_3 <- ifelse( (grepl("dopa",meds_icu_longitudinal_vasopressor$medication_generic_drug_name,ignore.case = T) & 
                                                              meds_icu_longitudinal_vasopressor$ug_kg_min>5) |
                                                             (grepl("epine",meds_icu_longitudinal_vasopressor$medication_generic_drug_name,ignore.case = T) & 
                                                                meds_icu_longitudinal_vasopressor$ug_kg_min<=0.1), 1, 0)
meds_icu_longitudinal_vasopressor$sofa_cardio_4 <- ifelse( (grepl("dopa",meds_icu_longitudinal_vasopressor$medication_generic_drug_name,ignore.case = T) & 
                                                              meds_icu_longitudinal_vasopressor$ug_kg_min>15) |
                                                             (grepl("epine",meds_icu_longitudinal_vasopressor$medication_generic_drug_name,ignore.case = T) & 
                                                                meds_icu_longitudinal_vasopressor$ug_kg_min>0.1), 1, 0)
meds_icu_longitudinal_vasopressor$vasopressor_use <- ifelse( (grepl("dopa",meds_icu_longitudinal_vasopressor$medication_generic_drug_name,ignore.case = T) & 
                                                                meds_icu_longitudinal_vasopressor$ug_kg_min>10) |
                                                               (grepl("epine",meds_icu_longitudinal_vasopressor$medication_generic_drug_name,ignore.case = T) & 
                                                                  meds_icu_longitudinal_vasopressor$ug_kg_min>0.1) |
                                                               (grepl("dobu",meds_icu_longitudinal_vasopressor$medication_generic_drug_name,ignore.case = T) & 
                                                                  meds_icu_longitudinal_vasopressor$ug_kg_min>1), 1, 0) #NED : DOPA : DOBU = 1:100:10 according to Fabrice for vaso
#kleine check
meds_icu_longitudinal_vasopressor_2 <- unique(meds_icu_longitudinal_vasopressor[meds_icu_longitudinal_vasopressor$sofa_cardio_2==1,c("ICU_ID_from_datasource", "diff_admin_meds", "sofa_cardio_2")] )
table(duplicated(meds_icu_longitudinal_vasopressor_2[,c("ICU_ID_from_datasource", "diff_admin_meds", "sofa_cardio_2")])) #klopt

meds_icu_longitudinal_vasopressor_3 <- unique(meds_icu_longitudinal_vasopressor[meds_icu_longitudinal_vasopressor$sofa_cardio_3==1,c("ICU_ID_from_datasource", "diff_admin_meds", "sofa_cardio_3")] )
table(duplicated(meds_icu_longitudinal_vasopressor_3[,c("ICU_ID_from_datasource", "diff_admin_meds", "sofa_cardio_3")])) #klopt

meds_icu_longitudinal_vasopressor_4 <- unique(meds_icu_longitudinal_vasopressor[meds_icu_longitudinal_vasopressor$sofa_cardio_4==1,c("ICU_ID_from_datasource", "diff_admin_meds", "sofa_cardio_4")] )
table(duplicated(meds_icu_longitudinal_vasopressor_4[,c("ICU_ID_from_datasource", "diff_admin_meds", "sofa_cardio_4")])) #klopt

meds_icu_longitudinal_vasopressor_use <- unique(meds_icu_longitudinal_vasopressor[meds_icu_longitudinal_vasopressor$vasopressor_use==1,c("ICU_ID_from_datasource", "diff_admin_meds", "vasopressor_use")] )
table(duplicated(meds_icu_longitudinal_vasopressor_use[,c("ICU_ID_from_datasource", "diff_admin_meds", "vasopressor_use")])) #klopt

sofa13_mis_sofa <- merge(sofa13_mis_sofa, meds_icu_longitudinal_vasopressor_2, by.x = c("ICU_ID_from_datasource", "diff_admin_daily"), 
                 by.y = c("ICU_ID_from_datasource","diff_admin_meds"), all.x = T)
sofa13_mis_sofa <- merge(sofa13_mis_sofa, meds_icu_longitudinal_vasopressor_3, by.x = c("ICU_ID_from_datasource", "diff_admin_daily"), 
                 by.y = c("ICU_ID_from_datasource","diff_admin_meds"), all.x = T)
sofa13_mis_sofa <- merge(sofa13_mis_sofa, meds_icu_longitudinal_vasopressor_4, by.x = c("ICU_ID_from_datasource", "diff_admin_daily"), 
                 by.y = c("ICU_ID_from_datasource","diff_admin_meds"), all.x = T)
sofa13_mis_sofa <- merge(sofa13_mis_sofa, meds_icu_longitudinal_vasopressor_use, by.x = c("ICU_ID_from_datasource", "diff_admin_daily"), 
                 by.y = c("ICU_ID_from_datasource","diff_admin_meds"), all.x = T)

sofa13_mis_sofa <- sofa13_mis_sofa[order( sofa13_mis_sofa$diff_admin_daily, sofa13_mis_sofa$ICU_ID_from_datasource),]
sofa13_mis_sofa$sofa_cardio_2 <- ifelse(is.na(sofa13_mis_sofa$sofa_cardio_2), 0, sofa13_mis_sofa$sofa_cardio_2)
sofa13_mis_sofa$sofa_cardio_3 <- ifelse(is.na(sofa13_mis_sofa$sofa_cardio_3), 0, sofa13_mis_sofa$sofa_cardio_3)
sofa13_mis_sofa$sofa_cardio_4 <- ifelse(is.na(sofa13_mis_sofa$sofa_cardio_4), 0, sofa13_mis_sofa$sofa_cardio_4)
sofa13_mis_sofa$vasopressor_use <- ifelse(is.na(sofa13_mis_sofa$vasopressor_use), 0, sofa13_mis_sofa$vasopressor_use)

#sofa cardiovascular
table(sofa13_mis_sofa$SOFA_Circulation, exclude = NULL)
sofa13_mis_sofa$min_map<- pmin(sofa13_mis_sofa$ABPm_Min, sofa13_mis_sofa$ABPm_Min_24h, na.rm=T  )
sofa13_mis_sofa$SOFA_Circulation<-ifelse(sofa13_mis_sofa$min_map>=70, 0, NA)
sofa13_mis_sofa$SOFA_Circulation<-ifelse(sofa13_mis_sofa$min_map<70, 1, sofa13_mis_sofa$SOFA_Circulation)
sofa13_mis_sofa$SOFA_Circulation<-ifelse(sofa13_mis_sofa$sofa_cardio_2==1, 2, sofa13_mis_sofa$SOFA_Circulation)
sofa13_mis_sofa$SOFA_Circulation<-ifelse(sofa13_mis_sofa$sofa_cardio_3==1, 3, sofa13_mis_sofa$SOFA_Circulation)
sofa13_mis_sofa$SOFA_Circulation<-ifelse(sofa13_mis_sofa$sofa_cardio_4==1, 4, sofa13_mis_sofa$SOFA_Circulation)
table(sofa13_mis_sofa$SOFA_Circulation, exclude = NULL)

sofa13_mis_sofa$SOFAtot<-(sofa13_mis_sofa$SOFA_Circulation + as.numeric(sofa13_mis_sofa$SOFA_Coagulation) +
                            as.numeric(sofa13_mis_sofa$SOFA_Liver) + sofa13_mis_sofa$SOFA_Renal +
                            as.numeric(sofa13_mis_sofa$SOFA_Respiration))
sofa13_mis_sofa$SOFAtot_imputed<-sofa13_mis_sofa$SOFAtot
sofa13_mis_sofa$SOFA_Circulation_imputed<-sofa13_mis_sofa$SOFA_Circulation
sofa13_mis_sofa$SOFA_Coagulation_imputed<-as.numeric(sofa13_mis_sofa$SOFA_Coagulation)
sofa13_mis_sofa$SOFA_Liver_imputed<-as.numeric(sofa13_mis_sofa$SOFA_Liver)
sofa13_mis_sofa$SOFA_Renal_imputed<-sofa13_mis_sofa$SOFA_Renal
sofa13_mis_sofa$SOFA_Respiration_imputed<-as.numeric(sofa13_mis_sofa$SOFA_Respiration)

ggplot(subset(sofa13_mis_sofa,is.na(sofa13_mis_sofa$SOFAtot_imputed) & sofa13_mis_sofa$diff_admin_daily<31), aes(diff_admin_daily))+
  geom_bar()+
  scale_y_continuous(name = "Missing SOFA score after imputation based on known lab values", limits=c(0, 700) )

sofas14 <- merge(sofas13, sofa13_mis_sofa[,c("ICU_ID_from_datasource",
                                                     "ICU_Admission_TK",
                                                     "diff_admin_daily",
                                                     "SOFAtot_imputed",           
                                                     "SOFA_Circulation_imputed",
                                                     "SOFA_Coagulation_imputed",
                                                     "SOFA_Liver_imputed",
                                                     "SOFA_Renal_imputed",        
                                                     "SOFA_Respiration_imputed")], by=c("ICU_ID_from_datasource",
                                                                                        "ICU_Admission_TK",
                                                                                        "diff_admin_daily"), all.x = T)

sofas14 <- unique (sofas14 )

summary(sofas14$SOFAtot)
sofas14$SOFAtot<-ifelse(is.na(sofas14$SOFAtot), sofas14$SOFAtot_imputed, sofas14$SOFAtot)
sofas14$SOFA_Circulation<-ifelse(is.na(sofas14$SOFA_Circulation), sofas14$SOFA_Circulation_imputed, sofas14$SOFA_Circulation)
sofas14$SOFA_Coagulation<-ifelse(is.na(sofas14$SOFA_Coagulation), sofas14$SOFA_Coagulation_imputed, sofas14$SOFA_Coagulation)
sofas14$SOFA_Liver<-ifelse(is.na(sofas14$SOFA_Liver), sofas14$SOFA_Liver_imputed, sofas14$SOFA_Liver)
sofas14$SOFA_Renal<-ifelse(is.na(sofas14$SOFA_Renal), sofas14$SOFA_Renal_imputed, sofas14$SOFA_Renal)
sofas14$SOFA_Respiration<-ifelse(is.na(sofas14$SOFA_Respiration), sofas14$SOFA_Respiration_imputed, sofas14$SOFA_Respiration)
summary(sofas14$SOFAtot)

sofas14$SOFAtot_narm<-mapply(sum,as.numeric(sofas14$SOFA_Circulation), 
                             as.numeric(sofas14$SOFA_Coagulation),
                             as.numeric(sofas14$SOFA_Liver), 
                             as.numeric(sofas14$SOFA_Renal),
                             as.numeric(sofas14$SOFA_Respiration), na.rm=T) 
sofas14$SOFAtot_narm <- ifelse(is.na(as.numeric(sofas14$SOFA_Circulation))& 
                                 is.na(as.numeric(sofas14$SOFA_Coagulation))&
                                 is.na(as.numeric(sofas14$SOFA_Liver))& 
                                 is.na(as.numeric(sofas14$SOFA_Renal))&
                                 is.na(as.numeric(sofas14$SOFA_Respiration)), NA, sofas14$SOFAtot_narm)        
summary(sofas14$SOFAtot_narm) 

#order df
sofas15 <- sofas14[order(  sofas14$ICU_ID_from_datasource, sofas14$diff_admin_daily),]

#save a workable file for Nerissa
setwd("D:/Nerissa")

sofas16 <- sofas15[,c("ICU_ID_from_datasource",
                      "Index_ICU_admittance_datetime",
                      "ICU_Admission_TK",
                      "Assessment_datetime",
                      "diff_admin_daily",
                      "SOFAtot",
                      "SOFA_Circulation",
                      "SOFA_Coagulation",
                      "SOFA_Liver",
                      "SOFA_Renal",
                      "SOFA_Respiration")]

save(sofas16, file="MARS_SOFAS.RData")




