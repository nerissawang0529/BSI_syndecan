rm(list = ls())

## install.packages("pacman")
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)
install.packages("haven")
library(haven)

## ready data
episode_of_sepsis <- read_sas("Documents/BSI/R_code/original_data/episode_of_sepsis.sas7bdat")
admission <- read_sas("Documents/BSI/R_code/original_data/icu_admission.sas7bdat")
patient <- read_sas("Documents/BSI/R_code/original_data/patient.sas7bdat")
hospital <- read_sas("Documents/BSI/R_code/original_data/hospital_admission (2).sas7bdat")
hospital$Batch_id <- NULL
hospital$Original_TK<- NULL

patient <- merge(patient, hospital[,c("Hospital_Admission_TK", "Patient_TK", "Index_ICU_admittance_datetime")], by=c("Patient_TK", "Index_ICU_admittance_datetime"), all.x = T)
patient <- merge(patient, admission[,c("Hospital_Admission_TK", "ICU_ID_from_datasource","Index_ICU_admittance_datetime")], by=c("Hospital_Admission_TK", "Index_ICU_admittance_datetime"), all.x = T)


################ episode of infection

#make infection at admission cohort
episode_of_sepsis  <- merge(episode_of_sepsis , admission[,c("ICU_ID_from_datasource", "ICU_Admission_TK")], by="ICU_Admission_TK", all.x = T)
episode_of_sepsis <- episode_of_sepsis[!(episode_of_sepsis$Infection_likelihood_primary=="none" &
                                           episode_of_sepsis$Infection_likelihood_secondary=="none"),] #waarbij rijen dubbel 'none' zijn, kunnen er sowieso uit.
#reshape to make every row an infection
episode_of_sepsis2<-reshape(episode_of_sepsis, 
                            direction='long', 
                            varying=c("Infection_likelihood_primary","Site_of_infection_primary", "Infection_likelihood_secondary", "Site_of_infection_secondary"), 
                            timevar='Group',
                            times=c('primary', 'secondary'),
                            v.names=c('likelihood', 'site'),
                            idvar='Episode_Of_Sepsis_TK')
episode_of_sepsis2 <- episode_of_sepsis2[!(episode_of_sepsis2$likelihood==""),]

row.names(episode_of_sepsis2) <- str_sub(row.names(episode_of_sepsis2),0, -9)
episode_of_sepsis2$Index_ICU_admittance_datetime2<-strptime(episode_of_sepsis2$Index_ICU_admittance_datetime,"%Y-%m-%d")

episode_of_sepsis2$startday <- round(difftime(as.POSIXct(episode_of_sepsis2$Sepsis_startdate,format="%Y-%m-%d"),
                                              as.POSIXct(episode_of_sepsis2$Index_ICU_admittance_datetime,format="%Y-%m-%d"),units="days"))
table(episode_of_sepsis2$startday)
episode_of_sepsis3 <- subset(episode_of_sepsis2, episode_of_sepsis2$startday<=2)
episode_of_sepsis4 <- subset(episode_of_sepsis3, !(episode_of_sepsis3$likelihood=="none"))
summary(duplicated(episode_of_sepsis4$ICU_ID_from_datasource)) #3017 admissions with sepsis in MARS (36.2%)
episode_of_sepsis4  <- merge(episode_of_sepsis4 , patient[,c("ICU_ID_from_datasource", "Patient_TK")], by="ICU_ID_from_datasource", all.x = T)
episode_of_sepsis4$Patient_TK <- as.numeric(episode_of_sepsis4$Patient_TK) 
summary(duplicated(episode_of_sepsis4$Patient_TK))


#############################################causative organism

causative_organism_primary <- read_sas("Documents/BSI/R_code/original_data/causative_organism_primary.sas7bdat", NULL)
causative_organism_secondary <- read_sas("Documents/BSI/R_code/original_data/causative_organism_secondary.sas7bdat", NULL)


############################################# microbiology result
event_microbiological_specimen <- read_sas("Documents/BSI/R_code/original_data/event_microbiological_specimen.sas7bdat", NULL)

#only blood
event_microbiological_specimen_blood <- event_microbiological_specimen[grepl("blood",event_microbiological_specimen$Material,ignore.case = T)==T,]

event_microbiological_specimen_blood <- event_microbiological_specimen_blood[!grepl("serum",event_microbiological_specimen_blood$Material,ignore.case = T)==T,]


#merge results
event_microbiological_result <- read_sas("Documents/BSI/R_code/original_data/event_microbiological_result.sas7bdat", NULL)
event_microbiological_result_blood <- event_microbiological_result[event_microbiological_result$Event_Microbiological_Specim_TK %in% 
                                                                     event_microbiological_specimen_blood$Event_Microbiological_Specim_TK &
                                                                     !(event_microbiological_result$Growth_density %in% "") &
                                                                     event_microbiological_result$Method=="culture",]

event_microbiological_result_blood <- event_microbiological_result_blood[,c("Event_Microbiological_Specim_TK", "Microbe")]

ss_c <- event_microbiological_result_blood
ss_c <- unique(ss_c)
ss_c <- ss_c[order(ss_c$Event_Microbiological_Specim_TK),]
ss_c <- within(ss_c, {
  microbe_result <- ave(Microbe,Event_Microbiological_Specim_TK, FUN=function(x) Reduce(paste, x)  )
})
event_microbiological_result_blood_flat<-unique(ss_c[,c("Event_Microbiological_Specim_TK", "microbe_result")] )
ss_c <- NULL

event_microbiological_specimen_blood <- merge(event_microbiological_specimen_blood,event_microbiological_result_blood_flat,
                                              by= "Event_Microbiological_Specim_TK", all.x = T)


event_microbiological_specimen_blood <- merge(event_microbiological_specimen_blood, admission[,c("ICU_Admission_TK",
                                                                                                 "ICU_ID_from_datasource")], all.x = T)

event_microbiological_specimen_blood$Index_ICU_admittance_datetime2 <-strptime(event_microbiological_specimen_blood$Index_ICU_admittance_datetime,"%Y-%m-%d")
event_microbiological_specimen_blood$Specimen_datetime <- strptime(event_microbiological_specimen_blood$Specimen_datetime,"%Y-%m-%d")

event_microbiological_specimen_blood$culture_day <- round(difftime(as.POSIXct(event_microbiological_specimen_blood$Specimen_datetime,format="%Y-%m-%d"),
                                              as.POSIXct(event_microbiological_specimen_blood$Index_ICU_admittance_datetime2,format="%Y-%m-%d"),units="days"))


#code made by Joe for pathogens#### 
#this is the ICU_admition_data
load("Documents/BSI/R_code/original_data/df_biomarker_BSI_second_time_from_Hessel.RData")

#exclude the ones are in other hospital
df_biomarker_wide_erik33 <- subset(df_biomarker_wide_erik3, !(MARSID %in% c(176, 290, 1234, 1954, 2166, 3865, 4067, 6025, 6689, 14526)))

subset_MARSID <- df_biomarker_wide_erik33$MARSID[df_biomarker_wide_erik33$Microbe_groups == "Bacteremia"]

ss <- event_microbiological_specimen_blood[event_microbiological_specimen_blood$ICU_ID_from_datasource %in% subset_MARSID, ]
length(unique(ss$ICU_ID_from_datasource))

ss01 <- ss[ss$culture_day == 0 | ss$culture_day == 1 | ss$culture_day == -1, ]

ss01 <- ss01[!is.na(ss01$microbe_result),]
length(unique(ss01$ICU_ID_from_datasource))

diff_marsid_to_icu <- setdiff(subset_MARSID, ss01$ICU_ID_from_datasource)
#1063  1906  2129  2405  2632  2914  3652  4596 10006 10369 10737 11688 13805 were in df_biomarker_wide_erik3 bacteria group but not in ss01

length(unique(ss01$ICU_ID_from_datasource))

ssx <- ss01[,c("microbe_result","ICU_ID_from_datasource")]

ssx <- unique(ssx)

result <- aggregate(microbe_result ~ ICU_ID_from_datasource, data = ssx, FUN = function(x) paste(x, collapse = ", "))

result$Klebsiella <-  as.numeric(grepl("Klebsi", result$microbe_result ))
result$S.aureus <-  as.numeric(grepl("aureus", result$microbe_result ))
result$CoNS <-  as.numeric(grepl("epidermidi|haemolytic|hominis|Coagulase", result$microbe_result ))
result$Staph_unknown <- as.numeric(grepl("gram pos. coc, wrsch. Staphilococcus spp.", result$microbe_result ))
result$Enterococcus <- as.numeric(grepl("nterococcus|Enteroccus", result$microbe_result ))
result$E.coli <- as.numeric(grepl("coli", result$microbe_result ))
result$Streptococcus <- as.numeric(grepl("trepto", result$microbe_result ))



result$microbenum <- rowSums(result[,3:9])

table(result$microbenum)

result0 <- result[result$microbenum==0,]
