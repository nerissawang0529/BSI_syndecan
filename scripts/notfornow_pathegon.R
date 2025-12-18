rm(list = ls())

## install.packages("pacman")
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)
install.packages("haven")
library(haven)

## ready data
admission <- read_sas("Documents/BSI/R_code/original_data/icu_admission.sas7bdat")
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
#this is the ICU_addmition_data
clinical_marker_scource <- import("Documents/BSI/R_code/original_data/clinical_marker_scource.csv")

subset_MARSID <- clinical_marker_scource$MARSID[clinical_marker_scource$Microbe_groups == "Bacteremia"]

ss <- event_microbiological_specimen_blood[event_microbiological_specimen_blood$ICU_ID_from_datasource %in% subset_MARSID, ]
length(unique(ss$ICU_ID_from_datasource))

ss01 <- ss[ss$culture_day == 0 | ss$culture_day == 1 | ss$culture_day == -1, ]

ss01 <- ss01[!is.na(ss01$microbe_result),]
length(unique(ss01$ICU_ID_from_datasource))

diff_marsid_to_icu <- setdiff(subset_MARSID, ss01$ICU_ID_from_datasource)
#1906 2129 2405 2632 2914 3652  4596 10006 10369 10737 11688 13805 were in clinical_marker_scource bacteria group but not in ss01
#these IDs already excluded in code 'clinical_marker'

length(unique(ss01$ICU_ID_from_datasource))

ssx <- ss01[,c("microbe_result","ICU_ID_from_datasource")]

ssx <- unique(ssx)
length(unique(ssx$ICU_ID_from_datasource))
result <- aggregate(microbe_result ~ ICU_ID_from_datasource, data = ssx, FUN = function(x) paste(x, collapse = ", "))

# i want the bacteria group which number >10
#Klebsiella, S.aureus, CoNS, Enterococcus, E.coli, Streptococcus,Staph_unknown chosen by Joe,Staph_unknown <10.
#Pseudomonas,Enterobacter, clostridium, bacteroides made by myself, but Clostridium and Bacteroides were < 10.

result$Klebsiella <-  as.numeric(grepl("Klebsi", result$microbe_result ))
result$S.aureus <-  as.numeric(grepl("aureus", result$microbe_result ))
result$CoNS <-  as.numeric(grepl("epidermidi|haemolytic|hominis|Coagulase", result$microbe_result ))
result$Enterococcus <- as.numeric(grepl("nterococcus|Enteroccus", result$microbe_result ))
result$E.coli <- as.numeric(grepl("coli", result$microbe_result ))
result$Streptococcus <- as.numeric(grepl("trepto", result$microbe_result ))
result$Pseudomonas <- as.numeric(grepl("Pseudomonas", result$microbe_result ))
result$Enterobacter <- as.numeric(grepl("Enterobacter", result$microbe_result ))
#result$Clostridium <- as.numeric(grepl("Clostridium", result$microbe_result )) #6
#result$Bacteroides <- as.numeric(grepl("Bacteroides", result$microbe_result )) #4
#result$Staph_unknown <- as.numeric(grepl("gram pos. coc, wrsch. Staphilococcus spp.", result$microbe_result ))

result$microbenum <- rowSums(result[,3:10])
table(result$microbenum)
result0 <- result[result$microbenum==0,]

result$other_bacteria <- ifelse(result$microbenum == 0, 1, 0)

result0$microbe_result
#export clinical_marker_unique data
destination_folder <- "Documents/BSI/R_code/original_data" 
export_file_name <- "result0.csv" 
write.csv(result0, file.path(destination_folder, export_file_name), row.names = FALSE)
cat("File exported to:", file.path(destination_folder, export_file_name), "\n")

#check if the viral and fungal are only have viral and fungal
#viral
viral_marsid <- clinical_marker_scource %>% filter(Microbe_groups == "Viral") %>%select(MARSID)
# 4463 is not in the event_microbiological_specimen_blood, event_microbiological_specimen_blood$microbe_result for 4534,4810,4828,4938,4951,5259,5433,12376,12549,12629,12715 are NA
#fungal
fungal_marsid <- clinical_marker_scource %>% filter(Microbe_groups == "Fungal") %>%select(MARSID)
# 999, 2841, 3915, 5338, 5421,12503 (only fungal within -1 to 1 day)


### GN
result$GN <- 0
result$GN <- as.numeric(grepl("coli|Pseudomonas|Klebs|Haemophil|Neisseria|Salmonella|Fusobact|Bacteroides|Serratia|Citrobact|Enterobacter|negatieve", result$microbe_result))
result$GP <- 0
result$GP <- as.numeric(grepl("trepto|taph|Enterococcus|Enteroccus|Clostrid|Bacillus", result$microbe_result))


result[result$GN==0 & result$GP==0,]

table(result$GN, result$GP)

#merge together
clinical_marker_scource_pathogen <- merge(clinical_marker_scource, result, 
                                          by.x = "MARSID", by.y = "ICU_ID_from_datasource", 
                                          all = TRUE)


#export clinical_marker_unique data
destination_folder <- "Documents/BSI/R_code/original_data" 
export_file_name <- "clinical_marker_scource_pathogen.csv" 
write.csv(clinical_marker_scource_pathogen, file.path(destination_folder, export_file_name), row.names = FALSE)
cat("File exported to:", file.path(destination_folder, export_file_name), "\n")





#Question
#1.the valuable for viral?
#2.how the non-infection was made?