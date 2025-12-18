rm(list = ls ())

## install.packages("pacman")
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)
install.packages("haven")
library(haven)

#Ready data of olink
ICU_olink <- read.csv("Original_data/ICU_olink.csv")
ICU_olink$ICU_ID <- sub("^LRM10", "", ICU_olink$SampleID)


## ready data for pathogen
admission <- read_sas("original_data/icu_admission.sas7bdat")


############################################# microbiology result
event_microbiological_specimen <- read_sas("original_data/event_microbiological_specimen.sas7bdat", NULL)
#only blood
event_microbiological_specimen_blood <- event_microbiological_specimen[grepl("blood",event_microbiological_specimen$Material,ignore.case = T)==T,]
event_microbiological_specimen_blood <- event_microbiological_specimen_blood[!grepl("serum",event_microbiological_specimen_blood$Material,ignore.case = T)==T,]

#merge results
event_microbiological_result <- read_sas("original_data/event_microbiological_result.sas7bdat", NULL)
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


ss <- event_microbiological_specimen_blood[event_microbiological_specimen_blood$ICU_ID_from_datasource %in% ICU_olink$ICU_ID, ]
length(unique(ss$ICU_ID_from_datasource))
ss01 <- ss[ss$culture_day == 0 | ss$culture_day == 1 | ss$culture_day == -1, ]
length(unique(ss01$ICU_ID_from_datasource))

# Write the data frame to a CSV file
destination_folder <- "Original_data" 
export_file_name <- "ICU_pathogen_Iron.csv" 
write.csv(ss01, file.path(destination_folder, export_file_name), row.names = FALSE)
cat("File exported to:", file.path(destination_folder, export_file_name), "\n")

##8702 CNS 732 CNS, others is NA.
