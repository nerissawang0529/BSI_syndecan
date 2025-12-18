rm(list = ls())

## install.packages("pacman")
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)
install.packages("dplyr")
library(dplyr)

#this is the ICU_admition_data
load("Documents/BSI/R_code/original_data/df_biomarker_BSI_second_time_from_Hessel.RData")
#this is the ICU_aquried_data
load("Documents/BSI/R_code/original_data/df_biomarker_BSI_ICUA.RData")

##Ready studydata(MARS database)
studydata <- import("Documents/ER_WARD_ICU/R_code_Syndecan_ER,WARD,ICU_after02052024/02052024/Original_Data/MARS database.xlsx")


# if df_biomarker_wide_erik3 are present in event_microbiological_specimen_blood####
all_in_icu <- all(df_biomarker_wide_erik3$MARSID %in% event_microbiological_specimen_blood$ICU_ID_from_datasource)
# Print the result
if (all_in_icu) {
  print("All MARSID values are present in ICU_ID_from_datasource.")
} else {
  print("Not all MARSID values are present in ICU_ID_from_datasource.")
}
# Find which MARSID values are not present in ICU_ID_from_datasource
missing_marsid <- df_biomarker_wide_erik3$MARSID[!(df_biomarker_wide_erik3$MARSID %in% event_microbiological_specimen_blood$ICU_ID_from_datasource)]
# Check if there are any missing values and print the result
if (length(missing_marsid) == 0) {
  print("All MARSID values are present in ICU_ID_from_datasource.")
} else {
  print("The following MARSID values are not present in ICU_ID_from_datasource:")
  print(missing_marsid)
}
#which df_biomarker_wide_erik3 are not present in event_microbiological_specimen_blood, they are non-infection
# 567   579  1378  1850  1967  2173  4463  4548  5752  6620  6689  6826  8209 10056 10062 10962 11815 12127 13484 13958####


# if df_biomarker_wide_erik3 are present in studydata####
all_included <- all(df_biomarker_wide_erik3$MARSID %in% studydata$ICU_ID_from_datasource)
# Output result
if (all_included) {
  print("All MARSID values are present in ICU_ID_from_datasource")
} else {
  print("Not all MARSID values are present in ICU_ID_from_datasource")
}
missing_marsid <- df_biomarker_wide_erik3$MARSID[!(df_biomarker_wide_erik3$MARSID %in% studydata$ICU_ID_from_datasource)]
# Check if there are any missing values and print the result
if (length(missing_marsid) == 0) {
  print("All MARSID values are present in ICU_ID_from_datasource.")
} else {
  print("The following MARSID values are not present in ICU_ID_from_datasource:")
  print(missing_marsid)
}
#which df_biomarker_wide_erik3 are not present in studydata
# 176   290  1234  1954  2166  3865  4067  6025  6689 14526####these patients are from other hospital, so we can exclude these for now.####

#code made by Joe for pathogens#### 
ss <- event_microbiological_specimen_blood[event_microbiological_specimen_blood$ICU_ID_from_datasource %in% df_biomarker_wide_erik3$MARSID, ]
length(unique(ss$ICU_ID_from_datasource))

ss01 <- ss[ss$culture_day == 0 | ss$culture_day == 1 | ss$culture_day == -1, ]

ss01 <- ss01[!is.na(ss01$microbe_result),]

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


#exclude the ones are in other hospital
result$ICU_ID_from_datasource
result <- subset(result, !(ICU_ID_from_datasource %in% c(176, 290, 1234, 1954, 2166, 3865, 4067, 6025, 6689, 14526)))


#export result

destination_folder <- "Documents/BSI/R_code/original_data" 
export_file_name <- "pathegon_result.csv" 
write.csv(result,  file.path(destination_folder, export_file_name), row.names = FALSE)
cat("File exported to:", file.path(destination_folder, export_file_name), "\n")
result$ICU_ID_from_datasource




####merge the pathogen with clinical data####
clinical_pathogen <- merge(df_biomarker_wide_erik3, result, by.x = "MARSID", by.y = "ICU_ID_from_datasource", all = TRUE)
otherhospital <- c(176, 290, 1234, 1954, 2166, 3865, 4067, 6025, 6689, 14526)
clinical_pathogen <- clinical_pathogen[!(clinical_pathogen$MARSID %in% otherhospital), ]
selected_data <- clinical_pathogen %>%
  select(MARSID, Microbe_groups, culture_day,49:58)

#1063,1906,2129,2405,2632,2914,3652,4596,10006,10369,10737,11688,13805