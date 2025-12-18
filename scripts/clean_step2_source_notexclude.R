rm(list = ls())

#the data is from Hessel

## install.packages("pacman")
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)
library(haven)


#data for BSI 
data <- read.csv("original_data/clinical_marker_unique.csv")

#combine the ICU_ID_from_datasource into episode_of_sepsis
episode_of_sepsis <- read_sas("original_data/episode_of_sepsis.sas7bdat")
admission <- read_sas("original_data/icu_admission.sas7bdat")
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

all(data$MARSID %in% episode_of_sepsis2$ICU_ID_from_datasource)
#missing_ids_2 <- data$ICU_ID_from_datasource[!data$ICU_ID_from_datasource %in% episode_of_sepsis2$ICU_ID_from_datasource]
#3286  1063  3613 14497 #these IDs already deleted from data (the code is clinical_marker)

episode_of_sepsis_data <- episode_of_sepsis2[episode_of_sepsis2$ICU_ID_from_datasource %in% data$MARSID, ]
length(unique(episode_of_sepsis_data$ICU_ID_from_datasource))
table(episode_of_sepsis_data$site)

episode_of_sepsis_data$Is_Admission_Event <- NULL
episode_of_sepsis_data$Pre_sepsis_SOFA_CNS <- NULL
episode_of_sepsis_data$Pre_sepsis_SOFA_Circulation <- NULL
episode_of_sepsis_data$Pre_sepsis_SOFA_Respiration <- NULL
episode_of_sepsis_data$Pre_sepsis_SOFA_Renal <- NULL
episode_of_sepsis_data$Pre_sepsis_SOFA_Liver <- NULL
episode_of_sepsis_data$Pre_sepsis_SOFA_coagulation <- NULL

#Define groups based on criteria from Hessel:
#categorized into “pure” sources without coinfection or “mixed” infections (more than one source); 
#sources were respiratory, abdominal, urinary tract, cardiovascular, CNS, skin, unknown source, and other. 
#“Other” consisted of several categories with few patients (eTable 1, electronic supplemental material). 
episode_of_sepsis_data$respiratory <- ifelse(grepl("CAP", episode_of_sepsis_data$site, ignore.case=T) |
                                               grepl("HAP", episode_of_sepsis_data$site, ignore.case=T) |
                                               grepl("VAP", episode_of_sepsis_data$site, ignore.case=T) |
                                               grepl("lung absceepisode_of_sepsis_data", episode_of_sepsis_data$site, ignore.case=T) |
                                               grepl("lung abscess", episode_of_sepsis_data$site, ignore.case=T) |
                                               grepl("Tracheobronchitis", episode_of_sepsis_data$site, ignore.case=T), 1, 0)
episode_of_sepsis_data$abdominal <- ifelse(grepl("abdominal", episode_of_sepsis_data$site, ignore.case=T) |
                                             grepl("gastro", episode_of_sepsis_data$site, ignore.case=T), 1, 0)
episode_of_sepsis_data$cardiovascular <- ifelse(grepl("Cardiovascular : catheter-related", episode_of_sepsis_data$site, ignore.case=T) |
                                                  grepl("endocarditis", episode_of_sepsis_data$site, ignore.case=T) |
                                                  grepl("Myocarditis", episode_of_sepsis_data$site, ignore.case=T), 1, 0)
episode_of_sepsis_data$urinary <- ifelse(grepl("UTI", episode_of_sepsis_data$site, ignore.case=T) |
                                           grepl("urinary tract", episode_of_sepsis_data$site, ignore.case=T), 1, 0)
episode_of_sepsis_data$CNS_source <- ifelse(grepl("brain absceepisode_of_sepsis_data", episode_of_sepsis_data$site, ignore.case=T) |
                                              grepl("primary meningitis", episode_of_sepsis_data$site, ignore.case=T) |
                                              grepl("secondary meningitis", episode_of_sepsis_data$site, ignore.case=T) |
                                              grepl("brain abscess", episode_of_sepsis_data$site, ignore.case=T) |
                                              grepl("CNS : spinal absceepisode_of_sepsis_data", episode_of_sepsis_data$site, ignore.case=T), 1, 0)
episode_of_sepsis_data$skin <- ifelse(grepl("Postoperative wound", episode_of_sepsis_data$site, ignore.case=T) |
                                        grepl("skin", episode_of_sepsis_data$site, ignore.case=T), 1, 0)
episode_of_sepsis_data$other_source <- ifelse(grepl("bone", episode_of_sepsis_data$site, ignore.case=T) |
                                                grepl("mediastinitis", episode_of_sepsis_data$site, ignore.case=T) |
                                                grepl("eye", episode_of_sepsis_data$site, ignore.case=T) |
                                                grepl("ear", episode_of_sepsis_data$site, ignore.case=T) |
                                                grepl("sinusitis", episode_of_sepsis_data$site, ignore.case=T) |
                                                grepl("pharyngitis", episode_of_sepsis_data$site, ignore.case=T) |
                                                grepl("oral", episode_of_sepsis_data$site, ignore.case=T) |
                                                grepl("endometritis", episode_of_sepsis_data$site, ignore.case=T), 1, 0)
episode_of_sepsis_data$unknown_source <- ifelse(grepl("Cardiovascular : primary", episode_of_sepsis_data$site, ignore.case=T) |
                                                  grepl("unknown", episode_of_sepsis_data$site, ignore.case=T) |
                                                  grepl("viral", episode_of_sepsis_data$site, ignore.case=T)|
                                                  grepl("Other \\(please specify\\)", episode_of_sepsis_data$site, ignore.case=T), 1, 0)

install.packages("dplyr")
library(dplyr)

# Combine the same IDs, if the same source showed more than one, then just keep one
if ("reshape2" %in% loadedNamespaces()) detach("package:DOSE", unload = TRUE)
if ("plyr" %in% loadedNamespaces()) detach("package:qvalue", unload = TRUE)

library(dplyr)
library(purrr)  

# 
exclude_cols <- c("ICU_ID_from_datasource", "Group", "respiratory", "abdominal", "cardiovascular",
                  "urinary", "CNS_source", "skin", "other_source", "unknown_source")

# 
episode_of_sepsis_data_2 <- episode_of_sepsis_data %>%
  dplyr::group_by(ICU_ID_from_datasource) %>%
  dplyr::summarise(
    Group = paste(unique(Group), collapse = " "),
    
    respiratory = if (all(is.na(respiratory))) NA else max(respiratory, na.rm = TRUE),
    abdominal = if (all(is.na(abdominal))) NA else max(abdominal, na.rm = TRUE),
    cardiovascular = if (all(is.na(cardiovascular))) NA else max(cardiovascular, na.rm = TRUE),
    urinary = if (all(is.na(urinary))) NA else max(urinary, na.rm = TRUE),
    CNS_source = if (all(is.na(CNS_source))) NA else max(CNS_source, na.rm = TRUE),
    skin = if (all(is.na(skin))) NA else max(skin, na.rm = TRUE),
    other_source = if (all(is.na(other_source))) NA else max(other_source, na.rm = TRUE),
    unknown_source = if (all(is.na(unknown_source))) NA else max(unknown_source, na.rm = TRUE),
    dplyr::across(
      .cols = setdiff(names(episode_of_sepsis_data), exclude_cols) %>%
        purrr::keep(~ is.numeric(episode_of_sepsis_data[[.x]])),
      .fns = ~ if (all(is.na(.x))) NA else max(.x, na.rm = TRUE),
      .names = "max_{.col}"
    )
  )

  

#merge together
clinical_marker_scource <- merge(data, episode_of_sepsis_data_2, 
                                 by.x = "MARSID", 
                                 by.y = "ICU_ID_from_datasource", 
                                 all = TRUE)


#export clinical_marker_unique data
destination_folder <- "original_data" 
export_file_name <- "clinical_marker_scource.csv" 
write.csv(clinical_marker_scource, file.path(destination_folder, export_file_name), row.names = FALSE)
cat("File exported to:", file.path(destination_folder, export_file_name), "\n")
