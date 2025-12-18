rm(list = ls())
install.packages("rms")
library(rms)
install.packages("FSA") 
library(FSA)
install.packages("VGAM")

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


#ready data
data <- read.csv("original_data/merged_data.csv")
#check subsofa score
##new sofa from hessel
load("original_data/MARS_SOFAS.RData")

library(dplyr)

#get new sofa from hessel within the IDs I need
sofas16_filtered <- sofas16 %>%
  filter(ICU_ID_from_datasource %in% data$MARSID)

# Define the columns to rename
cols_to_rename <- c("SOFAtot", "SOFA_Circulation","SOFA_Liver", "SOFA_Coagulation",
                    "SOFA_Live", "SOFA_Renal", "SOFA_Respiration")
# Find indices of these columns in the dataframe
indices <- which(names(sofas16_filtered) %in% cols_to_rename)
# Append "_new" to their names
names(sofas16_filtered)[indices] <- paste0(names(sofas16_filtered)[indices], "_new")

library(dplyr)

sofas16_filtered <- sofas16_filtered %>% 
  filter(diff_admin_daily %in% c(0, 1))

sofas16_filtered$Index_ICU_admittance_datetime <- NULL
sofas16_filtered$ICU_Admission_TK <- NULL
sofas16_filtered$Assessment_datetime <- NULL
sofas16_filtered$diff_admin_daily <- NULL

data <- merge(data, sofas16_filtered, 
              by.x = "MARSID", 
              by.y = "ICU_ID_from_datasource", 
              all.x = TRUE)

data <- data %>%
  dplyr::select(-matches("^SOFA(?!.*_new)", perl = TRUE))


# Write the data frame to a CSV file
destination_folder <- "original_data" 
export_file_name <- "final_data.csv" 
write.csv(data, file.path(destination_folder, export_file_name), row.names = FALSE)
cat("File exported to:", file.path(destination_folder, export_file_name), "\n")
