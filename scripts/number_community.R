rm(list = ls())

data <- read.csv("original_data/clinical_marker_scource_pathogen.csv")

selected_merged_data <- data[, c("ICU_ID_from_datasource","microbe_result","Klebsiella","S.aureus",
                                        "Enterococcus","E.coli","Streptococcus","Pseudomonas","Enterobacter","other_culture",
                                        "First_Hospital_admittance_date",
                                        "First_ICU_admittance_date","Index_ICU_discharge_datetime","Final_ICU_discharge_date","Final_Hospital_discharge_date","Index_hospital_discharge_date"
)]

library(dplyr)

# Assuming selected_merged_data is your dataframe
# Convert the columns to Date format
selected_merged_data <- selected_merged_data %>%
  mutate(
    First_Hospital_admittance_date = as.Date(First_Hospital_admittance_date),
    First_ICU_admittance_date = as.Date(First_ICU_admittance_date)
  )

# Calculate the difference in days
selected_merged_data <- selected_merged_data %>%
  mutate(diff_data = as.numeric(First_ICU_admittance_date - First_Hospital_admittance_date)) %>%
  # Filter out rows where the difference is greater than 2
  filter(diff_data <= 2) %>%
  # Optionally remove the diff_data column
  select(-diff_data)

