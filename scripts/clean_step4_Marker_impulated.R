rm(list = ls())

pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)
install.packages("mice")
install.packages("dplyr")
install.packages("tidyverse")
install.packages("VIM")

library(mice)
library(dplyr)
library(tidyverse)


# Read data
data <- read.csv("original_data/clinical_marker_scource_pathogen.csv")
## where is an outlier which did not be deleted by mahalobis distance
data <- data[which(data[["Syndecan.1.CD138..29."]] < 290000), ]#1038
#also need to exclude #715 because the pathogen is bacillus species
data <- data[data$ICU_ID_from_datasource != "715", ]
#change the #8821 into group enteroccus because the bacteria 2 is lactobacillus species
data[data$ICU_ID_from_datasource == 8821, "FinalGroup_pathegon"] <- "Enterococcus"
data[data$ICU_ID_from_datasource == 8821, "microbe_result"] <- "Enterococcus"

#rename the biomarkers
# Original column names to be shortened
original_names <- c("Ang.1..64.", "Ang.2..26.", "CD163..28.", "Coagulation.Factor.III.Tissue.Factor..62.", 
                    "CX3CL1.Fractalkine..46.", "Cystatin.C..75.", "D.dimer..43.", "ESM.1..22.", 
                    "IL.10..22.", "IL.18..78.", "IL.1ra.IL.1F3..30.", "IL.6..13.", "IL.8..18.", 
                    "MMP.8..27.", "NGAL..53.", "Procalcitonin..39.", "RAGE.AGER..45.", 
                    "Syndecan.1.CD138..29.", "Tenascin.C..35.", "Thrombomodulin.BDCA.3..47.", 
                    "TREM.1..65.")

# Shortened column names
short_names <- c("Ang1", "Ang2", "CD163", "CoagulationFactorIII", "CX3CL1", "CystatinC", "Ddimer", 
                 "ESM1", "IL10", "IL18", "IL1ra", "IL6", "IL8", "MMP8", "NGAL", "Procalcitonin", 
                 "RAGE", "Syndecan1", "TenascinC", "Thrombomodulin", "TREM1")

# Rename only the specified columns
colnames(data)[colnames(data) %in% original_names] <- short_names



#export clinical_marker_unique data
destination_folder <- "original_data" 
export_file_name <- "clinical_marker_scource_pathogen_short_markername.csv" 
write.csv(data, file.path(destination_folder, export_file_name), row.names = FALSE)
cat("File exported to:", file.path(destination_folder, export_file_name), "\n")





#Data from Hessel
data_moremarker <- import("original_data/LUMINEX_MARS_FINAL20150409_export_intern.csv")

#cleaning the data
# only choose the day which closed to the admission day
data_moremarker <- data_moremarker %>%
  mutate(admis.date = dmy(admis.date), 
         sample.date = dmy(sample.date))

data_moremarker <- data_moremarker %>%
  mutate(date_diff = abs(as.numeric(sample.date - admis.date)))

data_moremarker <- data_moremarker %>%
  group_by(MARSID) %>%
  slice_min(order_by = date_diff, n = 1) %>%
  ungroup()


#Data for IDs
IDs <- read.csv("original_data/clinical_marker_scource_pathogen_short_markername.csv")

#
# Check if all IDs from IDs$ICU_ID_from_datasource are in data$MARSID
all_present <- all(IDs$ICU_ID_from_datasource %in% data_moremarker$MARSID)

# Print the result
if (all_present) {
  print("All IDs in 'IDs$MARSID' are present in 'data_moremarker$MARSID'.")
} else {
  print("Not all IDs in 'IDs$MARSID' are present in 'data_moremarker$MARSID'.")
}


# Find IDs not present in data$MARSID
missing_ids <- IDs$ICU_ID_from_datasource[!IDs$ICU_ID_from_datasource %in% data_moremarker$MARSID]

# Print missing IDs
print("IDs not found in data_moremarker$MARSID:")
print(missing_ids) # missing data 37 

#the following is for the PCA sensitivity analysis####
missing_id_df <- data.frame(ICU_ID_from_datasource = missing_ids)

write.csv(
  missing_id_df,
  "original_data/AT_PROC_missing_ids.csv",
  row.names = FALSE
)
cat("Exported missing AT/PROC IDs to original_data/AT_PROC_missing_ids.csv\n")
#####


#merge the two dataframes
IDs <- IDs %>%
  mutate(MARSID = as.character(ICU_ID_from_datasource))
merged_data <- left_join(IDs, data_moremarker %>% dplyr::select(MARSID, Antitrombin, PROC), 
                         by = c("MARSID" = "MARSID"))


##impulation
marker_imputation <- merged_data
marker_imputation <- merged_data %>%
  dplyr::select(ICU_ID_from_datasource,
         Antitrombin, PROC, Ang1, Ang2, CD163, CoagulationFactorIII, CX3CL1, CystatinC, Ddimer, 
         ESM1, IL10, IL18, IL1ra, IL6, IL8, MMP8, NGAL, Procalcitonin, RAGE, TenascinC, Thrombomodulin, TREM1,Platelets_value_1,PT_Max_24h,
         APACHE_IV_Score, APACHE_IV_Acute_Physiology_Score)



str(marker_imputation)

# Convert commas to dots and then change to numeric
options(digits = 10)
# Convert commas to dots and then change to numeric
marker_imputation <- marker_imputation %>%
  mutate(
    Antitrombin = as.numeric(gsub(",", ".", Antitrombin)),
    PROC = as.numeric(gsub(",", ".", PROC))
  )

str(marker_imputation)

#Check Missing Data Pattern
colnames(marker_imputation)[colSums(is.na(marker_imputation)) > 0] #Antitrombin, PROC,CX3CL1,MMP8,Platelets_value_1,PT_Max_24h


##imputation
set.seed(0529)
### Create a new dataframe and set the ID column to 0
marker_imputation_check <- marker_imputation #for check if the final imputation is right
marker_imputation_ID_0 <- marker_imputation

# Remove ID column and convert character variables to factors
marker_imputation_ID_0 <- marker_imputation %>%
  dplyr::select(-ICU_ID_from_datasource) %>%
  mutate_if(is.character, as.factor)

### Initialize the mice object
imp <- mice(marker_imputation_ID_0, maxit=0)
predM <- imp$predictorMatrix
meth <- imp$method

### Specify the variables to impute
vars_to_impute <- c('Antitrombin','PROC','CX3CL1','MMP8','Platelets_value_1','PT_Max_24h' )

### Update the predictor matrix and imputation methods
predM[ , vars_to_impute] <- 1
meth[vars_to_impute] <- "cart"

### Perform the imputation process
marker_imputation_imp <- mice(data=marker_imputation_ID_0, m=20, maxit=20, meth=meth, pred=predM, printFlag=FALSE)
imputation_long <- complete(marker_imputation_imp, action="long", include = FALSE)

### Add the original ID column to the long format data
imputation_long$OriginalID <- rep(marker_imputation$ICU_ID_from_datasource, each = 20)

### Calculate the median imputed values for each original ID
median_imputed <- aggregate(imputation_long[, vars_to_impute], 
                            by = list(OriginalID = imputation_long$OriginalID), 
                            FUN = function(x) median(x, na.rm = TRUE))

### Sort by the original ID order
median_imputed <- median_imputed[order(match(median_imputed$OriginalID, marker_imputation$ICU_ID_from_datasource)), ]
### Ensure the ID order is consistent
stopifnot(all(median_imputed$OriginalID == marker_imputation$MARSID))

### Update the original dataset with the imputed median values
for (var in vars_to_impute) {
  is_na <- is.na(marker_imputation[[var]])
  marker_imputation[is_na, var] <- median_imputed[is_na, var]
}

### Restore the original ID column
marker_imputation$MARSID <- median_imputed$OriginalID


# Ensure ICU_ID_from_datasource is a character in both dataframes
marker_imputation$MARSID <- as.character(marker_imputation$MARSID)
merged_data$MARSID <- as.character(merged_data$MARSID)

# Select only the imputed columns along with the ID
imputed_values <- marker_imputation %>%
  dplyr::select(MARSID, Antitrombin, PROC, CX3CL1, MMP8, Platelets_value_1, PT_Max_24h)

# Update the missing values in merged_data
merged_data <- merged_data %>%
  dplyr::left_join(imputed_values, by = "MARSID", suffix = c("", "_imp"))

# Replace NAs in original columns with the imputed values
vars_to_impute <- c('Antitrombin','PROC','CX3CL1','MMP8','Platelets_value_1','PT_Max_24h')

for (var in vars_to_impute) {
  merged_data[[var]][is.na(merged_data[[var]])] <- merged_data[[paste0(var, "_imp")]][is.na(merged_data[[var]])]
}

# Drop the temporary imputed columns
merged_data <- merged_data %>%
  dplyr::select(-ends_with("_imp"))

# Check if all missing values were replaced
sum(is.na(merged_data[vars_to_impute])) # Should be 0 if all were imputed

# Convert commas to dots and then change to numeric
options(digits = 10)
# Convert commas to dots and then change to numeric
merged_data <- merged_data %>%
  mutate(
    Antitrombin = as.numeric(gsub(",", ".", Antitrombin)),
    PROC = as.numeric(gsub(",", ".", PROC))
  )



##export form
destination_folder <- "Original_Data/" 
export_file_name <- "merged_data.csv" 
write.csv(merged_data, file.path(destination_folder, export_file_name), row.names = FALSE)
cat("File exported to:", file.path(destination_folder, export_file_name), "\n")
