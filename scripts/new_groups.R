rm(list = ls())

# install.packages("pacman")
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)
install.packages("diptest")
library(diptest)

#Ready patients data
merged_data <- read.csv("original_data/merged_data.csv")
##Ready list_of_healthy
list_of_healthy <- read.csv2("Original_Data/EB1.0 43-plex cleaned.csv")
list_of_healthy <- list_of_healthy[list_of_healthy$ID >= 2000 & list_of_healthy$ID <= 3000, ]
##add New_ID
list_of_healthy <- data.frame(New_ID = paste("Ward_EB_healthy", list_of_healthy$ID, sep = "_"), list_of_healthy)
##delete unnecessary columns
list_of_healthy$Protein_C <-NULL
list_of_healthy$Plaat <- NULL
list_of_healthy$ICU <- NULL
list_of_healthy$FerGr<- NULL
list_of_healthy$Sepsis<- NULL
list_of_healthy$Study<- NULL
list_of_healthy$Pentraxin_3 <- NULL
list_of_healthy$ID <-NULL
list_of_healthy$Factor_XIV_protein_C <-NULL
##rename the colnames
colnames(list_of_healthy)[colnames(list_of_healthy) == "Lipocalin_2_NGAL"] <- "NGAL"
##as numeric
list_of_healthy <- list_of_healthy %>%
  mutate_at(vars(-New_ID), as.numeric)
##CRP* 1000000
list_of_healthy <- list_of_healthy %>%
  mutate(CRP = CRP * 1000000)
list_of_healthy$DcR3 <- NULL
list_of_healthy$Cystatin_C <- NULL
list_of_healthy$Gas6 <- NULL
list_of_healthy$PCSK9 <- NULL
str(list_of_healthy)


#log
list_of_healthy <- list_of_healthy %>%
  mutate_at(vars(2:36), ~ log(.))

#mean and standard raw data
mean_healthy <- mean(list_of_healthy$Syndecan, na.rm = TRUE)
sd_healthy <- sd(list_of_healthy$Syndecan, na.rm = TRUE)

#z score for patients within healthy
merged_data$Z_Syndecan_raw <- (merged_data$Syndecan1 - mean_healthy) / sd_healthy

#z score for patients
#mean_patients <- mean(merged_data$syndecan_1, na.rm = TRUE)
#sd_patients <- sd(merged_data$syndecan_1, na.rm = TRUE)
#merged_data$Z_Syndecan_pa <- (merged_data$syndecan_1 - mean_patients) / sd_patients



library(ggplot2)

ggplot(merged_data, aes(x = Z_Syndecan_raw, y = Syndecan1)) +
  geom_point(color = "darkgreen", size = 2, alpha = 0.6) +  # Smaller, transparent points
  labs(
    x = "Syndecan-1 Z-score",
    y = "Syndecan-1 Level",
    title = "Scatter Plot of Syndecan-1 Levels vs. Z-Score"
  ) +
  scale_x_continuous(
    breaks = seq(0, max(merged_data$Z_Syndecan_raw, na.rm = TRUE), by = 2),  # Detailed x-axis every 1 unit
    minor_breaks = seq(0, max(merged_data$Z_Syndecan_raw, na.rm = TRUE), by = 0.5)  # Even smaller ticks every 0.5
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Centered title
    axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5),  # Rotate labels for clarity
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold")
  )






###higher than 3*sd, logged#### 

#healthy####
#Find the number of patients with healthy syndecan_1 > 3 * sd_healthy
num_high <- sum(merged_data$syndecan_1 > (mean_healthy + 3 * sd_healthy), na.rm = TRUE) 
# Create a new column 'sd_group' with 'High' if syndecan_1 > 3*sd_healthy, otherwise 'Low'
merged_data$sd_group_based_healthy <- ifelse(merged_data$syndecan_1 > (mean_healthy + 3 * sd_healthy), "High", "Low")
#161 high, 28 low####


# Categorize patients into three groups based on syndecan-1 levels####
merged_data$sd_group_based_healthy_3 <- cut(
  merged_data$syndecan_1, 
  breaks = c(-Inf, mean_healthy + 3 * sd_healthy, mean_healthy + 5 * sd_healthy, Inf), 
  labels = c("1", "2", "3"),
  right = FALSE
)
#28 low, 41 moderate, 120 high####
dip.test(merged_data$syndecan_original)










# Count the number of patients in each group
table(merged_data$sd_group_based_healthy_3)

# Write the data frame to a CSV file
destination_folder <- "original_data" 
export_file_name <- "merged_data_new_group.csv" 
write.csv(merged_data, file.path(destination_folder, export_file_name), row.names = FALSE)
cat("File exported to:", file.path(destination_folder, export_file_name), "\n")





#patients, not do this way####
merged_data$sd_based_patients <- cut(
  merged_data$syndecan_1,
  breaks = c(-Inf, mean_patients, mean_patients + sd_patients, Inf),
  labels = c("Low", "Medium", "High"),
  include.lowest = TRUE
)
#87 low, 77 medium, 25 high####


##kmeans#### only for syndecan-1 or for all biomarkers?####
# Force 3 clusters
set.seed(123) # for reproducibility
km <- kmeans(merged_data$syndecan_1, centers = 3)

# Add cluster labels
merged_data$kmeans <- factor(km$cluster, levels = 1:3, 
                               labels = c("Group1", "Group2", "Group3")) #the number is 47,47,95, are they ordered by syndecan-1?

#higher then 3 standard 

#Model-Based Clustering (Gaussian Mixture Models)####


#z scores within healthy or patients ####
#z score with healthy in patients
merged_data$Z_Syndecan1 <- (merged_data$syndecan_1 - mean_healthy) / sd_healthy

#z score in patients
mean_patients <- mean(merged_data$syndecan_1, na.rm = TRUE)
sd_patients <- sd(merged_data$syndecan_1, na.rm = TRUE)
merged_data$Z_Syndecan1_patients <- (merged_data$syndecan_1 - mean_patients) / sd_patients
# Categorize z-scores
merged_data$Category <- cut(merged_data$Z_Syndecan1_patients, 
                            breaks = c(-Inf, -1.65, 1.65, Inf), 
                            labels = c("Low", "Normal", "High")) 
#low 13, normal 166, high 10####