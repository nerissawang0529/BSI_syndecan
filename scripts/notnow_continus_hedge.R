rm(list = ls())
# install.packages("pacman")
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)
library(rms)  # Load the package
# Load required libraries
library(rms)  # For restricted cubic splines
library(pheatmap)

# Read the dataset
merged_data <- read.csv("original_data/final_data.csv")
linear_non_p<- read.csv("original_data/linear_nonlinear_p_value.csv")

# List of biomarkers (excluding Syndecan1, which is the predictor)
markers <- c("Syndecan1","Antitrombin", "PROC", "Ang1", "Ang2", "CD163", "CoagulationFactorIII",
             "CX3CL1", "Ddimer", "ESM1", "IL10", "IL18", "IL1ra", "IL6", 
             "IL8", "MMP8", "NGAL", "Procalcitonin", "RAGE", "TenascinC", "Thrombomodulin",
             "TREM1", "Platelets_value_1", "PT_Max_24h")

markers_data <- merged_data[,c("Syndecan1","Antitrombin", "PROC", "Ang1", "Ang2", "CD163", "CoagulationFactorIII",
                               "CX3CL1", "Ddimer", "ESM1", "IL10", "IL18", "IL1ra", "IL6", 
                               "IL8", "MMP8", "NGAL", "Procalcitonin", "RAGE", "TenascinC", "Thrombomodulin",
                               "TREM1", "Platelets_value_1", "PT_Max_24h","SOFAtot_new","APACHE_IV_Score",
                               "gender","ckd","Cerebrovascular_disease","diabetes","age_yrs","Immune_deficiency","Past_myocardial_infarction")]
biomarker_domains <- list(
  "Systemic Inflammation and Organ Damage" = c("CD163", "MMP8", "NGAL", "Procalcitonin", "RAGE", "TenascinC", "TREM1"),
  "Cytokines" = c("IL10", "IL18", "IL1ra", "IL6", "IL8"),
  "Coagulation Activation" = c("CoagulationFactorIII", "Ddimer", "Platelets_value_1",  "PT_Max_24h", "Antitrombin", "PROC"),
  "Endothelial Cell Activation and Dysfunction" = c("Ang1", "Ang2", "CX3CL1", "Thrombomodulin", "ESM1")
)

# Log-transform  the biomarkers
markers_data <- markers_data %>%
  mutate(across(all_of(markers), ~ log(. + 1))) # Adding 1 to avoid log(0) issues

merged_markers_log_studydata <- markers_data

linear_model <- lm(RAGE ~ Syndecan1, data = merged_markers_log_studydata)


# Create a sequence of Syndecan1 values
new_data <- data.frame(Syndecan1 = seq(6.5,12,by=0.5))

# Predict values for the given Syndecan1 levels
predictions <- predict(linear_model, newdata = new_data)

# Calculate the difference
diff_predictions <- diff(predictions)
# Print the result
print(diff_predictions)

diff_predictions <- as.data.frame(t(diff_predictions))

# Define breaks and colors for the heatmap

# Load required libraries
library(pheatmap)

# Check the range of your data
range_values <- range(diff_predictions, na.rm = TRUE)

# Ensure breaks are unique
breaks <- pretty(range_values, n = 50) # Adjust 'n' as needed to prevent duplicate values

# Define color scale correctly
colors <- colorRampPalette(c("blue", "white", "red"))(length(breaks) - 1)

# Generate heatmap
pheatmap(diff_predictions,
         scale = "none",
         show_rownames = FALSE,
         show_colnames = FALSE,
         fontsize = 25,
         border_color = "NA",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colors,
         breaks = breaks,
         legend = FALSE)



# Generate predictions for specific Syndecan1 range
predictions_l <- data.frame(
  Syndecan1 = seq(6.9, 11.2, by = 0.01),
  RAGE = predict(linear_model, newdata = data.frame(Syndecan1 = seq(6.9, 11.2, by = 0.01)))
)

# Add points to an existing plot (assuming a base plot exists)
plot(predictions_l$Syndecan1, predictions_l$RAGE, 
     col = "blue", pch = 16, 
     xlab = "Syndecan1", ylab = "RAGE",
     main = "Restricted Cubic Spline Predictions")

# Plot original data
plot(merged_markers_log_studydata$Syndecan1, merged_markers_log_studydata$RAGE, 
     ylab = "RAGE", xlab = "Syndecan1", 
     main = paste())

# Overlay predicted values
points(predictions_l$Syndecan1, predictions_l$RAGE, col = "blue", pch = 16)




###spline
# Fit a restricted cubic spline model
spline_model <- lm(CX3CL1 ~ rcs(Syndecan1, 3), data = merged_markers_log_studydata)

# Generate new data with smaller intervals for Syndecan1
new_data <- data.frame(Syndecan1 = seq(min(merged_markers_log_studydata$Syndecan1),
                                       max(merged_markers_log_studydata$Syndecan1),
                                       length.out = 100))  # Increase points for smooth gradient

# Predict values using the fitted spline model
predictions <- predict(spline_model, newdata = new_data)

# Calculate the difference in predictions
diff_predictions <- diff(predictions)

# Convert to data frame and transpose
diff_predictions <- as.data.frame(t(diff_predictions))

# Define breaks and colors for the heatmap with more resolution
range_values <- range(diff_predictions)
breaks <- seq(range_values[1], range_values[2], length.out = 200)  # More breaks for finer gradient
colors <- colorRampPalette(c("blue", "white", "red"))(199)  # More colors for smoother transition

# Generate heatmap with a smooth gradient
pheatmap(diff_predictions,
         scale = "none",
         show_rownames = FALSE,  # Hide row names
         show_colnames = FALSE,  # Hide column names
         fontsize = 25,
         border_color = NA,  
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colors,
         breaks = breaks,
         legend = TRUE)

# Generate predictions for specific Syndecan1 range
predictions_nl <- data.frame(
  Syndecan1 = seq(6.9, 11.2, by = 0.01),
  CX3CL1 = predict(spline_model, newdata = data.frame(Syndecan1 = seq(6.9, 11.2, by = 0.01)))
)

# Add points to an existing plot (assuming a base plot exists)
plot(predictions_nl$Syndecan1, predictions_nl$CX3CL1, 
     col = "red", pch = 16, 
     xlab = "Syndecan1", ylab = "CX3CL1",
     main = "Restricted Cubic Spline Predictions")

# Plot original data
plot(merged_markers_log_studydata$Syndecan1, merged_markers_log_studydata$CX3CL1, 
     ylab = "CX3CL1", xlab = "Syndecan1", 
     main = paste())

# Overlay predicted values
points(predictions_nl$Syndecan1, predictions_nl$CX3CL1, col = "red", pch = 16)

