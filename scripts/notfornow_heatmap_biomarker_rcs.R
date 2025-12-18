# Load required libraries
library(rms)  # For restricted cubic splines
library(pheatmap)

# Fit a restricted cubic spline model
spline_model <- lm(Ang2 ~ rcs(Syndecan1, 3), data = merged_markers_log_studydata)

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
  Ang2 = predict(spline_model, newdata = data.frame(Syndecan1 = seq(6.9, 11.2, by = 0.01)))
)

# Add points to an existing plot (assuming a base plot exists)
plot(predictions_nl$Syndecan1, predictions_nl$Ang2, 
     col = "red", pch = 16, 
     xlab = "Syndecan1", ylab = "Ang2",
     main = "Restricted Cubic Spline Predictions")

# Plot original data
plot(merged_markers_log_studydata$Syndecan1, merged_markers_log_studydata$Ang2, 
     ylab = "Ang2", xlab = "Syndecan1", 
     main = paste())

# Overlay predicted values
points(predictions_nl$Syndecan1, predictions_nl$Ang2, col = "red", pch = 16)









#Ang2 0.01053214 0.04250369
#Ang2
#IL18
#NGAL
