linear_model <- lm(Antitrombin ~ Syndecan1, data = merged_markers_log_studydata)


# Create a sequence of Syndecan1 values
new_data <- data.frame(Syndecan1 = seq(min(merged_markers_log_studydata$Syndecan1),
                                       max(merged_markers_log_studydata$Syndecan1),
                                       length.out = 100))  # Increase points for smooth gradient


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
predictions_l <- data.frame(
  Syndecan1 = seq(6.9, 11.2, by = 0.01),
  Antitrombin = predict(linear_model, newdata = data.frame(Syndecan1 = seq(6.9, 11.2, by = 0.01)))
)

# Add points to an existing plot (assuming a base plot exists)
plot(predictions_l$Syndecan1, predictions_l$Antitrombin, 
     col = "blue", pch = 16, 
     xlab = "Syndecan1", ylab = "Antitrombin",
     main = "Restricted Cubic Spline Predictions")

# Plot original data
plot(merged_markers_log_studydata$Syndecan1, merged_markers_log_studydata$Antitrombin, 
     ylab = "Antitrombin", xlab = "Syndecan1", 
     main = paste())

# Overlay predicted values
points(predictions_l$Syndecan1, predictions_l$Antitrombin, col = "blue", pch = 16)

