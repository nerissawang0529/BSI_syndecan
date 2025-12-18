# Clear workspace
rm(list = ls())

# Load required libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, rms, pheatmap, gridExtra, grid, cowplot)

# Load the dataset
merged_data <- read.csv("original_data/merged_data.csv")
colnames(merged_data)[colnames(merged_data) == "CoagulationFactorIII"] <- "Tissue factor"
colnames(merged_data)[colnames(merged_data) == "Platelets_value_1"] <- "Platelets"
colnames(merged_data)[colnames(merged_data) == "PT_Max_24h"] <- "PT"

results_2 <- read.csv("original_data/linear_nonlinear_p_value.csv")

# List of biomarkers (excluding Syndecan1, which is the predictor)
#variables_Inflammation_combine
markers <- c("Syndecan1","CD163", "MMP8", "NGAL", "Procalcitonin", "RAGE", "TenascinC", "TREM1","IL10", "IL18", "IL1ra", "IL6", "IL8")
#variables_Endothelial
markers <- c("Syndecan1","Ang1", "Ang2", "CX3CL1", "Thrombomodulin", "ESM1")
#variables_Coagulation
markers <- c("Syndecan1","Tissue factor", "Ddimer", "Platelets",  "PT", "Antitrombin", "PROC")


# Select relevant columns
markers_data <- merged_data[,c("Syndecan1","Antitrombin", "PROC", "Ang1", "Ang2", "CD163", "Tissue factor",
                               "CX3CL1", "Ddimer", "ESM1", "IL10", "IL18", "IL1ra", "IL6", 
                               "IL8", "MMP8", "NGAL", "Procalcitonin", "RAGE", "TenascinC", "Thrombomodulin",
                               "TREM1", "Platelets", "PT","SOFAtot","APACHE_IV_Score",
                               "gender","ckd","Cerebrovascular_disease","diabetes","age_yrs","Immune_deficiency","Past_myocardial_infarction")]

#variables_Inflammation_combine
markers_data <- merged_data[,c("Syndecan1","CD163", "MMP8", "NGAL", "Procalcitonin", "RAGE", "TenascinC", "TREM1","IL10", "IL18", "IL1ra", "IL6", "IL8",
                               "SOFAtot","APACHE_IV_Score.x",
                               "gender","ckd","Cerebrovascular_disease","diabetes","age_yrs","Immune_deficiency","Past_myocardial_infarction")]

#variables_Endothelial
markers_data <- merged_data[,c("Syndecan1","Ang1", "Ang2", "CX3CL1", "Thrombomodulin", "ESM1",
                               "SOFAtot","APACHE_IV_Score.x",
                               "gender","ckd","Cerebrovascular_disease","diabetes","age_yrs","Immune_deficiency","Past_myocardial_infarction")]

#variables_Coagulation
markers_data <- merged_data[,c("Syndecan1","Tissue factor", "Ddimer", "Platelets","PT", "Antitrombin", "PROC",
                               "SOFAtot_new.x","APACHE_IV_Score.x",
                               "gender","ckd","Cerebrovascular_disease","diabetes","age_yrs","Immune_deficiency","Past_myocardial_infarction")]


# Log-transform the biomarkers (avoid log(0) issues)
markers_data <- markers_data %>%
  mutate(across(all_of(markers), ~ log(. + 1)))

merged_data <- markers_data  # Overwrite for easier reference

# Function to generate heatmaps without white gaps and with a SINGLE black frame around each bar
generate_heatmap <- function(marker, model) {
  range_values <- c(-0.05, 0.05)  # Ensure small values are visible
  breaks <- seq(range_values[1], range_values[2], length.out = 100)
  colors <- colorRampPalette(c("darkgreen", "yellow", "darkred"))(99)
  
  # Fit the model based on whether it's linear or nonlinear
  if (model == "Linear") {
    fit <- lm(reformulate("Syndecan1", marker), data = merged_data)
  } else if (model == "Nonlinear") {
    fit <- lm(reformulate(paste0("rcs(Syndecan1, 3)"), marker), data = merged_data)
  }
  
  # Generate predictions
  new_data <- data.frame(Syndecan1 = seq(min(merged_data$Syndecan1), max(merged_data$Syndecan1), length.out = 100))
  predictions <- predict(fit, newdata = new_data)
  diff_predictions <- diff(predictions)  # Compute differences
  diff_predictions <- as.data.frame(t(diff_predictions))  # Transpose
  
  # Prevent white gaps by ensuring no NA values
  diff_predictions[abs(diff_predictions) < 0.001] <- 0  # Set small values to 0 instead of NA
  
  # Generate heatmap without internal borders
  heatmap <- pheatmap(
    diff_predictions,
    scale = "none",
    show_rownames = FALSE,
    show_colnames = FALSE,
    fontsize = 10,
    border_color = NA,  # NO individual cell borders
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = colors,
    breaks = breaks,
    legend = FALSE,  
    width = 2, height = 15, dpi = 600
  )
  
  # Convert to grob for alignment
  heatmap_grob <- heatmap$gtable
  
  # Create a marker name grob with smaller font
  marker_label <- textGrob(marker, gp = gpar(fontsize = 10, fontface = "bold"))
  
  # Draw a single black rectangle around each heatmap bar
  frame_grob <- rectGrob(gp = gpar(col = "black", fill = NA, lwd = 2))  # Thick black border
  
  # Combine label, heatmap, and black frame
  combined_grob <- arrangeGrob(
    marker_label, 
    gTree(children = gList(heatmap_grob, frame_grob)), 
    ncol = 2, widths = c(1, 4)
  )
  
  return(combined_grob)
}

# Extract markers and their model types
filtered_results <- results_2[results_2$Marker %in% markers, ]
marker_list <- filtered_results$Marker
model_types <- filtered_results$Model

# Generate labeled heatmaps
heatmaps <- mapply(generate_heatmap, marker_list, model_types, SIMPLIFY = FALSE)

# Combine plots into a single column without gaps
heatmap_layout <- grid.arrange(grobs = heatmaps, ncol = 1, heights = rep(1, length(heatmaps)))  # Ensure equal spacing


png(file="figures/heatmap.png",
    width=600, height=700)
grid.arrange(grobs = heatmaps, ncol = 1, heights = rep(1, length(heatmaps)))
dev.off()


# Generate a unified legend
legend_matrix <- matrix(seq(-0.05, 0.05, length.out = 100), ncol = 1)  # Legend matrix
legend_plot <- pheatmap(
  legend_matrix,
  color = colorRampPalette(c("darkgreen", "yellow", "darkred"))(99),
  legend = TRUE,
  cluster_rows = FALSE,  
  cluster_cols = FALSE,
  border_color = NA,  # Remove grid lines
  show_rownames = FALSE,
  show_colnames = FALSE
)

# Display the legend
grid.draw(legend_plot$gtable)


#get the plot
ggplot(markers_data, aes(x = Syndecan1, y = IL8)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1.5) +
  theme_minimal()
ggplot(markers_data, aes(x = Syndecan1, y = NGAL)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE, color = "red", size = 1.5) +
  theme_minimal()
