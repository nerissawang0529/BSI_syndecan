# Load required libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, rms, pheatmap, gridExtra, grid, cowplot)

# Clear workspace
rm(list = ls())

# Load the dataset
merged_data <- read.csv("original_data/merged_data.csv")
colnames(merged_data)[colnames(merged_data) == "CoagulationFactorIII"] <- "Tissue factor"
colnames(merged_data)[colnames(merged_data) == "Platelets_value_1"] <- "Platelets"
colnames(merged_data)[colnames(merged_data) == "PT_Max_24h"] <- "PT"

results_2 <- read.csv("original_data/linear_nonlinear_p_value.csv")

# List of biomarkers (excluding Syndecan1, which is the predictor)
markers <- c("Syndecan1","Antitrombin", "PROC", "Ang1", "Ang2", "CD163", "Tissue factor",
             "CX3CL1", "Ddimer", "ESM1", "IL10", "IL18", "IL1ra", "IL6", 
             "IL8", "MMP8", "NGAL", "Procalcitonin", "RAGE", "TenascinC", "Thrombomodulin",
             "TREM1", "Platelets", "PT")

# Select relevant columns
markers_data <- merged_data[,c("Syndecan1","Antitrombin", "PROC", "Ang1", "Ang2", "CD163", "Tissue factor",
                               "CX3CL1", "Ddimer", "ESM1", "IL10", "IL18", "IL1ra", "IL6", 
                               "IL8", "MMP8", "NGAL", "Procalcitonin", "RAGE", "TenascinC", "Thrombomodulin",
                               "TREM1", "Platelets", "PT","SOFAtot_new.x","APACHE_IV_Score.x",
                               "gender","ckd","Cerebrovascular_disease","diabetes","age_yrs","Immune_deficiency","Past_myocardial_infarction")]

# Log-transform the biomarkers (avoid log(0) issues)
markers_data <- markers_data %>%
  mutate(across(all_of(markers), ~ log(. + 1)))

merged_data <- markers_data  # Overwrite for easier reference

# Function to generate heatmaps using predictions normalized as z-scores
generate_heatmap <- function(marker, model) {
  range_values <- c(-2, 2)  # Standardized range for z-scores
  breaks <- seq(range_values[1], range_values[2], length.out = 100)
  colors <- colorRampPalette(c("darkblue", "blue", "gray90", "red", "darkred"))(99)
  
  
  merged_data$level <- merged_data[,marker]
  merged_data$levelz <- scale(merged_data$level)
  # Fit the model based on whether it's linear or nonlinear
  if (model == "Linear") {
    fit <- lm(levelz ~ Syndecan1, data = merged_data)
  } else if (model == "Nonlinear") {
    fit <- lm(levelz ~ rcs(Syndecan1, 3), data = merged_data)
  }
  
  # Generate predictions
  max_synd <- mean(merged_data$Syndecan1) + 3*sd(merged_data$Syndecan1)
  min_synd <- mean(merged_data$Syndecan1) - 3*sd(merged_data$Syndecan1)
  
  new_data <- data.frame(Syndecan1 = seq(min_synd, max_synd, length.out = 100))
  predictions <- predict(fit, newdata = new_data)
  
  # Normalize predictions to z-score
  predictions_z <- predictions
  
  # Convert to dataframe format
  predictions_df <- as.data.frame(t(predictions_z))  # Transpose
  
  # Generate heatmap without internal borders
  heatmap <- pheatmap(
    predictions_df,
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
marker_list <- results_2$Marker
model_types <- results_2$Model

# Generate labeled heatmaps
heatmaps <- mapply(generate_heatmap, marker_list, model_types, SIMPLIFY = FALSE)

# Combine plots into a single column without gaps
heatmap_layout <- grid.arrange(grobs = heatmaps, ncol = 1, heights = rep(1, length(heatmaps)))  # Ensure equal spacing

# Generate a unified legend
legend_matrix <- matrix(seq(-2, 2, length.out = 100), ncol = 1)  # Legend matrix
legend_plot <- pheatmap(
  legend_matrix,
  color = colorRampPalette(c("blue", "lightgray", "red"))(99),
  legend = TRUE,
  cluster_rows = FALSE,  
  cluster_cols = FALSE,
  border_color = NA,  # Remove grid lines
  show_rownames = FALSE,
  show_colnames = FALSE
)

# Display the legend
grid.draw(legend_plot$gtable)
