#use the sofa score from Hessel
rm(list = ls())

## install.packages("pacman")
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)
install.packages("dplyr")
install.packages("FSA")
install.packages("rms")
# Install patchwork if not installed
if (!requireNamespace("patchwork", quietly = TRUE)) {
  install.packages("patchwork")
}
install.packages("ggrepel")  # If not already installed
library(ggrepel)
library(patchwork)
library(dplyr)
library(mice)
library(tidyr)
library(FSA)
library(rms)
library(ggplot2)
library(dplyr)  


#draw the lines one by one
# Read data
merged_data <- read.csv("original_data/merged_data.csv")
# Sort the dataframe by Syndecan1 from low to high
data_sorted <- merged_data %>% arrange(Syndecan1)

# List of markers (y-axis variables)
markers <- c("Antitrombin", "PROC", "Ang1", "Ang2", "CD163", "CoagulationFactorIII",
  "CX3CL1", "Ddimer", "ESM1", "IL10", "IL18", "IL1ra", "IL6", 
  "IL8", "MMP8", "NGAL", "Procalcitonin", "RAGE", "TenascinC", "Thrombomodulin",
  "TREM1", "Platelets_value_1", "PT_Max_24h")

# original####
for (marker in markers) {
  # Create the plot
  p <- ggplot(data_sorted, aes(x = Syndecan1, y = .data[[marker]])) +
    geom_point() +  # Add points
    geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), color = "blue") +  # Add a trend line
    labs(title = paste("Syndecan1 vs", marker),  # Add title
         x = "Syndecan1", 
         y = marker) +
    theme_minimal()  # Use a minimal theme
  
  # Display the plot
  print(p)
}


# log####
for (marker in markers) {
  # Create the plot
  p <- ggplot(data_sorted, aes(x = Syndecan1, y = .data[[marker]])) +
    geom_point() +  # Add points
    geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = FALSE, color = "blue") +  # Add a trend line
    scale_x_continuous(trans = "log", labels = scales::trans_format("log", math_format(e^.x))) +  # Natural log x-axis
    scale_y_continuous(trans = "log", labels = scales::trans_format("log", math_format(e^.x))) +  # Natural log y-axis
    labs(title = paste("Syndecan1 vs", marker),  # Add title
         x = "Syndecan1 (ln scale)", 
         y = paste(marker, "(ln scale)")) +
    theme_minimal()  # Use a minimal theme
  
  # Display the plot
  print(p)
}




#domian####
# Define biomarker domains
biomarker_domains <- list(
  "Systemic Inflammation and Organ Damage" = c("CD163", "MMP8", "NGAL", "Procalcitonin", "RAGE", "TenascinC", "TREM1"),#7
  "Cytokines" = c("IL10", "IL18", "IL1ra", "IL6", "IL8"),  #5
  "Coagulation Activation" = c("CoagulationFactorIII", "Ddimer", "Platelets_value_1",  "PT_Max_24h", "Antitrombin", "PROC"),#6
  "Endothelial Cell Activation and Dysfunction" = c("Ang1", "Ang2", "CX3CL1", "Thrombomodulin", "ESM1") #5
)

# with dots,log####
for (domain in names(biomarker_domains)) {
  
  markers <- biomarker_domains[[domain]]
  
  # Reshape data to long format for ggplot
  data_long <- data_sorted %>%
    select(Syndecan1, all_of(markers)) %>%
    pivot_longer(cols = -Syndecan1, names_to = "Marker", values_to = "Value")
  
  # Create a single plot with multiple trend lines for the domain
  p <- ggplot(data_long, aes(x = Syndecan1, y = Value, color = Marker)) +
    geom_point(alpha = 0.6) +  # Add points with slight transparency
    geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = FALSE) +  # Add trend lines per marker
    scale_x_continuous(trans = "log", labels = scales::trans_format("log", math_format(e^.x))) +  # Natural log x-axis
    scale_y_continuous(trans = "log", labels = scales::trans_format("log", math_format(e^.x))) +  # Natural log y-axis
    labs(title = paste(domain, "- Syndecan1 vs Biomarkers"),  # Add domain in title
         x = "Syndecan1 (ln scale)", 
         y = "Marker Values (ln scale)",
         color = "Markers") +
    theme_minimal()  # Use a minimal theme
  
  # Display the plot
  print(p)
}

# without dots,original####
for (domain in names(biomarker_domains)) {
  
  markers <- biomarker_domains[[domain]]
  
  # Reshape data to long format for ggplot
  data_long <- data_sorted %>%
    select(Syndecan1, all_of(markers)) %>%
    pivot_longer(cols = -Syndecan1, names_to = "Marker", values_to = "Value")
  
  # Create a single plot with only trend lines for the domain
  p <- ggplot(data_long, aes(x = Syndecan1, y = Value, color = Marker)) +
    geom_smooth(method = "loess", se = FALSE, linewidth = 1.2) +  # Add trend lines per marker
    #scale_x_continuous(trans = "log", labels = scales::trans_format("log", math_format(e^.x))) +  # Natural log x-axis
    #scale_y_continuous(trans = "log", labels = scales::trans_format("log", math_format(e^.x))) +  # Natural log y-axis
    labs(title = paste(domain, "- Syndecan1 vs Biomarkers"),  # Add domain in title
         x = "Syndecan1", 
         y = "Marker Values",
         color = "Markers") +
    theme_minimal()  # Use a minimal theme
  
  # Display the plot
  print(p)
}

#Normalize markers using Z-score####
# Normalize markers using Z-score, the log data for y and x axis#####
for (domain in names(biomarker_domains)) {
  
  markers <- biomarker_domains[[domain]]
  
  # Reshape data to long format for ggplot
  data_long <- data_sorted %>%
    select(Syndecan1, all_of(markers)) %>%
    pivot_longer(cols = -Syndecan1, names_to = "Marker", values_to = "Value") %>%
    mutate(Log_Value = log(Value + 1)) %>%  # Apply log transformation (log +1 to avoid log(0) issues)
    group_by(Marker) %>%
    mutate(Normalized_Log_Value = (Log_Value - mean(Log_Value, na.rm = TRUE)) / sd(Log_Value, na.rm = TRUE)) %>%
    ungroup()
  
  # Create a single plot with only trend lines for the domain
  p <- ggplot(data_long, aes(x = log(Syndecan1 + 1), y = Normalized_Log_Value, color = Marker)) +
    geom_smooth(method = "loess", se = FALSE, linewidth = 1.2) +  # Add trend lines per marker
    labs(title = paste(domain, "- Syndecan1 vs Biomarkers (Log-Normalized)"),  # Add domain in title
         x = "Log(Syndecan1 + 1)", 
         y = "Normalized Log Marker Values (Z-score)",
         color = "Markers") +
    theme_minimal()  # Use a minimal theme
  
  # Display the plot
  print(p)
}


# original data####
# Original data with marker labels alongside trend lines
for (domain in names(biomarker_domains)) {
  
  markers <- biomarker_domains[[domain]]
  
  # Reshape data to long format for ggplot
  data_long <- data_sorted %>%
    select(Syndecan1, all_of(markers)) %>%
    pivot_longer(cols = -Syndecan1, names_to = "Marker", values_to = "Value") %>%
    group_by(Marker) %>%
    mutate(Normalized_Value = (Value - mean(Value, na.rm = TRUE)) / sd(Value, na.rm = TRUE)) %>%
    ungroup()
  
  # Identify label positions (at max Syndecan1)
  label_positions <- data_long %>%
    group_by(Marker) %>%
    filter(Syndecan1 == max(Syndecan1, na.rm = TRUE)) %>%
    ungroup()
  
  # Create a single plot with only trend lines and marker labels
  p <- ggplot(data_long, aes(x = Syndecan1, y = Normalized_Value, color = Marker)) +
    geom_smooth(method = "loess", se = FALSE, linewidth = 1.2) +  # Add trend lines per marker
    geom_text_repel(data = label_positions, aes(label = Marker), hjust = 0, nudge_x = 500, size = 5, show.legend = FALSE) +  # Add marker labels at max Syndecan1
    labs(title = paste(domain, "- Syndecan1 vs Biomarkers (Normalized)"),  # Add domain in title
         x = "Syndecan1", 
         y = "Normalized Marker Values (Z-score)",
         color = "Markers") +
    theme_minimal()  # Use a minimal theme
  
  # Display the plot
  print(p)
}



#take syndecan-1 = 13249.52 as a reference ####
median(data_sorted$syndecan_original) #13249.52

#markers one by one 
# Convert Data to Matrix for XGBoost
X <- as.matrix(data_sorted$Syndecan1)  # Feature (independent variable)

# Reference Syndecan-1 value
ref_value <- 13249.52

# XGBoost parameters
params <- list(
  objective = "reg:squarederror",
  max_depth = 2,
  eta = 0.1,
  nthread = 2,
  subsample = 0.8
)

# Store models and predictions
models <- list()
predictions <- data.frame(Syndecan1 = data_sorted$Syndecan1)

# Train XGBoost models and predict values
for (marker in markers) {
  y_train <- as.numeric(data_sorted[[marker]])
  
  # Ensure no missing values
  non_na_idx <- !is.na(y_train)
  X_non_na <- X[non_na_idx, , drop = FALSE]
  y_train <- y_train[non_na_idx]
  
  # Convert to xgboost DMatrix
  dtrain <- xgb.DMatrix(data = X_non_na, label = y_train)
  
  # Train model
  models[[marker]] <- xgb.train(
    params = params,
    data = dtrain,
    nrounds = 100
  )
  
  # Predict for all values of Syndecan-1
  predictions[[marker]] <- predict(models[[marker]], xgb.DMatrix(data = X))
}

# Predict value at reference Syndecan-1 = 13249.52
ref_matrix <- matrix(ref_value, nrow = 1, ncol = 1)
ref_predictions <- sapply(models, function(model) {
  predict(model, xgb.DMatrix(data = ref_matrix))
})

# Generate and display plots one by one
for (marker in markers) {
  # Adjust y-values relative to the reference value
  predictions$Adjusted <- predictions[[marker]] - ref_predictions[marker]
  
  # Compute confidence interval using LOESS regression
  loess_model <- loess(predictions$Adjusted ~ predictions$Syndecan1, span = 0.5)
  loess_pred <- predict(loess_model, se = TRUE)
  
  # Store smooth predictions & standard errors
  predictions$Smooth <- loess_pred$fit
  predictions$Lower <- loess_pred$fit - 1.96 * loess_pred$se.fit
  predictions$Upper <- loess_pred$fit + 1.96 * loess_pred$se.fit
  
  # Generate plot
  plot <- ggplot(predictions, aes(x = Syndecan1, y = Smooth)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "gray", alpha = 0.3) + # Confidence interval
    geom_line(color = "black", size = 1) +  # Smoothed line
    xlab("Syndecan-1 (pg/ml)") +
    ylab(paste0(marker, " Level (Relative to ", ref_value, ")")) + 
    ggtitle(paste0("Association Between Syndecan-1 and ", marker)) + 
    theme_bw() +
    theme(
      aspect.ratio = 0.5/1,
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10)
    ) +
    geom_hline(yintercept= 0, linetype='solid', col = 'black') +  # Center Y at zero
    geom_vline(xintercept = ref_value, linetype = 'dashed', col = 'black')  # Reference line
  
  # Display each plot one by one
  print(plot)
  
  # Pause execution until user presses a key (optional)
  readline(prompt = "Press [Enter] to continue to the next plot...")
}
######


#markers one by domains
# Define biomarker domains
biomarker_domains <- list(
  "Systemic Inflammation and Organ Damage" = c("CD163", "MMP8", "NGAL", "Procalcitonin", "RAGE", "TenascinC", "TREM1"),
  "Cytokines" = c("IL10", "IL18", "IL1ra", "IL6", "IL8"),
  "Coagulation Activation" = c("CoagulationFactorIII", "Ddimer", "Platelets_value_1",  "PT_Max_24h", "Antitrombin", "PROC"),
  "Endothelial Cell Activation and Dysfunction" = c("Ang1", "Ang2", "CX3CL1", "Thrombomodulin", "ESM1")
)
# Convert Data to Matrix for XGBoost
X <- as.matrix(data_sorted$Syndecan1)  # Feature (independent variable)

# Reference Syndecan-1 value
ref_value <- 13249.52

# XGBoost parameters
params <- list(
  objective = "reg:squarederror",
  max_depth = 2,
  eta = 0.1,
  nthread = 2,
  subsample = 0.8
)

# Store models and predictions
models <- list()
predictions <- data.frame(Syndecan1 = data_sorted$Syndecan1)

# Train XGBoost models and predict values
for (domain in names(biomarker_domains)) {
  for (marker in biomarker_domains[[domain]]) {
    y_train <- as.numeric(data_sorted[[marker]])
    
    # Ensure no missing values
    non_na_idx <- !is.na(y_train)
    X_non_na <- X[non_na_idx, , drop = FALSE]
    y_train <- y_train[non_na_idx]
    
    # Convert to xgboost DMatrix
    dtrain <- xgb.DMatrix(data = X_non_na, label = y_train)
    
    # Train model
    models[[marker]] <- xgb.train(
      params = params,
      data = dtrain,
      nrounds = 100
    )
    
    # Predict for all values of Syndecan-1
    predictions[[marker]] <- predict(models[[marker]], xgb.DMatrix(data = X))
  }
}

# Predict value at reference Syndecan-1 = 13249.52
ref_matrix <- matrix(ref_value, nrow = 1, ncol = 1)
ref_predictions <- sapply(models, function(model) {
  predict(model, xgb.DMatrix(data = ref_matrix))
})

# Generate one figure per domain with smooth LOESS curves
for (domain in names(biomarker_domains)) {
  domain_data <- data.frame(Syndecan1 = data_sorted$Syndecan1)
  
  # Adjust markers relative to the reference
  for (marker in biomarker_domains[[domain]]) {
    domain_data[[marker]] <- predictions[[marker]] - ref_predictions[marker]
  }
  
  # Convert to long format for ggplot
  domain_data_long <- tidyr::pivot_longer(domain_data, cols = -Syndecan1, names_to = "Marker", values_to = "Value")
  
  # Apply LOESS smoothing for each marker
  plot <- ggplot(domain_data_long, aes(x = Syndecan1, y = Value, color = Marker)) +
    geom_smooth(se = FALSE, method = "loess", span = 0.3, size = 1) +  # Smooth lines (like first figure)
    xlab("Syndecan-1 (pg/ml)") +
    ylab("Biomarker Level (Relative to 13249.52)") + 
    ggtitle(paste0("Domain: ", domain)) + 
    theme_bw() +
    theme(
      aspect.ratio = 0.5/1,
      legend.position = "right",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10)
    ) +
    geom_hline(yintercept= 0, linetype='solid', col = 'black') +  # Center Y at zero
    geom_vline(xintercept = ref_value, linetype = 'dashed', col = 'black')  # Reference line
  
  # Display the plot
  print(plot)
}
#####

#Normalize markers using Z-score####
# Define biomarker domains
biomarker_domains <- list(
  "Systemic Inflammation and Organ Damage" = c("CD163", "MMP8", "NGAL", "Procalcitonin", "RAGE", "TenascinC", "TREM1"),
  "Cytokines" = c("IL10", "IL18", "IL1ra", "IL6", "IL8"),
  "Coagulation Activation" = c("CoagulationFactorIII", "Ddimer", "Platelets_value_1",  "PT_Max_24h", "Antitrombin", "PROC"),
  "Endothelial Cell Activation and Dysfunction" = c("Ang1", "Ang2", "CX3CL1", "Thrombomodulin", "ESM1")
)

# Convert Data to Matrix for XGBoost
X <- as.matrix(data_sorted$Syndecan1)  # Feature (independent variable)

# Reference Syndecan-1 value
ref_value <- 13249.52

# XGBoost parameters
params <- list(
  objective = "reg:squarederror",
  max_depth = 2,
  eta = 0.1,
  nthread = 2,
  subsample = 0.8
)

# Store models and predictions
models <- list()
predictions <- data.frame(Syndecan1 = data_sorted$Syndecan1)

# Train XGBoost models and predict values
for (domain in names(biomarker_domains)) {
  for (marker in biomarker_domains[[domain]]) {
    y_train <- as.numeric(data_sorted[[marker]])
    
    # Ensure no missing values
    non_na_idx <- !is.na(y_train)
    X_non_na <- X[non_na_idx, , drop = FALSE]
    y_train <- y_train[non_na_idx]
    
    # Convert to xgboost DMatrix
    dtrain <- xgb.DMatrix(data = X_non_na, label = y_train)
    
    # Train model
    models[[marker]] <- xgb.train(params = params, data = dtrain, nrounds = 100)
    
    # Predict for all values of Syndecan-1
    predictions[[marker]] <- predict(models[[marker]], xgb.DMatrix(data = X))
  }
}

# Predict value at reference Syndecan-1 = 13249.52
ref_matrix <- matrix(ref_value, nrow = 1, ncol = 1)
ref_predictions <- sapply(models, function(model) {
  predict(model, xgb.DMatrix(data = ref_matrix))
})

# Generate one figure per domain with Z-score normalization
for (domain in names(biomarker_domains)) {
  domain_data <- data.frame(Syndecan1 = data_sorted$Syndecan1)
  
  # Normalize markers using Z-score
  for (marker in biomarker_domains[[domain]]) {
    marker_values <- predictions[[marker]]
    marker_mean <- mean(marker_values, na.rm = TRUE)
    marker_sd <- sd(marker_values, na.rm = TRUE)
    
    # Z-score normalization
    domain_data[[marker]] <- (marker_values - marker_mean) / marker_sd
  }
  
  # Convert to long format for ggplot
  domain_data_long <- tidyr::pivot_longer(domain_data, cols = -Syndecan1, names_to = "Marker", values_to = "Z_Score")
  
  # Get label positions (at max Syndecan-1)
  label_positions <- domain_data_long %>%
    group_by(Marker) %>%
    filter(Syndecan1 == max(Syndecan1)) %>%
    ungroup()
  
  # Apply LOESS smoothing for each marker and plot
  plot <- ggplot(domain_data_long, aes(x = Syndecan1, y = Z_Score, color = Marker)) +
    geom_smooth(se = FALSE, method = "loess", span = 0.3, size = 1) +  # Smooth lines
    geom_text_repel(data = label_positions, aes(label = Marker), hjust = 0, nudge_x = 500, size = 5, show.legend = FALSE) +  # Label markers at end of curves
    xlab("Syndecan-1 (pg/ml)") +
    ylab("Z-Score Normalized Biomarker Level") +
    ggtitle(paste0("Domain: ", domain)) +
    theme_bw() +
    theme(
      aspect.ratio = 0.5/1,
      legend.position = "none",  # Hide legend since labels are on the plot
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10)
    ) +
    geom_hline(yintercept = 0, linetype = 'solid', col = 'black') +  # Center Y at zero
    geom_vline(xintercept = ref_value, linetype = 'dashed', col = 'black')  # Reference line
  
  # Display the plot
  print(plot)
}



#####
hist(merged_data$Syndecan1, 
     main = "Histogram of Syndecan-1 in BSI",
     xlab = "Syndecan-1 Levels",
     ylab = "Density",
     col = "lightblue", 
     border = "black", 
     probability = TRUE, 
     xlim = c(0, max(merged_data$Syndecan1, na.rm = TRUE) * 1.1),
     ylim = c(0, 4/100000 * 1.2))  # Extend x-axis

lines(density(merged_data$Syndecan1, na.rm = TRUE, bw = "nrd0"), 
      col = "red", 
      lwd = 2)  # Adjust bandwidth if needed
#####