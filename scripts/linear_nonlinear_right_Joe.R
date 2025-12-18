# Clear environment
rm(list = ls())

# Load required libraries
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)
p_load(rms)

# Read the dataset
merged_data <- read.csv("original_data/merged_data.csv")

# List of biomarkers (excluding Syndecan1, which is the predictor)
biomarkers <- c("Syndecan1","Antitrombin", "PROC", "Ang1", "Ang2", "CD163", "CoagulationFactorIII",
             "CX3CL1", "CystatinC", "Ddimer", "ESM1", "IL10", "IL18", "IL1ra", "IL6", 
             "IL8", "MMP8", "NGAL", "Procalcitonin", "RAGE", "TenascinC", "Thrombomodulin",
             "TREM1", "Platelets_value_1", "PT_Max_24h")

markers_data <- merged_data[,c("Syndecan1","Antitrombin", "PROC", "Ang1", "Ang2", "CD163", "CoagulationFactorIII",
                               "CX3CL1", "CystatinC", "Ddimer", "ESM1", "IL10", "IL18", "IL1ra", "IL6", 
                               "IL8", "MMP8", "NGAL", "Procalcitonin", "RAGE", "TenascinC", "Thrombomodulin",
                               "TREM1", "Platelets_value_1", "PT_Max_24h","SOFAtot_new.x","APACHE_IV_Score.x")]

# Log-transform  the biomarkers
markers_data <- markers_data %>%
  mutate(across(all_of(biomarkers), ~ log(. + 1))) # Adding 1 to avoid log(0) issues

# Standardize biomarkers after log transformation (only for complete cases)
markers_data[biomarkers] <- scale(markers_data[biomarkers])

merged_markers_log <- markers_data

# Set empty dataframe to store results
resultsdf_nool <- data.frame(biomarker=biomarkers, NonlinearP=NA)
row.names(resultsdf_nool) <- resultsdf_nool$biomarker

results_df <- data.frame(
  Biomarker = character(), 
  P_Value = numeric(), 
  AIC_Nonlinear = numeric(), 
  AIC_Linear = numeric(), 
  stringsAsFactors = FALSE
)

merged_markers_log <- merged_markers_log %>% 
  filter(!is.na(SOFAtot_new.x))

# **Fix: Set datadist for rms**
dd <- datadist(merged_markers_log)
options(datadist = "dd")

# Calculate p-values and store AIC
biomarkers <- biomarkers[biomarkers!="Syndecan1"]

for (g in biomarkers) {
  
  formula_nl <- as.formula(paste(g, " ~ SOFAtot_new.x + APACHE_IV_Score.x + rcs(Syndecan1, 3)"))
  model_nl <- ols(formula_nl, data = merged_markers_log)
  anres <- anova(model_nl)

  
  roundp <- anres[4,5]   # Check row position for nonlinearity p-value
  resultsdf_nool[g, "NonlinearP"] <- roundp
  
  # Store AIC for nonlinear model
  aic_nl <- AIC(model_nl)
  
  plot(merged_markers_log$Syndecan1, merged_markers_log[[g]], ylab=g, xlab="Syndecan1", 
       main=paste("Test for nonlinearity: P =", roundp), type="p", lwd=2, pch=19, cex=0.5)
  
  if (roundp < 0.05) {  # If nonlinear p < 0.05
    predictions_nl <- Predict(model_nl, Syndecan1=seq(-3.6,3.04,by=0.01))
    points(predictions_nl$Syndecan1, predictions_nl$yhat, col="red", pch=16, add=T)
    svglite::svglite(file=paste0("figures/", g, "_log.svg"), width=5, height=6)
    plot(merged_markers_log$Syndecan1, merged_markers_log[[g]], ylab=g, xlab="Syndecan1", 
         main=paste("Test for nonlinearity: P =", roundp))
    points(predictions_nl$Syndecan1, predictions_nl$yhat, col="red", pch=16, add=T)
    dev.off()
    
    aic_l <- NA  # No linear model AIC for nonlinear cases
    
  } else {  # If nonlinear p > 0.05
    formula_l <- as.formula(paste(g, " ~ SOFAtot_new.x + APACHE_IV_Score.x + Syndecan1"))
    model_l <- ols(formula_l, data = merged_markers_log)
    predictions_l <- Predict(model_l, Syndecan1=seq(-3.6,3.04,by=0.01))
    points(predictions_l$Syndecan1, predictions_l$yhat, col="blue", pch=16, add=T)
    svglite::svglite(file=paste0("figures/", g, "_log.svg"), width=5, height=6)
    plot(merged_markers_log$Syndecan1, merged_markers_log[[g]], ylab=g, xlab="Syndecan1", 
         main=paste("Test for nonlinearity: P =", roundp))
    points(predictions_l$Syndecan1, predictions_l$yhat, col="blue", pch=16, add=T)
    dev.off()
    
    # Store AIC for linear model
    aic_l <- AIC(model_l)
  }
  
  # Store results
  results_df <- rbind(results_df, data.frame(
    Biomarker = g, 
    P_Value = roundp, 
    AIC_Nonlinear = aic_nl, 
    AIC_Linear = aic_l
  ))
}

# Adjust p-values
results_df$P_ValueAdj <- p.adjust(results_df$P_Value, method = "BH")

# Reset options after analysis
options(datadist = NULL)
