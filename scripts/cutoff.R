rm(list = ls())

install.packages("mgcv")  # If not installed
install.packages("ggplot2")
install.packages("reshape2")
install.packages("gratia", dependencies=TRUE)


library(mgcv)
library(ggplot2)
library(reshape2)
library(gratia)


# Read data
merged_data <- read.csv("original_data/merged_data.csv")
# List of markers (y-axis variables)
marker_data <- merged_data[,c("Syndecan1","Antitrombin", "PROC", "Ang1", "Ang2", "CD163", "CoagulationFactorIII",
             "CX3CL1", "CystatinC", "Ddimer", "ESM1", "IL10", "IL18", "IL1ra", "IL6", 
             "IL8", "MMP8", "NGAL", "Procalcitonin", "RAGE", "TenascinC", "Thrombomodulin",
             "TREM1", "Platelets_value_1", "PT_Max_24h")]
# Reshape data: Convert all biomarkers into one column
long_data <- melt(marker_data, id.vars = "Syndecan1", variable.name = "Biomarker", value.name = "Value")
# Fit a GAM model with smooth spline for Syndecan-1 across all biomarkers
gam_model <- gam(Value ~ s(Syndecan1, k = 5) + Biomarker, data = long_data)

# Print model summary
summary(gam_model)

ggplot(long_data, aes(x = Syndecan1, y = Value, color = Biomarker)) +
  geom_point(alpha = 0.5) +  # Add points with slight transparency
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5), se = FALSE) +  # Add smooth GAM line
  labs(title = "Syndecan-1 vs Biomarkers (Smooth Fit)", x = "Syndecan-1", y = "Biomarker Value")



install.packages("gratia")  # If not installed
library(gratia)

# Compute the first derivative of the GAM model
deriv <- derivatives(gam_model, term = "s(Syndecan1)")

# Find where the derivative is significantly different from zero (potential breakpoints)
breakpoints <- deriv[abs(deriv$derivative) > 0.1, ]
print(breakpoints)

# Plot derivative to see where the biggest changes occur
draw(deriv)
