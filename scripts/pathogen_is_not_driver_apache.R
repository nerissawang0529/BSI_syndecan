# Load necessary libraries
library(ggplot2)
library(dplyr)

# Create a clean environment
rm(list = ls())

# Load data
data <- read.csv("original_data/clinical_marker_scource_pathogen.csv")

# Define infection source categories
data$FinalSource <- "unknown"
data$FinalSource[data$CNS_source == 1] <- "CNS"
data$FinalSource[data$abdominal == 1] <- "Abdominal"
data$FinalSource[data$respiratory == 1] <- "Respiratory"
data$FinalSource[data$urinary == 1] <- "Urinary"
data$FinalSource[data$cardiovascular == 1] <- "Cardiovascular"
data$FinalSource[data$skin == 1] <- "Skin"
data$FinalSource[data$other_source == 1] <- "Other"

# Fit linear model with APACHE IV Score
fit0 <- lm(data = data, log(data$Syndecan1) ~ 
             FinalGroup_pathegon + 
             FinalSource +
             APACHE_IV_Score.x)  # Using APACHE IV Score

# Extract variance explained
var_explained <- anova(fit0)
var_explained$PctExp <- (var_explained$"Sum Sq" / sum(var_explained$"Sum Sq")) * 100  # Convert to percentage

# Create dataframe for plotting
var_df <- data.frame(Variable = rownames(var_explained), Variance = var_explained$PctExp)

# Rename variables for clarity
var_df$Variable <- recode(var_df$Variable,
                          "APACHE_IV_Score.x" = "APACHE IV Score",  # Rename APACHE IV Score
                          "FinalGroup_pathegon" = "Pathogen Type",
                          "FinalSource" = "Infection Site")

# Define category groups as factors with custom order
var_df$Category <- factor(ifelse(var_df$Variable == "Residuals", "Residuals",
                                 ifelse(var_df$Variable == "APACHE IV Score", "Disease Severity (APACHE IV Score)",
                                        ifelse(var_df$Variable == "Pathogen Type", "Pathogen Type", "Infection Site"))),
                          levels = c("Residuals", "Disease Severity (APACHE IV Score)", "Pathogen Type", "Infection Site"))  # Custom order

# Sort bars in order of category levels
var_df <- var_df[order(match(var_df$Category, levels(var_df$Category))), ]

# Define **colorblind-friendly** color palette (Okabe-Ito palette)
category_colors <- c("Residuals" = "#999999",  # Gray
                     "Disease Severity (APACHE IV Score)" = "#D55E00",  # Distinct Orange for APACHE IV
                     "Pathogen Type" = "#56B4E9",  # Sky Blue
                     "Infection Site" = "#009E73")  # Green

# Generate high-quality publication-style bar chart with correct order
ggplot(var_df, aes(x = reorder(Variable, Variance), y = Variance, fill = Category)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +  # Thin black border for clarity
  coord_flip() +  
  labs(x = "", y = "Variance Explained (%)") +  # Remove x-axis label for cleaner look
  scale_fill_manual(name = "Category", values = category_colors, breaks = levels(var_df$Category)) +  # Ensure legend order
  guides(fill = "none") +  # Remove legend
  theme_minimal(base_size = 14) +  # Use minimal theme with readable font size
  theme(
    panel.grid.major.x = element_blank(),  # Remove major gridlines
    panel.grid.minor.x = element_blank(),  # Remove minor gridlines
    panel.grid.major.y = element_line(color = "grey80", linetype = "dashed"),  # Light dashed y-axis grid
    axis.text = element_text(face = "bold"),  # Make axis labels bold
    legend.position = "none",  # Remove legend
    plot.title = element_text(hjust = 0.5, face = "bold")  # Center and bold the title
  ) +
  #  *** Add significance stars for APACHE IV Score
  geom_text(data = var_df[var_df$Variable == "APACHE IV Score", ], 
            aes(label = "***", y = Variance + 2), size = 6, fontface = "bold")

