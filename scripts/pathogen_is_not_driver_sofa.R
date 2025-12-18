# 🧹 Clean environment
rm(list = ls())

# 📦 Load required libraries
library(ggplot2)
library(dplyr)

# 📂 Load your data
data <- read.csv("original_data/final_data.csv")

# ✅ Define infection source
data$FinalSource <- "unknown"
data$FinalSource[data$CNS_source == 1] <- "CNS"
data$FinalSource[data$abdominal == 1] <- "Abdominal"
data$FinalSource[data$respiratory == 1] <- "Respiratory"
data$FinalSource[data$urinary == 1] <- "Urinary"
data$FinalSource[data$cardiovascular == 1] <- "Cardiovascular"
data$FinalSource[data$skin == 1] <- "Skin"
data$FinalSource[data$other_source == 1] <- "Other"

# 🧮 Fit linear model with all specified covariates
fit <- lm(data = data, log(Syndecan1) ~ 
            FinalGroup_pathegon + 
            FinalSource +
            SOFAtot_new +
            age_yrs + gender + 
            ckd + Cerebrovascular_disease + 
            Past_myocardial_infarction + 
            diabetes + Immune_deficiency)

# 📊 Extract variance explained
var_explained <- anova(fit)
var_explained$PctExp <- (var_explained$"Sum Sq" / sum(var_explained$"Sum Sq")) * 100

# 🧱 Create dataframe for plotting
var_df <- data.frame(Variable = rownames(var_explained), Variance = var_explained$PctExp)

# ✍️ Rename variables for plotting
var_df$Variable <- recode(var_df$Variable,
                          "SOFAtot_new" = "SOFA Score",
                          "FinalGroup_pathegon" = "Pathogen Type",
                          "FinalSource" = "Infection Site",
                          "age_yrs" = "Age",
                          "gender" = "Gender",
                          "ckd" = "CKD",
                          "Cerebrovascular_disease" = "Cerebrovascular Disease",
                          "Past_myocardial_infarction" = "Past MI",
                          "diabetes" = "Diabetes",
                          "Immune_deficiency" = "Immune Deficiency")

# 🟡 Categorize variables for plotting
var_df$Category <- ifelse(var_df$Variable == "Residuals", "Residuals",
                          ifelse(var_df$Variable == "SOFA Score", "Disease Severity",
                                 ifelse(var_df$Variable %in% c("Pathogen Type", "Infection Site"), "Pathogen Info", 
                                        "Comorbidities")))

var_df$Category <- factor(var_df$Category,
                          levels = c("Residuals", "Disease Severity", "Comorbidities", "Pathogen Info"))


complication_vars <- c("CKD", "Cerebrovascular Disease", "Past MI", "Diabetes", "Immune Deficiency")
complication_sum <- var_df %>%
  filter(Variable %in% complication_vars) %>%
  summarise(Variance = sum(Variance)) %>%
  mutate(Variable = "Complication")


var_df_filtered <- var_df %>%
  filter(Variable %in% c("SOFA Score", "Pathogen Type", "Infection Site", "Age", "Gender", "Residuals")) %>%
  bind_rows(complication_sum)

library(ggplot2)
library(dplyr)

# 
var_df_filtered$Variable <- factor(var_df_filtered$Variable, levels = rev(c(
  "SOFA Score", "Pathogen Type", "Infection Site", "Complication", "Age", "Gender", "Residuals"
)))

# 
bar_colors <- c(
  "SOFA Score" = "#B22222",
  "Pathogen Type" = "#E67E22",
  "Infection Site" = "#F1C40F",
  "Complication" = "#7F8C8D",
  "Age" = "#3498DB",
  "Gender" = "#2980B9",
  "Residuals" = "#D3D3D3"
)

var_df_filtered$Color <- bar_colors[as.character(var_df_filtered$Variable)]


ggplot(var_df_filtered, aes(x = Variance, y = Variable, fill = Color)) +
  geom_col(width = 0.7) +  # ⬅️ 没有 color = "black"，去除黑边框
  scale_fill_identity() +
  scale_y_discrete(expand = c(0.01, 0.01)) +
  labs(x = "Variance Explained (%)", y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 13, face = "bold", margin = margin(r = 4)),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14, face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey90"),
    plot.margin = margin(t = 5, r = 10, b = 5, l = 10),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "none"
  ) +
  geom_text(data = var_df_filtered %>% filter(Variable == "SOFA Score"),
            aes(label = "***", x = Variance + 1.5),
            size = 5, fontface = "bold", color = "black")
