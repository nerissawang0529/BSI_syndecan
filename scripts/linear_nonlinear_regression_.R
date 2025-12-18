rm(list = ls())
install.packages("rms")
library(rms)
install.packages("FSA") 
library(FSA)
install.packages("VGAM")

## install.packages("pacman")
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)
install.packages("dplyr")
install.packages("FSA")
library(dplyr)
library(mice)
library(tidyr)
library(FSA)


#sensitivy analysis 1. ALL patients (N=190) 2. only those with known single pathogen (N=150)

#ready data
data <- read.csv("original_data/clinical_marker_scource_pathogen.csv")

#check subsofa score
##new sofa from hessel
load("original_data/MARS_SOFAS.RData")

library(dplyr)

#get new sofa from hessel within the IDs I need
sofas16_filtered <- sofas16 %>%
  filter(ICU_ID_from_datasource %in% data$ICU_ID_from_datasource)


# Define the columns to rename
cols_to_rename <- c("SOFAtot", "SOFA_Circulation","SOFA_Liver", "SOFA_Coagulation",
                    "SOFA_Live", "SOFA_Renal", "SOFA_Respiration")
# Find indices of these columns in the dataframe
indices <- which(names(sofas16_filtered) %in% cols_to_rename)
# Append "_new" to their names
names(sofas16_filtered)[indices] <- paste0(names(sofas16_filtered)[indices], "_new")

library(dplyr)

sofas16_filtered <- sofas16_filtered %>% 
  filter(diff_admin_daily %in% c(0, 1))

sofas16_filtered$Index_ICU_admittance_datetime <- NULL
sofas16_filtered$ICU_Admission_TK <- NULL
sofas16_filtered$Assessment_datetime <- NULL
sofas16_filtered$diff_admin_daily <- NULL

data <- merge(data, sofas16_filtered, by = "ICU_ID_from_datasource", all.x = TRUE)

sofa_data <- data[,c("ICU_ID_from_datasource","SOFAtot","SOFAtot_new","SOFA_Circulation","SOFA_Circulation_new","SOFA_Coagulation","SOFA_Coagulation_new",
                     "SOFA_Liver","SOFA_Liver_new","SOFA_Renal","SOFA_Renal_new","SOFA_Respiration","SOFA_Respiration_new")]


#data_DAG <- data[,c("ICU_ID_from_datasource", "FinalGroup_pathegon", "APACHE_IV_Acute_Physiology_Score.x", "SOFAtot_new", 
                    "mort365", "D.dimer..43.", "ANG2_1", "gender", "age_yrs", "Cerebrovascular_disease", 
                    "ckd", "Chronic_heart_failure",  "COPD", "Immune_deficiency","Metastatic_malignancy", "Past_myocardial_infarction", 
                    "Bilirubin","Creatinine_value_1","Platelets_value_1", "Blood_Urea_Nitrogen_value_1","Syndecan.1.CD138..29.")]

#destination_folder <- "original_data" 
#export_file_name <- "data_DAG.csv" 
#write.csv(data_DAG, file.path(destination_folder, export_file_name), row.names = FALSE)


#data <- data[!is.na(data$SOFAtot), ]
#data$syndecan_1 <- log(data$syndecan_original)
data$SOFA_Renal
data_subset <- data[, c("SOFA_Renal", "syndecan_1")]
# Remove rows with missing values in these columns
data_clean <- data_subset[complete.cases(data_subset), ]
data_clean <- data[complete.cases(data$SOFA_Circulation, data$syndecan_1), ]
data_clean <- data[complete.cases(data$SOFA_Circulation, data$syndecan_1), ]

library(rms)

# Set up the datadist object on the clean data
dd <- datadist(data_clean)
options(datadist = "dd")

#proportional odds model
# Convert SOFA_Circulation into an ordered factor with levels 0, 1, 2, 3, 4
data$SOFA_Renal <- factor(data$SOFA_Renal, 
                                levels = c(0, 1, 2, 3, 4), 
                                ordered = TRUE)

fit_rms <- orm(SOFA_Renal ~ syndecan_1, data = data_clean)
summary(fit_rms)


# Load the VGAM package
library(VGAM)

# Fit the cumulative logit (proportional odds) model
fit_vgam <- vglm(SOFA_Circulation ~ syndecan_1,
                 family = cumulative(parallel = TRUE),
                 data = data_clean)

# Display the model summary
summary(fit_vgam)

















#pathogen
data$FinalGroup_pathegon2 <- ifelse(data$FinalGroup_pathegon=="Mixed_pathegon", NA, data$FinalGroup_pathegon)
data$FinalGroup_pathegon2 <- factor(data$FinalGroup_pathegon2, levels = c("E.coli", "Enterobacter", "Enterococcus", "Klebsiella", "Other_pathogens", "Pseudomonas", "S.aureus","Streptococcus"))
#by this code, the E.coli is the reference group, normally take the largest number groups as the reference.

data <- data %>% 
  filter(!is.na(FinalGroup_pathegon2))


#lmer model 
install.packages("lme4")
library(lme4)
model2 <- lmer(syndecan_1 ~ SOFAtot + (1 | FinalGroup_pathegon2), data = data)

data$syndecan_1 <- as.numeric(as.character(data$syndecan_1))
data$SOFAtot <- as.numeric(as.character(data$SOFAtot))

model3<-lm(syndecan_1 ~ SOFAtot, data = data)
anova(model2, model3)

#make the source values
data$FinalSource <- "unknown"
data$FinalSource[data$CNS_source==1] <- "CNS"
data$FinalSource[data$abdominal==1] <- "abdominal"
data$FinalSource[data$respiratory==1] <- "respiratory"
data$FinalSource[data$urinary==1] <- "urinary"
data$FinalSource[data$cardiovascular==1] <- "cardiovascular"
data$FinalSource[data$skin==1] <- "skin"
data$FinalSource[data$other_source==1] <- "other"

#table(data$FinalSource)


fit <- lm(data=data, syndecan_1 ~ FinalGroup_pathegon2)
summary(fit)
anova(fit)
plot(fit)

fit1 <- lm(data=data, syndecan_1 ~ FinalGroup_pathegon2 + SOFAtot + APACHE_IV_Acute_Physiology_Score.x)
summary(fit1)
anova(fit1)


















install.packages("emmeans")
library(emmeans)
pairwise_comp <- emmeans(fit1, "FinalGroup_pathegon2")
pairwise_comp <- pairs(pairwise_comp) 
summary(pairwise_comp, adjust="tukey")

#sofa
fit_sofa <- lm(syndecan_1 ~ SOFAtot, data = data)
summary(fit_sofa)

data$SOFA_CNS <- NULL
fit_subsofa <- lm(data=data, log(data$Syndecan.1.CD138..29.) ~ SOFA_Circulation + SOFA_Coagulation + SOFA_Liver + SOFA_Renal + SOFA_Respiration)
summary(fit_subsofa)
anova(fit_subsofa)

#Respiratory System (PaO₂/FiO₂ ratio: PaO2_at_highest_A_aDO2_24h or PaO2_Min_24h/ FiO2_at_highest_A_aDO2_24h), not significant
#Coagulation (Platelet count: Platelets_value_1, Platelet_count_Max)
#Liver (Bilirubin level:Bilirubin half NA)
#Cardiovascular (Hypotension: Hypertension and vasopressor use: Use_of_vasoactive_medic_24h)
#Central Nervous System (Glasgow Coma Scale), which variable we used?
#Renal (Creatinine level:Creatinin_Min_24h, Creatinine_value_1; or urine output:Urine_output_Sum_24h), which variable we used?

install.packages("ggplot2")  # Install ggplot2 if not already installed
cor(log(data$Bilirubin),data$syndecan_1, use = "complete.obs")
library(ggplot2)


# Create a datadist object
data <- data[ , !(names(data) %in% "BNP_First_24h")]
# Remove constant variables
# Identify constant variables
constant_vars <- sapply(data, function(x) length(unique(x)) == 1)
# Display constant variables
names(data)[constant_vars]
# Remove constant variables
data <- data[ , !constant_vars]
d <- datadist(data)
options(datadist = "d")

# Set the reference value for Syndecan-1 concentration to log(2143.6)
log_ref <- log(2143.6)
d$limits["Adjust to", "syndecan_1"] <- log_ref

# Fit the regression model with restricted cubic splines for syndecan_1
unadj_syndecan <- ols(Bilirubin ~ rcs(syndecan_1, 3), data = data)
unadj_syndecan <- update(unadj_syndecan)

# Generate predictions with the specified reference
x1 <- Predict(unadj_syndecan, syndecan_1, ref.zero = TRUE)


# Plot the predictions
syndecan_plot <- ggplot(x1) +
  xlab("Log-transformed Syndecan-1 concentration (pg/ml)") +
  ylab("Bilirubin") +
  ggtitle("Relationship between Log-transformed Syndecan-1 and Bilirubin") +
  theme_bw() +
  theme(
    aspect.ratio = 0.5 / 1,
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  ) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_vline(xintercept = log_ref, linetype = "dashed", color = "black") +
  coord_cartesian()

# Display the plot
print(syndecan_plot)













# Create the plot with the correct column names
syndecan_plot <- ggplot(data, aes(x = syndecan_1, y = Bilirubin)) +
  geom_point(alpha = 0.7, color = 'blue', size = 2) +
  geom_smooth(method = 'loess', color = 'red', se = TRUE) +
  xlab("Syndecan-1 concentration (pg/ml)") +
  ylab("Bilirubin") +
  ggtitle("Relationship between Syndecan-1 Levels and Platelets") +
  theme_bw() +
  theme(
    aspect.ratio = 0.5/1,
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )

# Display the plot
print(syndecan_plot)





# Load ggplot2
library(ggplot2)

# Adjusted scatterplot with x-axis limit
# SOFA_Circulation
data <- data[!is.na(data$SOFA_Circulation), ]
ggplot(data, aes(x = as.factor(SOFA_Circulation), y = syndecan_original)) +
  geom_boxplot(fill = "lightblue", color = "black", outlier.color = "red", outlier.shape = 16) +
  stat_summary(fun = mean, geom = "point", shape = 4, size = 3, color = "darkred") + # Show mean as points
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "blue", linetype = "dashed") + # Add trendline
  labs(
    title = "Boxplot of Syndecan-1 Levels by Sub-SOFA Score with Trendline",
    x = "SOFA_Circulation",
    y = "Syndecan-1 Levels"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 13)
  )

kruskal_test <- kruskal.test(Syndecan.1.CD138..29. ~ factor(SOFA_Circulation), data = data)
dunn_test <- dunnTest(Syndecan.1.CD138..29. ~ factor(SOFA_Circulation), data = data, method = "none")


#SOFA_Coagulation
data <- data[!is.na(data$SOFA_Coagulation), ]
ggplot(data, aes(x = as.factor(SOFA_Circulation), y = Syndecan.1.CD138..29.)) +
  geom_boxplot(fill = "lightblue", color = "black", outlier.color = "red", outlier.shape = 16) +
  stat_summary(fun = mean, geom = "point", shape = 4, size = 3, color = "darkred") + # Show mean as points
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "blue", linetype = "dashed") + # Add trendline
  labs(
    title = "Boxplot of Syndecan-1 Levels by Sub-SOFA Score with Trendline",
    x = "SOFA_Coagulation",
    y = "Syndecan-1 Levels"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 13)
  )

kruskal_test <- kruskal.test(Syndecan.1.CD138..29. ~ factor(SOFA_Coagulation), data = data)
dunn_test <- dunnTest(Syndecan.1.CD138..29. ~ factor(SOFA_Coagulation), data = data, method = "none")


#SOFA_Liver
data <- data[!is.na(data$SOFA_Liver), ]
ggplot(data, aes(x = as.factor(SOFA_Liver), y = Syndecan.1.CD138..29.)) +
  geom_boxplot(fill = "lightblue", color = "black", outlier.color = "red", outlier.shape = 16) +
  stat_summary(fun = mean, geom = "point", shape = 4, size = 3, color = "darkred") + # Show mean as points
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "blue", linetype = "dashed") + # Add trendline
  labs(
    title = "Boxplot of Syndecan-1 Levels by Sub-SOFA Score with Trendline",
    x = "SOFA_Liver",
    y = "Syndecan-1 Levels"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 13)
  )

kruskal_test <- kruskal.test(Syndecan.1.CD138..29. ~ factor(SOFA_Liver), data = data)
dunn_test <- dunnTest(Syndecan.1.CD138..29. ~ factor(SOFA_Liver), data = data, method = "none")

#SOFA_Renal

data <- data[!is.na(data$SOFA_Renal), ]
ggplot(data, aes(x = as.factor(SOFA_Renal), y = Syndecan.1.CD138..29.)) +
  geom_boxplot(fill = "lightblue", color = "black", outlier.color = "red", outlier.shape = 16) +
  stat_summary(fun = mean, geom = "point", shape = 4, size = 3, color = "darkred") + # Show mean as points
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "blue", linetype = "dashed") + # Add trendline
  labs(
    title = "Boxplot of Syndecan-1 Levels by Sub-SOFA Score with Trendline",
    x = "SOFA_Renal",
    y = "Syndecan-1 Levels"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 13)
  )

kruskal_test <- kruskal.test(Syndecan.1.CD138..29. ~ factor(SOFA_Renal), data = data)
dunn_test <- dunnTest(Syndecan.1.CD138..29. ~ factor(SOFA_Renal), data = data, method = "none")

#SOFA_Respiration
data <- data[!is.na(data$SOFA_Respiration), ]
ggplot(data, aes(x = as.factor(SOFA_Respiration), y = Syndecan.1.CD138..29.)) +
  geom_boxplot(fill = "lightblue", color = "black", outlier.color = "red", outlier.shape = 16) +
  stat_summary(fun = mean, geom = "point", shape = 4, size = 3, color = "darkred") + # Show mean as points
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "blue", linetype = "dashed") + # Add trendline
  labs(
    title = "Boxplot of Syndecan-1 Levels by Sub-SOFA Score with Trendline",
    x = "SOFA_Respiration",
    y = "Syndecan-1 Levels"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 13)
  )

kruskal_test <- kruskal.test(Syndecan.1.CD138..29. ~ factor(SOFA_Respiration), data = data)
dunn_test <- dunnTest(Syndecan.1.CD138..29. ~ factor(SOFA_Respiration), data = data, method = "none")







##
library("variancePartition")


SynExpr <- t(log(data$Syndecan.1.CD138..29.))
colnames(SynExpr) <- data$ICU_ID_from_datasource
row.names(SynExpr) <- "Synd"
exp2 <- rbind(SynExpr, SynExpr)
row.names(exp2) <- c("Synd","S2")
bsinfo <- data[,c("age_yrs","FinalSource", "FinalGroup_pathegon","SOFAtot","APACHE_IV_Acute_Physiology_Score.x","APACHE_IV_Score.x","Enterobacter","Enterococcus","Klebsiella","Pseudomonas","S.aureus","Streptococcus")]
bsinfo <- data[,c("ICU_ID_from_datasource", "Klebsiella", "S.aureus", "Enterococcus", "E.coli", "Streptococcus", "Pseudomonas", "Enterobacter", "other_culture",
                  "FinalGroup_pathegon", "Shock_on_admis", "APACHE_IV_Score.x", "APACHE_IV_Acute_Physiology_Score.x", "SOFAtot",
                  "mort30","mort365",  "Ang.1..64.", "Ang.2..26.", "Coagulation.Factor.III.Tissue.Factor..62.", "CX3CL1.Fractalkine..46.",
                  "D.dimer..43.", "ESM.1..22.", "IL.10..22.", "IL.18..78.", "IL.1ra.IL.1F3..30.", "IL.6..13.", "IL.8..18.", "MMP.8..27.", "NGAL..53.",
                  "Procalcitonin..39.", "Syndecan.1.CD138..29.", "Thrombomodulin.BDCA.3..47.", "TREM.1..65.", 
                  "ANG2_1", "gender", "age_yrs", "Stroke", "Thrombolysis", "Cerebrovascular_disease", "Chronic_cardiovascular_insuf",
                  "ckd", "Chronic_heart_failure",  "COPD", "Immune_deficiency","Metastatic_malignancy", "Past_myocardial_infarction", 
                  "Bilirubin",  "Creatinine_value_1","Platelets_value_1", "Blood_Urea_Nitrogen_value_1", "WBC_2_1", "White_cell_count_Min_24h"
                  )]

row.names(bsinfo) <- data$ICU_ID_from_datasource
bsinfo <- as.data.frame(bsinfo)
form <- ~ (1| FinalGroup_pathegon) + SOFAtot + age_yrs + APACHE_IV_Acute_Physiology_Score.x + APACHE_IV_Score.x + Shock_on_admis + mort30 +
          mort365 + Ang.1..64. + D.dimer..43. + gender + age_yrs + Cerebrovascular_disease + ckd + Chronic_heart_failure + COPD + Immune_deficiency + Bilirubin +
          Creatinine_value_1 + Platelets_value_1 + Blood_Urea_Nitrogen_value_1 + WBC_2_1
colSums(is.na(bsinfo))
bsinfo <- na.omit(bsinfo)

varPart <- fitExtractVarPartModel(exp2, form, bsinfo)
vp <- sortCols(varPart)
plotPercentBars(vp)


### remove missing

bsinfo2 <- bsinfo[complete.cases(bsinfo),]
exp3 <- exp2[,row.names(bsinfo2)]
form <- ~ (1| FinalGroup_pathegon) + (1| FinalSource) + SOFAtot + age_yrs + APACHE_IV_Acute_Physiology_Score.x + APACHE_IV_Score.x
form <- ~ (1| FinalGroup_pathegon) + SOFAtot + age_yrs + APACHE_IV_Acute_Physiology_Score.x + APACHE_IV_Score.x + Shock_on_admis + mort30 +
  mort365 + Ang.1..64. + D.dimer..43. + gender + age_yrs + Cerebrovascular_disease + ckd + Chronic_heart_failure + COPD + Immune_deficiency + Bilirubin +
  Creatinine_value_1 + Platelets_value_1 + Blood_Urea_Nitrogen_value_1 + WBC_2_1

varPart <- fitExtractVarPartModel(exp3, form, bsinfo2)
vp <- sortCols(varPart)
plotPercentBars(vp)
vp



#non-linear model
# Create a datadist object
data <- data[ , !(names(data) %in% "BNP_First_24h")]
# Remove constant variables
# Identify constant variables
constant_vars <- sapply(data, function(x) length(unique(x)) == 1)
# Display constant variables
names(data)[constant_vars]
# Remove constant variables
data <- data[ , !constant_vars]
dd <- datadist(data)
options(datadist = "dd")

# Dynamically build the formula as a character string
formula_string <- paste(
  "syndecan_1 ~ FinalGroup_pathegon2",
  "rcs(SOFAtot, 3)",
  "rcs(APACHE_IV_Acute_Physiology_Score.x, 3)",
  sep = " + "
)

# Convert the character string to a formula object
fit_spline <- lm(
  formula(formula_string),
  data = data
)


summary(fit_spline)
anova(fit_spline)

anova_result <- anova(fit1, fit_spline)
print(anova_result)




###xgboost
rm(list = ls())

# Install and Load Required Libraries
if (!require("xgboost")) install.packages("xgboost")
if (!require("caret")) install.packages("caret")
if (!require("Matrix")) install.packages("Matrix")
if (!require("SHAPforxgboost")) install.packages("SHAPforxgboost")
library(xgboost)
library(caret)
library(Matrix)
library(SHAPforxgboost)

# Load Data
data <- read.csv("original_data/clinical_marker_scource_pathogen.csv")
data <- data[!(data$FinalGroup_pathegon %in% c("Mixed_pathegon", "Other_pathogens")), ]

# Cut Syndecan-1 into groups (if applicable for classification)
data$Syndecan_group <- cut(data$Syndecan.1.CD138..29.,
                           breaks = quantile(data$`Syndecan.1.CD138..29.`, probs = seq(0, 1, length = 4), na.rm = TRUE),
                           include.lowest = TRUE,
                           labels = c(1, 2, 3))

# Train-Test Split
set.seed(3456)
trainIndex <- createDataPartition(data$FinalGroup_pathegon, p = .8, list = FALSE)
train <- data[trainIndex, ]
test <- data[-trainIndex, ]

# One-hot encode 'FinalGroup_pathegon' and create sparse matrices
train_sparse <- sparse.model.matrix(~ FinalGroup_pathegon + APACHE_IV_Score.x - 1, data = train)
test_sparse <- sparse.model.matrix(~ FinalGroup_pathegon + APACHE_IV_Score.x - 1, data = test)

# Prepare the target variable for regression
train_label <- train$`Syndecan.1.CD138..29.`  # Syndecan-1 as the target
test_label <- test$`Syndecan.1.CD138..29.`

# Train XGBoost Model
xgb_model <- xgboost(data = train_sparse, 
                     label = train_label, 
                     max.depth = 3, 
                     eta = 0.1, 
                     nthread = 2, 
                     nrounds = 100, 
                     objective = "reg:squarederror", 
                     verbose = 0)


# Convert train_sparse to a dense matrix for SHAP
train_dense <- as.matrix(train_sparse)
# SHAP Values Calculation
shap_values <- shap.prep(xgb_model = xgb_model, X_train = train_dense)
# SHAP Summary Plot
shap.plot.summary(shap_values)



# Predictions on Test Data
predicted_values <- predict(xgb_model, test_sparse)

# Scatter Plot of Actual vs Predicted Values
plot(test_label, predicted_values,
     xlab = "Actual Syndecan-1 Levels", 
     ylab = "Predicted Syndecan-1 Levels", 
     main = "Scatter Plot: Actual vs Predicted",
     pch = 19, col = "blue")
abline(0, 1, col = "red", lwd = 2)  # Perfect Prediction Line


# Create Bins for Calibration (Fixing Duplicate Breaks)
bins <- seq(min(predicted_values, na.rm = TRUE), max(predicted_values, na.rm = TRUE), length.out = 10)

# Bin the predicted values
binned_data <- data.frame(
  Actual = test_label,
  Predicted = predicted_values,
  Bin = cut(predicted_values, bins, include.lowest = TRUE)
)

# Calculate mean actual and predicted values per bin
calibration_data <- aggregate(cbind(Actual, Predicted) ~ Bin, data = binned_data, mean)

# Calibration Plot
plot(as.numeric(calibration_data$Predicted), calibration_data$Actual, 
     xlab = "Mean Predicted Syndecan-1", 
     ylab = "Mean Actual Syndecan-1", 
     main = "Calibration Plot",
     pch = 19, col = "blue")
abline(0, 1, col = "red", lwd = 2)  # Perfect Calibration Line

