rm(list = ls())

install.packages("drat", repos="https://cran.rstudio.com")
drat:::addRepo("dmlc")
install.packages("xgboost", repos="http://dmlc.ml/drat/", type = "source")
install.packages("caret")
require(xgboost)
library(xgboost)
library(Matrix)
library(dplyr)

#ready data
data <- read.csv("original_data/clinical_marker_scource_pathogen.csv")
data <- data[!(data$FinalGroup_pathegon %in% c("Mixed_pathegon", "Other_pathogens")),]
data <- data[!is.na(data$SOFAtot), ]


# Cutting the sorted data into three equal parts (1, 2, 3)
data$Syndecan_group <- cut(data$Syndecan.1.CD138..29.,
                           breaks = quantile(data$`Syndecan.1.CD138..29.`, probs = seq(0, 1, length = 4), na.rm = TRUE),
                           include.lowest = TRUE,
                           labels = c(1, 2, 3))


#split the data into train and test
library(caret)
set.seed(3456)
trainIndex <- createDataPartition(data$FinalGroup_pathegon, p = .8, 
                                  list = FALSE, 
                                  times = 1)
train <- data[trainIndex, ]   # Rows selected by trainIndex (Training set)
train <- train[!is.na(train$SOFAtot), ]
test  <- data[-trainIndex, ]  # Rows NOT in trainIndex (Testing set)

# Encode 'FinalGroup_pathegon' as a numeric feature using one-hot encoding
train_one_hot <- model.matrix(~ FinalGroup_pathegon + SOFAtot + APACHE_IV_Acute_Physiology_Score.x - 1, data = train)

# Adjust the labels to start from 0
label <- as.numeric(train$Syndecan_group) - 1  # Convert to 0, 1, 2


# Convert the one-hot encoded matrix to a sparse matrix for XGBoost
data_sparse <- Matrix::sparse.model.matrix(~.-1, data = as.data.frame(train_one_hot))

# Train the XGBoost model
bstSparse <- xgboost(data = data_sparse, 
                     label = label, 
                     max.depth = 2, 
                     eta = 1, 
                     nthread = 2, 
                     nrounds = 2, 
                     objective = "multi:softmax", 
                     num_class = length(unique(label)))


#Test

#Prepare the test data by encoding 'FinalGroup_pathegon' using one-hot encoding
test_one_hot <- model.matrix(~ FinalGroup_pathegon + SOFAtot + APACHE_IV_Acute_Physiology_Score.x - 1, data = test)
# Convert the one-hot encoded matrix to a sparse matrix for prediction
test_sparse <- Matrix::sparse.model.matrix(~.-1, data = as.data.frame(test_one_hot))
# Make predictions
pred <- predict(bstSparse, test_sparse)
# Size of the prediction vector
print(length(pred))
# Print the predicted classes
print(pred)

prediction <- as.numeric(pred > 0.5)
print(head(prediction))

err <- mean(as.numeric(pred > 0.5) != test$label)
print(paste("test-error=", err))


# Get feature names from the sparse matrix
feature_names <- colnames(data_sparse)

# Compute the importance of features
importance <- xgb.importance(feature_names = feature_names, model = bstSparse)

# Display the importance
print(importance)

# Plot the importance for better visualization
xgb.plot.importance(importance)


# Update the XGBoost training to use "multi:softprob" for probabilities
bstSparse <- xgboost(data = data_sparse, 
                     label = label, 
                     max.depth = 2, 
                     eta = 1, 
                     nthread = 2, 
                     nrounds = 2, 
                     objective = "multi:softprob", 
                     num_class = length(unique(label)))

# Predict probabilities for each class on the test set
pred_prob <- predict(bstSparse, test_sparse)

# Convert predictions into a matrix (each column corresponds to a class probability)
pred_prob_matrix <- matrix(pred_prob, ncol = length(unique(label)), byrow = TRUE)

# Load the pROC library
library(pROC)

# Initialize a list to store AUCs for each class
auc_list <- list()

# Calculate AUC for each class
for (i in 1:length(unique(label))) {
  # One-vs-all approach for AUC calculation
  roc_curve <- roc(as.numeric(test$Syndecan_group == (i)), pred_prob_matrix[, i])
  auc_value <- auc(roc_curve)
  auc_list[[i]] <- auc_value
  
  # Print AUC for each class
  print(paste("AUC for class", i, ":", auc_value))
  
  # Plot the ROC curve
  plot(roc_curve, main = paste("ROC Curve for Class", i), col = i, lwd = 2)
}

# Average AUC across all classes
avg_auc <- mean(unlist(auc_list))
print(paste("Average AUC:", avg_auc))

