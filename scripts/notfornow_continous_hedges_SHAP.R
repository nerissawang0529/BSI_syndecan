rm(list = ls())

#ready packages
library(dbarts)
library(shapviz)
library(fastshap)

#ready data
merged_data <- read.csv("original_data/merged_data.csv")

# Define the target variable and predictor markers
target <- "Syndecan1"
markers <- c("Antitrombin", "PROC", "Ang1", "Ang2", "CD163", "CoagulationFactorIII",
             "CX3CL1", "Ddimer", "ESM1", "IL10", "IL18", "IL1ra", "IL6", 
             "IL8", "MMP8", "NGAL", "Procalcitonin", "RAGE", "TenascinC", "Thrombomodulin",
             "TREM1", "Platelets_value_1", "PT_Max_24h")

# Prepare data: Select relevant columns and remove NA values
merged_data <- merged_data[, c(target, markers)]


# Train-test split
set.seed(123)
n <- nrow(merged_data)
ntrain <- round(0.6 * n)
train_indices <- sample(1:n, ntrain)

train_data <- merged_data[train_indices, ]
test_data <- merged_data[-train_indices, ]

X_train <- train_data[, markers]
Y_train <- train_data[, target]
X_test <- test_data[, markers]
Y_test <- test_data[, target]

# Fit Linear Model
linmod <- lm(Y_train ~ ., data = train_data)

# Fit BART Model
fitb <- bart(x.train = X_train, y.train = Y_train, keeptrees = TRUE, ndpost = 1000)

# Predictions
pred_lm <- predict(linmod, newdata = test_data)
pred_bart <- apply(predict(fitb, newdata = X_test, type = "response", combineChains = TRUE), 2, mean)

# Compute RMSE
rmse_lm <- sqrt(mean((pred_lm - Y_test)^2))
rmse_bart <- sqrt(mean((pred_bart - Y_test)^2))

# Print RMSE results
cat("RMSE Linear Model:", rmse_lm, "\n")
cat("RMSE BART Model:", rmse_bart, "\n")

# Compute SHAP values for BART
pfun <- function(object, newdata) {
  colMeans(predict(object, newdata = newdata))
}

shap_values <- explain(fitb, X = X_test, pred_wrapper = pfun, nsim = 3, adjust = TRUE, shap_only = FALSE)

# Convert to shapviz object for visualization
sv <- shapviz(shap_values)

# Plot SHAP values for the top 3 markers
par(mfrow = c(1, 3))
top_markers <- names(sort(apply(abs(sv$S), 2, mean), decreasing = TRUE))[1:3]

for (marker in top_markers) {
  marker_index <- which(names(X_test) == marker)
  plot(X_test[, marker_index], sv$S[, marker_index],
       main = paste("SHAP for", marker),
       xlab = marker, ylab = "SHAP Value")
  abline(h = 0, col = "red")
}