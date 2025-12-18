rm(list = ls())
library(DALEX)
library(purrr)

# Fit a Linear Model
data(mtcars)
lm_model <- lm(mpg ~ wt + hp + cyl, data = mtcars)

# Define a prediction function wrapper
pfun <- function(model, newdata) predict(model, newdata)

# Create an explainer
explainer_lm <- explain(
  model = lm_model,
  X = mtcars[, c("wt", "hp", "cyl")],  # Features
  y = mtcars$mpg, 
  label = "Linear Model",
  predict_function = pfun
)

# Select multiple new observations (e.g., first 5 rows)
new_obs <- mtcars[1:5, c("wt", "hp", "cyl")]

# Compute SHAP values for multiple observations
shap_results <- map_dfr(1:nrow(new_obs), function(i) {
  shap_values <- predict_parts(
    explainer = explainer_lm,
    new_observation = new_obs[i, , drop = FALSE],  # Select row i
    type = "shap"
  )
  shap_values$observation_id <- i  # Add ID column
  return(shap_values)
})

# Print results
print(shap_results)