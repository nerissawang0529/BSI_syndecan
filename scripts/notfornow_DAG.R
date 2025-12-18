# Install dagitty if not installed
if (!requireNamespace("dagitty", quietly = TRUE)) {
  install.packages("dagitty")
}

# Load the dagitty library
library(dagitty)

dag <- dagitty("
dag {
  Syndecan1 [outcome]
  Coagulation [exposure]
  Thrombosis
  Platelets
  Diabetes
  CKD

  Coagulation -> Syndecan1
  Thrombosis -> Coagulation
  Platelets -> Coagulation
  Diabetes -> Syndecan1
  CKD -> Syndecan1
  Diabetes -> CKD
}
")

plot(dag)
# Plot the DAG
plot(dag)

# Optional: Display the adjustment sets for valid causal inference
adjustmentSets(dag)
