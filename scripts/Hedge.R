rm(list = ls())
# Load required libraries
library(esvis)  # Load library for effect size visualization
library(pheatmap)  # Load library for creating heatmaps
library(tidyverse)  # Load library for data manipulation and visualization

# Read data
data <- read.csv("Original_data/final_data.csv")
# List of markers (y-axis variables)
markers <- c("Antitrombin", "PROC", "Ang1", "Ang2", "CD163", "CoagulationFactorIII",
             "CX3CL1", "Ddimer", "ESM1", "IL10", "IL18", "IL1ra", "IL6", 
             "IL8", "MMP8", "NGAL", "Procalcitonin", "RAGE", "TenascinC", "Thrombomodulin",
             "TREM1", "Platelets_value_1", "PT_Max_24h")
# the following is for delete the antithrombin and protein C
markers <- c("Ang1", "Ang2", "CD163", "CoagulationFactorIII",
             "CX3CL1", "Ddimer", "ESM1", "IL10", "IL18", "IL1ra", "IL6", 
             "IL8", "MMP8", "NGAL", "Procalcitonin", "RAGE", "TenascinC", "Thrombomodulin",
             "TREM1", "Platelets_value_1", "PT_Max_24h")
#log the markers data
data[markers] <- log10(data[markers])

#Grouped the patients based on syndecan
data$Syndecan_group <- cut(data$Syndecan1,breaks = quantile(data$Syndecan1, probs = seq(0, 1, length = 4), na.rm = TRUE),
                           include.lowest = TRUE,labels = c(1, 2, 3))


# List of domains
#variables_Inflammation_combine <- c("CD163", "MMP8", "NGAL", "Procalcitonin", "RAGE", "TenascinC", "TREM1","IL10", "IL18", "IL1ra", "IL6", "IL8")
#variables_Endothelial <- c("Ang1", "Ang2", "CX3CL1", "Thrombomodulin", "ESM1")
variables_Coagulation <- c("CoagulationFactorIII", "Ddimer", "Platelets_value_1",  "PT_Max_24h")

# List to store the results
results_list <- list()

# Perform hedg_g for each dependent variable
####if we want to get the Hedges'g for different domains, we need to change the variables_domains##############################################################
for (var in variables_Coagulation) {
  # Conduct hedg_g analysis
  result <- hedg_g(data, 
                   as.formula(paste(var, "~ Syndecan_group")), # column with group information
                   ref_group = "1") # specify the reference group
  result <- mutate(result, Dependent_Variable = var)
  results_list[[var]] <- result
}

# Combine the results into a single DataFrame
combined_results <- bind_rows(results_list)

# Print the combined DataFrame
print(combined_results)

# Prepare data for the heatmap
data_heatmap_hedg_g <- select(combined_results, Dependent_Variable, hedg_g, Syndecan_group_foc)  # Select relevant columns
data_heatmap_hedg_g <- data_heatmap_hedg_g %>%
  pivot_wider(names_from = Syndecan_group_foc, values_from = hedg_g)  # Reshape data for heatmap
data_heatmap_hedg_g <- column_to_rownames(data_heatmap_hedg_g, var = "Dependent_Variable")  # Set variable as row names


# Print the heatmap data
print(data_heatmap_hedg_g)

# Define breaks and colors for the heatmap
breaks <- c(-Inf, -1.5, -0.8, -0.5, -0.2, 0.2, 0.5, 0.8, 1.5, Inf)

#for the figure
##Question:Is the color I set correct?
colors <- c("#f70d1a","#ff533f","#ff9e88","#ffbfae", "white", "#9F9FFF","#4444FF","#0000E7","#00008B")

####if we want to get the Hedges'g for different domains, we need to change the variables_domains##############################################################
#output the figure
output_pdf_file <- "figures/hegde_coagulation_groups.pdf"
pdf(file = output_pdf_file, width = 7, height = 7)


# Create the heatmap using pheatmap with custom breaks and colors
pheatmap(data_heatmap_hedg_g,
         scale = "none",
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize = 25,
         border_color = "black",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         width = 2, height = 5, dpi = 600,
         color = colors,
         breaks = breaks,
         legend = FALSE
)

dev.off()
