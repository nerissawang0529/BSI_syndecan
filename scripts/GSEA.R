set.seed(123)
library("fgsea")
library(RColorBrewer)
library(fgsea)
library(enrichplot)
library(ggplot2)

# Define pathway IDs or descriptions for "Hemostasis" and "Extracellular Matrix Organization"
hemostasis_id <- "R-HSA-109582"  
ecm_id <- "R-HSA-1474244"  
platelet_id <- "R-HSA-76002"

# Filter the GSEA results for Hemostasis and ECM
hemostasis_results <- pathway_results[pathway_results$ID == hemostasis_id, ]
ecm_results <- pathway_results[pathway_results$ID == ecm_id, ]
platelet_results <- pathway_results[pathway_results$ID == platelet_id, ]
# Print the results to check
print(hemostasis_results)
print(ecm_results)
print(platelet_results)

# Plot GSEA for Hemostasis
p_hemostasis <- gseaplot2(y, geneSetID = hemostasis_id,  pvalue_table = TRUE, title = "GSEA for Hemostasis")

# Plot GSEA for Extracellular Matrix Organization
p_ecm <- gseaplot2(y, geneSetID = ecm_id, pvalue_table = TRUE, title = "GSEA for Extracellular Matrix Organization")

#Plot GSEA for platelet
p_platelet <- gseaplot2(y, geneSetID = platelet_id, pvalue_table = TRUE, title = "GSEA for Platelet activation, signaling and aggregation")



# Print the plots
print(p_hemostasis)
print(p_ecm)
print(p_platelet)
file_path <- "Documents/ER_WARD_ICU/DGE AND Reactome/hemostasis_GSEA_after_Mews.svg" 
ggsave(file_path, plot = p_hemostasis, width = 8, height = 8)

file_path <- "Documents/ER_WARD_ICU/DGE AND Reactome/ecm_GSEA_after_Mews.svg" 
ggsave(file_path, plot = p_ecm, width = 8, height = 8)
