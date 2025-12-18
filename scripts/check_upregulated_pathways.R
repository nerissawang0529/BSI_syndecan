#need to run the transcriptome script first




# Load necessary libraries
library(dplyr)
library(ggplot2)
library(pheatmap)

#Check logFC for Hemostasis ("R-HSA-109582") & ECM Pathways (, "R-HSA-1474244")================

# Convert row names to a column named "gene_symbol"
results_limma <- results_limma %>%
  rownames_to_column(var = "gene_symbol")

# Extract genes from significant pathways (Example: Hemostasis and ECM pathways)
significant_pathways <- c("R-HSA-1474244")  # Hemostasis & ECM Pathways
significant_genes <- pathway_results %>%
  filter(ID %in% significant_pathways) %>%
  select(core_enrichment) %>%
  unlist() %>%
  unique()

library(org.Hs.eg.db)
# Split the character string into a vector of Entrez IDs
significant_genes <- unlist(strsplit(significant_genes, "/"))

# Ensure it's a unique list of valid numeric Entrez IDs
significant_genes <- unique(significant_genes)

# Convert Entrez IDs to Gene Symbols
gene_symbols <- mapIds(org.Hs.eg.db, 
                       keys = significant_genes, 
                       column = "SYMBOL", 
                       keytype = "ENTREZID", 
                       multiVals = "first")

# Remove NA values
gene_symbols <- na.omit(gene_symbols)

# Convert to a character vector for filtering
gene_symbols <- as.character(gene_symbols)


# Filter results_limma to keep only rows where gene_symbol is in gene_symbols
hemostasis_results <- results_limma %>%
  filter(gene_symbol %in% gene_symbols)


# Plot LogFC distribution
ggplot(hemostasis_results, aes(x = logFC)) +
  geom_histogram(binwidth = 0.2, fill = "skyblue", color = "black", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "LogFC Distribution for Significant Extracellular matrix organization Genes",
       x = "log Fold Change (logFC)", y = "Frequency") +
  theme_minimal()


