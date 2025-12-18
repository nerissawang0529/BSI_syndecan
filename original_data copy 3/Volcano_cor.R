rm(list = ls())

# Load required libraries
library(ggplot2)
library(dplyr)
library(ggrepel)

# Load required libraries
library(ggplot2)
library(dplyr)
library(ggrepel)

# Prepare data for volcano plot
results_limma_cor <- results_limma_cor %>%
  mutate(
    log10_padj = -log10(adjusted_p_value),  # Convert adjusted p-values to -log10 scale
    Significance = case_when(
      adjusted_p_value < 0.05 & estimate > 0 ~ "Positively Correlated",
      adjusted_p_value < 0.05 & estimate < 0 ~ "Negatively Correlated",
      TRUE ~ "Not Significant"
    )
  )

# Select top 5 most significant (lowest adjusted p-value) genes on both sides
top_neg <- results_limma_cor %>%
  filter(Significance == "Negatively Correlated") %>%
  arrange(adjusted_p_value) %>%
  head(5)  # Select top 5 negatively correlated genes

top_pos <- results_limma_cor %>%
  filter(Significance == "Positively Correlated") %>%
  arrange(adjusted_p_value) %>%
  head(5)  # Select top 5 positively correlated genes

significant_genes <- bind_rows(top_neg, top_pos)  # Combine both sets

# Volcano plot
volcano_gene_plot <- ggplot(results_limma_cor, aes(x = estimate, y = log10_padj, color = Significance)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = c("Positively Correlated" = "red", "Negatively Correlated" = "blue", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-0.4, 0.4), linetype = "dashed", color = "black") +  # Threshold for correlation
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  # Significance threshold
  geom_text_repel(data = significant_genes, aes(label = GeneID), size = 3, box.padding = 0.4, max.overlaps = 10) +  # Label top genes
  theme_minimal() +
  labs(x = "Spearman Correlation Estimate",
       y = "-log10 Adjusted p-value") +
  theme(
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

# Define significance thresholds
p_value_threshold <- 0.05
correlation_threshold <- 0.2

# Count significantly upregulated genes (positively correlated)
upregulated_count <- sum(results_limma_cor$adjusted_p_value < p_value_threshold & 
                           results_limma_cor$estimate > correlation_threshold)

# Count significantly downregulated genes (negatively correlated)
downregulated_count <- sum(results_limma_cor$adjusted_p_value < p_value_threshold & 
                             results_limma_cor$estimate < -correlation_threshold)

# Print results
cat("Upregulated genes:", upregulated_count, "\n")
cat("Downregulated genes:", downregulated_count, "\n")

print(volcano_gene_plot)