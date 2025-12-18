library(ggrepel)
library(biomaRt)
library(RSQLite)
library(tidyverse)
library(dplyr)
library(rstatix)
library(ggpubr)
library(dunn.test)


# Specify the Lancet-style theme
theme_lancet <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(family = "Gill Sans"),
      plot.title = element_text(hjust = 0.5, size = 16),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white")
    )
}


# Calculate the number of significantly upregulated and downregulated genes
significant_upregulated <- results_limma %>% filter(logFC > 0 & adj.P.Val < 0.05) %>%  arrange(desc(logFC))
significant_downregulated <- results_limma %>% filter(logFC < 0 & adj.P.Val < 0.05) %>% arrange(desc(abs(logFC)))

num_upregulated <- nrow(significant_upregulated)
num_downregulated <- nrow(significant_downregulated)


# Create a new column for point color based on significance and logFC
results_limma$Color <- ifelse(results_limma$logFC < 0 & results_limma$adj.P.Val < 0.05, "Lower in high Syndecan", 
                              ifelse(results_limma$logFC > 0 & results_limma$adj.P.Val < 0.05, "Higher in high Syndecan", "Not Significant"))

results_limma$gene_name <- rownames(results_limma)

# Select top 5 upregulated and top 5 downregulated genes based on adjusted p-value
top5_upregulated <- results_limma %>%
  filter(logFC > 0) %>%
  arrange(adj.P.Val) %>%
  head(10)

top5_downregulated <- results_limma %>%
  filter(logFC < 0) %>%
  arrange(adj.P.Val) %>%
  head(10)

top_genes <- bind_rows(top5_upregulated, top5_downregulated)
rownames(top_genes) <- top_genes$gene_name

# Highlight top 5 genes
results_limma$label <- ifelse(row.names(results_limma) %in% row.names(top_genes), row.names(results_limma), NA)


# Apply the Lancet-style theme and customize colors and labels
dev.off()
library(ggrepel)
p <- ggplot(results_limma, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = Color), size = 1) +
  scale_color_manual(values = c("Lower in high Syndecan" = "blue", "Higher in high Syndecan" = "red", "Not Significant" = "gray")) +
  labs(x = "logFC", y = "-log10(adj.P.Val)", title = "Volcano Plot") +
  theme_lancet() +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray") +
  geom_text_repel(data = subset(results_limma, !is.na(label)), aes(label = label), size = 4) +
  theme_minimal(base_size = 11) +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white")
  )
p
