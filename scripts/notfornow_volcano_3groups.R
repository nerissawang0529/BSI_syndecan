rm(list = ls())

## install.packages("pacman")
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("Biobase")
library(Biobase)
BiocManager::install("limma")
if (!requireNamespace("biomaRt", quietly = TRUE)) {
  install.packages("biomaRt")
}
install.packages("RSQLite")
install.packages("ggrepel")

library(ggrepel)
library(biomaRt)
library(RSQLite)
library(tidyverse)
library(dplyr)
library(rstatix)
library(ggpubr)
library(dunn.test)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("limma", quietly = TRUE)) {
  install.packages("limma")
}
BiocManager::install("limma")
library(limma)

# Read data
BSIonly3platforms <-  readRDS("original_data/BSIonly3platforms.Rdata")
BSIpatientdata_280 <- readRDS("original_data/BSIpatientdata_280.Rdata")

BSIonly3platforms$platformID <- rownames(BSIonly3platforms)
BSIpatientdata_280_2 <- BSIpatientdata_280[,c("platformID","MARSID")]

merged_BSI <- merge(BSIonly3platforms, BSIpatientdata_280_2, by = "platformID", all.x = TRUE)
merged_BSI$Case <- NULL

#data for syndecan group
data <- read.csv("original_data/clinical_marker_scource_pathogen.csv")

# Cutting the sorted data into three equal parts (1, 2, 3)
data$Syndecan_group <- cut(data$Syndecan.1.CD138..29.,
breaks = quantile(data$`Syndecan.1.CD138..29.`, probs = seq(0, 1, length = 4), na.rm = TRUE),
include.lowest = TRUE,
labels = c(1, 2, 3))




data2 <- data[,c("ICU_ID_from_datasource","Syndecan_group","APACHE_IV_Score.x","SOFAtot_new")]
data2 <- data2[!is.na(data2$SOFAtot_new), ]
#data2 <- data[,c("ICU_ID_from_datasource","Syndecan_group")]

BSI_transcriptome <- merge(data2, merged_BSI, by.x = "ICU_ID_from_datasource",by.y = "MARSID", all = FALSE)
BSI_transcriptome$platformID <- NULL

##pick Group1 and Group3
BSI_transcriptome_2 <-  BSI_transcriptome %>% filter(Syndecan_group != 2)
BSI_transcriptome_2$Group <- factor(BSI_transcriptome_2$Syndecan_group, levels = c(1, 3))

BSI_transcriptome_2$Syndecan_group <- NULL

# Load your data and create the design matrix###################################################################
#design <- model.matrix(~ 0+ Group, data = BSI_transcriptome_2)
design <- model.matrix(~ 0 + Group + APACHE_IV_Score.x + SOFAtot_new, data = BSI_transcriptome_2)

# Check the design matrix
small <- BSI_transcriptome_2[ , c(4:13230)] %>% as.data.frame()
row.names(small) <- BSI_transcriptome_2$geneID
small$geneID <- NULL

small <- data.frame(lapply(small, function(x) as.numeric(as.character(x))))
small <- as.matrix(small)

# Fit linear model
fit <- lmFit(t(small), design)

# Create contrasts
cont.matrix <- makeContrasts(
  logFC = Group3 - Group1,
  levels = design
)

# Fit the contrasts
fit2 <- contrasts.fit(fit, cont.matrix)
fit3 <- eBayes(fit2)

# Get differential expression results, 'by performing BH adjusted moderated t-statistics using limma'
fit3 <- topTable(fit3, number = Inf, adjust.method = "BH")
results_limma <- fit3

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
file_path <- "Documents/ER_WARD_ICU/DGE AND Reactome/DGE_vocanol_after_mews_20241015.svg" 
ggsave(file_path, plot = p, width = 10, height = 10)
