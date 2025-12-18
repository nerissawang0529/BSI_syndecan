rm(list = ls())
library(dplyr)
library(limma)
library(AnnotationDbi)
library(org.Hs.eg.db)  # For human gene annotations
install.packages("pacman")
library(pacman)
library(ggplot2)

# Read data
BSIonly3platforms <-  readRDS("original_data/BSIonly3platforms.Rdata")
BSIpatientdata_280 <- readRDS("original_data/BSIpatientdata_280.Rdata")

BSIonly3platforms$platformID <- rownames(BSIonly3platforms)
BSIpatientdata_280_2 <- BSIpatientdata_280[,c("platformID","MARSID")]

merged_BSI <- merge(BSIonly3platforms, BSIpatientdata_280_2, by = "platformID", all.x = TRUE)
merged_BSI$Case <- NULL

#data for group and severity
data <- read.csv("original_data/merged_data_new_group.csv")
data2 <- data[,c("ICU_ID_from_datasource","sd_group_based_healthy_3","APACHE_IV_Score.x","SOFAtot_new.x")]
data2 <- data2[!is.na(data2$SOFAtot_new), ]

BSI_transcriptome <- merge(data2, merged_BSI, by.x = "ICU_ID_from_datasource",by.y = "MARSID", all = FALSE)
BSI_transcriptome$platformID <- NULL

##pick Group1 and Group3
BSI_transcriptome_2 <- BSI_transcriptome %>% filter(sd_group_based_healthy_3 != 2)
BSI_transcriptome_2$Group <- factor(BSI_transcriptome_2$sd_group_based_healthy_3, levels = c(1, 3))
BSI_transcriptome_2$sd_group_based_healthy_3 <- NULL

# Load your data and create the design matrix###################################################################
design <- model.matrix(~ 0+ Group, data = BSI_transcriptome_2)
#design <- model.matrix(~ 0 + Group + APACHE_IV_Score.x + SOFAtot_new.x, data = BSI_transcriptome_2)


# Check the design matrix
small <- BSI_transcriptome_2[ , c(4:13230)] %>% as.data.frame()
row.names(small) <- BSI_transcriptome_2$ICU_ID_from_datasource

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

#
gene_list <- results_limma$t
names(gene_list) <- row.names(results_limma)
gene_list <- sort(gene_list, decreasing = TRUE)

# Convert gene symbols to Entrez IDs
gene_df <- data.frame(gene_symbol = names(gene_list), logFC = gene_list)
library(org.Hs.eg.db)
ens_keys <- keys(org.Hs.eg.db, keytype = "ENSEMBL")


rownames(gene_df) <- names(gene_list)

# Assuming the rownames of your data contain gene symbols, use keytype = "SYMBOL"
gene_df$entrez <- mapIds(
  org.Hs.eg.db,
  keys = rownames(gene_df),     # Assuming the rownames contain gene symbols
  keytype = "SYMBOL",           # Use SYMBOL since your data contains gene names
  column = "ENTREZID",          # The column you want to map to
  multiVals = "first"           # Handle duplicates by taking the first match
)

# Merge converted IDs back to gene listgene_list_entrez <- gene_df$logFC

gene_list_entrez <- gene_df$logFC
names(gene_list_entrez) <- gene_df$entrez
gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)



#####
set.seed(11)
p_load(ReactomePA, enrichplot)
y <- gsePathway(geneList = gene_list_entrez,
                pvalueCutoff=1, 
                pAdjustMethod="BH", verbose = F, maxGSSize = 10000, minGSSize = 0, seed = T, eps = 0, nPermSimple = 10000)

pathway_results <- y %>% as.data.frame()

## last updated pathways => 18-03-2024
### All plots together
### only host response pathways (adaptive, innate, cytokine signaling, programmed cell death, hemostasis)
RPR <- read.table("original_data/reactomepathways.txt", sep="", stringsAsFactors = F)
names(RPR) <- c("V1", "V2")

selected_pathway_id1 <- c("R-HSA-109582")
selected_pathway_id2 <- c("R-HSA-1474244")

gen0 <- data.frame(ID=selected_pathway_id1, gen=0)
gen1 <- data.frame(ID=subset(RPR, V1 %in% gen0$ID)$V2, gen=1)
gen2 <- data.frame(ID=subset(RPR, V1 %in% gen1$ID)$V2, gen=2)
gen3 <- data.frame(ID=subset(RPR, V1 %in% gen2$ID)$V2, gen=3)
hemo<- rbind(gen0,gen1,gen2,gen3)

gen0 <- data.frame(ID=selected_pathway_id2, gen=0)
gen1 <- data.frame(ID=subset(RPR, V1 %in% gen0$ID)$V2, gen=1)
gen2 <- data.frame(ID=subset(RPR, V1 %in% gen1$ID)$V2, gen=2)
gen3 <- data.frame(ID=subset(RPR, V1 %in% gen2$ID)$V2, gen=3)
extr <- rbind(gen0,gen1,gen2,gen3)

yi <- pathway_results[pathway_results$ID %in% c(hemo$ID, extr$ID) ,]
yi_1 <- yi
yi_1$new_adjusted_pvalue <- yi_1$p.adjust  
library(dplyr)

# Define colors
yi_1 <- yi_1 %>% mutate(
  Color = case_when(
    new_adjusted_pvalue < 0.05 & NES < 0 ~ ifelse(ID == 'gen0', '#1f78b4', '#a6cee3'), # blue, bright blue
    new_adjusted_pvalue < 0.05 & NES > 0 ~ ifelse(ID == 'gen0', '#e31a1c', '#fb9a99'), # red, bright red
    TRUE ~ '#bdbdbd' # grey
  )
)

# Create the plot
p <- ggplot(yi_1, aes(x = Description, y = NES, fill = Color)) +
  geom_bar(stat = 'identity', position = position_dodge()) + # Ensure consistent bar heights
  scale_fill_identity() +
  coord_flip() +
  scale_x_discrete(limits = rev) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(),  # Remove background
    plot.background = element_blank(),   # Remove plot background
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    strip.background = element_blank(),  # Remove strip background
    strip.text = element_blank()         # Remove strip text
  ) +
  labs(title = "Pathway NES Plot",
       x = "Pathway Description",
       y = "NES") +
  geom_hline(yintercept = 0, color = "black")  # Add a black line at y = 0

p


##volcanal
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
p_vol <- ggplot(results_limma, aes(x = logFC, y = -log10(adj.P.Val))) +
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
p_vol









#checked the syndecan groups
# 清除环境
rm(list = ls())

# 加载包
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, limma, AnnotationDbi, org.Hs.eg.db, ggrepel, ReactomePA, enrichplot, ggplot2, patchwork)

# 读取数据
BSIonly3platforms <- readRDS("original_data/BSIonly3platforms.Rdata")
BSIpatientdata_280 <- readRDS("original_data/BSIpatientdata_280.Rdata")
data <- read.csv("original_data/final_data.csv")

# 合并表达矩阵与 MARSID
BSIonly3platforms$platformID <- rownames(BSIonly3platforms)
BSIpatientdata_280_2 <- BSIpatientdata_280[, c("platformID", "MARSID")]
merged_BSI <- merge(BSIonly3platforms, BSIpatientdata_280_2, by = "platformID", all.x = TRUE)
merged_BSI$Case <- NULL

# 创建Syndecan_group（三分位分组）
data$Syndecan_group <- cut(data$Syndecan1,
                           breaks = quantile(data$Syndecan1, probs = seq(0, 1, length = 4), na.rm = TRUE),
                           include.lowest = TRUE,
                           labels = c(1, 2, 3))

# 提取分析所需字段并去除缺失
data2 <- data[, c("ICU_ID_from_datasource", "Syndecan_group", "APACHE_IV_Score", "SOFAtot_new")]
data2 <- data2[!is.na(data2$SOFAtot_new), ]

# 合并表达与分组数据
BSI_transcriptome <- merge(data2, merged_BSI, by.x = "ICU_ID_from_datasource", by.y = "MARSID", all = FALSE)

# 保留组1和组3
BSI_transcriptome_2 <- BSI_transcriptome %>% filter(Syndecan_group != 2)
BSI_transcriptome_2$Group <- factor(BSI_transcriptome_2$Syndecan_group, levels = c(1, 3))

# 构建设计矩阵 + 表达矩阵
design <- model.matrix(~ 0 + Group, data = BSI_transcriptome_2)
expr <- BSI_transcriptome_2[, 4:ncol(BSI_transcriptome_2)] %>% mutate(across(everything(), as.numeric)) %>% as.matrix()
rownames(expr) <- BSI_transcriptome_2$ICU_ID_from_datasource

# 差异表达分析
fit <- lmFit(t(expr), design)
cont.matrix <- makeContrasts(logFC = Group3 - Group1, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit3 <- eBayes(fit2)
results_limma <- topTable(fit3, number = Inf, adjust.method = "BH")

# 转换 SYMBOL -> ENTREZ
gene_list <- results_limma$t
names(gene_list) <- rownames(results_limma)
gene_list <- sort(gene_list, decreasing = TRUE)
gene_df <- data.frame(gene_symbol = names(gene_list), logFC = gene_list)
gene_df$entrez <- mapIds(org.Hs.eg.db, keys = gene_df$gene_symbol, keytype = "SYMBOL", column = "ENTREZID", multiVals = "first")
gene_list_entrez <- gene_df$logFC
names(gene_list_entrez) <- gene_df$entrez
gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)

# GSEA Reactome
set.seed(11)
gsea_result <- gsePathway(geneList = gene_list_entrez,
                          pvalueCutoff = 1,
                          pAdjustMethod = "BH",
                          verbose = FALSE,
                          maxGSSize = 10000,
                          minGSSize = 0,
                          eps = 0,
                          seed = TRUE,
                          nPermSimple = 10000)
pathway_results <- as.data.frame(gsea_result)

# 读取通路层级信息
RPR <- read.table("original_data/reactomepathways.txt", sep = "", stringsAsFactors = FALSE)
colnames(RPR) <- c("parent", "child")

extract_pathways <- function(parent_id, label, depth = 3) {
  result <- data.frame(ID = parent_id, gen = 0)
  current <- parent_id
  for (i in 1:depth) {
    next_gen <- RPR[RPR$parent %in% current, "child", drop = TRUE]
    if (length(next_gen) == 0) break
    result <- rbind(result, data.frame(ID = next_gen, gen = i))
    current <- next_gen
  }
  result$parent_cat <- label
  result$is_mother <- result$gen == 0
  return(result)
}

hemo <- extract_pathways("R-HSA-109582", "Hemostasis")
ecm  <- extract_pathways("R-HSA-1474244", "Extracellular Matrix")
all_paths <- bind_rows(hemo, ecm)
all_paths$parent_cat <- factor(all_paths$parent_cat, levels = c("Hemostasis", "Extracellular Matrix"))

# 合并路径富集结果
df <- merge(pathway_results, all_paths, by = "ID")
df$sig_label <- ifelse(df$p.adjust < 0.05, "FDR < 0.05", "Not significant")

# 保留母通路 + Top10 子通路
df_top <- df %>% filter(!is_mother) %>%
  group_by(parent_cat) %>% slice_max(order_by = abs(NES), n = 10, with_ties = FALSE) %>%
  ungroup()
df_mothers <- df %>% filter(is_mother)
df <- bind_rows(df_mothers, df_top)

df <- df %>%
  mutate(
    is_mother = as.numeric(as.logical(is_mother)),
    Description_fmt = ifelse(is_mother == 1, paste0("**", Description, "**"), Description)
  ) %>%
  group_by(parent_cat) %>%
  arrange(desc(is_mother), desc(abs(NES))) %>%
  mutate(Description_factor = factor(Description_fmt, levels = rev(unique(Description_fmt)))) %>%
  ungroup()

# dot plot
p_gsea <- ggplot(df, aes(x = NES, y = Description_factor, fill = NES)) +
  geom_point(aes(size = sig_label, shape = sig_label), color = "black") +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", name = "NES") +
  scale_size_manual(values = c("FDR < 0.05" = 4, "Not significant" = 2)) +
  scale_shape_manual(values = c("FDR < 0.05" = 21, "Not significant" = 23)) +
  facet_wrap(~ parent_cat, scales = "free_y", ncol = 1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = ggtext::element_markdown(size = 12),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    strip.text = element_blank(),
    legend.position = "none"
  )

# 火山图
results_limma <- results_limma %>%
  mutate(
    log10_padj = -log10(adj.P.Val),
    Color = case_when(
      logFC > 0 & adj.P.Val < 0.05 ~ "Higher in high Syndecan",
      logFC < 0 & adj.P.Val < 0.05 ~ "Lower in high Syndecan",
      TRUE ~ "Not Significant"
    ),
    gene_name = rownames(.)
  )

top_genes <- bind_rows(
  results_limma %>% filter(Color == "Higher in high Syndecan") %>% arrange(adj.P.Val) %>% head(5),
  results_limma %>% filter(Color == "Lower in high Syndecan") %>% arrange(adj.P.Val) %>% head(5)
)

results_limma$label <- ifelse(rownames(results_limma) %in% rownames(top_genes), rownames(results_limma), NA)

p_volcano <- ggplot(results_limma, aes(x = logFC, y = log10_padj, color = Color)) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("Higher in high Syndecan" = "red", "Lower in high Syndecan" = "blue", "Not Significant" = "grey")) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "black") +
  geom_text_repel(data = subset(results_limma, !is.na(label)), aes(label = label), size = 3) +
  theme_minimal(base_size = 14) +
  labs(title = "Volcano Plot", x = "logFC", y = "-log10(adj.P.Val)")

# 合并
combined_plot <- p_gsea + p_volcano + plot_layout(ncol = 2, widths = c(1, 1.5))
combined_plot
