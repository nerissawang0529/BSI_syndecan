rm(list = ls())

#
if (!require("pacman")) install.packages("pacman")
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               kmed, survival, survminer, skimr, grDevices, ggrepel, limma, ppcor,
               AnnotationDbi, org.Hs.eg.db)

# 
BSIonly3platforms <-  readRDS("original_data/BSIonly3platforms.Rdata")
BSIpatientdata_280 <- readRDS("original_data/BSIpatientdata_280.Rdata")

BSIonly3platforms$platformID <- rownames(BSIonly3platforms)
BSIpatientdata_280_2 <- BSIpatientdata_280[,c("platformID","MARSID")]

merged_BSI <- merge(BSIonly3platforms, BSIpatientdata_280_2, by = "platformID", all.x = TRUE)
merged_BSI$Case <- NULL

# 
data <- read.csv("original_data/final_data.csv")
data2 <- data[,c("ICU_ID_from_datasource","Syndecan1","APACHE_IV_Score","SOFAtot_new")]
data2 <- data2[!is.na(data2$SOFAtot_new), ]

BSI_transcriptome <- merge(data2, merged_BSI, by.x = "ICU_ID_from_datasource", by.y = "MARSID", all = FALSE)
BSI_transcriptome$platformID <- NULL

# 
BSI_transcriptome[, -c(1, 2, 3, 4)] <- lapply(BSI_transcriptome[, -c(1, 2, 3, 4)], as.numeric)

# 
gene_ids <- character()
p_values <- numeric()
estimates <- numeric()

# 
gene_columns <- setdiff(colnames(BSI_transcriptome), c("ICU_ID_from_datasource", "Syndecan1", "APACHE_IV_Score", "SOFAtot_new"))

# 
for (gene in gene_columns) {
  x <- BSI_transcriptome$Syndecan1
  y <- BSI_transcriptome[[gene]]
  
  if (sum(!is.na(x) & !is.na(y)) > 3) {  
    cor_test_result <- cor.test(x, y, method = "spearman")
    gene_ids <- c(gene_ids, gene)
    p_values <- c(p_values, cor_test_result$p.value)
    estimates <- c(estimates, unname(cor_test_result$estimate))
  }
}

# 
adjusted_p_values <- p.adjust(p_values, method = "BH")

# 
cor_results_df <- data.frame(
  GeneID = gene_ids,
  p_value = p_values,
  adjusted_p_value = adjusted_p_values,
  estimate = estimates
)

# 
results_limma_cor <- cor_results_df
rownames(results_limma_cor) <- results_limma_cor$GeneID


#Create a named vector of correlation estimates sorted in decreasing order
#gene_list_cor <- results_limma_cor$estimate
#names(gene_list_cor) <- row.names(results_limma_cor)
#gene_list_cor <- sort(gene_list_cor, decreasing = TRUE)

#Create a dataframe with gene symbols and estimate values
#gene_df_cor <- data.frame(gene_symbol = names(gene_list_cor), logFC = gene_list_cor)
#rownames(gene_df_cor) <- names(gene_list_cor)

# Prepare gene list for pathway analysis
#set.seed(11)  # Ensure reproducibility
#gene_list_cor_symbols <- gene_df_cor$logFC
#names(gene_list_cor_symbols) <- gene_df_cor$gene_symbol
#gene_list_cor_symbols <- sort(gene_list_cor_symbols, decreasing = TRUE)

# Convert gene symbols to Entrez IDs
#gene_df <- data.frame(gene_symbol = names(gene_list_cor_symbols), logFC = gene_list_cor_symbols)


library(org.Hs.eg.db)
ens_keys <- keys(org.Hs.eg.db, keytype = "ENSEMBL")


rownames(cor_results_df) <- cor_results_df$GeneID

# Assuming the rownames of your data contain gene symbols, use keytype = "SYMBOL"
cor_results_df$entrez <- mapIds(
  org.Hs.eg.db,
  keys = rownames(cor_results_df),     # Assuming the rownames contain gene symbols
  keytype = "SYMBOL",           # Use SYMBOL since your data contains gene names
  column = "ENTREZID",          # The column you want to map to
  multiVals = "first"           # Handle duplicates by taking the first match
)

# Merge converted IDs back to gene listgene_list_entrez <- gene_df$logFC
#set.seed(11)
gene_list_entrez <- cor_results_df$estimate
names(gene_list_entrez) <- cor_results_df$entrez
gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)
gene_list_entrez <- gene_list_entrez[!is.na(names(gene_list_entrez))]

#gene_list_entrez_df <- gene_list_entrez %>% as.data.frame()
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


# ==== 准备数据 ====
library(dplyr)
library(ggplot2)
library(forcats)
library(tibble)
library(ggtext)

# 假设你的 GSEA 结果保存在 `pathway_results`
# 并且 reactome parent-children 关系保存在 reactomepathways.txt
pathway_results <- y %>% as.data.frame()
RPR <- read.table("original_data/reactomepathways.txt", sep="", stringsAsFactors = FALSE)
colnames(RPR) <- c("parent", "child")

# 目标母通路
selected_parents <- c("R-HSA-109582", "R-HSA-1474244")  # Hemostasis & ECM


# ==== 获取母通路及子通路 ====
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

# ==== 合并富集结果 ====
df <- merge(pathway_results, all_paths, by = "ID")
df$sig_label <- ifelse(df$p.adjust < 0.05, "FDR < 0.05", "Not significant")

# ==== 每组保留 Top10 子通路 + 母通路 ====
df_top <- df %>%
  filter(!is_mother) %>%
  group_by(parent_cat) %>%
  slice_max(order_by = abs(NES), n = 10, with_ties = FALSE) %>%
  ungroup()

df_mothers <- df %>% filter(is_mother)
df <- bind_rows(df_mothers, df_top)

df <- df %>%
  mutate(
    Description = case_when(
      Description == "Factors involved in megakaryocyte development and platelet production" ~ 
        "Megakaryocyte development and platelet production",
      TRUE ~ Description
    )
  )


# ==== 设置排序顺序：母通路放在组内最上面 ====
df <- df %>%
  mutate(Description_fmt = ifelse(is_mother, paste0("**", Description, "**"), Description)) %>%
  group_by(parent_cat) %>%
  arrange(dplyr::desc(as.numeric(is_mother)), dplyr::desc(abs(NES))) %>%
  mutate(Description_factor = factor(Description_fmt, levels = rev(unique(Description_fmt)))) %>%
  ungroup()

# ==== 绘图 ====
# Define font sizes (you can modify these values as needed)
font_size <- 14
axis_font_size <- 14
axis_text_size <- 12

# ==== GSEA plot (p1) with updated font styling ====
p1 <- ggplot(df, aes(x = NES, y = Description_factor, fill = NES)) +
  geom_point(aes(size = sig_label, shape = sig_label), color = "black") +
  scale_fill_gradient2(
    low = "#4575b4", mid = "white", high = "#d73027",
    name = "Enrichment Score (NES)"
  ) +
  scale_size_manual(
    name = "Significance",
    values = c("FDR < 0.05" = 4, "Not significant" = 2)
  ) +
  scale_shape_manual(
    name = "Significance",
    values = c("FDR < 0.05" = 21, "Not significant" = 23)
  ) +
  facet_wrap(~ parent_cat, scales = "free_y", ncol = 1, strip.position = "top") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal(base_size = font_size) +
  theme(
    text = element_text(family = "sans", size = font_size),
    axis.text.y = ggtext::element_markdown(family = "sans", size = axis_text_size),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(family = "sans", size = axis_text_size),
    strip.text = element_blank(),
    legend.position = "none",
    plot.title = element_blank()
  )


p1

# ==== Save to file ====
ggsave("figures/GSEA_cor.svg", p1, width = 9, height = 6, dpi = 300)













# 添加横向 jitter 但保留原始 y 坐标
set.seed(42)
results_limma_cor_jitterx <- results_limma_cor %>%
  mutate(
    log10_padj = -log10(adjusted_p_value),
    estimate_jitter = estimate + runif(n(), min = -0.01, max = 0.01),
    Significance = case_when(
      adjusted_p_value < 0.05 & estimate > 0 ~ "Positively Correlated",
      adjusted_p_value < 0.05 & estimate < 0 ~ "Negatively Correlated",
      TRUE ~ "Not Significant"
    )
  )

# Top gene 标注 jitter 也同步
top_neg <- results_limma_cor_jitterx %>%
  filter(Significance == "Negatively Correlated") %>%
  arrange(adjusted_p_value) %>%
  head(5)

top_pos <- results_limma_cor_jitterx %>%
  filter(Significance == "Positively Correlated") %>%
  arrange(adjusted_p_value) %>%
  head(5)

significant_genes <- bind_rows(top_neg, top_pos)

# 画图（x轴加抖动，y轴逻辑一致）
p2 <- ggplot(results_limma_cor_jitterx, aes(x = estimate_jitter, y = log10_padj, color = Significance)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("Positively Correlated" = "red", "Negatively Correlated" = "blue", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-0.4, 0.4), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(data = significant_genes,
                  aes(x = estimate_jitter, y = log10_padj, label = GeneID),
                  size = 3, color = "black", box.padding = 0.4, max.overlaps = 10) +
  theme_minimal(base_size = font_size) +
  labs(
    y = "-log10 Adjusted p-value"  # removed x label here
  ) +
  theme(
    text = element_text(family = "sans", size = font_size),
    axis.title = element_text(family = "sans", size = axis_font_size),
    axis.text = element_text(family = "sans", size = axis_text_size),
    axis.title.x = element_blank(),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.title = element_blank()
  )

p2

# ==== Combine both with A/B tags ====
library(patchwork)
combined_plot <- p2 + p1 +
  plot_layout(ncol = 2, widths = c(1.2, 1)) +
  plot_annotation(
    tag_levels = 'A',
    theme = theme(
      plot.tag = element_text(size = 16, face = "bold", family = "sans"),
      plot.title = element_blank()
    )
  )

# ==== Save to file ====
ggsave("figures/combined_GSEA_Volcano_horizontal.svg", combined_plot, width = 14, height = 6, dpi = 300)

summary_counts <- results_limma_cor_jitterx %>%
  dplyr::filter(Significance != "Not Significant") %>%
  dplyr::count(Significance)

print(summary_counts)
#1 Negatively Correlated 135
#2 Positively Correlated 203









##method1 adjust severity####
# Initialize vectors to store the geneID, p-value, and correlation estimate
gene_ids <- character()
p_values <- numeric()
estimates <- numeric()

# Get the column names for gene expression
gene_columns <- setdiff(colnames(BSI_transcriptome), c("ICU_ID_from_datasource", "syndecan_original", "APACHE_IV_Score.x", "SOFAtot_new.x"))

# Loop through each gene and perform regression
for (gene in gene_columns) {
  model <- lm(BSI_transcriptome[[gene]] ~ BSI_transcriptome$syndecan_original + 
                BSI_transcriptome$APACHE_IV_Score.x + 
                BSI_transcriptome$SOFAtot_new.x, data = BSI_transcriptome)
  
  # Extract the coefficient, p-value, and gene ID
  summary_model <- summary(model)
  gene_ids <- c(gene_ids, gene)
  p_values <- c(p_values, summary_model$coefficients["BSI_transcriptome$syndecan_original", "Pr(>|t|)"])
  estimates <- c(estimates, summary_model$coefficients["BSI_transcriptome$syndecan_original", "Estimate"])
}

# Adjust the p-values using the Benjamini-Hochberg method
adjusted_p_values <- p.adjust(p_values, method = "BH")

# Create a dataframe with the results
cor_results_df <- data.frame(
  GeneID = gene_ids,
  p_value = p_values,
  adjusted_p_value = adjusted_p_values,
  estimate = estimates
)

# Optionally, store this in a variable for further analysis
results_limma_cor <- cor_results_df
rownames(results_limma_cor) <- results_limma_cor$GeneID

# Create a sorted vector for pathway analysis
gene_list_cor <- results_limma_cor$estimate
names(gene_list_cor) <- rownames(results_limma_cor)
gene_list_cor <- sort(gene_list_cor, decreasing = TRUE)

# Create a dataframe with gene symbols and estimate values
gene_df_cor <- data.frame(gene_symbol = names(gene_list_cor), logFC = gene_list_cor)
rownames(gene_df_cor) <- names(gene_list_cor)

# Prepare gene list for pathway analysis
set.seed(11)  # Ensure reproducibility
gene_list_cor_symbols <- gene_df_cor$logFC
names(gene_list_cor_symbols) <- gene_df_cor$gene_symbol
gene_list_cor_symbols <- sort(gene_list_cor_symbols, decreasing = TRUE)

# Convert gene symbols to Entrez IDs
gene_df <- data.frame(gene_symbol = names(gene_list_cor_symbols), logFC = gene_list_cor_symbols)


library(org.Hs.eg.db)
ens_keys <- keys(org.Hs.eg.db, keytype = "ENSEMBL")


rownames(gene_df) <- names(gene_list_cor_symbols)

# Assuming the rownames of your data contain gene symbols, use keytype = "SYMBOL"
gene_df$entrez <- mapIds(
  org.Hs.eg.db,
  keys = rownames(gene_df),     # Assuming the rownames contain gene symbols
  keytype = "SYMBOL",           # Use SYMBOL since your data contains gene names
  column = "ENTREZID",          # The column you want to map to
  multiVals = "first"           # Handle duplicates by taking the first match
)

# Merge converted IDs back to gene listgene_list_entrez <- gene_df$logFC
set.seed(11)
gene_list_entrez <- gene_df$logFC
names(gene_list_entrez) <- gene_df$entrez
gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)
gene_list_entrez_df <- gene_list_entrez %>% as.data.frame()

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
yi_1$new_adjusted_pvalue <- p.adjust(yi$pvalue, method = "BH")

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
#####

##method two:Partial Correlation####
library(ppcor)
# Initialize vectors to store gene information
gene_ids <- character()
p_values <- numeric()
estimates <- numeric()
# Get the column names except for syndecan_original, and metadata columns
gene_columns <- setdiff(colnames(BSI_transcriptome), c("ICU_ID_from_datasource", "syndecan_original", "APACHE_IV_Score.x", "SOFAtot_new.x"))

# Loop through each gene and calculate partial correlation
for (gene in gene_columns) {
  pc <- pcor.test(BSI_transcriptome[[gene]], BSI_transcriptome$syndecan_original,
                  cbind(BSI_transcriptome$APACHE_IV_Score.x, BSI_transcriptome$SOFAtot_new.x), method = "spearman")
  
  gene_ids <- c(gene_ids, gene)
  p_values <- c(p_values, pc$p.value)
  estimates <- c(estimates, pc$estimate)
}

# Adjust the p-values
adjusted_p_values <- p.adjust(p_values, method = "BH")

# Create results dataframe
cor_results_df <- data.frame(
  GeneID = gene_ids,
  p_value = p_values,
  adjusted_p_value = adjusted_p_values,
  estimate = estimates
)

# Optionally, store this in a variable for further analysis
results_limma_cor <- cor_results_df
rownames(results_limma_cor) <- results_limma_cor$GeneID

# Sort genes for pathway analysis
#gene_list_cor <- cor_results_df$estimate
#names(gene_list_cor) <- cor_results_df$GeneID
#gene_list_cor <- sort(gene_list_cor, decreasing = TRUE)


#Create a dataframe with gene symbols and estimate values
#gene_df_cor <- data.frame(gene_symbol = names(gene_list_cor), logFC = gene_list_cor)
#rownames(gene_df_cor) <- names(gene_list_cor)

# Prepare gene list for pathway analysis
#set.seed(11)  # Ensure reproducibility
#gene_list_cor_symbols <- gene_df_cor$logFC
#names(gene_list_cor_symbols) <- gene_df_cor$gene_symbol
#gene_list_cor_symbols <- sort(gene_list_cor_symbols, decreasing = TRUE)

# Convert gene symbols to Entrez IDs
#gene_df <- data.frame(gene_symbol = names(gene_list_cor_symbols), logFC = gene_list_cor_symbols)


library(org.Hs.eg.db)
ens_keys <- keys(org.Hs.eg.db, keytype = "ENSEMBL")


rownames(cor_results_df) <- cor_results_df$GeneID

# Assuming the rownames of your data contain gene symbols, use keytype = "SYMBOL"
cor_results_df$entrez <- mapIds(
  org.Hs.eg.db,
  keys = rownames(cor_results_df),     # Assuming the rownames contain gene symbols
  keytype = "SYMBOL",           # Use SYMBOL since your data contains gene names
  column = "ENTREZID",          # The column you want to map to
  multiVals = "first"           # Handle duplicates by taking the first match
)

# Merge converted IDs back to gene listgene_list_entrez <- gene_df$logFC

gene_list_entrez <- cor_results_df$estimate
names(gene_list_entrez) <- cor_results_df$entrez
gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)
gene_list_entrez_df <- gene_list_entrez %>% as.data.frame()


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

p#####



# ==== GSEA plot (p1) with proper legends ====

p1 <- ggplot(df, aes(x = NES, y = Description_factor)) +
  
  # ---- 主点图：颜色 = NES, shape + size = significance ----
geom_point(aes(fill = NES, size = sig_label, shape = sig_label),
           color = "black", stroke = 0.3) +
  
  # ---- 颜色 legend：NES ----
scale_fill_gradient2(
  low = "#4575b4", mid = "white", high = "#d73027",
  name = "NES"
) +
  
  # ---- shape legend ----
scale_shape_manual(
  name = "Significance",
  values = c("FDR < 0.05" = 21, "Not significant" = 23)
) +
  
  # ---- size legend ----
scale_size_manual(
  name = "Significance",
  values = c("FDR < 0.05" = 4, "Not significant" = 2)
) +
  
  facet_wrap(~ parent_cat, scales = "free_y", ncol = 1, strip.position = "top") +
  
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  theme_minimal(base_size = font_size) +
  theme(
    text = element_text(family = "sans", size = font_size),
    axis.text.y = ggtext::element_markdown(family = "sans", size = axis_text_size),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(family = "sans", size = axis_text_size),
    
    strip.text = element_blank(),
    
    # ---- 关键：显示 legend ----
    legend.position = "right",
    legend.title = element_text(size = 12, family = "sans"),
    legend.text  = element_text(size = 11, family = "sans")
  )

p1
# 把当前的 GSEA NES 点图 p1 导出为 SVG
ggsave(
  filename = "figures/GSEA_NES_plot.svg",  # 保存文件名和路径
  plot     = p1,                           # 要导出的图对象
  width    = 9,                            # 宽（inch，可按需要改）
  height   = 6,                            # 高（inch，可按需要改）
  dpi      = 300                           # 分辨率，对 SVG 影响不大，但可以保留
)
