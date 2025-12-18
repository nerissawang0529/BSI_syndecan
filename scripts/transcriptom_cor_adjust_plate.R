rm(list = ls())

if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
  nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
  kmed, survival, survminer, skimr, grDevices, ggrepel, limma, ppcor,
  AnnotationDbi, org.Hs.eg.db, ReactomePA, enrichplot
)

## 1. 读入表达矩阵 & 临床数据 -----------------------------
BSIonly3platforms   <- readRDS("original_data/BSIonly3platforms.Rdata")
BSIpatientdata_280  <- readRDS("original_data/BSIpatientdata_280.Rdata")

BSIonly3platforms$platformID <- rownames(BSIonly3platforms)
BSIpatientdata_280_2 <- BSIpatientdata_280[, c("platformID", "MARSID")]

merged_BSI <- merge(
  BSIonly3platforms,
  BSIpatientdata_280_2,
  by = "platformID",
  all.x = TRUE
)
merged_BSI$Case <- NULL

# 临床数据：多加一个 Platelets_value_1
data  <- read.csv("original_data/final_data.csv")
data2 <- data[, c("ICU_ID_from_datasource",
                  "Syndecan1",
                  "APACHE_IV_Score",
                  "SOFAtot_new",
                  "Platelets_value_1")]
data2 <- data2[!is.na(data2$SOFAtot_new), ]

# 合并表达 + 临床
BSI_transcriptome <- merge(
  data2,
  merged_BSI,
  by.x = "ICU_ID_from_datasource",
  by.y = "MARSID",
  all = FALSE
)
BSI_transcriptome$platformID <- NULL

# 基因表达列转 numeric
BSI_transcriptome[, -c(1, 2, 3, 4, 5)] <-
  lapply(BSI_transcriptome[, -c(1, 2, 3, 4, 5)], as.numeric)

## 2. 计算“校正血小板后的”偏相关 --------------------------
gene_ids <- character()
p_values <- numeric()
estimates <- numeric()

# 基因列 = 除去 ID、Syndecan1、APACHE、SOFA、Platelets
gene_columns <- setdiff(
  colnames(BSI_transcriptome),
  c("ICU_ID_from_datasource",
    "Syndecan1",
    "APACHE_IV_Score",
    "SOFAtot_new",
    "Platelets_value_1")
)

for (gene in gene_columns) {
  x <- BSI_transcriptome$Syndecan1
  y <- BSI_transcriptome[[gene]]
  z <- BSI_transcriptome$Platelets_value_1   # 要调整的血小板
  
  idx <- complete.cases(x, y, z)
  if (sum(idx) > 3) {
    # Spearman 偏相关：Syndecan1 ~ gene | Platelets_value_1
    pc <- ppcor::pcor.test(
      x[idx],
      y[idx],
      z[idx],
      method = "spearman"
    )
    
    gene_ids  <- c(gene_ids, gene)
    p_values  <- c(p_values, pc$p.value)
    estimates <- c(estimates, unname(pc$estimate))  # 偏相关系数
  }
}

# 多重校正（FDR）
adjusted_p_values <- p.adjust(p_values, method = "BH")

cor_results_df <- data.frame(
  GeneID          = gene_ids,
  p_value         = p_values,
  adjusted_p_value= adjusted_p_values,
  estimate        = estimates
)
rownames(cor_results_df) <- cor_results_df$GeneID

## 3. 基因 ID 转 Entrez，并生成 GSEA 用的向量 ----------------
cor_results_df$entrez <- mapIds(
  org.Hs.eg.db,
  keys     = rownames(cor_results_df),  # 这里假设 rownames 是 SYMBOL
  keytype  = "SYMBOL",
  column   = "ENTREZID",
  multiVals= "first"
)

# 以“偏相关系数”作为 rank 值，名字是 Entrez ID
gene_list_entrez <- cor_results_df$estimate
names(gene_list_entrez) <- cor_results_df$entrez
gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)
gene_list_entrez <- gene_list_entrez[!is.na(names(gene_list_entrez))]

## 4. 运行 Reactome GSEA（和你之前相同） ---------------------
set.seed(11)
y <- gsePathway(
  geneList     = gene_list_entrez,
  pvalueCutoff = 1,
  pAdjustMethod= "BH",
  verbose      = FALSE,
  maxGSSize    = 10000,
  minGSSize    = 0,
  seed         = TRUE,
  eps          = 0,
  nPermSimple  = 10000
)

pathway_results <- y %>% as.data.frame()

# 之后的 hemostasis / ECM 选路、分组、画图部分类似你原来那段继续往下接就可以

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
ggsave("figures/GSEA_cor_adjust_platelet.svg", p1, width = 9, height = 6, dpi = 300)
