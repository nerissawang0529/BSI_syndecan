# 清除环境
rm(list = ls())

# 加载必要包
if (!require("pacman")) install.packages("pacman")
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, ggtext,
               limma, ggrepel, ReactomePA, enrichplot, AnnotationDbi,
               org.Hs.eg.db, patchwork)

# === 数据准备 ===
# 读取数据
BSIonly3platforms <- readRDS("original_data/BSIonly3platforms.Rdata")
BSIpatientdata_280 <- readRDS("original_data/BSIpatientdata_280.Rdata")
data <- read.csv("original_data/final_data.csv")

# merge transcriptome with clinical
BSIonly3platforms$platformID <- rownames(BSIonly3platforms)
BSIpatientdata_280_2 <- BSIpatientdata_280[, c("platformID", "MARSID")]
merged_BSI <- merge(BSIonly3platforms, BSIpatientdata_280_2, by = "platformID", all.x = TRUE)
merged_BSI$Case <- NULL

# 创建Syndecan三分组
data$Syndecan_group <- cut(data$Syndecan1,
                           breaks = quantile(data$Syndecan1, probs = seq(0, 1, length = 4), na.rm = TRUE),
                           include.lowest = TRUE,
                           labels = c(1, 2, 3))

data2 <- data[, c("ICU_ID_from_datasource", "Syndecan_group", "APACHE_IV_Score", "SOFAtot_new")]
data2 <- data2[!is.na(data2$SOFAtot_new), ]

BSI_transcriptome <- merge(data2, merged_BSI, by.x = "ICU_ID_from_datasource", by.y = "MARSID", all = FALSE)
BSI_transcriptome_2 <- BSI_transcriptome %>% filter(Syndecan_group != 2)
BSI_transcriptome_2$Group <- factor(BSI_transcriptome_2$Syndecan_group, levels = c(1, 3))

# 建模
design <- model.matrix(~ 0 + Group, data = BSI_transcriptome_2)
small <- BSI_transcriptome_2[, 4:ncol(BSI_transcriptome_2)] %>% mutate(across(everything(), as.numeric)) %>% as.matrix()
fit <- lmFit(t(small), design)
cont.matrix <- makeContrasts(logFC = Group3 - Group1, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit3 <- eBayes(fit2)
results_limma <- topTable(fit3, number = Inf, adjust.method = "BH")

# GSEA 准备
gene_list <- results_limma$t
names(gene_list) <- rownames(results_limma)
gene_list <- sort(gene_list, decreasing = TRUE)

gene_df <- data.frame(gene_symbol = names(gene_list), logFC = gene_list)
gene_df$entrez <- mapIds(org.Hs.eg.db, keys = rownames(gene_df), keytype = "SYMBOL", column = "ENTREZID", multiVals = "first")

gene_list_entrez <- gene_df$logFC
names(gene_list_entrez) <- gene_df$entrez
gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)

# 运行 Reactome GSEA
set.seed(11)
y <- gsePathway(geneList = gene_list_entrez,
                pvalueCutoff = 1,
                pAdjustMethod = "BH",
                verbose = FALSE,
                maxGSSize = 10000,
                minGSSize = 0,
                seed = TRUE,
                eps = 0,
                nPermSimple = 10000)

pathway_results <- as.data.frame(y)

# 读取 parent-child 关系
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

# 合并富集结果
df <- merge(pathway_results, all_paths, by = "ID")
df$sig_label <- ifelse(df$p.adjust < 0.05, "FDR < 0.05", "Not significant")

# 保留母通路 + 每组 Top10
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


df <- df %>%
  mutate(
    is_mother = as.numeric(as.logical(is_mother)),  # 确保是 numeric 类型
    Description_fmt = ifelse(is_mother == 1, paste0("**", Description, "**"), Description)
  ) %>%
  group_by(parent_cat) %>%
  arrange(dplyr::desc(is_mother), dplyr::desc(abs(NES))) %>%
  mutate(Description_factor = factor(Description_fmt, levels = rev(unique(Description_fmt)))) %>%
  ungroup()


# GSEA dot plot
font_size <- 14
axis_text_size <- 12
axis_font_size <- 14

p1 <- ggplot(df, aes(x = NES, y = Description_factor, fill = NES)) +
  geom_point(aes(size = sig_label, shape = sig_label), color = "black") +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", name = "NES") +
  scale_size_manual(name = "Significance", values = c("FDR < 0.05" = 4, "Not significant" = 2)) +
  scale_shape_manual(name = "Significance", values = c("FDR < 0.05" = 21, "Not significant" = 23)) +
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
ggsave("figures/GSEA_horizontal_groups.svg", p1, width = 9, height = 6, dpi = 300)









# ========== Volcano Plot ==========
# 
results_limma <- results_limma %>%
  filter(!rownames(.) %in% c("Group", "SOFAtot_new"))

results_limma <- results_limma %>%
  mutate(
    log10_padj = -log10(adj.P.Val),
    estimate_jitter = logFC + runif(n(), min = -0.01, max = 0.01),
    Significance = case_when(
      adj.P.Val < 0.05 & logFC > 0 ~ "Higher in high Syndecan",
      adj.P.Val < 0.05 & logFC < 0 ~ "Lower in high Syndecan",
      TRUE ~ "Not Significant"
    ),
    gene_name = rownames(.)
  )

top_genes <- bind_rows(
  results_limma %>% filter(Significance == "Higher in high Syndecan") %>% arrange(adj.P.Val) %>% head(5),
  results_limma %>% filter(Significance == "Lower in high Syndecan") %>% arrange(adj.P.Val) %>% head(5)
)

p2 <- ggplot(results_limma, aes(x = estimate_jitter, y = log10_padj, color = Significance)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("Higher in high Syndecan" = "red", "Lower in high Syndecan" = "blue", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-0.4, 0.4), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(data = top_genes, aes(label = gene_name), size = 3, color = "black", box.padding = 0.4, max.overlaps = 10) +
  theme_minimal(base_size = font_size) +
  theme(
    text = element_text(family = "sans", size = font_size),
    axis.title = element_text(size = axis_font_size),
    axis.text = element_text(size = axis_text_size),
    axis.title.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  ) +
  labs(y = "-log10 Adjusted p-value")
p2

# 合并
combined_plot <- p2 + p1 +
  plot_layout(ncol = 2, widths = c(1.2, 1)) +
  plot_annotation(tag_levels = 'A',
                  theme = theme(plot.tag = element_text(size = 16, face = "bold", family = "sans")))

# 保存图
ggsave("figures/figureS_combined_GSEA_Volcano_GroupComparison.svg", combined_plot, width = 14, height = 6, dpi = 300)

# 展示
combined_plot




#heatmap####

# 加载包
# 加载包
if (!require("pacman")) install.packages("pacman")
pacman::p_load(biomaRt, dplyr, pheatmap, grid, ggplot2)

# 创建保存路径
dir.create("figures", showWarnings = FALSE)

# 连接 Ensembl
hsmart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# 获取 Reactome 路径中 top5 和 bottom5 基因函数
get_top_bottom_genes <- function(results, pathway_id) {
  pathway_genes <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "reactome",
    values = pathway_id,
    mart = hsmart
  )
  results$gene_name <- rownames(results)
  pathway_results <- results %>% filter(gene_name %in% pathway_genes$hgnc_symbol)
  top_genes <- pathway_results %>% filter(logFC > 0) %>% arrange(adj.P.Val) %>% head(5)
  bottom_genes <- pathway_results %>% filter(logFC < 0) %>% arrange(adj.P.Val) %>% head(5)
  return(unique(c(top_genes$gene_name, bottom_genes$gene_name)))
}

# 热图绘制函数（按Group排序样本）
generate_heatmap <- function(gene_list, file_name) {
  genes_available <- intersect(gene_list, colnames(BSI_transcriptome_2))
  if (length(genes_available) == 0) {
    warning("No valid genes found in expression matrix.")
    return(NULL)
  }
  
  # 提取表达与Group信息
  expr <- BSI_transcriptome_2[, genes_available, drop = FALSE]
  rownames(expr) <- BSI_transcriptome_2$ICU_ID_from_datasource
  group_info <- BSI_transcriptome_2[, c("ICU_ID_from_datasource", "Group")]
  rownames(group_info) <- group_info$ICU_ID_from_datasource
  
  # 按Group排序样本
  group_info_sorted <- group_info[order(group_info$Group), , drop = FALSE]
  expr <- expr[rownames(group_info_sorted), , drop = FALSE]
  annotation <- data.frame(Group = group_info_sorted$Group)
  rownames(annotation) <- rownames(group_info_sorted)
  
  # Z-score标准化
  expr_scaled <- scale(expr)
  expr_matrix <- t(as.matrix(expr_scaled))
  
  # 颜色参数
  ann_colors <- list(Group = c("1" = "#17becf", "3" = "#ff7f0e"))
  colors <- colorRampPalette(c("blue", "white", "red"))(200)
  breaks <- c(seq(-2, 0, length.out = 100), seq(0, 3, length.out = 100)[-1])
  
  # 绘图
  p <- pheatmap(
    expr_matrix,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    annotation_col = annotation,
    annotation_colors = ann_colors,
    color = colors,
    breaks = breaks,
    show_colnames = FALSE,
    show_rownames = TRUE,
    fontsize = 10,                     # 行名字体大小
    fontsize_row = 10,                # 行名（基因）字体
    fontsize_col = 10,                # 样本名（不显示也保留大小）
    fontfamily = "sans"               # 字体统一为 sans
  )
  
  # 保存为 PDF
  pdf_filename <- file.path("figures", file_name)
  pdf(pdf_filename, width = 8, height = 4, family = "sans")
  grid::grid.draw(p$gtable)
  dev.off()
  
  # 显示热图
  grid.newpage()
  grid.draw(p$gtable)
}


# 示例使用
platelet_degranulation_pathway_id <- "R-HSA-76005" #this is Response to elevated platelet cytosolic Ca2+
genes_platelet <- get_top_bottom_genes(results_limma, platelet_degranulation_pathway_id)
generate_heatmap(genes_platelet, "heatmap_platelet.pdf")


##three plates pathway genes analysis
path_ids <- c(
  Kinesins                      = "R-HSA-983189",
  Platelet_degranulation       = "R-HSA-114608",
  Response_to_platelet_Ca2plus = "R-HSA-76005"
)

get_pathway_symbols <- function(pathway_id) {
  getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters    = "reactome",
    values     = pathway_id,
    mart       = hsmart
  ) %>% filter(hgnc_symbol != "")
}

path_symbols_list <- lapply(unname(path_ids), get_pathway_symbols)
platelet_union_symbols <- unique(unlist(lapply(path_symbols_list, `[[`, "hgnc_symbol")))

results_limma$gene_name <- rownames(results_limma)
in_union <- results_limma %>% filter(gene_name %in% platelet_union_symbols)

up5 <- in_union %>%
  filter(logFC > 0) %>%
  arrange(adj.P.Val, desc(logFC)) %>%
  slice_head(n = 5)

down5 <- in_union %>%
  filter(logFC < 0) %>%
  arrange(adj.P.Val, logFC) %>%      
  slice_head(n = 5)

selected_genes <- unique(c(up5$gene_name, down5$gene_name))

if (length(selected_genes) < 10) {
  remainder <- in_union %>%
    filter(!gene_name %in% selected_genes) %>%
    arrange(adj.P.Val)
  need <- 10 - length(selected_genes)
  selected_genes <- unique(c(selected_genes, head(remainder$gene_name, need)))
}

#
ordered_genes <- unique(c(up5$gene_name, down5$gene_name))
# 
ordered_genes <- intersect(ordered_genes, selected_genes)
# 
if (length(ordered_genes) < length(selected_genes)) {
  ordered_genes <- c(ordered_genes, setdiff(selected_genes, ordered_genes))
}

generate_heatmap(
  gene_list = ordered_genes,
  file_name = "FigS7_platelet_union_Top5Up_Top5Down.pdf"
)
