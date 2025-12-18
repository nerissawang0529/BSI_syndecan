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

#data for syndecan group
data <- read.csv("original_data/clinical_marker_scource_pathogen.csv")

# Cutting the sorted data into three equal parts (1, 2, 3)
data$Syndecan_group <- cut(data$Syndecan.1.CD138..29.,
                           breaks = quantile(data$`Syndecan.1.CD138..29.`, probs = seq(0, 1, length = 4), na.rm = TRUE),
                           include.lowest = TRUE,
                           labels = c(1, 2, 3))  # Now using 4 groups


data2 <- data[,c("ICU_ID_from_datasource","Syndecan_group","APACHE_IV_Score.x","SOFAtot_new")]
data2 <- data2[!is.na(data2$SOFAtot_new), ]

BSI_transcriptome <- merge(data2, merged_BSI, by.x = "ICU_ID_from_datasource",by.y = "MARSID", all = FALSE)
BSI_transcriptome$platformID <- NULL

##pick Group1 and Group3
BSI_transcriptome_2 <- BSI_transcriptome %>% filter(Syndecan_group != 2)

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

#
gene_list <- results_limma$t
names(gene_list) <- row.names(results_limma)
gene_list <- sort(gene_list, decreasing = TRUE)

# Convert gene symbols to Entrez IDs
gene_df <- data.frame(gene_symbol = names(gene_list), logFC = gene_list)

#gene_df <- gene_df[gene_df$gene_symbol != "Group", ]
#gene_df <- gene_df[gene_df$gene_symbol != "APACHE_IV_Score.x",]
#gene_df <- gene_df[gene_df$gene_symbol != "SOFAtot",]


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
set.seed(11)
gene_list_entrez <- gene_df$logFC
names(gene_list_entrez) <- gene_df$entrez
gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)
gene_list_entrez_df <- gene_list_entrez %>% as.data.frame()
#####

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



















### only host response pathways (adaptive, innate, cytokine signaling, programmed cell death, hemostasis)
RPR <- read.table("original_data/reactomepathways.txt", sep="", stringsAsFactors = F)
#if the code use this, then if the followling code
names(RPR) <- c("parent", "child")
## select wanted pathways
#use reactome order####
# List of parent and child pathways

#all of the pathways within Hemostasis and Extracellular Matrix Pathway
pathways_bar <- list(
  #Extracellular Matrix Pathway
  'R-HSA-1474244' = c('R-HSA-3000170','R-HSA-1474290','R-HSA-1566977','R-HSA-1566948','R-HSA-3000157','R-HSA-3000171','R-HSA-3000178','R-HSA-1474228','R-HSA-216083','R-HSA-8941237'),
  #Hemostasis
  'R-HSA-109582'  = c('R-HSA-418346','R-HSA-75892','R-HSA-76002','R-HSA-140877','R-HSA-75205','R-HSA-202733','R-HSA-983231')
)


# Function to filter and order data within each parent pathway with parent at the top
filter_and_order <- function(df, parent, children) {
  df_parent <- df %>% filter(ID == parent)
  df_children <- df %>% filter(ID %in% children) %>%
    arrange(desc(NES))
  df_filtered <- bind_rows(df_parent, df_children) %>%
    mutate(Type = ifelse(ID == parent, 'Parent', 'Child'),
           Group = parent)
  df_filtered
}

# Create a new dataframe for plotting
plot_data <- data.frame()
for (parent in names(pathways_bar)) {
  children <- pathways_bar[[parent]]
  filtered_data <- filter_and_order(pathway_results, parent, children)
  plot_data <- rbind(plot_data, filtered_data)
}


# Define colors
plot_data <- plot_data %>% mutate(
  Color = case_when(
    p.adjust < 0.05 & NES < 0 ~ ifelse(Type == 'Parent', '#1f78b4', '#a6cee3'), # blue, bright blue
    p.adjust < 0.05 & NES > 0 ~ ifelse(Type == 'Parent', '#e31a1c', '#fb9a99'), # red, bright red
    TRUE ~ '#bdbdbd' # grey
  )
)

# Reorder the plot_data
ordered_data <- data.frame()
for (parent in names(pathways_bar)) {
  temp <- plot_data %>% filter(Group == parent | is.na(Group))
  ordered_data <- rbind(ordered_data, temp)
}
# Removing NA Group rows added due to spacer placement

# Ensure the correct order for Description: parents first, then children, within each group
ordered_data <- ordered_data %>%
  arrange(Group, Type == "Child")

# Set Description as a factor and explicitly set the levels in the order of appearance in ordered_data
ordered_data$Description <- factor(ordered_data$Description, levels = unique(ordered_data$Description))


#reorder the description for all Hemostasis，Extracellular Matrix Pathways
ordered_data$Description <- factor(ordered_data$Description, levels = c('Syndecan interactions','Extracellular matrix organization','Elastic fibre formation','ECM proteoglycans','Degradation of the extracellular matrix','Collagen formation','Non-integrin membrane-ECM interactions','Fibronectin matrix formation','Laminin interactions','Integrin cell surface interactions','Invadopodia formation','Hemostasis','Formation of Fibrin Clot (Clotting Cascade)','Platelet activation, signaling and aggregation','Platelet Adhesion to exposed collagen','Platelet homeostasis','Factors involved in megakaryocyte development and platelet production','Dissolution of Fibrin Clot','Cell surface interactions at the vascular wall'))
#reorder the description for directly related with syndecan-1 within Hemostasis，Extracellular Matrix Pathways
#ordered_data$Description <- factor(ordered_data$Description, levels = c('Syndecan interactions', 'Extracellular matrix organization', 'Degradation of the extracellular matrix', 'Fibronectin matrix formation', 'Laminin interactions', 'Hemostasis', 'Formation of Fibrin Clot (Clotting Cascade)', 'Platelet activation, signaling and aggregation', 'Platelet Adhesion to exposed collagen', 'Factors involved in megakaryocyte development and platelet production', 'Cell surface interactions at the vascular wall'))



# Create the plot
p <- ggplot(ordered_data, aes(x = Description, y = NES, fill = Color, shape = Type)) +
  geom_bar(stat = 'identity', position = position_dodge()) + # Ensure consistent bar heights
  scale_fill_identity() +
  coord_flip() +
  scale_x_discrete(limits = rev(levels(ordered_data$Description))) + # Correct the x-axis order
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(),  # Remove background
    plot.background = element_blank(),   # Remove plot background
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14, family = "Arial"),
    axis.text.y = element_text(size = 14, family = "Arial"),  # Adjust font size and family for y-axis labels
    axis.title.x = element_text(size = 14, family = "Arial"), # Adjust font size and family for x-axis title
    axis.title.y = element_text(size = 14, family = "Arial"), # Adjust font size and family for y-axis title
    plot.title = element_text(size = 14, family = "Arial"),   # Adjust font size and family for title
    strip.background = element_blank(),  # Remove strip background
    strip.text = element_blank()         # Remove strip text
  ) +
  labs(title = "Pathway NES Plot",
       x = "Pathway Description",
       y = "NES") +
  geom_hline(yintercept = 0, color = "black") + # Add a black line at y = 0
  scale_shape_manual(values = c('Parent' = 24, 'Child' = 21)) + # Custom shapes
  scale_color_identity()

p

#adjusted the p value####
ordered_data$new_adjusted_pvalue <- p.adjust(ordered_data$pvalue, method = "BH")
ordered_data$qvalue <- NULL
ordered_data$rank  <- NULL
ordered_data$Group  <- NULL
ordered_data$Type <- NULL
ordered_data$core_enrichment <- NULL
ordered_data$leading_edge <- NULL
ordered_data_1 <- ordered_data[,c("ID","Description","pvalue","p.adjust","new_adjusted_pvalue")]

# Create a new dataframe for plotting
plot_data <- data.frame()
for (parent in names(pathways_bar)) {
  children <- pathways_bar[[parent]]
  filtered_data <- filter_and_order(pathway_results, parent, children)
  plot_data <- rbind(plot_data, filtered_data)
}

plot_data <- plot_data %>%
  left_join(ordered_data %>% dplyr::select(ID, new_adjusted_pvalue), by = "ID")

# Define colors
plot_data <- plot_data %>% mutate(
  Color = case_when(
    new_adjusted_pvalue < 0.05 & NES < 0 ~ ifelse(Type == 'Parent', '#1f78b4', '#a6cee3'), # blue, bright blue
    new_adjusted_pvalue < 0.05 & NES > 0 ~ ifelse(Type == 'Parent', '#e31a1c', '#fb9a99'), # red, bright red
    TRUE ~ '#bdbdbd' # grey
  )
)

# Reorder the plot_data
ordered_data <- data.frame()
for (parent in names(pathways_bar)) {
  temp <- plot_data %>% filter(Group == parent | is.na(Group))
  ordered_data <- rbind(ordered_data, temp)
}
# Removing NA Group rows added due to spacer placement

# Ensure the correct order for Description: parents first, then children, within each group
ordered_data <- ordered_data %>%
  arrange(Group, Type == "Child")

# Set Description as a factor and explicitly set the levels in the order of appearance in ordered_data
ordered_data$Description <- factor(ordered_data$Description, levels = unique(ordered_data$Description))

#reorder the description for all Hemostasis，Extracellular Matrix Pathways
ordered_data$Description <- factor(ordered_data$Description, levels = c('Syndecan interactions','Extracellular matrix organization','Elastic fibre formation','ECM proteoglycans','Degradation of the extracellular matrix','Collagen formation','Non-integrin membrane-ECM interactions','Fibronectin matrix formation','Laminin interactions','Integrin cell surface interactions','Invadopodia formation','Hemostasis','Formation of Fibrin Clot (Clotting Cascade)','Platelet activation, signaling and aggregation','Platelet Adhesion to exposed collagen','Platelet homeostasis','Factors involved in megakaryocyte development and platelet production','Dissolution of Fibrin Clot','Cell surface interactions at the vascular wall'))

#reorder the description for directly related with syndecan-1 within Hemostasis，Extracellular Matrix Pathways
#ordered_data$Description <- factor(ordered_data$Description, levels = c('Syndecan interactions', 'Extracellular matrix organization', 'Degradation of the extracellular matrix', 'Fibronectin matrix formation', 'Laminin interactions', 'Hemostasis', 'Formation of Fibrin Clot (Clotting Cascade)', 'Platelet activation, signaling and aggregation', 'Platelet Adhesion to exposed collagen', 'Factors involved in megakaryocyte development and platelet production', 'Cell surface interactions at the vascular wall'))


# Create the plot
p <- ggplot(ordered_data, aes(x = Description, y = NES, fill = Color)) +
  geom_bar(stat = 'identity', position = position_dodge(), width = 0.8) + # Ensure consistent bar heights and adjust bar width
  scale_fill_identity() +
  coord_flip() +
  scale_x_discrete(limits = rev) +
  scale_y_continuous(expand = c(0, 0), breaks = c(-0.5, 0, 0.5, 1,1.5,2)) + # Set y-axis breaks to -1, 0, 1, 2
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(),  # Remove background
    plot.background = element_blank(),   # Remove plot background
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14, family = "Arial"),
    axis.text.y = element_text(size = 14, family = "Arial"),  # Adjust font size and family for y-axis labels
    axis.title.x = element_text(size = 14, family = "Arial"), # Adjust font size and family for x-axis title
    axis.title.y = element_text(size = 14, family = "Arial"), # Adjust font size and family for y-axis title
    plot.title = element_text(size = 14, family = "Arial"),   # Adjust font size and family for title
    strip.background = element_blank(),  # Remove strip background
    strip.text = element_blank()         # Remove strip text
  ) +
  labs(title = "Pathway NES Plot",
       x = "Pathway Description",
       y = "NES") +
  geom_hline(yintercept = 0, color = "black") # Add a black line at y = 0

p

