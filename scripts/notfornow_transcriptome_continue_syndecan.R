#transcriptome for continus 
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

#log syndecan-1
data$syndecan_1 <- log(data$Syndecan.1.CD138..29.)

data2 <- data[,c("ICU_ID_from_datasource","syndecan_1","APACHE_IV_Score.x","SOFAtot")]
data2 <- data2[!is.na(data2$SOFAtot), ]


BSI_transcriptome <- merge(data2, merged_BSI, by.x = "ICU_ID_from_datasource",by.y = "MARSID", all = FALSE)
BSI_transcriptome$platformID <- NULL
BSI_transcriptome_2 <- BSI_transcriptome

#design <- model.matrix(~ syndecan_1, data = BSI_transcriptome_2)
design <- model.matrix(~ syndecan_1 + APACHE_IV_Score.x + SOFAtot - 1, data = BSI_transcriptome_2)

# Check the design matrix
small <- BSI_transcriptome_2[ , c(4:13231)] %>% as.data.frame()
row.names(small) <- BSI_transcriptome_2$geneID
small$geneID <- NULL

small <- data.frame(lapply(small, function(x) as.numeric(as.character(x))))
small <- as.matrix(small)

# Fit linear model
fit2 <- lmFit(t(small), design)

fit3 <- eBayes(fit2)

# Get differential expression results, 'by performing BH adjusted moderated t-statistics using limma'
fit3 <- topTable(fit3, number = Inf, adjust.method = "BH")
results_limma <- fit3

#
gene_list <- results_limma$F
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

