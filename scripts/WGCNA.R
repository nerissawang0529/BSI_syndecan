#
rm(list = ls())

#ready packages
install.packages("WGCNA", repos = "https://cloud.r-project.org")
BiocManager::install("impute")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("preprocessCore")
BiocManager::install("clusterProfiler")
BiocManager::install("ReactomePA")
BiocManager::install("org.Hs.eg.db")
install.packages("ggplot2")
install.packages("igraph")
install.packages("ggraph")
install.packages("pheatmap")
install.packages("dplyr") 

# Load Libraries
library(WGCNA)
library(dplyr)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(ggplot2)
library(igraph)
library(ggraph)
library(pheatmap)

options(stringsAsFactors = FALSE)

#ready data
BSIonly3platforms <-  readRDS("original_data/BSIonly3platforms.Rdata")
BSIpatientdata_280 <- readRDS("original_data/BSIpatientdata_280.Rdata")

BSIonly3platforms$platformID <- rownames(BSIonly3platforms)
BSIpatientdata_280_2 <- BSIpatientdata_280[,c("platformID","MARSID")]

merged_BSI <- merge(BSIonly3platforms, BSIpatientdata_280_2, by = "platformID", all.x = TRUE)
merged_BSI$Case <- NULL

#data for syndecan group
data <- read.csv("original_data/clinical_marker_scource_pathogen.csv")
data$syndecan_1 <- log(data$Syndecan.1.CD138..29.)
# Split into 4 Quantile-Based Groups
data$Syndecan_group <- cut(data$Syndecan.1.CD138..29.,
                           breaks = quantile(data$`Syndecan.1.CD138..29.`, probs = seq(0, 1, length = 5), na.rm = TRUE),
                           include.lowest = TRUE,
                           labels = c(1, 2, 3, 4))  # Now using 4 groups

data2 <- data[,c("ICU_ID_from_datasource","syndecan_1", "Syndecan_group")]

BSI_transcriptome <- merge(data2, merged_BSI, by.x = "ICU_ID_from_datasource",by.y = "MARSID", all = FALSE)
BSI_transcriptome$platformID <- NULL

# Extract Expression Data
rownames(BSI_transcriptome) <- BSI_transcriptome$ICU_ID_from_datasource
library(dplyr)
datExpr <- BSI_transcriptome %>% select(-c(ICU_ID_from_datasource, syndecan_1, Syndecan_group))
#if the rows can be matched?

# Remove Lowly Expressed Genes (Top 10% Most Variable)
geneVariance <- apply(datExpr, 2, var)
topGenes <- names(sort(geneVariance, decreasing = TRUE))[1:(0.10 * length(geneVariance))]
datExpr <- datExpr[, topGenes]


# Check for missing values and outliers
gsg <- WGCNA:::goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}
# Cluster samples to detect outliers
sampleTree <- hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", xlab="", sub="")


# Construct Scale-Free Network (Signed)
powers <- c(1:20)
library(WGCNA)
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed", verbose = 5)

# Choose Optimal Power
softPower <- 15  # Based on `SFT.R.sq` > 0.9
adjacency <- adjacency(datExpr, power = softPower, type = "signed")

# Convert adjacency matrix into a Topological Overlap Matrix (TOM)
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM

# Module Detection with Dynamic Tree Cut
geneTree <- hclust(as.dist(dissTOM), method = "average")
minModuleSize <- 40  # Adjusted to match paper's methodology
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize)
dynamicColors <- labels2colors(dynamicMods)

# Merge Similar Modules (correlation > 0.75)
mergedColors <- mergeCloseModules(datExpr, dynamicColors, cutHeight = 0.25, verbose = 3)$colors


# Calculate Module Eigengenes & Filter by |KME| > 0.7
MEList <- moduleEigengenes(datExpr, colors = mergedColors)
MEs <- MEList$eigengenes
MEs <- data.frame(lapply(MEs, as.numeric))
kME <- signedKME(datExpr, MEs)
filteredGenes <- colnames(datExpr)[apply(kME, 2, max) >= 0.7]

# ============================================
# Merge `MEs` with `BSI_transcriptome`
# ============================================
rownames(MEs) <-rownames(datExpr) 
MEs$SampleID <- rownames(MEs)
BSI_transcriptome_withMEs <- merge(BSI_transcriptome, MEs, by.x = "ICU_ID_from_datasource", by.y = "SampleID")



# Create module-trait relationships
traitData <- BSI_transcriptome_withMEs[, c("ICU_ID_from_datasource", "Syndecan_group")]
rownames(traitData) <- traitData$ICU_ID_from_datasource
traitData$ICU_ID_from_datasource <- NULL
moduleTraitCor <- cor(MEs, traitData, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# Convert to dataframe for ggplot
moduleTraitCor_df <- as.data.frame(moduleTraitCor)
moduleTraitCor_df$Module <- rownames(moduleTraitCor_df)
moduleTraitCor_df_melted <- reshape2::melt(moduleTraitCor_df, id.vars = "Module")

# Heatmap of module-trait relationships
ggplot(moduleTraitCor_df_melted, aes(x = variable, y = Module, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  labs(title = "Module - Age Decade Relationship", x = "Trait", y = "Module") +
  theme_minimal()
# Bar Plot of Gene Significance
geneSignificance <- apply(datExpr, 2, function(x) abs(cor(x, traitData$Syndecan_group, use = "p")))
moduleGeneSignificance <- tapply(geneSignificance, mergedColors, mean)

geneSig_df <- data.frame(Module = names(moduleGeneSignificance), GeneSignificance = moduleGeneSignificance)
geneSig_df$Module <- factor(geneSig_df$Module, levels = names(moduleGeneSignificance))

ggplot(geneSig_df, aes(x = Module, y = GeneSignificance, fill = Module)) +
  geom_bar(stat = "identity") +
  labs(title = "Gene Significance across Modules", x = "Module", y = "Gene Significance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Hierarchical clustering of eigengenes
dissimilarity <- 1 - cor(MEs)
eigengeneTree <- hclust(as.dist(dissimilarity), method = "average")

# Plot dendrogram
plot(eigengeneTree, main = "Hierarchical Clustering of Module Eigengenes",
     sub = "", xlab = "", ylab = "Height")


# Eigengene adjacency heatmap
pheatmap(dissimilarity, clustering_method = "average", 
         main = "Eigengene Adjacency Heatmap",
         color = colorRampPalette(c("blue", "white", "red"))(50))


library(igraph)
library(ggraph)

# Select module (change color as needed)
names(mergedColors) <- colnames(datExpr)
selectedModule <- "yellow"
moduleGenes <- names(mergedColors)[mergedColors == selectedModule]

geneCorMatrix <- cor(datExpr[, moduleGenes], use = "p")

# Convert to adjacency matrix
adjMatrix <- ifelse(abs(geneCorMatrix) > 0.3, 1, 0)  # Threshold correlation at 0.3

# Create graph object
networkGraph <- graph_from_adjacency_matrix(adjMatrix, mode = "undirected", diag = FALSE)

# Plot network
ggraph(networkGraph, layout = "stress") +  # Use "stress" or "kk" for better performance
  geom_edge_link(aes(alpha = 0.5), show.legend = FALSE) +
  geom_node_point(size = 5, color = "blue") +
  geom_node_text(aes(label = V(networkGraph)$name), repel = TRUE) +
  labs(title = paste(selectedModule, "Module Network")) +
  theme_void()
