## Gene expression analyses on raw U219 data 
## Start script: 14-03-2025 


# raw U219 data (already log2 converted)
# selects only admission samples 
# includes comparison with healthies 
# no coconut 


###### loading packages and data ###############################################

#### import basic packages
#pacman::p_load(pacman, dplyr, GGally, ggplot2, ggthemes, ggpubr, httr, lubridate, plotly, rio, rmarkdown, stringr, tidyr, tableone, tibble, plyr, dunn.test, lattice, MatchIt, survival, survminer, cobalt)
pacman::p_load(pacman, dplyr, ggplot2, httr, rio, tidyverse, tableone, Biobase)
dselect <- dplyr::select

##load raw U216 data 
load("L:/basic/divg/CEMM-Infectious disease/Susanne Doeleman/1_projecten/1_aspiratie_mars/1_data/U219_MARS/MARS_U219_nofilt.eset")

## load ascap data
ascap <- import("L:/basic/divg/CEMM-Infectious disease/Susanne Doeleman/1_projecten/1_aspiratie_mars/1_data/ASCAP_final_20250304.xlsx")

###### prep expression data ###############################################
exData <- exprs(mars_pax_nofilt)
dim(exData)

## Filter out patients that were sampled on day of ICU acquired complication, 
## only select those registered as "opname"
ALL_PAX_diamond <- read.csv("L:/basic/divg/CEMM-Infectious disease/Susanne Doeleman/1_projecten/1_aspiratie_mars/1_data/U219_MARS/platinum_pax_guide.csv",header=T, na.strings="",stringsAsFactors = F)
colnames(exData) <- gsub(",", ".", colnames(exData))
table(colnames(exData) %in% ALL_PAX_diamond$chipID)
head(colnames(exData)[!colnames(exData) %in% ALL_PAX_diamond$chipID ])
u2pdat <- ALL_PAX_diamond[ALL_PAX_diamond$chipID %in% colnames(exData), ]
table(u2pdat$Event)

### select only Admission "Opname"
u2pdat <- u2pdat[u2pdat$Event=="Opname",]
table(ascap$ID %in% u2pdat$MARSID)

ascap_u2 <- merge(ascap, u2pdat, by.x="ID", by.y="MARSID")
table(ascap_u2$group)
n_distinct(ascap_u2$ID)

# Remove one duplicate row (paxID 2803200 for patient 10413)
ascap_u2 <- ascap_u2[ascap_u2$PAXID!=2803200,]

### select expression data for healthies and CAP groups
exData_hvcap <- exData[, colnames(exData) %in% ascap_u2$chipID| grepl("hv", colnames(exData)) ]
dim(exData_hvcap ) # correct, 85 AP/CAP and 42 healthies 

## save as expression 
expression <- exData_hvcap
tail(colnames(expression))

dim(expression)

###### prep p data #############################################################
pData <- phenoData(mars_pax_nofilt)@data
table(pData$health) 
pData_hvcap <- pData[rownames(pData) %in% ascap_u2$chipID| grepl("hv", rownames(pData)),]

## add group
x <- dselect(ascap_u2, ID, group)
pData_hvcap_g <- merge(pData_hvcap, x, by.x="MARSID", by.y="ID", all.x = TRUE)
table((pData_hvcap_g$group))
pData_hvcap_g[is.na(pData_hvcap_g$group), "group"] <- "HV"
pData_hvcap_g <- dselect(pData_hvcap_g,chipID, MARSID, group, age, Patient_gender)

# save as pData 
pData <- pData_hvcap_g
rownames(pData) <- pData$chipID

rm(pData_hvcap_g)
rm(pData_hvcap)


pData <- pData[colnames(expression),]
pData$Patient_gender <- factor(pData$Patient_gender, levels= c("female", "male"))
str(pData$Patient_gender)

## Transfer from probes to genes
# BiocManager::install("GEOquery")
library(GEOquery)
gpl13667 <- getGEO('GPL13667',dest="L:/basic/divg/CEMM-Infectious disease/Susanne Doeleman/1_projecten/1_aspiratie_mars/1_data/U219_MARS")

f  <- as.data.frame(rownames(expression),stringsAsFactors=F)
row.names(f) <- f[,1]
fdt <- Table(gpl13667)[,c(1,15,21,26)]
row.names(fdt) <- fdt[,1]
fdt2 <- merge(f, fdt, by="row.names")

row.names(fdt2) <- fdt2[,2]
fdt2$Row.names <- NULL
fdt2$ID <- NULL
featureData <- new("AnnotatedDataFrame", data= fdt2) # sort of eset 

dim(featureData)
head(fdt2)


###### limma // adjusted for age and sex #######################################

eset <- ExpressionSet(assayData=expression,
                      phenoData = AnnotatedDataFrame(pData),
                      featureData= AnnotatedDataFrame(fdt2))

# BiocManager::install("limma")
library(limma)
str(pData)

# main model => linear model for each gene / probe accounting for age/sex 
design_cm <- model.matrix(~0 + group + age + Patient_gender, data=pData(eset)) # ~0 removes the intercept
head(design_cm) 

# check: right number of patients
colSums(design_cm) 
table(pData(eset)[,"group"]) 

# make contrasts matrix to define which groups to compare
cm <- makeContrasts(status = groupAP - groupHV, levels = design_cm)
cm 

# fit linear model with design matrix 
fit_cm <- lmFit(eset, design_cm)
head(fit_cm)

# apply contrast matrix to compare between groups AP vs health
fit_cm2 <- contrasts.fit(fit_cm, contrasts=cm)
head(fit_cm2$coefficients, 3)

# adding bayesian moderation (key feature of limma, moderate variance estimates, stabilizing them)
fit_cm2 <- eBayes(fit_cm2) 
results <- decideTests(fit_cm2)
summary(results)

# create top table and order by probe name (then same order as fit_cm2)
ttAP <- topTable(fit_cm2, n=Inf) 
ttAP <- ttAP[order(ttAP$rownames.expression.),] # order 
table(ttAP$rownames.expression. == fit_cm2$genes$`rownames(expression)`)  

# calculate cohens d and add to top table 
fit_cm2$cohensd <- fit_cm2$coefficients / sqrt(fit_cm2$s2.post)
ttAP$cohensD <- fit_cm2$cohensd 

# calculate hedges g and add to toptable 
n_ap <- sum(pData$group == "AP")
n_cap <- sum(pData$group == "CAP")
n_hv <- sum(pData$group == "HV")
ttAP$hedgesG <- ttAP$cohensD * (1 - (3/(4*(n_ap+n_hv)-9)))
plot(ttAP$cohensD, ttAP$hedgesG)

# select those significant based on p value 
ttAP$sig <- ttAP$adj.P.Val<0.05
table(ttAP$sig)

# select those significant based on hedges g 
ttAP$sig_hg <- abs(ttAP$hedgesG) > 0.8
table(ttAP$sig_hg)
table(ttAP$sig)
 
# checks 
#plot(ttAP$logFC, ttAP$cohensD)
#plot(ttAP$logFC, -log10(ttAP$adj.P.Val))
#plot(ttAP$cohensD, -log10(ttAP$adj.P.Val))
#plot(ttAP$t , -log10(ttAP$adj.P.Val))
head(ttAP[order(-ttAP$logFC),])
head(fit_cm2$coefficients[order(-fit_cm2$coefficients),])

## interpretation log2 FC => 11730765_at has log2FC 6.203975, meaning mean is 2^6.203975 = 73 times higher in AP than healthy





## CAP vs HV ################################################################### 

# main model => linear model for each gene / probe accounting for age/sex 
design_cm <- model.matrix(~0 + group + age + Patient_gender, data=pData(eset)) # ~0 removes the intercept
head(design_cm) 

# make contrasts matrix to define which groups to compare
cm <- makeContrasts(status = groupCAP - groupHV, levels = design_cm)

# fit linear model with design matrix 
fit_cm <- lmFit(eset, design_cm)
head(fit_cm)

# apply contrast matrix to compare between groups AP vs health
fit_cm2 <- contrasts.fit(fit_cm, contrasts=cm)
head(fit_cm2$coefficients, 3)

# adding bayesian moderation (key feature of limma, moderate variance estimates, stabilizing them)
fit_cm2 <- eBayes(fit_cm2) 
results <- decideTests(fit_cm2)
summary(results)

# create top table and order by probe name (then same order as fit_cm2)
ttCAP <- topTable(fit_cm2, n=Inf) 
ttCAP <- ttCAP[order(ttCAP$rownames.expression.),] # order 
table(ttCAP$rownames.expression. == fit_cm2$genes$`rownames(expression)`)  

# calculate cohens d and add to top table 
fit_cm2$cohensd <- fit_cm2$coefficients / sqrt(fit_cm2$s2.post)
ttCAP$cohensD <- fit_cm2$cohensd 

# calculate hedges g and add to toptable 
n_ap <- sum(pData$group == "AP")
n_cap <- sum(pData$group == "CAP")
n_hv <- sum(pData$group == "HV")
ttCAP$hedgesG <- ttCAP$cohensD * (1 - (3/(4*(n_cap+n_hv)-9)))
plot(ttCAP$cohensD, ttCAP$hedgesG)

# select those significant based on p value 
ttCAP$sig <- ttCAP$adj.P.Val<0.05
table(ttCAP$sig)

# select those significant based on hedges g 
ttCAP$sig_hg <- abs(ttCAP$hedgesG) > 0.8
table(ttCAP$sig_hg)
table(ttCAP$sig)

# checks 
#plot(ttCAP$logFC, ttCAP$cohensD)
#plot(ttCAP$logFC, -log10(ttCAP$adj.P.Val))
#plot(ttCAP$cohensD, -log10(ttCAP$adj.P.Val))
#plot(ttCAP$t , -log10(ttCAP$adj.P.Val))
head(ttCAP[order(-ttCAP$logFC),])
head(fit_cm2$coefficients[order(-fit_cm2$coefficients),])

### merge AP and CAP toptables

ttm <- merge(ttCAP, ttAP, by=0)
table(ttm$sig.x, ttm$sig.y)
colnames(ttm)[colnames(ttm) == "sig_hg.x"] <- "CAP_significant"
colnames(ttm)[colnames(ttm) == "sig_hg.y"] <- "AP_significant"
ttm$CAP_significant <- as.logical(ttm$CAP_significant)
ttm$AP_significant <- as.logical(ttm$AP_significant)

table(ttm$CAP_significant, ttm$AP_significant)

# plot VENN diagram  
library(eulerr)
mat <- cbind(ttm$CAP_significant, ttm$AP_significant)
colnames(mat) <- c("non-AP","AP")
fit2 <- euler(mat)

fill_colors <- c("grey95",                 # Very light grey for CAP
                 rgb(112/255, 137/255, 185/255, alpha = 0.25))  # Light blue for AP
(fill_colors2 <- c("#EBE8E0",                # Light grey for CAP
                  rgb(255/255, 37/255, 17/255)))  # Red with alpha 0.25 for AP
plot(fit2, 
     fills = fill_colors2,
     main = "Probes differentially expressed versus healthy volunteers (Hedges' g > 0.8)",    
     quantities = TRUE,
     edges = "transparent",
     labels = list(
       col = "black",  # Set label color to dark grey
       font = 2,        # Bold font
       cex = 1.2        # Font size
     ),
     main_args = list(cex = 0.6)
)

fill_colors2 <- c("#FADADD",                        # Light pink (Economist-style) for CAP
                    rgb(255/255, 37/255, 17/255))  # Soft red for AP

# Corresponding outlines (darker versions)
outline_colors <- c(
  "#D36C6C",   # Darker pink for CAP
  "#990000"    # Dark red for AP
)

plot(fit2, 
     fills = fill_colors2,
     borders = outline_colors,
     main = "Probes differentially expressed versus healthy volunteers (Hedges' g > 0.8)",    
     quantities = TRUE,
     labels = list(
       col = "black",  
       font = 2,
       cex = 1.2
     ),
     main_args = list(cex = 0.6)
)



## volcano of AP vs CAP genes ##################################################

# redo contrast matrix for AP vs CAP directly 

# main model => linear model for each gene / probe accounting for age/sex 
design_cm <- model.matrix(~0 + group + age + Patient_gender, data=pData(eset)) # ~0 removes the intercept
head(design_cm) 

# make contrasts matrix to define which groups to compare
cm <- makeContrasts(status = groupAP - groupCAP, levels = design_cm)

# fit linear model with design matrix 
fit_cm <- lmFit(eset, design_cm)
head(fit_cm)

# apply contrast matrix to compare between groups AP vs health
fit_cm2 <- contrasts.fit(fit_cm, contrasts=cm)
head(fit_cm2$coefficients, 3)

# adding bayesian moderation (key feature of limma, moderate variance estimates, stabilizing them)
fit_cm2 <- eBayes(fit_cm2) 
results <- decideTests(fit_cm2)
summary(results)

# create top table and order by probe name (then same order as fit_cm2)
ttboth <- topTable(fit_cm2, n=Inf) 
ttboth <- ttboth[order(ttboth$rownames.expression.),] # order 
table(ttboth$rownames.expression. == fit_cm2$genes$`rownames(expression)`)  

# add colors for volcano 
ttboth$col <- "A" 
ttboth$col <- ifelse(ttboth$logFC < 0 & ttboth$adj.P.Val<0.05, "B", ttboth$col) # lower 
ttboth$col <- ifelse(ttboth$logFC>0 & ttboth$adj.P.Val<0.05, "C", ttboth$col) # higher 

# add row for labeling  
ttboth$label2 <- NA 
ttboth$label2[ttboth$adj.P.Val<0.05] <- ttboth$Gene.Symbol[ttboth$adj.P.Val<0.05]

# ggplot volcano
library(ggplot2)
#install.packages("ggrepel")
library(ggrepel)

table(ttboth$col)

volcano_rna <- ggplot(ttboth, aes(x=logFC, y= -log10(adj.P.Val), label=label2)) +
  geom_point(aes(col=col), size=1) +
  scale_color_manual(values = c("A" = "#EBE8E0", "B" =	"#2F2440", "C" = 	"#FF2511")) +
  scale_shape_manual(values=c(16,16,20)) + 
  theme_minimal() + 
  geom_text_repel(max.overlaps = 33,
                  #nudge_x = 0.05,   # adjust these values as needed
                  #nudge_y = 0.1,    # adjust these values as needed
                  point.padding = unit(0.4, "lines"),   # increase space around points
                  box.padding = unit(0.8, "lines"),    # increase the repelling area
                  segment.color = 'grey50',            # adjust line color
                  segment.size = 0.2) + 
  theme(legend.position = "none",
        plot.title = element_text(face="bold", size = 15),
        axis.title.y = element_text(face="bold", size=15),
        axis.title.x = element_text(size=15),
        axis.text.x = element_text(size=12, color = "grey50"),
        axis.text.y = element_text(size=12, color = "grey50"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = 'grey50', fill=NA, size=1),
        plot.background = element_rect(fill = "white", color = NA),  # Keep the overall plot background white
        panel.background = element_rect(fill = "gray99", color = NA))  +
  coord_cartesian(xlim = c(-3, 3), ylim=c(0, -log10(min(ttboth$adj.P.Val))+1 )) + 
  xlab("Log2(FC)") + 
  ylab(expression('-log'[10]*'(BH adjusted P)')) + 
  ggtitle("Gene expression in AP versus non-AP") + 
  geom_hline(yintercept=-log10(0.05), linetype='dotted', col = "#2F2440") +
  geom_vline(xintercept=0, linetype='dotted', col = "#2F2440") +
  annotate("text", x = -1.75, y = 4.5, label = "Number of downregulated \n probes in AP = 220", size = 4, col = "#2F2440") +
  annotate("text", x = 1.75, y = 4.5, label = "Number of upregulated \n probes in AP = 531", size = 4, col = "#FF2511")

print(volcano_rna)


################################################################################
## make volcano on gene level - select the most significant probes only 
gene_level <- ttboth

# select most significant probe per gene (pathway only runs on single entrez gene)
gene_level_select <- gene_level %>% arrange(Entrez.Gene, P.Value)
gene_level_select <- filter(gene_level_select, gene_level_select$Entrez.Gene != "---")
gene_level_select <- filter(gene_level_select, !duplicated(gene_level_select$Entrez.Gene))

table(gene_level_select$col)

# Remove all dashed gene symbols 
gene_level_select$Gene.Symbol2 <- sub(" ///.*", "", gene_level_select$Gene.Symbol)

# visually inspect changed rows => 993 rows 
x <- select(gene_level_select, Gene.Symbol, Gene.Symbol2, P.Value)
x <- filter(x, Gene.Symbol != Gene.Symbol2)

# remove duplicate values, removed ~ 500 
gene_level_select <- gene_level_select %>% arrange(Gene.Symbol2, P.Value)
gene_level_select <- filter(gene_level_select, !duplicated(gene_level_select$Gene.Symbol2))

table(gene_level_select$col)


# rename label and make italic
gene_level_select$label2[gene_level_select$adj.P.Val<0.05] <- gene_level_select$Gene.Symbol2[gene_level_select$adj.P.Val<0.05]
gene_level_select$label2_expr <- ifelse(is.na(gene_level_select$label2),NA,paste0("italic('", gene_level_select$label2, "')"))

# install aptos as a font
install.packages("showtext")
library(showtext)
font_add(family = "Aptos", regular = "Aptos.ttf")  # Or "Aptos" if already installed
showtext_auto()

# final volcano 
volcano_rna <- ggplot(gene_level_select, aes(x=logFC, y= -log10(adj.P.Val), label=label2_expr)) +
  geom_point(aes(col=col), size=1) +
  scale_color_manual(values = c("A" = "#EBE8E0", "B" =	"#2F2440", "C" = 	"#FF2511")) +
  scale_shape_manual(values=c(16,16,20)) + 
  theme_minimal() + 
  geom_label_repel(
    parse = TRUE,
    max.overlaps = 19,
    point.padding = unit(0.4, "lines"),
    box.padding = unit(0.8, "lines"),
    segment.color = 'grey50',
    segment.size = 0.2,
    size = 3.5,  # size of text
    label.size = 0.25,  # thickness of label border
    label.r = unit(0.15, "lines"),  # corner radius of the box
    label.padding = unit(0.2, "lines"),  # space between text and box edge
    fill = "white",  # background color of the label box
    color = 'grey30'  # text and border color
  ) + 
  theme(legend.position = "none",
        text = element_text(family = "Aptos"),
        plot.title = element_text(face = "bold", size = 15, family = "Aptos"),
        axis.title.y = element_text(face="bold", size=15, family = "Aptos"),
        axis.title.x = element_text(size=15, family = "Aptos"),
        axis.text.x = element_text(size=12, color = "grey50",family = "Aptos"),
        axis.text.y = element_text(size=12, color = "grey50",family = "Aptos"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = 'grey50', fill=NA, size=1),
        plot.background = element_rect(fill = "white", color = NA),  # Keep the overall plot background white
        panel.background = element_rect(fill = "gray99", color = NA))  +
  coord_cartesian(xlim = c(-3, 3), ylim=c(0, -log10(min(ttboth$adj.P.Val))+1 )) + 
  xlab("Log2(FC)") + 
  ylab(expression('-log'[10]*'(BH adjusted P)')) + 
  ggtitle("Gene expression in AP versus non-AP (gene level)") + 
  geom_hline(yintercept=-log10(0.05), linetype='dotted', col = "#2F2440") +
  geom_vline(xintercept=0, linetype='dotted', col = "#2F2440") +
  annotate("text", x = -1.75, y = 4.5, label = "Downregulated genes in AP \n 163 / 19469 (0.8%)", size = 4, col = "#2F2440") +
  annotate("text", x = 1.75, y = 4.5, label = "Upregulated genes in AP \n 379 / 19469 (2.0%)", size = 4, col = "#FF2511")

print(volcano_rna)



######### pathway analysis of V2 genes #############################################
 
#install.packages("BiocManager")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")
#BiocManager::install("org.Hs.eg.db")

# install.packages("yulab.utils")
library("clusterProfiler")
library(ggplot2)
library(enrichplot)
library(org.Hs.eg.db)
library(dplyr)
library(ReactomePA) # for gseapathway

# get log2fold change or t statistic 
AP_pathway <- ttboth

# prep df / select most significant probe per entrez (pathway only runs on single entrez gene)
AP_pathway_sorted <- AP_pathway %>% arrange(Entrez.Gene, P.Value)
AP_pathway_sorted$Entrez.Gene <- sub(" ///.*", "", AP_pathway_sorted$Entrez.Gene)
AP_pathway_sorted <- filter(AP_pathway_sorted, AP_pathway_sorted$Entrez.Gene != "---")
AP_pathway_sorted <- filter(AP_pathway_sorted, !duplicated(AP_pathway_sorted$Entrez.Gene))

str(AP_pathway_sorted)

# name the vector according to entrez gene name 
AP_pathway <- AP_pathway_sorted %>% dplyr::arrange(desc(t))
gene_list <- AP_pathway$t
colnames(AP_pathway)
names(gene_list) <- AP_pathway$Entrez.Gene
head(gene_list)
 
# select parent pathways 
parents <- import("https://reactome.org/download/current/ReactomePathwaysRelation.txt", header = FALSE)
colnames(parents) <- c("parent","child") 

# get gene ID for adaptive immunity, innate immunity, hemostasis, cytokine signalling, programmed cell death and ECM
relevant_pathways <- parents %>% filter(parent %in% c("R-HSA-1280218", # adaptive
                                                      "R-HSA-168249", # innate
                                                      "R-HSA-1280215", # cytokine
                                                      "R-HSA-109582",  # hemostasis
                                                      "R-HSA-5357801")) # programmed cell death
#"R-HSA-1474244"))  # extracellular matrix organisation

# get all values (p value cut off at 1)
enrichment <- gsePathway(gene_list, 
                         pvalueCutoff = 1,
                         pAdjustMethod = "BH", 
                         verbose = FALSE)

enrichment_df <- as.data.frame(enrichment)
head(enrichment_df)

# make a new enrichment object "select_enrich" with only relevant_pathways
select_enrich <- enrichment
select_enrich@result <- enrichment[enrichment@result[["ID"]] %in% relevant_pathways$child]
head(select_enrich)

viewPathway("E2F mediated regulation of DNA replication", readable = TRUE, 
            foldChange = gene_list)

## NES: normalized enrichment score / higher 
# display top 15 enriched (child) pathways of all relevant (parent) pathways in tree plot 
tree_relevant <- pairwise_termsim(select_enrich)
pr <- treeplot(tree_relevant, 
               color = "NES", 
               showCategory = 10,
               cluster.params = list(label_words_n = "ward.D2")) + 
               scale_color_gradient(low = "#7089B9", high = 	"#F6423C") + 
               labs(color="NES") # change legend title to "NES"

pr


# for extra information on treeplot: 
?treeplot




set.seed(106678)
fgy <- gsePathway(gene_list, 
                  pvalueCutoff = 1,
                  pAdjustMethod = "BH", 
                  nPerm=10000,
                  minGSSize = 5, # minimum gene set size 
                  maxGSSize = 350, # to prevent huge sets like "immune system", reduce number of tests
                  verbose = FALSE,
                  seed=T)
head(fgy)


fgy2 <- as.data.frame(fgy)




library(reactome.db)
library(org.Hs.eg.db)


selected_pathway_ids <- c("R-HSA-1280218", # adaptive
                          "R-HSA-168249", # innate
                          "R-HSA-1280215", # cytokine
                          "R-HSA-109582") # Hemostasis
                          #"R-HSA-1430728",  # Metabolism
                          #"R-HSA-5357801") # Programmed Cell Death

# Get gene IDs for selected Reactome pathway IDs
# xx <- as.list(reactomePATHID2EXTID )
# 
# xxss <- lapply(selected_pathway_ids, function(id) if (id %in% names(xx)) xx[[id]] else NULL)
# 
# ### select only immune pathways
# setwd("J:/My Drive/Pathogen")

RPR <- import("https://reactome.org/download/current/ReactomePathwaysRelation.txt", header = FALSE)
# mydatHSA <- RPR[ grepl("HSA", RPR$V1)  ,]
# 
# test <- unique(RPR$V1[ !RPR$V1 %in% RPR$V2 ])
# parents <- test[grepl("HSA", test)]
# length(unique(parents))

gen0 <- data.frame(ID=selected_pathway_ids, gen=0)
gen1 <- data.frame(ID=subset(RPR, V1 %in% gen0$ID)$V2, gen=1)
gen2 <- data.frame(ID=subset(RPR, V1 %in% gen1$ID)$V2, gen=2)
gen3 <- data.frame(ID=subset(RPR, V1 %in% gen2$ID)$V2, gen=3)
gen4 <- data.frame(ID=subset(RPR, V1 %in% gen3$ID)$V2, gen=4)
gen5 <- data.frame(ID=subset(RPR, V1 %in% gen4$ID)$V2, gen=5)
gen6 <- data.frame(ID=subset(RPR, V1 %in% gen5$ID)$V2, gen=6)
gen7 <- data.frame(ID=subset(RPR, V1 %in% gen6$ID)$V2, gen=7)

allgens <- rbind(gen0,gen1,gen2) #gen3 #,gen4,gen5,gen6,gen7)

table(allgens$gen)


fgy2[grepl("TLR", fgy2$Description ),]
fgy2[grepl("granul", fgy2$Description ),]

fgys <- fgy2[fgy2$ID %in% allgens$ID, ]


head(fgys[order(fgys$pvalue),], 40)
head(fgys[order(fgys$NES),], 30)
head(fgys[order(-fgys$NES),], 30)

# load("J:/My Drive/Axes/XPN.Rdata")


# upsig <- merge( psig_up, xpn, by.x=0, by.y="DB_ID" )
# 
# head(upsig[order(upsig$`sapply(biggens$ID, pathsig)`),],20)
# 
# boxplot( psout_up[, "R-HSA-1632852"] ~ clinex$mort7  )


#### merge to parents

gen01 <- RPR[RPR$V2 %in% gen1$ID,]

names(gen01) <- c("gen0","gen1")

gen012 <- merge(gen01, RPR, by.x="gen1", by.y="V1")  

names(gen012)[3] <- "gen2"

gen0123 <- merge(gen012, RPR, by.x="gen2", by.y="V1")  

names(gen0123)[4] <- "gen3"

gen01234 <- merge(gen0123, RPR, by.x="gen3", by.y="V1")  

names(gen01234)[5] <- "gen4"


gen012345 <- merge(gen01234, RPR, by.x="gen4", by.y="V1")  

names(gen012345)[6] <- "gen5"


gen0123456 <- merge(gen012345, RPR, by.x="gen5", by.y="V1")  

names(gen0123456)[7] <- "gen6"


fgys <- fgy2[fgy2$ID %in% allgens$ID, ]


table(fgys$ID %in% gen01$gen1 )
table(fgys$ID %in% gen012$gen2 )
table(fgys$ID %in% gen0123$gen3 )
table(fgys$ID %in% gen01234$gen4 )
table(fgys$ID %in% gen012345$gen5 )
table(fgys$ID %in% gen0123456$gen6 )

fgys$parent <- NA
fgys$parent[fgys$ID %in% gen01$gen1] <- gen01$gen0[match(fgys$ID[fgys$ID %in% gen01$gen1], gen01$gen1)]
fgys$parent[fgys$ID %in% gen012$gen2] <- gen012$gen0[match(fgys$ID[fgys$ID %in% gen012$gen2], gen012$gen2)]
fgys$parent[fgys$ID %in% gen0123$gen3] <- gen0123$gen0[match(fgys$ID[fgys$ID %in% gen0123$gen3], gen0123$gen3)]
fgys$parent[fgys$ID %in% gen01234$gen4] <- gen01234$gen0[match(fgys$ID[fgys$ID %in% gen01234$gen4], gen01234$gen4)]
fgys$parent[fgys$ID %in% gen012345$gen5] <- gen012345$gen0[match(fgys$ID[fgys$ID %in% gen012345$gen5], gen012345$gen5)]
fgys$parent[fgys$ID %in% gen0123456$gen6] <- gen0123456$gen0[match(fgys$ID[fgys$ID %in% gen0123456$gen6], gen0123456$gen6)]

head(fgys)

table(fgys$p.adjust<0.05)

fgyss <- fgys[fgys$p.adjust < 0.05,]
fgyss <- fgys[fgys$p.adjust < 0.05,]


# doesnt order to NES yet 
library(ggplot2)
ggplot(data=fgyss, aes(x=NES, y=Description, group=parent)  )+
  geom_point()+
  facet_grid(parent ~ ., scales = "free", space = "free") +
  theme(strip.text.y = element_text(angle = 0))

fgyss$parent <- factor(
  fgyss$parent,
  levels = c("Hemostasis", "Innate immunity", "Cytokine signaling", "Adaptive immunity")
)

# order by NES
fgyss <- fgyss %>%
  group_by(parent) %>%
  mutate(Description = factor(Description, levels = Description[order(NES)])) %>%
  ungroup()

# add bins for p values
colnames(fgyss)
fgyss$p_bin <- cut(
  fgyss$p.adjust,
  breaks = c(-Inf, 0.01, 0.05),
  labels = c("< 0.01", "< 0.05")
)

table(fgyss$p_bin)

# adjust names for parents manually 
fgyss$parent[fgyss$parent == "R-HSA-1280218"] <- "Adaptive immunity"
fgyss$parent[fgyss$parent == "R-HSA-168249"]  <- "Innate immunity"
fgyss$parent[fgyss$parent == "R-HSA-1280215"] <- "Cytokine signaling"
fgyss$parent[fgyss$parent == "R-HSA-109582"]  <- "Hemostasis"
fgyss$parent[fgyss$parent == "R-HSA-1430728"] <- "Metabolism"
fgyss$parent[fgyss$parent == "R-HSA-5357801"] <- "Programmed cell death"


# with aptos?
ggplot(data = fgyss, aes(x = NES, y = Description, group = parent)) +
  geom_point(aes(fill = NES, size = p_bin), shape = 21, color = "black", stroke = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  facet_grid(parent ~ ., scales = "free", space = "free") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "#F6423C") +
  scale_size_manual(
    values = c("< 0.01" = 4, "< 0.05" = 2),
    name = "p-value"
  ) +
  scale_x_continuous(limits = c(-3, 3)) +
  theme_minimal(base_family = "Aptos") +  # apply Aptos as base font
  theme(
    strip.text.y = element_text(angle = 0, family = "Aptos"),
    axis.text = element_text(family = "Aptos"),
    axis.title = element_text(family = "Aptos"),
    legend.text = element_text(family = "Aptos"),
    legend.title = element_text(family = "Aptos")
  )+
  ylab(NULL)



