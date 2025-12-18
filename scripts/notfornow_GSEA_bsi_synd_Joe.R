rm(list = ls())

#package

library(readxl)
install.packages("BiocManager")
BiocManager::install("GEOquery")
library(GEOquery)
install.packages("ggplot2")
library(ggplot2)

#Ready data
bsi_synd <- as.data.frame( read_excel("original_data/Kopie van BSI pathogen list TVDPOLL2.xlsx") )
bsi_synd <- bsi_synd[bsi_synd$ID != 715, ] #omit this because 'bacillus species'

load("original_data/u2pdat.Rdata")
load("original_data/exg.Rdata")
load("original_data/u219_expression_data.Rdata")

table(bsi_synd$ID %in% u2pdat$MARSID)
bsi_syn_u2 <- bsi_synd[bsi_synd$ID %in% u2pdat$MARSID,]
table(bsi_syn_u2$`Pathogen group`)
u2pdat_bsi <- u2pdat[u2pdat$MARSID %in% bsi_synd$ID , ]
bsi_m <- merge(u2pdat, bsi_synd, by.x="MARSID", by.y="ID" )
exg <- exg[!(row.names(exg) == '13_03_2013_C02_715.CEL'), ]

#dfmcb <- dfmc[,colnames(dfmc) %in% u2pdat_bsi$chipID  ]
#dim(dfmcb)
#tex <- as.data.frame((dfmcb))
#save(tex, file="J:/My Drive/Nerissa/BSI_syndecan/u219_expression_data.Rdata")



### select most variable probe per gene
#probevars <- apply(tex,1, var)
#probevars <- as.data.frame(probevars)
#head(probevars)

#gpl13667 <- getGEO('GPL13667')
#gpl13667 <- getGEO(filename='GPL13667.soft')

#f <- as.data.frame(rownames(dfmc), stringsAsFactors = FALSE)
#row.names(f) <- f[,1]
#fdt <- Table(gpl13667)[,c(1,15,21,26)]
#row.names(fdt) <- fdt[,1]
#fdt2 <- merge(f, fdt, by="row.names")
#row.names(fdt2) <- fdt2[,2]
#fdt2$Row.names <- NULL
#fdt2$ID <- NULL
#featureData = new("AnnotatedDataFrame", data= fdt2)

#featannot <- featureData@data
#featannot$EntrezU <- sapply( strsplit( featannot$`Entrez Gene`, " /// "  ), `[`, 1 )

#library(splitstackshape)

#featannot2 <- cSplit(featannot, "Gene Symbol", sep = " /// ", direction = "long")
#featannot2 <- featannot2[featannot2$`Gene Symbol`!="/",]


#library(data.table)
#Probevarm <- merge(probevars, featannot2[,c("rownames(dfmc)","Gene Symbol","EntrezU")] , by.x=0 ,by.y="rownames(dfmc)"  )

#Probevarm <- as.data.table(Probevarm)

#GenesU219 <- Probevarm[Probevarm[, .I[which.max(probevars)], by='EntrezU']$V1]
#GenesU219 <- as.data.frame(GenesU219)

#dim(GenesU219)


#texg <- tex[row.names(tex) %in% GenesU219$Row.names, ]
#texg <- merge(GenesU219[,c("Row.names","Gene Symbol","EntrezU")],texg, by.x="Row.names", by.y=0)



### flip transpose
#row.names(texg) <- texg$EntrezU
#texg$EntrezU <- NULL
#texg$`Gene Symbol` <- NULL
#texg$Row.names <- NULL

#exg <- as.data.frame(t(texg))

### 1.E.coli 

EC <- bsi_m$chipID[bsi_m$`Pathogen group`=="E.coli"]

exg$group <- row.names(exg) %in% EC
table(exg$group)

dfq2 <- data.frame(row.names = colnames(exg)[!grepl("group",colnames(exg))], 
                    gene= colnames(exg)[!grepl("group",colnames(exg))] ,  
                    t=NA )

for( g in dfq2$gene ){
  
  dfq2[g,"t"] <- t.test(exg[,g] ~ exg$group)$statistic
}

head(dfq2[order(-dfq2$t),])

### GSEA
der <- dfq2[order(-dfq2$t),]
derv <- der$t
names(derv) <- der$gene

set.seed(55513)

library(ReactomePA )

y2 <- gsePathway(derv, 
                pvalueCutoff = 0.2,
                pAdjustMethod = "BH", 
                verbose = FALSE,
                seed=T)
head(y2)


y2df <- as.data.frame(y2)

selected_pathway_id1 <- c("R-HSA-109582")
selected_pathway_id2 <- c("R-HSA-1474244")

RPR <- read.table("original_data/reactomepathways.txt")

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

yi <- y2df[y2df$ID %in% c(hemo$ID, extr$ID) ,]
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


### 2.Strep

STR <- bsi_m$chipID[bsi_m$`Pathogen group`=="Streptococcus"]

exg$group <- row.names(exg) %in% STR
table(exg$group )

dfq2 <- data.frame( row.names = colnames(exg)[!grepl("group",colnames(exg))], 
                    gene= colnames(exg)[!grepl("group",colnames(exg))] ,  
                    t=NA )

for( g in dfq2$gene ){
  
  dfq2[g,"t"] <- t.test(exg[,g] ~ exg$group  )$statistic
}


head(dfq2[order(-dfq2$t),])


### GSEA
der <- dfq2[order(-dfq2$t),]
derv <- der$t
names(derv) <- der$gene


set.seed(55513)

library(ReactomePA )

y2 <- gsePathway(derv, 
                 pvalueCutoff = 0.2,
                 pAdjustMethod = "BH", 
                 verbose = FALSE,
                 seed=T)
head(y2)


y2df <- as.data.frame(y2)

selected_pathway_id1 <- c("R-HSA-109582")
selected_pathway_id2 <- c("R-HSA-1474244")

RPR <- read.table("original_data/reactomepathways.txt")

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



yi <- y2df[y2df$ID %in% c(hemo$ID, extr$ID) ,]
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





yi[yi$p.adjust<0.05,]


#


####3.S.aureus##
S.aureus <- bsi_m$chipID[bsi_m$`Pathogen group`=="S.aureus"]

exg$group <- row.names(exg) %in% S.aureus
table(exg$group )

dfq2 <- data.frame( row.names = colnames(exg)[!grepl("group",colnames(exg))], 
                    gene= colnames(exg)[!grepl("group",colnames(exg))] ,  
                    t=NA )

for( g in dfq2$gene ){
  
  dfq2[g,"t"] <- t.test(exg[,g] ~ exg$group  )$statistic
}


head(dfq2[order(-dfq2$t),])


### GSEA
der <- dfq2[order(-dfq2$t),]
derv <- der$t
names(derv) <- der$gene


set.seed(55513)

library(ReactomePA )

y2 <- gsePathway(derv, 
                 pvalueCutoff = 0.2,
                 pAdjustMethod = "BH", 
                 verbose = FALSE,
                 seed=T)
head(y2)


y2df <- as.data.frame(y2)

selected_pathway_id1 <- c("R-HSA-109582")
selected_pathway_id2 <- c("R-HSA-1474244")

RPR <- read.table("original_data/reactomepathways.txt")

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



yi <- y2df[y2df$ID %in% c(hemo$ID, extr$ID) ,]
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



###4.Klebsiella####
Klebsiella <- bsi_m$chipID[bsi_m$`Pathogen group`=="Klebsiella"]

exg$group <- row.names(exg) %in% Klebsiella
table(exg$group )




dfq2 <- data.frame( row.names = colnames(exg)[!grepl("group",colnames(exg))], 
                    gene= colnames(exg)[!grepl("group",colnames(exg))] ,  
                    t=NA )

for( g in dfq2$gene ){
  
  dfq2[g,"t"] <- t.test(exg[,g] ~ exg$group  )$statistic
}


head(dfq2[order(-dfq2$t),])


### GSEA
der <- dfq2[order(-dfq2$t),]
derv <- der$t
names(derv) <- der$gene


set.seed(55513)

library(ReactomePA )

y2 <- gsePathway(derv, 
                 pvalueCutoff = 0.2,
                 pAdjustMethod = "BH", 
                 verbose = FALSE,
                 seed=T)
head(y2)


y2df <- as.data.frame(y2)

selected_pathway_id1 <- c("R-HSA-109582")
selected_pathway_id2 <- c("R-HSA-1474244")

RPR <- read.table("original_data/reactomepathways.txt")

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



yi <- y2df[y2df$ID %in% c(hemo$ID, extr$ID) ,]

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


####5.Enterococcus####
Enterococcus <- bsi_m$chipID[bsi_m$`Pathogen group`=="Enterococcus"]

exg$group <- row.names(exg) %in% Enterococcus
table(exg$group )




dfq2 <- data.frame( row.names = colnames(exg)[!grepl("group",colnames(exg))], 
                    gene= colnames(exg)[!grepl("group",colnames(exg))] ,  
                    t=NA )

for( g in dfq2$gene ){
  
  dfq2[g,"t"] <- t.test(exg[,g] ~ exg$group  )$statistic
}


head(dfq2[order(-dfq2$t),])


### GSEA
der <- dfq2[order(-dfq2$t),]
derv <- der$t
names(derv) <- der$gene


set.seed(55513)

library(ReactomePA )

y2 <- gsePathway(derv, 
                 pvalueCutoff = 0.2,
                 pAdjustMethod = "BH", 
                 verbose = FALSE,
                 seed=T)
head(y2)


y2df <- as.data.frame(y2)

selected_pathway_id1 <- c("R-HSA-109582")
selected_pathway_id2 <- c("R-HSA-1474244")

RPR <- read.table("original_data/reactomepathways.txt")

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



yi <- y2df[y2df$ID %in% c(hemo$ID, extr$ID) ,]

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

####6.Pseudomonas###
Pseudomonas <- bsi_m$chipID[bsi_m$`Pathogen group`=="Pseudomonas"]

exg$group <- row.names(exg) %in% Pseudomonas
table(exg$group )

dfq2 <- data.frame( row.names = colnames(exg)[!grepl("group",colnames(exg))], 
                    gene= colnames(exg)[!grepl("group",colnames(exg))] ,  
                    t=NA )

for( g in dfq2$gene ){
  
  dfq2[g,"t"] <- t.test(exg[,g] ~ exg$group  )$statistic
}


head(dfq2[order(-dfq2$t),])


### GSEA
der <- dfq2[order(-dfq2$t),]
derv <- der$t
names(derv) <- der$gene


set.seed(55513)

library(ReactomePA )

y2 <- gsePathway(derv, 
                 pvalueCutoff = 0.5,
                 pAdjustMethod = "BH", 
                 verbose = FALSE,
                 seed=T)
head(y2)


y2df <- as.data.frame(y2)

selected_pathway_id1 <- c("R-HSA-109582")
selected_pathway_id2 <- c("R-HSA-1474244")

RPR <- read.table("original_data/reactomepathways.txt")

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



yi <- y2df[y2df$ID %in% c(hemo$ID, extr$ID) ,]

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


