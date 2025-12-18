rm(list = ls())

# install.packages("pacman")
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)


#outliers
biom <- read.csv("Documents/BSI/R_code/original_data/clinical_marker_unique_log.csv")

##log
#all marker
#biom_x <- biom[, 51:71] 
#not include Cystatin.C..75.
biom_x <- biom[, c(51:55, 57:71)]
#endothelial markers
#biom_x <- biom[, c(51:56, 68)]
#only for syndecan
#biom_x <- biom[, 68, drop = FALSE]
str(biom_x)
rownames(biom_x) <- biom$ICU_ID_from_datasource


#deal with NA
#in "CX3CL1.Fractalkine..46." "MMP.8..27." there is NA
biom_x$CX3CL1.Fractalkine..46.[is.na(biom_x$CX3CL1.Fractalkine..46.)] <- median(biom_x$CX3CL1.Fractalkine..46., na.rm = TRUE)
biom_x$MMP.8..27.[is.na(biom_x$MMP.8..27.)] <- median(biom_x$MMP.8..27., na.rm = TRUE)


##mahalobis distance
m.dist <- mahalanobis(biom_x, colMeans(biom_x), cov(biom_x))
hist(m.dist, breaks = 10)
boxplot(m.dist, main="Boxplot of Mahalanobis Distance", ylab="Mahalanobis Distance")
head(m.dist[order(-m.dist)])

### remove 1% most extreme  points (outliers)
percentage.to.remove <- 1 # Remove 1% of points
number.to.remove     <- trunc(nrow(biom_x) * percentage.to.remove / 100)
m.dist               <- mahalanobis(biom_x, colMeans(biom_x), cov(biom_x))
m.dist.order         <- order(m.dist, decreasing=TRUE)
rows.to.keep.index   <- m.dist.order[(number.to.remove+1):nrow(biom_x)]
biom_x               <- biom_x[rows.to.keep.index,]


### see which outliers removed
(1:307)[ !(1:307) %in% rows.to.keep.index ]
### outlier patient IDs
biom$ICU_ID_from_datasource[(1:307)[ !(1:307) %in% rows.to.keep.index ]]

#The outliers for all of the markers are 6611 1967 2791
#The outliers for marker except for Cystatin.C..75. 6611 1967 2791
#The outliers for endothelial are 4298 14507  1535
#The outliers for syndecan-1 are 273 4570 1038

