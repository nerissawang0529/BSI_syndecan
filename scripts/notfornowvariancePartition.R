rm(list = ls())

#sensitivy analysis 1. ALL patients (N=190) 2. only those with known single pathogen (N=150)

#ready data
data <- read.csv("original_data/clinical_marker_scource_pathogen.csv")
#make the source values
data$FinalSource <- "unknown"
data$FinalSource[data$CNS_source==1] <- "CNS"
data$FinalSource[data$abdominal==1] <- "abdominal"
data$FinalSource[data$respiratory==1] <- "respiratory"
data$FinalSource[data$urinary==1] <- "urinary"
data$FinalSource[data$cardiovascular==1] <- "cardiovascular"
data$FinalSource[data$skin==1] <- "skin"
data$FinalSource[data$other_source==1] <- "other"



fit0 <- lm(data=data, log(data$Syndecan1) ~ 
             E.coli + Enterobacter + Enterococcus + Klebsiella + Pseudomonas + S.aureus + Streptococcus + 
             unknown_source + CNS_source + abdominal+ respiratory + urinary + cardiovascular + skin + other_source +
             age_yrs + gender + 
             Past_myocardial_infarction + Cerebrovascular_disease + ckd + Immune_deficiency+
             SOFAtot + MEWS_score + APACHE_IV_Score.x + Shock_on_admis)

summary(fit0)

# Extract standardized coefficients
coefficients_std <- summary(fit0)$coefficients[, 1]

# Calculate percentage contribution
contributions <- abs(coefficients_std) / sum(abs(coefficients_std)) * 100

# Combine coefficients with variable names
contribution_df <- data.frame(Variable = names(coefficients_std), Contribution = contributions)

# Sort contributions in descending order
contribution_df <- contribution_df[order(-contribution_df$Contribution), ]

# Print the percentage contributions
print(contribution_df)
#notes: the values need to check pathogen (mix_group, Other_pathogens) ; cardiovascular, site (unknown); Diabetes 



### variance partition
library("variancePartition")


SynExpr <- t(log(data$Syndecan1))
colnames(SynExpr) <- data$ICU_ID_from_datasource
row.names(SynExpr) <- "Synd"


exp2 <- rbind(SynExpr, SynExpr)
row.names(exp2) <- c("Synd","S2")



bsinfo <- data[,c("unknown_source","CNS_source","abdominal","respiratory","urinary","cardiovascular","skin","other_source", "E.coli","Enterobacter","Enterococcus","Klebsiella","Pseudomonas","S.aureus","Streptococcus","SOFAtot","Platelets_value_1",
                  "FinalGroup_pathegon", "FinalSource",
                  "APACHE_IV_Score.x", "Shock_on_admis",
                  "age_yrs", "gender",
                    "Past_myocardial_infarction","Cerebrovascular_disease","ckd","Immune_deficiency",
                    "SOFAtot","MEWS_score")]
row.names(bsinfo) <- data$ICU_ID_from_datasource

bsinfo$logplatelets <- log(bsinfo$Platelets_value_1)

bsinfo <- as.data.frame(bsinfo)

form <- ~ (1| FinalSource) + (1| FinalGroup_pathegon)
form <- ~ (1| FinalGroup_pathegon) + (1| FinalSource) + SOFAtot 
form <- ~ (1| FinalGroup_pathegon) + (1| FinalSource) + SOFAtot + logplatelets
form <- ~ (1| FinalGroup_pathegon) + (1| FinalSource) + SOFAtot + logplatelets + APACHE_IV_Score.x
form <- ~ (1|E.coli) + (1|Enterobacter) + (1|Enterococcus) + (1|Klebsiella) + (1|Pseudomonas) + (1|S.aureus) + (1|Streptococcus) + 
  (1|unknown_source) + (1|CNS_source) + (1|abdominal)+ (1|respiratory) + (1|urinary) + (1|cardiovascular) + (1|skin) + (1|other_source) +
  age_yrs + (1|gender) + 
  (1|Past_myocardial_infarction) + (1|Cerebrovascular_disease) + (1|ckd) + (1|Immune_deficiency) +
  SOFAtot + MEWS_score + APACHE_IV_Score.x + (1|Shock_on_admis)


varPart <- fitExtractVarPartModel(exp2, form, bsinfo)

vp <- sortCols(varPart)

plotPercentBars(vp)



### remove missing

bsinfo2 <- bsinfo[complete.cases(bsinfo),]

exp3 <- exp2[,row.names(bsinfo2)]


form <- ~ (1| FinalGroup_pathegon) + (1| FinalSource) 
form <- ~ (1| FinalGroup_pathegon) + (1| FinalSource) + SOFAtot 
form <- ~ (1| FinalGroup_pathegon) + (1| FinalSource) + SOFAtot + logplatelets


library(dplyr)

bsinfo2 <- bsinfo2 %>%
  mutate(across(c(E.coli,Enterobacter,Enterococcus,Klebsiella,Pseudomonas,S.aureus,Streptococcus,
                  unknown_source,CNS_source,abdominal,respiratory,urinary,cardiovascular,skin,other_source,
                  gender,Past_myocardial_infarction,Cerebrovascular_disease,ckd,Immune_deficiency,Shock_on_admis), as.factor))


form <- ~ (1|E.coli) + (1|Enterobacter) + (1|Enterococcus) + (1|Klebsiella) + (1|Pseudomonas) + (1|S.aureus) + (1|Streptococcus) + 
  (1|unknown_source) + (1|CNS_source) + (1|abdominal)+ (1|respiratory) + (1|urinary) + (1|cardiovascular) + (1|skin) + (1|other_source) +
  age_yrs + (1|gender) + 
  (1|Past_myocardial_infarction) + (1|Cerebrovascular_disease) + (1|ckd) + (1|Immune_deficiency) +
  SOFAtot + MEWS_score + APACHE_IV_Score.x + (1|Shock_on_admis)


varPart <- fitExtractVarPartModel(exp3, form, bsinfo2)

vp <- sortCols(varPart)

plotPercentBars(vp)

vp



# Initialize an empty list to store the results
results_list <- list()

# Repeat the random sampling 10 times
for (i in 1:10) {
  # Randomly sample 34 rows from the non-Streptococcus group
  random_rows <- bsinfo2 %>%
    filter(Streptococcus != 1) %>%
    sample_n(34)
  
  # Combine the Streptococcus == 1 rows with the random sample
  combined_data <- bsinfo2 %>%
    filter(Streptococcus == 1) %>%
    bind_rows(random_rows)
  
  form <- ~ (1|E.coli) + (1|Enterobacter) + (1|Enterococcus) + (1|Klebsiella) + (1|Pseudomonas) + (1|S.aureus) + (1|Streptococcus) + 
    (1|unknown_source) + (1|CNS_source) + (1|abdominal)+ (1|respiratory) + (1|urinary) + (1|cardiovascular) + (1|skin) + (1|other_source) +
    age_yrs + (1|gender) + 
    (1|Past_myocardial_infarction) + (1|Cerebrovascular_disease) + (1|ckd) + (1|Immune_deficiency) +
    SOFAtot + MEWS_score + APACHE_IV_Score.x + (1|Shock_on_admis)
  
  
  varPart <- fitExtractVarPartModel(exp3, form, bsinfo2)
  
  vp <- sortCols(varPart)
  
  # Append the combined data to the results list
  results_list[[i]] <- vp
}
