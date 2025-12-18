### repeat code pathegone but without CoNS

rm(list = ls())

## install.packages("pacman")
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)
install.packages("haven")
library(haven)

## ready data
admission <- read_sas("original_data/icu_admission.sas7bdat")
#causative_organism_primary <- read_sas("Documents/BSI/R_code/original_data/causative_organism_primary.sas7bdat", NULL)
#causative_organism_secondary <- read_sas("Documents/BSI/R_code/original_data/causative_organism_secondary.sas7bdat", NULL)


############################################# microbiology result
event_microbiological_specimen <- read_sas("original_data/event_microbiological_specimen.sas7bdat", NULL)
#only blood
event_microbiological_specimen_blood <- event_microbiological_specimen[grepl("blood",event_microbiological_specimen$Material,ignore.case = T)==T,]
event_microbiological_specimen_blood <- event_microbiological_specimen_blood[!grepl("serum",event_microbiological_specimen_blood$Material,ignore.case = T)==T,]

#merge results
event_microbiological_result <- read_sas("original_data/event_microbiological_result.sas7bdat", NULL)
event_microbiological_result_blood <- event_microbiological_result[event_microbiological_result$Event_Microbiological_Specim_TK %in% 
                                                                     event_microbiological_specimen_blood$Event_Microbiological_Specim_TK &
                                                                     !(event_microbiological_result$Growth_density %in% "") &
                                                                     event_microbiological_result$Method=="culture",]

event_microbiological_result_blood <- event_microbiological_result_blood[,c("Event_Microbiological_Specim_TK", "Microbe")]

ss_c <- event_microbiological_result_blood
ss_c <- unique(ss_c)
ss_c <- ss_c[order(ss_c$Event_Microbiological_Specim_TK),]
ss_c <- within(ss_c, {
  microbe_result <- ave(Microbe,Event_Microbiological_Specim_TK, FUN=function(x) Reduce(paste, x)  )
})
event_microbiological_result_blood_flat<-unique(ss_c[,c("Event_Microbiological_Specim_TK", "microbe_result")] )
ss_c <- NULL

event_microbiological_specimen_blood <- merge(event_microbiological_specimen_blood,event_microbiological_result_blood_flat,
                                              by= "Event_Microbiological_Specim_TK", all.x = T)


event_microbiological_specimen_blood <- merge(event_microbiological_specimen_blood, admission[,c("ICU_Admission_TK",
                                                                                                 "ICU_ID_from_datasource")], all.x = T)

event_microbiological_specimen_blood$Index_ICU_admittance_datetime2 <-strptime(event_microbiological_specimen_blood$Index_ICU_admittance_datetime,"%Y-%m-%d")
event_microbiological_specimen_blood$Specimen_datetime <- strptime(event_microbiological_specimen_blood$Specimen_datetime,"%Y-%m-%d")


#ready hospital
studydata <- import("original_data/MARS database.xlsx")
studydata <- studydata %>%
  distinct(ICU_ID_from_datasource, .keep_all = TRUE)

library(dplyr)

studydata <- studydata %>%
  mutate(ICU_ID_from_datasource = as.character(ICU_ID_from_datasource))

event_microbiological_specimen_blood <- event_microbiological_specimen_blood %>%
  mutate(ICU_ID_from_datasource = as.character(ICU_ID_from_datasource))

event_microbiological_specimen_blood <- studydata %>%
  inner_join(event_microbiological_specimen_blood, by = "ICU_ID_from_datasource")



event_microbiological_specimen_blood$culture_day <- round(difftime(as.POSIXct(event_microbiological_specimen_blood$Specimen_datetime,format="%Y-%m-%d"),
                                                                   as.POSIXct(event_microbiological_specimen_blood$First_Hospital_admittance_date,format="%Y-%m-%d"),units="days"))



#this is the ICU_addmition_data
clinical_marker_scource <- import("original_data/clinical_marker_scource.csv")


#exclude the IDs
# 3286  1063  3613 14497, there IDs have no source information, excluded them
exclude_values <- c(3286, 1063, 3613, 14497)
clinical_marker_scource <- clinical_marker_scource[!clinical_marker_scource$MARSID %in% exclude_values, ]

# 6611 1967 2791 were outliers for markers data for both all marker and also for marker except Cystatin.C..75.
exclude_values_2 <- c(6611, 1967, 2791)
clinical_marker_scource <- clinical_marker_scource[!clinical_marker_scource$MARSID %in% exclude_values_2, ]
#the group of microgroups is not right

# 5433, the syndecan-1 value is 'NA', delete it for now
clinical_marker_scource <- clinical_marker_scource[!clinical_marker_scource$MARSID == 5433, ]

# 4440 is fungal and contamination, delete it for now
clinical_marker_scource <- clinical_marker_scource[!clinical_marker_scource$MARSID == 4440, ]


# 12328, 13323, 13700, 8493, 8690 are ambiguous (vague) bacteria decriptions (n=5) in the code pathegone_exclude_CNS
exclude_values_4 <- c(12328, 13323, 13700, 8493, 8690)
clinical_marker_scource <- clinical_marker_scource[!clinical_marker_scource$MARSID %in% exclude_values_4, ]
#Gram negatieve staven, nadere uitslag volgt. 8493,8690
#gram pos. coc, wrsch. Staphilococcus spp.13323,12328,13700



# Microbe_groups == "Bacteremia" & likelihood == "none"
filtered_data <- clinical_marker_scource[clinical_marker_scource$Microbe_groups == "Bacteremia" & 
                                           clinical_marker_scource$likelihood == "none", ]
exclude_values_5 <- c(1056,2724,3652,4596,4908,6386,6675,8495,8991,11399,13744) #n=11
clinical_marker_scource <- clinical_marker_scource[!clinical_marker_scource$MARSID %in% exclude_values_5, ]


##
subset_MARSID <- clinical_marker_scource$MARSID[clinical_marker_scource$Microbe_groups == "Bacteremia"]
ss <- event_microbiological_specimen_blood[event_microbiological_specimen_blood$ICU_ID_from_datasource %in% subset_MARSID, ]
length(unique(ss$ICU_ID_from_datasource))
ss01 <- ss[ss$culture_day == 0 | ss$culture_day == 1 | ss$culture_day == -1, ]
ss01 <- ss01[!is.na(ss01$microbe_result),]
length(unique(ss01$ICU_ID_from_datasource))
length(unique(ss01$ICU_ID_from_datasource))
ssx <- ss01[,c("microbe_result","ICU_ID_from_datasource")]
ssx <- unique(ssx)
length(unique(ssx$ICU_ID_from_datasource))
result <- aggregate(microbe_result ~ ICU_ID_from_datasource, data = ssx, FUN = function(x) paste(x, collapse = ", "))


### check CNS
CNS <- result[grepl("epidermidi|haemolytic|hominis|oagulase|simulans|schleiferi|capitis" ,result$microbe_result),]
#totally 33
#I checked the list, only exclude the IDs which only have CNS. for the ones also have other bactaria keep them in the corhort for now.
exclude_values_6 <- c(12034,13104,13501,13906,1535,189,	2877,	3471,364,5305,5408,5764,5841,6312,6437,8563,8758,8811) #18

##
result <- result[!result$ICU_ID_from_datasource%in% exclude_values_6,]

##
clinical_marker_scource <- clinical_marker_scource[!clinical_marker_scource$MARSID %in% exclude_values_6,]

# 1906 2129 2405 2632 2914 10006 10369 10737 11688 13805 are delete from 'ss01 <- ss[ss$culture_day == 0 | ss$culture_day == 1 | ss$culture_day == -1, ]'
exclude_values_3 <- c(1906,2129,2405,2632,2914,10006,10369,10737,11688,13805) #10
clinical_marker_scource <- clinical_marker_scource[!clinical_marker_scource$MARSID %in% exclude_values_3, ]

#export clinical_marker_unique data
destination_folder <- "original_data" 
export_file_name <- "clinical_marker_scource_pathogen.csv" 
write.csv(clinical_marker_scource, file.path(destination_folder, export_file_name), row.names = FALSE)
cat("File exported to:", file.path(destination_folder, export_file_name), "\n")



# i want the bacteria group which number >10
result$Klebsiella <-  as.numeric(grepl("Klebsi", result$microbe_result )) #22
result$S.aureus <-  as.numeric(grepl("aureus", result$microbe_result )) #25
#result$CoNS <-  as.numeric(grepl("epidermidi|haemolytic|hominis|oagulase|simulans|schleiferi|capitis", result$microbe_result )) #15
result$Enterococcus <- as.numeric(grepl("nterococcus|Enteroccus", result$microbe_result )) #37
result$E.coli <- as.numeric(grepl("coli", result$microbe_result )) #50
result$Streptococcus <- as.numeric(grepl("trepto", result$microbe_result )) #38
result$Pseudomonas <- as.numeric(grepl("Pseudomonas", result$microbe_result )) #11
result$Enterobacter <- as.numeric(grepl("Enterobacter", result$microbe_result )) #12

result$other_culture <- as.numeric(grepl("Proteus|Haemophil|Neisseria|Salmonella|Fusobact|Bacteroides|Serratia|Citrobact|Clostrid|Bacillus|Morganella|Abiotrophia|Aeromonas|Lactobacillus|Eggerthella", result$microbe_result))
result$microbenum2 <- rowSums(result[,c("Klebsiella","S.aureus","Enterococcus","E.coli","Streptococcus","Pseudomonas","Enterobacter","other_culture")])

result2 <- result[result$microbenum2!=0,]

#6116  Morganella morganii ss morganii,Proteus mirabilis,Streptococcus haemolyticus groep G?     'two others' keep this one for the mix-group
#5881  Abiotrophia species, Fusobacterium necrophorum ?                                          'two others' change this one for the mix-group
#1073  Escherichia coli, Grampositieve staaf nader gedetermineerd ?                              Escherichia coli is negative
#12926, 4755,8371 Amoxicilline resistente Enteroccus faecium Enterococcus faecalis Pseudomonas aeruginosa    two kinds of 'Enterococcus'
#CNS

#### GN
result2$GN <- 0
result2$GN <- as.numeric(grepl("coli|Pseudomonas|Klebs|Haemophil|Neisseria|Salmonella|Fusobact|Bacteroides|Serratia|Citrobact|Enterobacter|Proteus|Morganella|Aeromonas", result2$microbe_result))
result2$GP <- 0
result2$GP <- as.numeric(grepl("trepto|aureus|Enterococcus|Enteroccus|Clostrid|Bacillus|Abiotroph|Lactobacillus|Eggerthella", result2$microbe_result))


table(result2$GN, result2$GP)
barplot(table(result2$microbenum2))
table(result2$microbenum2)
result2$FinalGroup_gram <- with(result2, ifelse(GN == 0 & GP == 1, "Gram_positive",
                                                ifelse(GN == 1 & GP == 0, "Gram_negative",
                                                       
                                                       ifelse(GN == 1 & GP == 1, "Both", "Not"))))

table(result2$FinalGroup_gram)


### Final categorical variable

result2$FinalGroup_pathegon <- "Other_pathogens"
result2$FinalGroup_pathegon[result2$microbenum2>=2] <- "Mixed_pathegon"
result2$FinalGroup_pathegon[result2$S.aureus==1 & result2$microbenum2==1] <- "S.aureus"
result2$FinalGroup_pathegon[result2$Klebsiella==1 & result2$microbenum2==1 ] <- "Klebsiella"
result2$FinalGroup_pathegon[result2$Enterococcus==1 & result2$microbenum2==1] <- "Enterococcus"
result2$FinalGroup_pathegon[result2$E.coli==1 & result2$microbenum2==1] <- "E.coli"
result2$FinalGroup_pathegon[result2$Streptococcus==1 & result2$microbenum2==1] <- "Streptococcus"
result2$FinalGroup_pathegon[result2$Pseudomonas==1 & result2$microbenum2==1] <- "Pseudomonas"
result2$FinalGroup_pathegon[result2$Enterobacter==1 & result2$microbenum2==1] <- "Enterobacter"

### change    #5881  Abiotrophia species, Fusobacterium necrophorum  from other to mixed
result2$FinalGroup_pathegon[result2$ICU_ID_from_datasource==5881] <- "Mixed_pathegon"
result2$microbenum2[result2$ICU_ID_from_datasource==5881] <- "2"

### change 13679# Gram negatieve staven, nadere uitslag volgt gram pos. coc, wrsch. Streptococcus spp., Gram negatieve staven, nadere uitslag volgt
result2$microbenum2[result2$ICU_ID_from_datasource==13679] <- "2"
result2$FinalGroup_pathegon[result2$ICU_ID_from_datasource==13679] <- "Mixed_pathegon"
result2$FinalGroup_gram[result2$ICU_ID_from_datasource==13679] <- "Both"

### change   #6116 Morganella morganii ss morganii,Proteus mirabilis,Streptococcus haemolyticus groep G
result2$microbenum2[result2$ICU_ID_from_datasource==6116] <- "3"


table(result2$FinalGroup_pathegon)
FinalGroup_pathegon_gram <- table(result2$FinalGroup_pathegon, result2$FinalGroup_gram)


destination_folder <- "original_data"
export_file_name <- "result2_pathegon.csv" 
write.csv(result2, file.path(destination_folder, export_file_name), row.names = TRUE)


clinical_marker_scource_pathogen <- merge(result2, clinical_marker_scource, 
                                          by.x = "ICU_ID_from_datasource", 
                                          by.y = "MARSID", 
                                          all.x = TRUE)

destination_folder <- "original_data"
export_file_name <- "clinical_marker_scource_pathogen.csv" 
write.csv(clinical_marker_scource_pathogen, file.path(destination_folder, export_file_name), row.names = TRUE)



destination_folder <- "original_data" 
export_file_name <- "FinalGroup_pathegon_gram.csv" 
write.csv(FinalGroup_pathegon_gram, file.path(destination_folder, export_file_name), row.names = TRUE)



table(result2$microbenum2)


### table patients with one bacteria (n=)
bact1 <- result2[result2$microbenum2==1,]
bact1 <- bact1 %>%
  arrange(FinalGroup_pathegon)

destination_folder <- "original_data" 
export_file_name <- "bact_one.csv" 
write.csv(bact1, file.path(destination_folder, export_file_name), row.names = FALSE)





### table patients with two bacteria (n=16)
bact2 <- result2[result2$microbenum2==2,]
destination_folder <- "original_data" 
export_file_name <- "bact_two.csv" 
write.csv(bact2, file.path(destination_folder, export_file_name), row.names = FALSE)

### table patients with three bacteria (n=3)
bact3 <- result2[result2$microbenum2==3,]
destination_folder <- "original_data" 
export_file_name <- "bact_three.csv" 
write.csv(bact3, file.path(destination_folder, export_file_name), row.names = FALSE)

### table patients with four bacteria (n=4)
bact4 <- result2[result2$microbenum2==4,]
destination_folder <- "original_data" 
export_file_name <- "bact_four.csv" 
write.csv(bact4, file.path(destination_folder, export_file_name), row.names = FALSE)


### Streptococcus mono (N=32)
strep_mono <- result2[result2$microbenum2==1 & result2$Streptococcus==1,]
table(strep_mono$microbe_result)

strep_mono$strep_group <- strep_mono$microbe_result
strep_mono$strep_group[grepl("Streptococcus pneumoniae",strep_mono$microbe_result)] <- "Streptococcus pneumoniae"
strep_mono$strep_group[grepl("groep A",strep_mono$microbe_result)] <- "Group A"
strep_mono$strep_group[grepl("Hem. streptococ groep B",strep_mono$microbe_result)] <- "Group B"
strep_mono$strep_group[grepl("groep G",strep_mono$microbe_result)] <- "groep G"
strep_mono$strep_group[grepl("vergroenend|Vergroenende|milleri|anginosus",strep_mono$microbe_result)] <- "Streptococcus viridans"
strep_mono$strep_group[grepl("wrsch",strep_mono$microbe_result)] <- "Streptococcus species (unspecified)"

table(strep_mono$strep_group)

str_tab <- as.data.frame(table(strep_mono$strep_group))
str_tab <- str_tab[order(-str_tab$Freq),]
str_tab

destination_folder <- "original_data" 
export_file_name <- "str_tab.csv" 
write.csv(str_tab, file.path(destination_folder, export_file_name), row.names = FALSE)

### Enterococcus mono (N=24)
entero_mono <- result2[result2$microbenum2==1 & result2$Enterococcus==1,]
table(entero_mono$microbe_result)

entero_mono$entero_group <- NA
entero_mono$entero_group[grepl("faecium",entero_mono$microbe_result)] <- "Enterotococcus faecium"
entero_mono$entero_group[grepl("faecalis",entero_mono$microbe_result)] <- "Enterococcus faecalis"
entero_mono$entero_group[grepl("Enterococcus faecalis Enterococcus faecium",entero_mono$microbe_result)] <- "Enterotococcus faecium and faecalis"
entero_mono$entero_group[grepl("species",entero_mono$microbe_result)] <- "Enterococcus species (unspecified)"

ent_tab <- as.data.frame(table(entero_mono$entero_group))
ent_tab <- ent_tab[order(-ent_tab$Freq),]
ent_tab

destination_folder <- "original_data" 
export_file_name <- "ent_tab.csv" 
write.csv(ent_tab, file.path(destination_folder, export_file_name), row.names = FALSE)

# Klebsiella (N=14)
klebsiella_mono <- result2[result2$microbenum2==1 & result2$Klebsiella==1,]
table(klebsiella_mono$microbe_result)

klebsiella_mono$klebsiella_group <- NA
klebsiella_mono$klebsiella_group[grepl("oxytoca",klebsiella_mono$microbe_result)] <- "Klebsiella oxytoca"
klebsiella_mono$klebsiella_group[grepl("pneumoniae",klebsiella_mono$microbe_result)] <- "Klebsiella pneumoniae"

kle_tab <- as.data.frame(table(klebsiella_mono$klebsiella_group))
kle_tab <- kle_tab[order(-kle_tab$Freq),]

destination_folder <- "original_data" 
export_file_name <- "kle_tab.csv" 
write.csv(kle_tab, file.path(destination_folder, export_file_name), row.names = FALSE)

#S.aureus (N=23)
S_aureus_mono <- result2[result2$microbenum2==1 & result2$S.aureus==1,]
table(S_aureus_mono$microbe_result)

S_aureus_mono$S_aureus_group <- NA
S_aureus_mono$S_aureus_group[grepl("aureus",S_aureus_mono$microbe_result)] <- "Staphylococcus aureus"

S_aureus_tab <- as.data.frame(table(S_aureus_mono$S_aureus_group))
S_aureus_tab <- S_aureus_tab[order(-S_aureus_tab$Freq),]

destination_folder <- "original_data" 
export_file_name <- "S_aureus_tab.csv" 
write.csv(S_aureus_tab, file.path(destination_folder, export_file_name), row.names = FALSE)


#
E_coli_mono <- result2[result2$microbenum2==1 & result2$E.coli==1,]
table(E_coli_mono$microbe_result)
#1073 Escherichia coli Grampositieve staaf nader gedetermineerd
E_coli_mono$E_coli_group <- NA
E_coli_mono$E_coli_group[grepl("coli",E_coli_mono$microbe_result)] <- "Escherichia coli"

E_coli_tab <- as.data.frame(table(E_coli_mono$E_coli_group))

destination_folder <- "original_data" 
export_file_name <- "E_coli_tab.csv" 
write.csv(E_coli_tab, file.path(destination_folder, export_file_name), row.names = FALSE)

#
Pseudomonas_mono <- result2[result2$microbenum2==1 & result2$Pseudomonas==1,]
table(Pseudomonas_mono$microbe_result)

Pseudomonas_mono$Pseudomonas_group <- NA
Pseudomonas_mono$Pseudomonas_group[grepl("Pseudomonas",Pseudomonas_mono$microbe_result)] <- "Pseudomonas aeruginosa"

Pseudomonas_tab <- as.data.frame(table(Pseudomonas_mono$Pseudomonas_group))


destination_folder <- "original_data" 
export_file_name <- "Pseudomonas_tab.csv" 
write.csv(Pseudomonas_tab, file.path(destination_folder, export_file_name), row.names = FALSE)


#
Enterobacter_mono <- result2[result2$microbenum2==1 & result2$Enterobacter==1,]
table(Enterobacter_mono$microbe_result)

Enterobacter_mono$Enterobacter_group <- NA
Enterobacter_mono$Enterobacter_group[grepl("aerogenes",Enterobacter_mono$microbe_result)] <- "Enterobacter aerogenes"
Enterobacter_mono$Enterobacter_group[grepl("cloacae",Enterobacter_mono$microbe_result)] <- "Enterobacter cloacae"
Enterobacter_mono$Enterobacter_group[grepl("species",Enterobacter_mono$microbe_result)] <- "Enterobacter species"

Enterobacter_tab <- as.data.frame(table(Enterobacter_mono$Enterobacter_group))
Enterobacter_tab <- Enterobacter_tab[order(-Enterobacter_tab$Freq),]

destination_folder <- "original_data" 
export_file_name <- "Enterobacter_tab.csv" 
write.csv(Enterobacter_tab, file.path(destination_folder, export_file_name), row.names = FALSE)

#need to check the outliers

