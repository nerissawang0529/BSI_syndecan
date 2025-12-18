#### repeat code pathegone but without CoNS

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

event_microbiological_specimen_blood$culture_day <- round(difftime(as.POSIXct(event_microbiological_specimen_blood$Specimen_datetime,format="%Y-%m-%d"),
                                                                   as.POSIXct(event_microbiological_specimen_blood$Index_ICU_admittance_datetime2,format="%Y-%m-%d"),units="days"))

#code made by Joe for pathogens#### 
#this is the ICU_addmition_data
clinical_marker_scource <- import("original_data/clinical_marker_scource.csv")

subset_MARSID <- clinical_marker_scource$MARSID[clinical_marker_scource$Microbe_groups == "Bacteremia"]

ss <- event_microbiological_specimen_blood[event_microbiological_specimen_blood$ICU_ID_from_datasource %in% subset_MARSID, ]
length(unique(ss$ICU_ID_from_datasource))

ss01 <- ss[ss$culture_day == 0 | ss$culture_day == 1 | ss$culture_day == -1, ]

ss01 <- ss01[!is.na(ss01$microbe_result.x),]
length(unique(ss01$ICU_ID_from_datasource))


### check CNS
CNS <- ss01[grepl("epidermidi|haemolytic|hominis|oagulase|simulans|schleiferi|capitis" ,ss01$microbe_result),]



#diff_marsid_to_icu <- setdiff(subset_MARSID, ss01$ICU_ID_from_datasource)
#1906 2129 2405 2632 2914 3652  4596 10006 10369 10737 11688 13805 were in clinical_marker_scource bacteria group but not in ss01
#these IDs already excluded in code 'clinical_marker'

length(unique(ss01$ICU_ID_from_datasource))

ssx <- ss01[,c("microbe_result.x","ICU_ID_from_datasource")]

ssx <- unique(ssx)
length(unique(ssx$ICU_ID_from_datasource))
result <- aggregate(microbe_result.x ~ ICU_ID_from_datasource, data = ssx, FUN = function(x) paste(x, collapse = ", "))

#export clinical_marker_unique data
destination_folder <- "Documents/BSI/R_code/original_data" 
export_file_name <- "result_pathogen.csv" 
write.csv(result, file.path(destination_folder, export_file_name), row.names = FALSE)
cat("File exported to:", file.path(destination_folder, export_file_name), "\n")




# i want the bacteria group which number >10
#Klebsiella, S.aureus, CoNS, Enterococcus, E.coli, Streptococcus,Staph_unknown chosen by Joe,Staph_unknown <10.
#Pseudomonas,Enterobacter, clostridium, bacteroides made by myself, but Clostridium and Bacteroides were < 10.

result$Klebsiella <-  as.numeric(grepl("Klebsi", result$microbe_result )) #22
result$S.aureus <-  as.numeric(grepl("aureus", result$microbe_result )) #25
result$CoNS <-  as.numeric(grepl("epidermidi|haemolytic|hominis|oagulase|simulans|schleiferi|capitis", result$microbe_result )) #17
result$Enterococcus <- as.numeric(grepl("nterococcus|Enteroccus", result$microbe_result )) #37
result$E.coli <- as.numeric(grepl("coli", result$microbe_result )) #50
result$Streptococcus <- as.numeric(grepl("trepto", result$microbe_result )) #38
result$Pseudomonas <- as.numeric(grepl("Pseudomonas", result$microbe_result )) #11
result$Enterobacter <- as.numeric(grepl("Enterobacter", result$microbe_result )) #13

#result$Clostridium <- as.numeric(grepl("Clostridium", result$microbe_result )) #6
#result$Bacteroides <- as.numeric(grepl("Bacteroides", result$microbe_result )) #4
#result$Staph_unknown <- as.numeric(grepl("gram pos. coc, wrsch. Staphilococcus spp.", result$microbe_result ))

#result$microbenum <- rowSums(result[,3:10])
#table(result$microbenum)
#result0 <- result[result$microbenum==0,]

#result$other_bacteria <- ifelse(result$microbenum == 0, 1, 0)
### we need to remove Enterobacter
result$other_culture <- as.numeric(grepl("Proteus|Candida|Haemophil|Neisseria|Salmonella|Fusobact|Bacteroides|Serratia|Citrobact|Clostrid|Bacillus", result$microbe_result))

result$microbenum2 <- rowSums(result[,c("Klebsiella","S.aureus","CoNS","Enterococcus","E.coli","Streptococcus","Pseudomonas","Enterobacter","other_culture")])

### Exclude ambiguous (vague) bacteria decriptions (n=5) #IDs are 12328, 13323, 13700, 8493, 8690
#ambiguous_ids <- result$ICU_ID_from_datasource[result$microbenum2 == 0]
#already exclude in the clinical_marker_source

result2 <- result[result$microbenum2!=0,]

### Exclude CoNS (n=19) 
#IDs are 12034, 13104, 13501, 13744, 13906, 1535, 189, 2877, 3471, 364, 5305, 5408, 5764, 5841, 6312, 6437, 8563, 8758, 8811
### patients with only CNS
#CNS <- result[result$CoNS==1 & result$microbenum2==1,]
#already exclude in the clinical_marker_source

#result2 <- result2[!result2$ICU_ID_from_datasource %in% CNS$ICU_ID_from_datasource,]


#result0$microbe_result
#export clinical_marker_unique data

#check if the viral and fungal are only have viral and fungal
#viral
viral_marsid <- clinical_marker_scource %>% filter(Microbe_groups == "Viral") %>%select(MARSID)
# 4463 is not in the event_microbiological_specimen_blood, event_microbiological_specimen_blood$microbe_result for 4534,4810,4828,4938,4951,5259,5433,12376,12549,12629,12715 are NA
#fungal
fungal_marsid <- clinical_marker_scource %>% filter(Microbe_groups == "Fungal") %>%select(MARSID)
# 999, 2841, 3915, 5338, 5421,12503 (only fungal within -1 to 1 day)


### GN
result2$GN <- 0
result2$GN <- as.numeric(grepl("coli|Pseudomonas|Klebs|Haemophil|Neisseria|Salmonella|Fusobact|Bacteroides|Serratia|Citrobact|Enterobacter", result2$microbe_result))
result2$GP <- 0
result2$GP <- as.numeric(grepl("trepto|taph|Enterococcus|Enteroccus|Clostrid|Bacillus", result2$microbe_result))


table(result2$GN, result2$GP)
barplot(table(result2$microbenum2))
table(result2$microbenum2)
result2$FinalGroup_gram <- with(result2, ifelse(GN == 0 & GP == 1, "Gram_positive",
                                           ifelse(GN == 1 & GP == 0, "Gram_negative",
                                                  ifelse(GN == 1 & GP == 1, "Both", "Not"))))
table(result2$FinalGroup_gram)
result2[result2$microbenum2==4,]

#table(result2$Pseudomonas)

### Final categorical variable

result2$FinalGroup_pathegon <- "Other_pathegon"
result2$FinalGroup_pathegon[result2$microbenum2>=2] <- "Mixed_pathegon"
result2$FinalGroup_pathegon[result2$S.aureus==1 & result2$microbenum2==1] <- "S.aureus"
result2$FinalGroup_pathegon[result2$Klebsiella==1 & result2$microbenum2==1 ] <- "Klebsiella"
result2$FinalGroup_pathegon[result2$Enterococcus==1 & result2$microbenum2==1] <- "Enterococcus"
result2$FinalGroup_pathegon[result2$E.coli==1 & result2$microbenum2==1] <- "E.coli"
result2$FinalGroup_pathegon[result2$Streptococcus==1 & result2$microbenum2==1] <- "Streptococcus"
result2$FinalGroup_pathegon[result2$Pseudomonas==1 & result2$microbenum2==1] <- "Pseudomonas"
result2$FinalGroup_pathegon[result2$Enterobacter==1 & result2$microbenum2==1] <- "Enterobacter"

table(result2$FinalGroup_pathegon)
table(result2$FinalGroup_pathegon, result2$FinalGroup_gram)


table(result2$microbenum2)

### table patients with two bacteria (n=29)
bact2 <- result2[result2$microbenum2==2,]

### table patients with three bacteria (n=4)
bact3 <- result2[result2$microbenum2==3,]

### table patients with four bacteria (n=2)
bact4 <- result2[result2$microbenum2==4,]

destination_folder <- "Documents/BSI/R_code/original_data" 
export_file_name <- "bact_two.csv" 
write.csv(bact2, file.path(destination_folder, export_file_name), row.names = FALSE)


### Streptococcus mono (N=34)
strep_mono <- result2[result2$microbenum2==1 & result2$Streptococcus==1,]
table(strep_mono$microbe_result)

strep_mono$strep_group <- strep_mono$microbe_result
strep_mono$strep_group[grepl("Streptococcus pneumoniae",strep_mono$microbe_result)] <- "Streptococcus pneumoniae"
strep_mono$strep_group[grepl("groep A",strep_mono$microbe_result)] <- "Group A"
strep_mono$strep_group[grepl("Hem. streptococ groep B",strep_mono$microbe_result)] <- "Group B"
strep_mono$strep_group[grepl("vergroenend|Vergroenende|milleri|anginosus",strep_mono$microbe_result)] <- "Streptococcus viridans"
strep_mono$strep_group[grepl("wrsch",strep_mono$microbe_result)] <- "Streptococcus species (unspecified)"

table(strep_mono$strep_group)

str_tab <- as.data.frame(table(strep_mono$strep_group))
str_tab <- str_tab[order(-str_tab$Freq),]
str_tab

### Enterococcus mono (N=21)
entero_mono <- result2[result2$microbenum2==1 & result2$Enterococcus==1,]
table(entero_mono$microbe_result)

entero_mono$entero_group <- NA
entero_mono$entero_group[grepl("faecium",entero_mono$microbe_result)] <- "Enterotococcus faecium"
entero_mono$entero_group[grepl("faecalis",entero_mono$microbe_result)] <- "Enterococcus faecalis"
entero_mono$entero_group[grepl("Enterococcus faecalis Enterococcus faecium",entero_mono$microbe_result)] <- "Enterotococcus faecium and faecalis"
entero_mono$entero_group[grepl("species",entero_mono$microbe_result)] <- "Enterococcus species (unspecified)"

ent_tab <- as.data.frame(table(entero_mono$entero_group))
ent_tab <- ent_tab[order(-ent_tab$Freq),]

### Other mono (N=16)
oth_mono <- result2[ result2$FinalGroup_pathegon=="Other_pathegon",]
oth_mono <- result2[result2$microbenum2==1 & result2$other_culture==1,]

oth_mono$microbe_result[oth_mono$microbe_result=="Abiotrophia species Fusobacterium necrophorum, Fusobacterium necrophorum"] <-
  "Fusobacterium necrophorum"

table(oth_mono$microbe_result)

oth_tab <- as.data.frame(table(oth_mono$microbe_result))
oth_tab <- oth_tab[order(-oth_tab$Freq),]

oth_tab

#export clinical_marker_unique data
destination_folder <- "original_data" 
export_file_name <- "result2.csv" 
write.csv(result2, file.path(destination_folder, export_file_name), row.names = FALSE)
cat("File exported to:", file.path(destination_folder, export_file_name), "\n")



#merge together
clinical_marker_scource_pathogen <- merge(clinical_marker_scource, result2, 
                                          by.x = "MARSID", by.y = "ICU_ID_from_datasource")

table(clinical_marker_scource_pathogen$likelihood)


# #### 8 likelihood none
# liknon <- clinical_marker_scource_pathogen[clinical_marker_scource_pathogen$likelihood=="none",]
# table(liknon$FinalGroup_pathegon)


#### 8 likelihood none
likposs <- clinical_marker_scource_pathogen[clinical_marker_scource_pathogen$likelihood=="possible",]
table(likposs$FinalGroup_pathegon)

table(likposs$site)

#export clinical_marker_unique data
destination_folder <- "Documents/BSI/R_code/original_data" 
export_file_name <- "clinical_marker_scource_pathogen.csv" 
write.csv(clinical_marker_scource_pathogen, file.path(destination_folder, export_file_name), row.names = FALSE)
cat("File exported to:", file.path(destination_folder, export_file_name), "\n")





#Question
#1.the valuable for viral?
#2.how the non-infection was made?