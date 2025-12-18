### plasma analysis swimmer
rm(list = ls())

library("readxl")

plasma <- read_excel("original_data/plasma data dates.xlsx")

class(plasma$ICU_admit_date)


### swim code

main1bcs <- merge(main1, bloodcultsBOTH, by.x="MARSID", 
                  by.y="ICU_ID_from_datasource", all.x=T)

main1bcs$PAX_DAY <- as.character(main1bcs$pax_date)

main1bcs$PAXhour0 <- paste0(main1bcs$PAX_DAY," ","00:00:00")
main1bcs$PAXhour0 <- as.POSIXct(main1bcs$PAXhour0, tz ='UTC')


main1bcs$ADMtime <- as.numeric(difftime(main1bcs$Index_ICU_admittance_datetime, 
                                        main1bcs$PAXhour0,  units="hours"))

main1bcs$PAXtime <- as.numeric(difftime(main1bcs$pax_datetime2, 
                                        main1bcs$PAXhour0,  units="hours"))


hist(main1bcs$PAXtime)
summary(main1bcs$PAXtime)



main1bcs$BCdatetime <- as.POSIXct(main1bcs$Specimen_datetime, format= "%d%b%Y:%H:%M:%OS",  tz ='UTC')

## time between pax day time0 and BC
main1bcs$BCtime <- as.numeric(difftime(main1bcs$BCdatetime, 
                                       main1bcs$PAXhour0,  units="hours"))

hist(main1bcs$BCtime)

## time between pax day time0 and ADMission
main1bcs$ADMtime <- as.numeric(difftime(main1bcs$Index_ICU_admittance_datetime, 
                                        main1bcs$PAXhour0,  units="hours"))
hist(main1bcs$ADMtime)


###
BCPOSx <- main1bcs[,c("Patient_TK.x","MARSID","BCtime","Microbe")] 

BCPOSx <- BCPOSx[BCPOSx$Microbe!="see description",]

BCPOSx$Patient_TK <- as.factor(BCPOSx$Patient_TK.x )

colnames(BCPOSx )[3] <- "time"

BCPOSx$Patient_TK.x <- NULL

BCPOSx$type <- "BC"




ADMx <- main1bcs[,c("Patient_TK.x","MARSID","ADMtime")] 

ADMx  <- ADMx[!duplicated(ADMx ),]

ADMx$Patient_TK <- as.factor(ADMx$Patient_TK.x )

colnames(ADMx  )[3] <- "time"

ADMx$Patient_TK.x <- NULL

ADMx$type <- "ADM"






PAXx <- main1bcs[,c("Patient_TK.x","MARSID","PAXtime")] 

PAXx  <- PAXx[!duplicated(PAXx ),]

PAXx$Patient_TK <- as.factor(PAXx$Patient_TK.x )

colnames(PAXx  )[3] <- "time"

PAXx$Patient_TK.x <- NULL

PAXx$type <- "PAXdraw"




ALLr <- rbind(PAXx, ADMx)

ALLr$Microbe <- as.character(NA)




PAXx2 <- main1bcs[,c("Patient_TK.x","MARSID","PAXtime", "pax_date")] 

PAXx2  <- PAXx2[!duplicated(PAXx2 ),]

PAXx2$Patient_TK <- as.factor(PAXx2$Patient_TK.x )

#colnames(PAXx  )[3] <- "time"

PAXx2$Patient_TK.x <- NULL

#PAXx$type <- "PAXdraw"

# pax_ALL <- ALL2r[ALL2r$type=="PAXdraw",]
# 
# table(pax_ALL$MARSID %in% plasma$MARSID)


pax_m <- merge(PAXx2, plasma, by="MARSID")

pax_m$dd <- difftime(pax_m$event_date, pax_m$pax_date)


####

PAXx2 <- main1bcs[,c("Patient_TK.x","MARSID", "PAXhour0")] 

PAXx2  <- PAXx2[!duplicated(PAXx2 ),]

PAXx2$Patient_TK <- as.factor(PAXx2$Patient_TK.x )

#colnames(ADMx2  )[3] <- "time"

PAXx2$Patient_TK.x <- NULL

# adm_m <- merge(ADMx2, plasma, by="MARSID")
# 
# adm_m$dd <- difftime(adm_m$event_date, adm_m$ICU_admittance_date, units="days")
# 
# table(adm_m$dd)



###
table(plasma$MARSID %in% ALLr$MARSID)


plasma2 <- plasma[plasma$MARSID %in% ALLr$MARSID,]

pax_m2 <- merge(PAXx2, plasma2, by="MARSID")

pax_m2$PLASMAtime <- as.numeric(difftime(pax_m2$sample_date, 
                                         pax_m2$PAXhour0,  units="hours"))


hist(pax_m2$PLASMAtime)
table(pax_m2$PLASMAtime)




PLASMAtoADM <- as.numeric(difftime(pax_m2$sample_date, 
                                         pax_m2$ICU_admit_date,  units="hours"))

table(PLASMAtoADM )

PLASx <- pax_m2[,c("Patient_TK.x","MARSID","PLASMAtime")] 

PLASx  <- PLASx[!duplicated(PLASx ),]

PLASx$Patient_TK <- as.factor(PLASx$Patient_TK.x )

colnames(PLASx  )[3] <- "time"

PLASx$Patient_TK.x <- NULL

PLASx$type <- "Plasma"
PLASx$Microbe <- NA

table(PLASx$time)

ALLrp <- rbind(ALLr[ALLr$Patient_TK %in% PLASx$Patient_TK,], PLASx)

length(unique(ALLrp$MARSID))

plotall <- ggplot( data = ALLrp, aes(x=time, y=Patient_TK, col=type))+
  geom_point()



pathtofigs <- paste0("J:/Shared drives/Joe B and HPS/Bacteremia paper/Plasma protein analyses/all.svg")

svg(pathtofigs, width=8, height=22)

plotall

dev.off()

##### SWIMMING PLOTS
######################################################################
#### STREP

GRPid <- ma_ngs_admx2$Patient_TK.x

GRPid <- GRPid[GRPid %in% PLASx$Patient_TK]

#GRPid <- ma_ngs_admx$Patient_TK[ma_ngs_admx$Streptococcus=="Case"]

BCPOSx1 <- BCPOSx[BCPOSx$Patient_TK %in% GRPid,]
#BCPOSx1 <- BCPOSx1[BCPOSx1$Microbe!="Staphylococcus epidermidi",]
BCPOSx1 <- BCPOSx1[BCPOSx1$time < 80 & BCPOSx1$time > -60,]


#BCPOSx1[BCPOSx1$Microbe=="Staphylococcus epidermidi",]
#BCPOSx1[BCPOSx1$Patient_TK==38834,]

#BCPOSx1 <- BCPOSx1[grepl("trep", BCPOSx1$Microbe),]


ALL2r <- rbind(ALLrp[ALLrp$Patient_TK %in% GRPid,],  BCPOSx1  )  

ALL2r <- ALL2r[ALL2r$MARSID != 5311,]

### relevel
paxt <- ALL2r[ALL2r$type=="PAXdraw",]
paxt <- paxt[order(-paxt$time),]

ALL2r$Patient_TK <- factor(ALL2r$Patient_TK, levels = paxt$Patient_TK )

plot <- ggplot( )+
  
  annotate("rect", xmin=0, xmax=24, ymin=-Inf, ymax=Inf, alpha=0.3, fill="gray") +
  annotate("rect", xmin=0, xmax=-24, ymin=-Inf, ymax=Inf, alpha=0.2, fill="gray") +
  annotate("rect", xmin=24, xmax=48, ymin=-Inf, ymax=Inf, alpha=0.2, fill="gray") +
  annotate("rect", xmin=-24, xmax=-48, ymin=-Inf, ymax=Inf, alpha=0.1, fill="gray") +
  annotate("rect", xmin=72, xmax=48, ymin=-Inf, ymax=Inf, alpha=0.1, fill="gray") +
  
  geom_segment( data = ALL2r[ALL2r$type=="ADM",] , col="gray",
                aes(x=time, y=Patient_TK, xend=100, yend=Patient_TK )) +
  
  geom_point(data = ALL2r[ALL2r$type=="PAXdraw",], size=3, alpha=0.5,
             aes(x=time, y=Patient_TK, shape=type, col=Microbe))+
  
  geom_point(data = ALL2r[ALL2r$type=="BC",], size=2,
             aes(x=time, y=Patient_TK, shape=type, col=Microbe))+
  
  geom_point(data = ALL2r[ALL2r$type=="Plasma",], size=2,
             aes(x=time, y=Patient_TK, shape=type, col=Microbe))+
  
  theme_classic()+
  coord_cartesian(xlim=c(-60,80)) +
  scale_x_continuous(name="Hours", 
                     breaks=c(-48,-24,0,24,48,72))+
  scale_shape_manual(values = c(16, 17, 4))+
  ggtitle("All Plasma samples")

plot

pathtofigs <- paste0("J:/Shared drives/Joe B and HPS/Bacteremia paper/Plasma protein analyses/all2.svg")

svg(pathtofigs, width=18, height=30)

plot

dev.off()


###

head(plasma[!plasma$MARSID %in% ALLr$MARSID,])


load("J:/My Drive/LORNA/RawData/COUNTSMX.Rdata")

load("J:/My Drive/SPRINTAR/FINAL2/Rscripts/Flowchart/MA_NGS_ADMX.Rdata")



in243 <- plasma$MARSID[plasma$MARSID %in% ALLr$MARSID]
