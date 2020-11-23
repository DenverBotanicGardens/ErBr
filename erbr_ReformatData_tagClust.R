## April Goebl
## Script modified 07-08-20
## Collaboration with Denver Botanic Gardens
## Re-format Eriogonum brandegeei data for using with IPMs
## Cluster individuals by related (same truncated number) tag ID



rm(list=ls())
dev.off()



## SET WD -----------------------------------------------------------------------------------------
setwd("C:/Users/april/Dropbox/CU_Boulder_PhD/DBG_Internship")
## ------------------------------------------------------------------------------------------------



## LOAD PACKAGES AND FUNCTIONS --------------------------------------------------------------------
library(dplyr)
## ------------------------------------------------------------------------------------------------



## LOAD DATA --------------------------------------------------------------------------------------
erbr <- read.csv("reportsdataformatprocedures/rawdata_2018_1.csv", header=TRUE)
clim3seas <- read.csv("erbr_climData3seasons.csv", header=TRUE)
## ------------------------------------------------------------------------------------------------





## MODIFY FORM OF DATA AND CLUSTER BY RELATED/ TRUNCATED TAG ---------------------------------------------------
erbr <- erbr[!duplicated(erbr),]              #Remove duplicate rows


## Remove Cleora sites since only 1 year of data collected 
erbr.1 <- erbr[erbr$Site!="Cleora",]


## Remove years 2013 onwards (for analyses using single year transitions only) 
erbr.1 <- erbr.1[which(erbr.1$Year <= 2013),]


## Change NAs in Rosettes and Infl columns (no data entered) to zero to indicate dead, missing, or subsumed
erbr.1$Rosettes[is.na(erbr.1$Rosettes)] <- 0
erbr.1$Infl[is.na(erbr.1$Infl)] <- 0




## Add a modified Tag column to hold truncated tag values only (for clustered analysis)
erbr.1$TagClust <- NA
erbr.1$TagClust[which(erbr.1$Site=="Garden Park East")] <- paste("E.", trunc(erbr.1$Tag[which(erbr.1$Site=="Garden Park East")], sep=''))
erbr.1$TagClust[which(erbr.1$Site=="Garden Park West")] <- paste("W.", trunc(erbr.1$Tag[which(erbr.1$Site=="Garden Park West")], sep=''))


## Combine size (Rosettes) and repro (Infl) for clusters of plts with same truncated tag number
erbr.2 <- erbr.1 %>% group_by(TagClust, Year) %>% mutate(RosClust=sum(Rosettes), InflClust=sum(Infl)) %>% ungroup()


## Remove rows that are duplicates in terms of TagClust and Year values
erbr.2 <- erbr.2[!duplicated(erbr.2[,c("TagClust","Year")]),]




## Add in climate variables
clim3seas.names <- c("PptFall","PptWinter","PptSummer","TempFall","TempWinter","TempSummer",
                     "TmaxFall","TmaxWinter","TmaxSummer","TminFall","TminWinter","TminSummer")
erbr.2[, clim3seas.names] <- NA

for (cc in 1:nrow(clim3seas)) {
  erbr.2$PptFall[which(erbr.2$Year==clim3seas$Year[cc])] <- clim3seas$Tot_fall_ppt[cc]
  erbr.2$PptWinter[which(erbr.2$Year==clim3seas$Year[cc])] <- clim3seas$Tot_winter_ppt[cc]
  erbr.2$PptSummer[which(erbr.2$Year==clim3seas$Year[cc])] <- clim3seas$Tot_summer_ppt[cc]
  
  erbr.2$TempFall[which(erbr.2$Year==clim3seas$Year[cc])] <- clim3seas$Mean_fall_temp[cc]
  erbr.2$TempWinter[which(erbr.2$Year==clim3seas$Year[cc])] <- clim3seas$Mean_winter_temp[cc]
  erbr.2$TempSummer[which(erbr.2$Year==clim3seas$Year[cc])] <- clim3seas$Mean_summer_temp[cc]
  
  erbr.2$TmaxFall[which(erbr.2$Year==clim3seas$Year[cc])] <- clim3seas$Mean_fall_tmax[cc]
  erbr.2$TmaxWinter[which(erbr.2$Year==clim3seas$Year[cc])] <- clim3seas$Mean_winter_tmax[cc]
  erbr.2$TmaxSummer[which(erbr.2$Year==clim3seas$Year[cc])] <- clim3seas$Mean_summer_tmax[cc]
  
  erbr.2$TminFall[which(erbr.2$Year==clim3seas$Year[cc])] <- clim3seas$Mean_fall_tmin[cc]
  erbr.2$TminWinter[which(erbr.2$Year==clim3seas$Year[cc])] <- clim3seas$Mean_winter_tmin[cc]
  erbr.2$TminSummer[which(erbr.2$Year==clim3seas$Year[cc])] <- clim3seas$Mean_summer_tmin[cc]
}



## Copy relevant columns so they are offset to serve as t+1 data
erbr.2 <- erbr.2 %>% mutate(TagClust1=lead(TagClust))  #For tag
erbr.2 <- erbr.2 %>% mutate(RosClust1=lead(RosClust))  #For ros number

## For climate
clim3seas.names.t1 <- c("PptFall1","PptWinter1","PptSummer1","TempFall1","TempWinter1","TempSummer1",
                        "TmaxFall1","TmaxWinter1","TmaxSummer1","TminFall1","TminWinter1","TminSummer1")
erbr.2[, clim3seas.names.t1] <- NA

erbr.2 <- erbr.2 %>% mutate(PptFall1=lead(PptFall), PptWinter1=lead(PptWinter), PptSummer1=lead(PptSummer),
                            TempFall1=lead(TempFall), TempWinter1=lead(TempWinter), TempSummer1=lead(TempSummer),
                            TmaxFall1=lead(TmaxFall), TmaxWinter1=lead(TmaxWinter), TmaxSummer1=lead(TmaxSummer),
                            TminFall1=lead(TminFall), TminWinter1=lead(TminWinter), TminSummer1=lead(TminSummer))
  



## Remove rows with mismatch truncated tags in t and t+1 (ie. transition b/w diff indivs)
## Note this removes row with last year of data for each plt
erbr.3 <- erbr.2[which(erbr.2$TagClust == erbr.2$TagClust1),] 



## Remove rows with zero as size data (ie. plt not present (dead or missing); empty row).
erbr.4 <- erbr.3[which(erbr.3$RosClust > 0),]


## Change size in t+1 to NA instead of zero if missing and not dead
erbr.4$RosClust1.new <- erbr.4$RosClust1
erbr.4$RosClust1.temp <- NA
## Set temporary column to zero if consecutive rows have same tag (no death)
erbr.4$RosClust1.temp[erbr.4$TagClust[1:(nrow(erbr.4)-1)] == erbr.4$TagClust[2:nrow(erbr.4)]] <- 0
## Set new Ros1 col to NA if zero size but not dead (ie. plt was missing)
erbr.4$RosClust1.new[erbr.4$RosClust1[1:nrow(erbr.4)] == erbr.4$RosClust1.temp[1:nrow(erbr.4)]] <- NA




## Make transect unique (ie. transect, site combo); to be used as predictor variable
erbr.4$TransectNew <- NA
erbr.4$TransectNew[which(erbr.4$Site=="Garden Park East")] <- paste("E.",erbr.4$Transect[which(erbr.4$Site=="Garden Park East")], sep='')
erbr.4$TransectNew[which(erbr.4$Site=="Garden Park West")] <- paste("W.",erbr.4$Transect[which(erbr.4$Site=="Garden Park West")], sep='')




## Add in necessary response variables
## Determine whether or not reproduction occurred (to be used for prob of repro response)
erbr.4$InflClustYesNo <- NA
erbr.4$InflClustYesNo[erbr.4$InflClust > 0] <- 1
erbr.4$InflClustYesNo[erbr.4$InflClust == 0] <- 0


## Calculate survival (if you have non-zero ros in t+1, then you survived)
erbr.4$SurvClust1 <- 1 
erbr.4$SurvClust1[erbr.4$RosClust1.new == 0] <- 0
## ------------------------------------------------------------------------------------------------





## KEEP RELEVANT COLUMNS ONLY----------------------------------------------------------------------
#erbr.5 <- erbr.4 %>% select(c(Year, TransectNew, TagClust, RosClust, RosClust1.new, InflClust, InflClustYesNo, SurvClust1,
#                              Ppt, Temp, Tmax, Tmin, PptFall, PptWinter, PptSummer, TempFall, TempWinter, TempSummer,
#                              TmaxFall, TmaxWinter, TmaxSummer, TminFall, TminWinter, TminSummer,PptFall1, PptWinter1, 
#                              PptSummer1, TempFall1, TempWinter1, TempSummer1,TmaxFall1, TmaxWinter1, TmaxSummer1, TminFall1, 
#                              TminWinter1, TminSummer1))

erbr.5 <- erbr.4 %>% dplyr::select(!c(Transect, Tag, X, Y, DiameterX, DiameterY, Rosettes, Infl, Rust, BrType, InflBr, Comments, TagClust1, 
                               RosClust1, RosClust1.temp))
## ------------------------------------------------------------------------------------------------




## SAVE FORMATTED RAW DATA
write.csv(erbr.5, "erbr_tagClust_200708.csv", row.names=FALSE)
## ------------------------------------------------------------------------------------------------

