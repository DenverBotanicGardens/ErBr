## April Goebl & Dan Doak
## Script modified 20-05-03
## Collaboration with Denver Botanic Gardens
## Re-format Eriogonum brandegeei data for use with JAGS and IPMs
## Either CLUSTER individuals by related (same truncated number) tag ID
## Or keep individuals SEPARATE by unique tag ID (i.e. truncated and decimal tag IDs)
## Either keep all data years and rows with missing/ no data 
## Or make pruned year datasets
## Add survival for current year
## Add lag columns for modeling with jags
## Add in column with probability of reproducing (yes or no)
## Michelle DePrenger-Levin 2021-06-04 using 2020 data repeat April's steps

rm(list=ls())
graphics.off()


## SET WD -----------------------------------------------------------------------------------------
# setwd("C:/Users/april/Dropbox/CU_Boulder_PhD/DBG_Internship")
setwd("C:/Users/DePrengm/Denver Botanic Gardens/Conservation - General/AllProjectsBySpecies/Eriogonum brandegeei/2020_Eriogonum-brandegeei_AprilGoebl_PVA/2021_Eriogonum-brandegeei_integratedpopulationmodels")
## ------------------------------------------------------------------------------------------------


## LOAD PACKAGES AND FUNCTIONS --------------------------------------------------------------------
library(dplyr)
library(bazar)
## ------------------------------------------------------------------------------------------------


## LOAD DATA --------------------------------------------------------------------------------------
# erbr <- read.csv("Files_from_Michelle/rawdata_2018_1.csv", header=TRUE)
erbr <- read.csv("Eriogonum-brandegeei_rawdata_2004-2020.csv")
clim3seas <- read.csv("erbr_climData3seasons_201022.csv", header=TRUE)

## ------------------------------------------------------------------------------------------------

## MODIFY FORM OF DATA ----------------------------------------------------------------------------
# nrow(erbr)
# nrow(erbr[!duplicated(erbr),])
# duptags <- erbr$Tag[duplicated(erbr)] # 259 is not ErBr, is ErJa so remove
# erbr[erbr$Tag == 128,] # repeated
erbr <- erbr[erbr$Tag != 259,] # E. jamesii
erbr <- erbr[!duplicated(erbr),]    #Remove duplicate rows



## Remove Cleora sites since only 1 year of data collected 
erbr.1 <- erbr[erbr$Site!="Cleora",]
table(erbr.1$Site, erbr.1$Year)

## Change zeros in Rosettes (no data entered) to NAs in Ros and Infl columns to indicate dead, missing, or subsumed
erbr.1$Rosettes[erbr.1$Rosettes==0] <- NA
erbr.1$Infl[erbr.1$Rosettes==0] <- NA
erbr.1$Infl[is.na(erbr.1$Rosettes)] <- NA
## ------------------------------------------------------------------------------------------------

# It appears that in 2020, if no new information was added, the row isn't output. 
erbr.1[erbr.1$Year == 2020,]



## OPTIONAL *******************************************************************
## FOR MAKING CONSECUTIVE-ONLY OR PRUNED YEAR DATASETS 
## Remove years 2013 onwards 
#erbr.1 <- erbr.1[which(erbr.1$Year <= 2013),]


## PRUNE DATA TO HAVE EVERY-OTHER YEARS ONLY 
#yrs <- unique(erbr.1$Year)

## To keep even years
#yrs.even <- yrs[lapply(yrs, "%%",2)==0] #Keep values where the remainder after division by 2 is zero
#erbr.1 <- erbr.1[which(erbr.1$Year %in% yrs.even),]

## To keep odd years
#yrs.odd <- yrs[lapply(yrs, "%%",2)!=0] 
#erbr.1 <- erbr.1[which(erbr.1$Year %in% yrs.odd),]


## PRUNE DATA TO HAVE SHORT CONSECUTIVE YEARS ONLY
## To keep 1st 5 yeras
#erbr.1 <- erbr.1[which(erbr.1$Year <=2008),]

## To keep 2nd 5 years
#erbr.1 <- erbr.1[which(erbr.1$Year >=2009),]
### ***************************************************************************




## CLUSTER BY RELATED TAG OR KEEP UNIQUE TAGS SEPARATE -----------------------------------

## RUN EITHER THIS ... 
## For tag CLUST
## Add a modified Tag column (site, transect, tag) to hold truncated tag values only (for clustered analysis)
erbr.1$TagNew <- NA
erbr.1$TagNew[which(erbr.1$Site=="Garden Park East")] <- paste("E",erbr.1$Transect[which(erbr.1$Site=="Garden Park East")],sep=".",trunc(erbr.1$Tag[which(erbr.1$Site=="Garden Park East")]))
erbr.1$TagNew[which(erbr.1$Site=="Garden Park West")] <- paste("W",erbr.1$Transect[which(erbr.1$Site=="Garden Park West")],sep=".",trunc(erbr.1$Tag[which(erbr.1$Site=="Garden Park West")]))

## Combine size (Rosettes) and repro (Infl) for clusters of plts with same truncated tag number
erbr.1 <- erbr.1 %>% group_by(TagNew, Year) %>% mutate(RosNew=sumNA(Rosettes,na.rm=TRUE), InflNew=sumNA(Infl,na.rm=TRUE)) %>% ungroup()

## Remove rows that are duplicates in terms of TagNew and Year values
erbr.1 <- erbr.1[!duplicated(erbr.1[,c("TagNew","Year")]),]


## Skip this 2021-08-24 to test if that is the difference and problem with the new data - to make this df match 
# "C:/Users/DePrengm/Denver Botanic Gardens/Conservation - General/AllProjectsBySpecies/Eriogonum brandegeei/+
# 2020_Eriogonum-brandegeei_AprilGoebl_PVA/erbr_TagClust_210510.csv"

#######################################################
## ... OR THIS
## For tag SEP
## Add a new Tag column with unique tag ID (i.e. incorporates transect and tag) 
erbr.1$TagNewAll <- NA
erbr.1$TagNewAll[which(erbr.1$Site=="Garden Park East")] <- paste("E",erbr.1$Transect[which(erbr.1$Site=="Garden Park East")],sep=".",erbr.1$Tag[which(erbr.1$Site=="Garden Park East")])
erbr.1$TagNewAll[which(erbr.1$Site=="Garden Park West")] <- paste("W",erbr.1$Transect[which(erbr.1$Site=="Garden Park West")],sep=".",erbr.1$Tag[which(erbr.1$Site=="Garden Park West")])

## Rename variables for consistency with above
erbr.1 <- rename(erbr.1, RosNew=Rosettes, InflNew=Infl)
## -----------------------------------------------------------------------------------
######################################################

erbr.1[erbr.1$Rosettes == erbr.1$RosNew & !is.na(erbr.1$Rosettes),]


## START OF CHANGES NEEDED FOR JAGs --------------------------------------------------
dats <- erbr.1

## Make a col that indicates if a line should be kept: removes lines w NAs before a plt appeared or after died. Keeps lines that are NAs but bracketed by sz data yrs
dats$save <- 0   #Start with save = 0, or no
dats$surv <- NA  #Column to show plt survival/ if plant is alive in current year

# Testing not truncated ones
# tags <- unique(dats$TagNewAll) 
tags <- unique(dats$TagNew)
for (tt in tags){
  szlines <- which(dats$TagNew==tt)
  # szlines <- which(dats$TagNewAll==tt) 
  
  szs <- dats$RosNew[szlines]   
  goodszslines <- szlines[is.na(szs)==FALSE] 
  

  badszlines <- szlines[is.na(szs)==TRUE] #For determining row representing 1st year dead
  badyrs <- dats$Year[badszlines]         
  goodyrs <- dats$Year[goodszslines]  
  

  if (length(goodszslines)>0){
    mingoodlines <- min(goodszslines)
    maxgoodlines <- max(goodszslines)
    dats$save[mingoodlines:maxgoodlines]=1}
  
  ## If statement that keeps row of data representing 1st year dead
  if (length(badyrs)>0 && length(goodyrs)>0 && max(badyrs,na.rm=TRUE) > max(goodyrs,na.rm=TRUE)) {
    dats$save[maxgoodlines+1] <- 1
    dats$surv[maxgoodlines+1] <- 0  #Change survival to zero for 1st year dead
  }
}

dats <- dats[dats$save==1,]        #Remove NA rows that are not in middle of the data
dats$surv[dats$RosNew>0] <- 1      #Change survival/ alive to 1 if plant is non-zero size

head(dats)
dats[is.na(dats$Rosettes),]

## Now add in rows for missing yrs in the set of data for each plt, & make a variable that indicates if the row's sz is a dependent variable & how far back to get to the last sz
dats$lagsrtsz <- -1    #This is a variable that will indicate if the row sz is a dependent value (0=no) or what the lag in time is back to the last observed sz (1,2,3 etc)
dats$lagforsurv <- -1  #Another lag variable that givens values out to the final yr of non-survival, for plts that died (i.e. if died, when was most recent sz measure?)
tags <- unique(dats$TagNew)
# tags <- unique(dats$TagNewAll) 
dats2 <- NULL         #Placeholder for the new data

# Test
# tt <- tags[1]

for (tt in tags){
  dds <- dats[which(dats$TagNew==tt),] #Temporary data
  # dds <- dats[which(dats$TagNewAll==tt),] #Temporary data
  if (length(dds$Year)>1){
   for (yy in 2:length(dds$Year)) {
      pastyrs <- dds$Year[1:(yy-1)]
      goodpastyrs <- pastyrs[is.na(dds$RosNew[1:(yy-1)])==FALSE] 
      if (is.na(dds$RosNew[yy])==FALSE) {  
        dds$lagsrtsz[yy] <- min(dds$Year[yy] - goodpastyrs)
        dds$lagforsurv[yy] <- min(dds$Year[yy] - goodpastyrs)  #lagforsurv has the same values as lagsrtsz for non-death years
      }
      ## if statement to add years since last measure for death years
      if (!is.na(dds$surv[yy]) && dds$surv[yy]==0) {
        dds$lagforsurv[yy] <- min(dds$Year[yy] - goodpastyrs) 
      }
   } # end yr loop

        ## Find and add in the missing year rows:
        allyrs <- min(dds$Year):max(dds$Year)
        yrs <- c(dds$Year) # years with a census
        missingyrs <- allyrs[which(allyrs%in%yrs ==FALSE)]
        ddsmissing <- do.call('rbind',replicate(length(missingyrs),dds[1,],simplify=FALSE))
        ddsmissing$Year <- missingyrs
        ddsmissing$X <- ddsmissing$Y <- ddsmissing$DiameterX <- ddsmissing$DiameterY <- 
          ddsmissing$Rosettes <- ddsmissing$Infl <- ddsmissing$RosNew <- ddsmissing$InflNew <- ddsmissing$Rust <- 
          ddsmissing$InflBr <- ddsmissing$Comments <- ddsmissing$surv <- ddsmissing$BrType <- NA
        
        ddsmissing$lagforsurv <- 0
        ddsmissing$lagsrtsz <- 0
        
        dds <- rbind(dds,ddsmissing)
        dds <- dds[order(dds$Year),] #Reordered, full record for this plt
      
      } #End if the plt was observed more than once
    dats2 <- rbind(dats2,dds)
  } #End going through each plt



# Now has NA with other info in missing years 
# dats2[dats2$TagNewAll == tags[2],]
dats2[dats2$TagNew == tags[2],]
erbr.2 <- dats2 
## -----------------------------------------------------------------------------------



## ADD IN CLIMATE VARIABLES ----------------------------------------------------------
clim3seas.names <- c("PptFall","PptWinter","PptSummer","TempFall","TempWinter","TempSummer")
                   
erbr.2[, clim3seas.names] <- NA

for (cc in 1:nrow(clim3seas)) {
  erbr.2$PptFall[which(erbr.2$Year==clim3seas$Year[cc])] <- clim3seas$Tot_fall_ppt[cc]
  erbr.2$PptWinter[which(erbr.2$Year==clim3seas$Year[cc])] <- clim3seas$Tot_winter_ppt[cc]
  erbr.2$PptSummer[which(erbr.2$Year==clim3seas$Year[cc])] <- clim3seas$Tot_summer_ppt[cc]
  
  erbr.2$TempFall[which(erbr.2$Year==clim3seas$Year[cc])] <- clim3seas$Mean_fall_temp[cc]
  erbr.2$TempWinter[which(erbr.2$Year==clim3seas$Year[cc])] <- clim3seas$Mean_winter_temp[cc]
  erbr.2$TempSummer[which(erbr.2$Year==clim3seas$Year[cc])] <- clim3seas$Mean_summer_temp[cc]
}
## -----------------------------------------------------------------------------------


## MAKE TRANSECT UNIQUE (ie. transect, site combo) FOR USE AS PREDICTOR VARIABLE -----
erbr.2$TransectNew <- NA
erbr.2$TransectNew[which(erbr.2$Site=="Garden Park East")] <- paste("E.",erbr.2$Transect[which(erbr.2$Site=="Garden Park East")], sep='')
erbr.2$TransectNew[which(erbr.2$Site=="Garden Park West")] <- paste("W.",erbr.2$Transect[which(erbr.2$Site=="Garden Park West")], sep='')
## -----------------------------------------------------------------------------------


## ADD PROBABILITY OF REPRO RESPONSE VARIABLE ----------------------------------------
## Determine whether or not reproduction occurred 
erbr.2$InflYesNo <- NA 
erbr.2$InflYesNo[erbr.2$InflNew > 0] <- 1
erbr.2$InflYesNo[erbr.2$InflNew == 0] <- 0
## -----------------------------------------------------------------------------------


## KEEP RELEVANT COLUMNS ONLY --------------------------------------------------------
erbr.3 <- erbr.2 %>% dplyr::select(!c(Transect, Tag, X, Y, DiameterX, DiameterY, Rust, BrType, InflBr, Comments, save)) 
## -----------------------------------------------------------------------------------


## SAVE FORMATTED DATA ---------------------------------------------------------------
# date <- as.character(210510)          #Enter date to be added to file name
date <- as.character(Sys.Date())
name <- as.character("TagClust_") #Enter name of file, e.g. Tagclust, 4to13, 4to13odd, 4to13even, 4to8, 9to13
# name <- as.character("TagNotClust_")

write.csv(erbr.3, file=paste("erbr_", name, date, ".csv", sep=""), row.names=FALSE)
## -----------------------------------------------------------------------------------



