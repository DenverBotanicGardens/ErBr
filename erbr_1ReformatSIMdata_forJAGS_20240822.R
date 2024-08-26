## April Goebl & Dan Doak
## Script modified 2024-08-24
## Eriogonum brandegeei collaboration with Denver Botanic Gardens
## Re-format simulated data for use with JAGS 
## Add lag columns for modeling with jags


rm(list=ls())
graphics.off()


## SET WD -----------------------------------------------------------------------------------------
setwd("C:/Users/april/Dropbox/CU_Boulder_PhD/DBG_Internship")
## ------------------------------------------------------------------------------------------------


## LOAD PACKAGES AND FUNCTIONS --------------------------------------------------------------------
library(dplyr)
library(bazar)
library(stringr) 
library(tidyr)   
## ------------------------------------------------------------------------------------------------


## LOAD DATA --------------------------------------------------------------------------------------
erbr <- read.csv("20240826_SimData.csv", header=TRUE)
clim3seas <- read.csv("erbr_climData3seas32yr_221114.csv", header=TRUE)
## ------------------------------------------------------------------------------------------------


## MODIFY FORM OF DATA ----------------------------------------------------------------------------
## Change zeros in Rosettes & Infls to NAs in RosNew & InflNew columns to indicate dead
erbr$RosNew[erbr$RosNew==0] <- NA
erbr[erbr$RosNew==0 & !is.na(erbr$RosNew),] #Confirm that no RosNew=0
erbr$InflNew[erbr$InflNew==0] <- NA
erbr$InflNew[is.na(erbr$RosNew) & !is.na(erbr$InflNew)] #Confirm
## Change Infl to zero from NA if Rosettes has data 
erbr$InflNew[erbr$RosNew>0 & !is.na(erbr$RosNew) & is.na(erbr$InflNew)] <- 0
## ------------------------------------------------------------------------------------------------




## Confirm no rows are duplicates in terms of TagNew and Year values
erbr[duplicated(erbr[,c("TagNew","Year")]),]

## ********************** ## 

# table(erbr.1$RosNew, useNA = "always") #Checks
# nrow(erbr.1) #Checks




## How many observed individuals (tag clusters) were there each year
#indivXyear <- erbr.1 %>%
#  dplyr::group_by(Year) %>%
#  dplyr::summarise(Indivs = n_distinct(TagNew[RosNew > 0]))
#as.data.frame(indivXyear) ## Only 4 transects were read at GPE instead of 7 in 2020. Totals are wrong. One year of NA?

### Save data for number and size of clusters in predicted years (after 2013)
# save(erbr.1, file = paste("C:/Users/deprengm/OneDrive - Denver Botanic Gardens/P drive/hackathon/ErBr/erbr.1",
#                          Sys.Date(),".Rdata", sep=""))


## MEDL ADDED CODE ENDS --------------------------------------------------------------





## ADD IN NA ROWS FOR 2020 AND 2022 IF MISSING --------------------------------------- 
tags <- unique(erbr.1$TagNew)
prevYr <- 2018    #All data collection years have rows in data sheet up until this year     
latestYr <- 2022  #Most recent year of data collection
dats <- NULL      #Placeholder for the new data

for (tt in tags){
    dds <- erbr.1[which(erbr.1$TagNew==tt),] #Temporary data
    
    if (length(dds$Year)>1){
      
      ## Find and add in the missing recent year (no data collected) rows:
      recentyrs <- prevYr:latestYr
      yrs <- c(dds$Year)
      missingyrs <- recentyrs[which(recentyrs%in%yrs ==FALSE)]
      missingyrs <- missingyrs[missingyrs==2020 | missingyrs==2022] #Keep only if 2020 and 2022 missing
      ddsmissing <- do.call('rbind',replicate(length(missingyrs),dds[1,],simplify=FALSE))
      ddsmissing$Year <- missingyrs
      ddsmissing$DiameterX=ddsmissing$DiameterY=
      ddsmissing$RosNew=ddsmissing$InflNew=ddsmissing$Rosettes=ddsmissing$Infl=
      ddsmissing$Rust=ddsmissing$BrType=ddsmissing$InflBr=ddsmissing$Comments=NA
      dds <- rbind(dds,ddsmissing)
      dds <- dds[order(dds$Year),] #Reordered, full record for this plt
        } #End if the plt was observed more than once
    
    dats <- rbind(dats,dds)
    } #End going through each plt


## Check if desired rows added
nrow(erbr.1[erbr.1$Year==2004,])
nrow(erbr.1[erbr.1$Year==2010,])
nrow(erbr.1[erbr.1$Year==2012,])
nrow(erbr.1[erbr.1$Year==2013,])
nrow(erbr.1[erbr.1$Year==2016,])
nrow(erbr.1[erbr.1$Year==2018,])
nrow(erbr.1[erbr.1$Year==2020,])
nrow(erbr.1[erbr.1$Year==2022,])

nrow(dats[dats$Year==2013,])
nrow(dats[dats$Year==2016,])
nrow(dats[dats$Year==2018,])
nrow(dats[dats$Year==2020,])
nrow(dats[dats$Year==2022,])
## -----------------------------------------------------------------------------------






## START OF CHANGES NEEDED FOR JAGS --------------------------------------------------

## Make a column that indicates if a line should be kept: removes lines w NAs before a plt appeared or after died. Keeps lines that are NAs but bracketed by sz data yrs
dats$save <- 0   #Start with save = 0, or no
dats$surv <- NA  #Column to show plt survival/ if plant is alive in current year

tags <- unique(dats$TagNew)
for (tt in tags){
  szlines <- which(dats$TagNew==tt) ## index of a tag
  szs <- dats$RosNew[szlines]       ## sizes as number of rosettes each year
  goodszslines <- szlines[is.na(szs)==FALSE] ## index of years with number of rosettes counted

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


## Now add in rows for missing yrs in the set of data for each plt, & make a variable that indicates if the row's sz is a dependent variable & how far back to get to the last sz
dats$lagsrtsz <- 0    #This is a variable that will indicate if the row sz is a dependent value (0=no) or what the lag in time is back to the last observed sz (1,2,3 etc)
dats$lagforsurv <- 0  #Another lag variable that givens values out to the final yr of non-survival, for plts that died (i.e. if died, when was most recent sz measure?)
tags <- unique(dats$TagNew)
dats2 <- NULL         #Placeholder for the new data



for (tt in tags){
  dds <- dats[which(dats$TagNew==tt),] #Temporary data
 # print(dds$Year)
  if (length(dds$Year)>1){
     for (yy in 2:length(dds$Year)) {
        pastyrs <- dds$Year[1:(yy-1)]
        goodpastyrs <- pastyrs[is.na(dds$RosNew[1:(yy-1)])==FALSE]
        if (is.na(dds$RosNew[yy])==FALSE) {
          dds$lagsrtsz[yy] <- min(dds$Year[yy] - goodpastyrs)    #MEDL: when there have been no good past years, returns Inf
          dds$lagforsurv[yy] <- min(dds$Year[yy] - goodpastyrs)  #lagforsurv has the same values as lagsrtsz for non-death years
        if(is.infinite(dds$lagforsurv[[yy]])) print(tt)
          }
        ## if statement to add years since last measure for death years
        if (!is.na(dds$surv[yy]) && dds$surv[yy]==0) {
          dds$lagforsurv[yy] <- min(dds$Year[yy] - goodpastyrs)
        }
        # print(dds$lagforsurv)    ##MEDL: e.g. test shows W.1.1 has Inf lag in year 4
     } # end yr loop

    ## Find and add in the missing year rows:
    allyrs <- min(dds$Year):max(dds$Year)
    yrs <- c(dds$Year)
    missingyrs <- allyrs[which(allyrs%in%yrs ==FALSE)]
    ddsmissing <- do.call('rbind',replicate(length(missingyrs),dds[1,],simplify=FALSE))
    ddsmissing$Year <- missingyrs
    ddsmissing$X=ddsmissing$Y=ddsmissing$DiameterX=ddsmissing$DiameterY=
    ddsmissing$RosNew=ddsmissing$InflNew=
      ddsmissing$Rust=ddsmissing$BrType=
      NA
    ddsmissing$InflBr=ddsmissing$Comments=ddsmissing$surv=NA
    ddsmissing$lagsrtsz <- 0
    dds <- rbind(dds,ddsmissing)
    dds <- dds[order(dds$Year),] #Reordered, full record for this plt

  } #End if the plt was observed more than once
  
  dats2 <- rbind(dats2,dds)
} #End going through each plt

## Check lag values to ensure biologically reasonable
table(dats2$lagforsurv) 
table(dats2$lagsrtsz) 
lagCheck <- dats2[dats2$lagsrtsz>5,]

erbr.2 <- dats2
## -----------------------------------------------------------------------------------



## ADD IN CLIMATE VARIABLES ----------------------------------------------------------
erbr.2 <- erbr.2 %>%
  left_join(clim3seas, by = c("Year" = "Year"))
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
erbr.3 <- erbr.2 %>% dplyr::select(!c(Transect, Tag, save)) 

### 2023.01.03 changes
# erbr.3 <- erbr.2 %>% dplyr::select(!c(Transect, # Tag, X, Y, DiameterX, DiameterY, Rust, BrType, InflBr, Comments,
#                                   save))
## -----------------------------------------------------------------------------------


## SAVE FORMATTED DATA ---------------------------------------------------------------
date <- Sys.Date() #as.character(210617)        #Enter date to be added to file name
date <- str_replace_all(date, "-", "")
name <- as.character("TagClust2022_")     #Enter name of file, e.g. Tagclust, 4to13, 4to13odd, 4to13even, 4to8, 9to13

write.csv(erbr.3, file=paste("erbr_", name, date, ".csv", sep=""), row.names=FALSE)
## -----------------------------------------------------------------------------------



