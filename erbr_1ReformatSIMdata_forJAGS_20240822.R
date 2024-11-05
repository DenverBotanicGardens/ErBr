## April Goebl & Dan Doak
## Script modified 2024-08-24
## Eriogonum brandegeei collaboration with Denver Botanic Gardens
## Re-format simulated data for use with JAGS 
## Add lag columns for modeling with jags


rm(list=ls())
graphics.off()


## SET WD -----------------------------------------------------------------------------------------
#setwd("C:/Users/april/Dropbox/CU_Boulder_PhD/DBG_Internship")
## ------------------------------------------------------------------------------------------------


## LOAD PACKAGES AND FUNCTIONS --------------------------------------------------------------------
library(dplyr)
library(bazar)
library(stringr) 
library(tidyr)   
## ------------------------------------------------------------------------------------------------


## LOAD DATA --------------------------------------------------------------------------------------
erbr <- read.csv("20240904_erbr_SimDat20yrR3.csv", header=TRUE)
clim3seas <- read.csv("erbr_climData3seas32yr_221114.csv", header=TRUE)
## ------------------------------------------------------------------------------------------------


## MODIFY FORM OF DATA ----------------------------------------------------------------------------
## Change zeros in Rosettes & Infls to NAs in RosNew & InflNew columns to indicate dead
erbr$RosNew[erbr$RosNew==0] <- NA
erbr[erbr$RosNew==0 & !is.na(erbr$RosNew),] #Confirm that no RosNew=0
erbr$InflNew[erbr$InflNew==0] <- NA
erbr$InflNew[is.na(erbr$RosNew) & !is.na(erbr$InflNew)] #Confirm
erbr$InflNew[is.na(erbr$RosNew)] <- NA
## Change Infl to zero from NA if Rosettes has data 
erbr$InflNew[erbr$RosNew>0 & !is.na(erbr$RosNew) & is.na(erbr$InflNew)] <- 0
## ------------------------------------------------------------------------------------------------





## Confirm no rows are duplicates in terms of TagNew and Year values
erbr[duplicated(erbr[,c("TagNew","Year")]),]




## How many observed individuals (tag clusters) were there each year
indivXyear <- erbr %>%
  dplyr::group_by(Year) %>%
  dplyr::summarise(Indivs = n_distinct(TagNew[RosNew > 0]))
## -----------------------------------------------------------------------------------






## START OF CHANGES NEEDED FOR JAGS --------------------------------------------------

## Make a column that indicates if a line should be kept: removes lines w NAs before a plt appeared or after died. Keeps lines that are NAs but bracketed by sz data yrs
dats <- erbr      #Placeholder for the new data
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
lagCheck <- dats2[dats2$lagsrtsz>3,]

erbr.1 <- dats2
## -----------------------------------------------------------------------------------



## ADD IN CLIMATE VARIABLES ----------------------------------------------------------
erbr.1 <- erbr.1 %>%
  left_join(clim3seas, by = c("ClimYr" = "Year"))
## -----------------------------------------------------------------------------------



## ADD PROBABILITY OF REPRO RESPONSE VARIABLE ----------------------------------------
## Determine whether or not reproduction occurred
erbr.1$InflYesNo <- NA
erbr.1$InflYesNo[erbr.1$InflNew > 0] <- 1
erbr.1$InflYesNo[erbr.1$InflNew == 0] <- 0
## -----------------------------------------------------------------------------------



## SAVE FORMATTED DATA ---------------------------------------------------------------
date <- Sys.Date()                             #Enter date to be added to file name
date <- str_replace_all(date, "-", "")
name <- as.character("SimDat20yrR3")           #Enter name of file, e.g. Tagclust, 4to13, simulated data...

write.csv(erbr.1, file=paste(date, "_erbr_", name, "_Format4JAGS", ".csv", sep=""), row.names=FALSE)
## -----------------------------------------------------------------------------------



