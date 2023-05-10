## April Goebl & Dan Doak
## Script modified 20-05-03
## Collaboration with Denver Botanic Gardens
## Re-format Eriogonum brandegeei data for use with JAGS and IPMs
## Either CLUSTER individuals by related (same truncated number) tag ID
## Or keep individuals SEPARATE by unique tag ID (i.e. truncated and demical tag IDs)
## Either keep all data years and rows with missing/ no data
## Or make pruned year datasets
## Add survival for current year
## Add lag columns for modeling with jags
## Add in column with probability of reproducing (yes or no)


rm(list=ls())
graphics.off()


## SET WD -----------------------------------------------------------------------------------------
setwd("C:/Users/april/Dropbox/CU_Boulder_PhD/DBG_Internship")
#username <- "deprengm" 
#setwd(paste("C:/Users/",username, 
#            "/Denver Botanic Gardens/Conservation - General/AllProjectsBySpecies/Eriogonum-brandegeei/Eriogonum-brandegeei_Data",
#            sep=""))
## ------------------------------------------------------------------------------------------------


## LOAD PACKAGES AND FUNCTIONS --------------------------------------------------------------------
library(dplyr)
library(bazar)
## 2023.01.03 additions
library(stringr) 
library(tidyr)   
## ------------------------------------------------------------------------------------------------


## LOAD DATA --------------------------------------------------------------------------------------
erbr <- read.csv("Files_from_Michelle/rawdata_2022.csv", header=TRUE)
#erbr.18 <- read.csv("Files_from_Michelle/rawdata_2018_1.csv", header=TRUE)

clim3seas <- read.csv("erbr_climData3seas32yr_221114.csv", header=TRUE) #Check that my clim data is same as Michelle's
#load("C:/Users/deprengm/OneDrive - Denver Botanic Gardens/P drive/hackathon/ErBr/erbr_climData3seas2022-10-11.Rdata")
#load("ErBr/erbr_climData3seas2022-10-11.Rdata") #This is not the same as my clim3seas. I think Michelle and I worked out that her 2022-10-11 needs updating.
## ------------------------------------------------------------------------------------------------


## MODIFY FORM OF DATA ----------------------------------------------------------------------------
erbr <- erbr[!duplicated(erbr),]    #Remove 10 duplicate rows
erbr <- erbr[erbr$Tag != 259,]      #MEDL: E. jamesii # added 2023.01.03 from 2020_erbr_1ReformatData_forJAGS_20210823.R

## Remove Cleora sites since only 1 year of data collected
erbr.1 <- erbr[erbr$Site!="Cleora",]

## Change zeros in Rosettes (no data entered) to NAs in Ros and Infl columns to indicate dead, missing, or subsumed
erbr.1$Rosettes[erbr.1$Rosettes==0] <- NA
erbr.1$Infl[erbr.1$Rosettes==0] <- NA
erbr.1$Infl[is.na(erbr.1$Rosettes)] <- NA
## Change Infl to zero from NA if Rosettes has data ***

## ------------------------------------------------------------------------------------------------


## MEDL Note: 2022-10-11 get num of clusters in recent yrs for ability of model to predict counts- skip this portion
## OPTIONAL *******************************************************************
## FOR MAKING CONSECUTIVE-ONLY OR PRUNED YEAR DATASETS
## Remove years 2013 onwards
# erbr.1 <- erbr.1[which(erbr.1$Year <= 2013),] ### commented out to get total individuals each year

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
# erbr.1 <- erbr.1[which(erbr.1$Year >=2009),] ## Commented out to get total individuals each year
### ***************************************************************************



## MDL Note: 2022-10-11 to get the number of clusters in 2020 and 2022 for fit of predictions skipped portion above
## CLUSTER BY RELATED TAG OR KEEP UNIQUE TAGS SEPARATE -----------------------------------

## RUN EITHER THIS ...
## For tag CLUST
## Add a modified Tag column (site, transect, tag) to hold truncated tag values only (for clustered analysis)
erbr.1$TagNew <- NA
erbr.1$TagNew[which(erbr.1$Site=="Garden Park East")] <- paste("E",erbr.1$Transect[which(erbr.1$Site=="Garden Park East")],sep=".",trunc(erbr.1$Tag[which(erbr.1$Site=="Garden Park East")]))
erbr.1$TagNew[which(erbr.1$Site=="Garden Park West")] <- paste("W",erbr.1$Transect[which(erbr.1$Site=="Garden Park West")],sep=".",trunc(erbr.1$Tag[which(erbr.1$Site=="Garden Park West")]))



## MEDL ADDED JAN 2023 TO ONLY CLUSTER PLTS W/IN 30cm OF TAG ----------------------------------------------------
### RUN THIS 2023.01.03 ...
### truncate only if within 30cm of original tag. Tags generally placed 10cm downhill from cluster. Due to the easily erroded soil, tags that did not wash away are sometimes used to mark new plants. These new plants can be
### MEDL 2022.12.15 realizing that some tags are used when closer tag is lost, need to limit truncation to a certain distance
## From the PHPmyadmin database, the tags, transects, and sites tables to join to get the location information
# tagsPHP <- read.csv("C:/Users/deprengm/Denver Botanic Gardens/Conservation - General/AllProjectsBySpecies/Eriogonum-brandegeei/Eriogonum-brandegeei_Projects/2020_Eriogonum-brandegeei_AprilGoebl_PVA/2021_Eriogonum-brandegeei_integratedpopulationmodels/_erbr_tags_PHP2022.csv")
# sitesPHP <- read.csv("C:/Users/deprengm/Denver Botanic Gardens/Conservation - General/AllProjectsBySpecies/Eriogonum-brandegeei/Eriogonum-brandegeei_Projects/2020_Eriogonum-brandegeei_AprilGoebl_PVA/2021_Eriogonum-brandegeei_integratedpopulationmodels/_erbr_sites_PHP2022.csv")
# transectsPHP <- read.csv("C:/Users/deprengm/Denver Botanic Gardens/Conservation - General/AllProjectsBySpecies/Eriogonum-brandegeei/Eriogonum-brandegeei_Projects/2020_Eriogonum-brandegeei_AprilGoebl_PVA/2021_Eriogonum-brandegeei_integratedpopulationmodels/_erbr_transects_PHP2022.csv")
tagsPHP <- read.csv("Files_from_Michelle/_erbr_tags_PHP2022.csv")
sitesPHP <- read.csv("Files_from_Michelle/_erbr_sites_PHP2022.csv")
transectsPHP <- read.csv("Files_from_Michelle/_erbr_transects_PHP2022.csv")

#names(tagsPHP)
#names(sitesPHP)
#names(transectsPHP)
#names(erbr.1)

## recorded as either "XXcm Dir of tag" or "Tag XXcm dir" sometimes XX cm sometimes XXcm
## soil erodes easily, tags are lost and used sparingly. Some plants are added to tags even at great distances.
## setting an arbitrary cutoff to keep what are most likely different plants separate
erbr.loc <- tagsPHP %>%
left_join(transectsPHP, by = "ErBr_transect_id") %>%
left_join(sitesPHP, by = "ErBr_site_id") %>%
left_join(erbr.1, by = c("tag_number" = "Tag", "transect_number" = "Transect", "site_name" = "Site")) %>%
dplyr::select(c(site_name ,Year,transect_number,tag_number ,X:Comments, TagNew, location_comment )) %>% ## keep all matching erbr.1 plus the location comment which has distance to tag
filter(site_name != "Cleora") %>%
filter(!is.na(Year))

## Missing some years, filter out, some are tags that no longer exist or were missing tag numbers
# head(erbr.loc[is.na(erbr.loc$Year),],100)

str(erbr.loc)
erbr.loc$location_comment <- as.character(erbr.loc$location_comment)
## by hand delete all the extra numbers that aren't distances
table(erbr.loc$location_comment, useNA = "always") ## There are no NAs; there are 3619 NULL
erbr.loc$location_comment[erbr.loc$location_comment == "x=24.9, y=0.46, location of plant (no tag), 140cm E of #698"] <- "140cm E of"
erbr.loc$location_comment[erbr.loc$location_comment == "22cm NW of tag 706"] <- "22cm NW of tag"
erbr.loc$location_comment[erbr.loc$location_comment == "downhill of nail, (plant 9cm S of tag) plant 10cm N of tag, tag 10N"] <- "downhill of nail, tag 10N"
erbr.loc$location_comment[erbr.loc$location_comment == "plant out of transect, plant at nail; plant within 50 cm of transect and at tag"] <- "plant out of transect, plant at nail; plant within"

## Add a column with the distance to tag information from the location comments 
erbr.loc1 <- erbr.loc %>%
mutate(TagDist = if_else(grepl("tag", location_comment),                 #If reference to distance of plant from tag
                           if_else(nchar(gsub("[^0-9.-]", "", location_comment)) > 0,
                           gsub("[^0-9.-]", "", location_comment), "0"), #Then pull numbers, some cases where there's the word "tag" but no distance given
                           "0")) %>%                                     #If not, the numbers do not refer to distance and are zero, NA would be fine too
       mutate(TagDist = if_else(is.na(TagDist), as.numeric(0), as.numeric(TagDist))) #Make Tag dist numeric


## Check that numbers that don't refer to distance from Tag were removed  
#erbr.loc1[which(nchar(erbr.loc1$TagDist)>4),] ## any that have more than one number or a number greater than 4 digits should catch all the locations with multiple numbers
#table(erbr.loc1$TagDist, useNA = "always") ## Distribution of distances of tags to plants (clusters) ## 290 with NA
#erbr.loc1[which(erbr.loc1$TagDist > 30),]  ##More checks
#erbr.loc1[erbr.loc1$TagDist > 30,] ## Woo hoo! (looks good)
#table(erbr.loc1$location_comment[which(erbr.loc1$TagDist > 30)]) #Looking to see that all looks good
names(erbr.loc1) <- c(names(erbr.1),"location_comment","TagDist")
erbr.1 <- erbr.loc1

### truncate only if distance is within 30cm of tag
erbr.1 <- erbr.loc1 %>%
  mutate(TagNew = case_when(Site == "Garden Park East" & TagDist < 30 ~ paste("E",Transect,trunc(Tag),sep="."),
                            Site == "Garden Park East" & TagDist >= 30 ~ paste("E",Transect,Tag,sep="."),
                            Site == "Garden Park West" & TagDist < 30 ~ paste("W",Transect,trunc(Tag),sep="."),
                            Site == "Garden Park West" & TagDist >= 30 ~ paste("W",Transect,Tag,sep=".") ))

## Checks 
head(erbr.1[which(erbr.1$TagDist > 30),],100)
head(erbr.1[erbr.1$TagDist > 30,], 100)
head(erbr.1[which(erbr.1$TagDist > 20 & erbr.1$TagDist < 35),],100)


## Combine size (Rosettes) and repro (Infl) for clusters of plts with same truncated tag number
## MEDL Note: TagNew has the Site and Transect info in it
erbr.1 <- erbr.1 %>% dplyr::group_by(TagNew, Year) %>% dplyr::mutate(RosNew=sumNA(Rosettes,na.rm=TRUE), InflNew=sumNA(Infl,na.rm=TRUE)) %>% 
          ungroup() 
## MEDL Note: bazar::sumNA returns NA instead of 0 when input contains only missing values

## Remove rows that are duplicates in terms of TagNew and Year values
erbr.1 <- erbr.1[!duplicated(erbr.1[,c("TagNew","Year")]),]
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


## ... OR RUN THIS
## For tag SEP
## Add a new Tag column with unique tag ID (i.e. incorporates transect and tag)
#erbr.1$TagNew <- NA
#erbr.1$TagNew[which(erbr.1$Site=="Garden Park East")] <- paste("E",erbr.1$Transect[which(erbr.1$Site=="Garden Park East")],sep=".",erbr.1$Tag[which(erbr.1$Site=="Garden Park East")])
#erbr.1$TagNew[which(erbr.1$Site=="Garden Park West")] <- paste("W",erbr.1$Transect[which(erbr.1$Site=="Garden Park West")],sep=".",erbr.1$Tag[which(erbr.1$Site=="Garden Park West")])

## Rename variables for consistency with above
#erbr.1 <- rename(erbr.1, RosNew=Rosettes, InflNew=Infl)
## -----------------------------------------------------------------------------------



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



## MANUALLY ADJUST OR REMOVE TAG NUMS FOR SOME PLTS THAT NEED UNIQUE ID OR ARE OUT OF PLOTS
dats$TagNew[dats$TagNew=="E.7.691" & dats$Year>=2018] <- "E.7.691.1" #This plt was 3 pieces from 2007-2010 ~20 ros, died, then new 1 ros in 2018
dats$TagNew[dats$TagNew=="E.3.110" & dats$Year>=2020] <- "E.3.110.1" #This plt was 2 pieces from 2004-2006 ~60 ros, died, then new 2 ros in 2020
dats <- dats[dats$TagNew !="W.4.304.1",]
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



