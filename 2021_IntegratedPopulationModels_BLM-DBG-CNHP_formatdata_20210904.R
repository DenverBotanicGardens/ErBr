# BLM data ingestion
rm(list=ls())

## LOAD PACKAGES AND FUNCTIONS --------------------------------------------------------------------
# format data
library(dplyr)
library(bazar)
library(tidyr)
library(readr)
library(prism)
library(raster)
library(sp)

# Bayesian analysis
library(lme4)
library(ggplot2)
library(rjags)
library(runjags)
library(coda)
library(corrplot)
## ---
## ------------------------------------------------------------------------------------------------


#I've attached the raw demographic data for 4 of the 5 BLM study sites. The fifth (Oil Well Flats) doesn't have tagged plants. 
#  The data is structured similarly to the BLM Sclerocactus data 0 = veg, 1 = repro, 
#  trailing blanks indicate before detection or post-mortality. One of the sites (Big Bend) 
#  had tags added to a couple transects retroactively so there is missing data for several transects in 2018. 
####################################################################################################################################
# 2016-2020
CastleGarden <- read.csv("C:/Users/DePrengm/Denver Botanic Gardens/Conservation - General/AllProjectsBySpecies/Eriogonum brandegeei/2021_Eriogonum-brandegeei_SSA/Eriogonum-brandegeei_BLM_demographicdata/2021_BLM_demographicdata_CastleGardenSalida.csv")
DroneyGulch_BLM <- read.csv("C:/Users/DePrengm/Denver Botanic Gardens/Conservation - General/AllProjectsBySpecies/Eriogonum brandegeei/2021_Eriogonum-brandegeei_SSA/Eriogonum-brandegeei_BLM_demographicdata/2021_BLM_demographicdata_DroneyGulch.csv")
GardenParkQuarry <- read.csv("C:/Users/DePrengm/Denver Botanic Gardens/Conservation - General/AllProjectsBySpecies/Eriogonum brandegeei/2021_Eriogonum-brandegeei_SSA/Eriogonum-brandegeei_BLM_demographicdata/2021_BLM_demographicdata_GardenParkQuarry.csv")
# 2018-2020
BigBend <- read.csv("C:/Users/DePrengm/Denver Botanic Gardens/Conservation - General/AllProjectsBySpecies/Eriogonum brandegeei/2021_Eriogonum-brandegeei_SSA/Eriogonum-brandegeei_BLM_demographicdata/2021_BLM_demographicdata_BigBend.csv")

# which columns are NA 
allNAcol <- apply(CastleGarden, 1, function(r) which(is.na(r)))
CastleGarden[388,] # in this example, no data in 2016, new in 2017, missing in 2018, dead in 2020
allRECORDcol <- apply(CastleGarden, 1, function(r) which(!is.na(r)))
# remove all dead in col[4] (first year) through to first real data (col > 4) and remove all NA in 2020 through previous missing
# keep all NA between a non-NA within cols 4 - length



head(BigBend)
head(CastleGarden)
# Pivot data from wide to long

seasonAllsites <- read.csv(file = "C:/Users/DePrengm/Denver Botanic Gardens/Conservation - General/AllProjectsBySpecies/Eriogonum brandegeei/2020_Eriogonum-brandegeei_AprilGoebl_PVA/ErBr_scripts_May2021/2021_Eriogonum-brandegeei_climate.csv")
seasonAllsites <- seasonAllsites[,-1]

## Modified from April Goebl & Dan Doak "erbr_1ReformatData_forJAGS_210503.R"
## Script modified 2021-05-27
## Re-format Eriogonum brandegeei data from BLM and CNHP for use with JAGS and IPMs and MPMs
## Add survival for current year
## Add lag columns for modeling with jags for missing data
## Add in column with probability of reproducing (yes or no)


## -------------------------- BLM demography ----------------------

CastleGarden <- pivot_longer(CastleGarden, cols = starts_with("X"), names_to = "Year", values_to = "Stage")
CastleGarden$Year <- extract_numeric(CastleGarden$Year)
CastleGarden$Site <- paste("CG",CastleGarden$Transect,sep="_")
CastleGarden[!is.na(CastleGarden$Stage),]
CastleGarden[is.na(CastleGarden$Stage) ==  TRUE & CastleGarden$Year == 2016,]

# add identifiers to the No Tags
DroneyGulch_BLM$Tag..[grep("No Tag", DroneyGulch_BLM$Tag..)] <- paste(DroneyGulch_BLM$Tag..[grep("No Tag", DroneyGulch_BLM$Tag..)],
      1:length(DroneyGulch_BLM$Tag..[grep("No Tag", DroneyGulch_BLM$Tag..)]), sep = ".")

DroneyGulch_BLM <- pivot_longer(DroneyGulch_BLM, cols = starts_with("X"), names_to = "Year", values_to = "Stage")
DroneyGulch_BLM$Year <- extract_numeric(DroneyGulch_BLM$Year)
DroneyGulch_BLM$Site <- paste("DGBLM",DroneyGulch_BLM$Transect, sep="_")

GardenParkQuarry <- pivot_longer(GardenParkQuarry, cols = starts_with("X"), names_to = "Year", values_to = "Stage")
GardenParkQuarry$Year <- extract_numeric(GardenParkQuarry$Year)
GardenParkQuarry$Site <- paste("GPQ",GardenParkQuarry$Transect, sep="_")

BigBend <- pivot_longer(BigBend, cols = starts_with("X"), names_to = "Year", values_to = "Stage")
BigBend$Year <- extract_numeric(BigBend$Year)
BigBend$Site <- paste("BB",BigBend$Transect,sep="_")
BigBend[!is.na(BigBend$Stage),]
###########################################################################################################


BLM_erbr <- do.call(rbind, list(CastleGarden, DroneyGulch_BLM, GardenParkQuarry, BigBend))
BLM_erbr$Data_id <- 1:nrow(BLM_erbr)
head(as.data.frame(BLM_erbr))

# Make a unique tag per plant; Site is the transect and Site name
BLM_erbr$Tag <- paste(BLM_erbr$Site, BLM_erbr$Tag.., sep="_")
BLM_erbr[is.na(BLM_erbr$Stage),] # leading NAs to a new plant, tail NAs mean dead, middle are missing

# write.csv(BLM_erbr, file = "C:/Users/DePrengm/Denver Botanic Gardens/Conservation - General/AllProjectsBySpecies/Eriogonum brandegeei/2020_Eriogonum-brandegeei_AprilGoebl_PVA/ErBr_scripts_May2021/erbr_BLM_demography_20210603.csv",
#           row.names = FALSE)

## ------------------------------------------------------------------------------------------------

##########################################################################################################
##        Prepare for Impute missing data
## -----------------------------------------


## START OF CHANGES NEEDED FOR JAGs --------------------------------------------------
dats_BLM <- BLM_erbr
# dats_BLM$Data_id <- 1:nrow(dats_BLM) ## did this first thing

## Make a col that indicates if a line should be kept: removes lines w NAs before a plt appeared or after died. Keeps lines that are NAs but bracketed by sz data yrs
dats_BLM$save <- 0   #Start with save = 0, or no
dats_BLM$surv <- NA  #Column to show plt survival/ if plant is alive in current year
## Now add in rows for missing yrs in the set of data for each plt, & make a variable that indicates if the row's stg is a dependent variable & how far back to get to the last stg
dats_BLM$lagsrtstg <- 0    #This is a variable that will indicate if the row lagsrtstg is a dependent value (0=no) or what the lag in time is back to the last observed stg (1,2,3 etc)
dats_BLM$lagforsurv <- 0  #Another lag variable that givens values out to the final yr of non-survival, for plts that died (i.e. if died, when was most recent stg measure?)

# Unlike DBG data, there are rows with NA in the middle already, no need to add those in

tags <- unique(dats_BLM$Tag) 
for (tt in tags){
  stglines <- which(dats_BLM$Tag==tt) 
  stgs <- dats_BLM$Stage[stglines]   
  goodstgslines <- stglines[is.na(stgs)==FALSE] 
  
  badstglines <- stglines[is.na(stgs)==TRUE] #For determining row representing 1st year dead with DBG data, in BLM is missing years
  badyrs <- dats_BLM$Year[badstglines] # years with missing data     
  goodyrs <- dats_BLM$Year[goodstgslines] # years with data
  
  if (length(goodstgslines)>0){
    mingoodlines <- min(goodstgslines)
    maxgoodlines <- max(goodstgslines)
    dats_BLM$save[mingoodlines:maxgoodlines] <- 1 # keep all rows, even NA between first and last good year
    dats_BLM$surv[mingoodlines:maxgoodlines] <- 1 # 
    }
  
  ## If statement that keeps row of data representing 1st year dead, where max badyrs > goodyrs
  if (length(badyrs)>0 && length(goodyrs)>0 && max(badyrs,na.rm=TRUE) > max(goodyrs,na.rm=TRUE)) {
    dats_BLM$save[maxgoodlines+1] <- 1
    dats_BLM$surv[maxgoodlines+1] <- 0  #Change survival to zero for 1st year dead
  }
  
  # Check it:
  # dats_BLM[stglines,]
  ## set lag years for NAs of survival and stage within good years
  if (length(badyrs)>0 && length(goodyrs)>0 && any(badyrs < max(goodyrs, na.rm = TRUE) &
                                                   badyrs > min(goodyrs, na.rm = TRUE))) { # max(badyrs,na.rm=TRUE) < max(goodyrs,na.rm=TRUE)) {
    # need to start at second good year, not first year
    for(yy in ((which(stglines == min(goodstgslines)))+1):length(stglines)){
      # pastyrs <- dats_BLM$Year[stglines[1]:stglines[(yy-1)]]
      firstgood <- stglines[which(stglines == min(goodstgslines))]
      pastyrs <- dats_BLM$Year[firstgood:stglines[(yy-1)]]
      goodpastyrs <- pastyrs[is.na(dats_BLM$Stage[firstgood:stglines[(yy-1)]])==FALSE]
      if(length(goodpastyrs) == 0){
        dats_BLM$lagsrtstg[stglines[yy]] <- 0
        dats_BLM$lagforsurv[stglines[yy]] <- 0
      }
                            
      if(is.na(dats_BLM$Stage[stglines[yy]])== TRUE){ # FALSE){
        dats_BLM$lagsrtstg[stglines[yy]] <- min(dats_BLM$Year[stglines[yy]] - goodpastyrs)
        dats_BLM$lagforsurv[stglines[yy]] <- min(dats_BLM$Year[stglines[yy]] - goodpastyrs)
        }
        ## for years of death, not missing
      if(!is.na(dats_BLM$Stage[stglines[yy]]) && dats_BLM$surv[stglines[yy]] == 0) {
        dats_BLM$lagforsurv[stglines[yy]] <- min(dats_BLM$Year[stglines[yy]] - goodpastyrs)
      }
      }
    } # End setting lag
  }
  
  
# Separate either need pull years of stage and next one... to get into JAGS. 
## scale to zero mean and var = 1 or not? needing to convert back to original scale, didn't to show/see if it matters
table(dats_BLM$Year[dats_BLM$lagsrtstg > 0]) # where there is a missing year of data for a plant
table(dats_BLM$lagsrtstg, dats_BLM$save)

# dats_BLM_nosave  <- dats_BLM[dats_BLM$save != 0,]
# dats_BLM_nosave[dats_BLM_nosave$Tag ==  "BB_43_2078",]
# dats_BLM[dats_BLM$Tag == "BB_43_2078",]

tagsnewdead <- dats_BLM$Tag[dats_BLM$save==0]  
dats_BLM[dats_BLM$Tag == "CG_1_2228",] # Inf for lag should be 0, first year seen? --> now fixed
dats_BLM[dats_BLM$Tag %in% tagsnewdead[-c(1:10)],] 

dats_BLM <- dats_BLM[dats_BLM$save==1,]        #Remove NA rows that are not in middle of the data
table(dats_BLM$surv, dats_BLM$Year) 
dats_BLM[dats_BLM$Year == 2020 & dats_BLM$surv == 0,]


table(dats_BLM$Site)
dats_BLM$SiteName <- "CastleGarden_BLM"
dats_BLM$SiteName[grep("DGBLM", dats_BLM$Site)] <- "Droney_BLM"
dats_BLM$SiteName[grep("GPQ", dats_BLM$Site)] <- "GardenParkQuarry_BLM"
dats_BLM$SiteName[grep("BB", dats_BLM$Site)] <- "BigBend_BLM"

## Number of tagged indiviudals per site
colSums(table(dats_BLM$Tag[!duplicated(dats_BLM[,c("Tag","SiteName")])],
      dats_BLM$SiteName[!duplicated(dats_BLM[,c("Tag","SiteName")])]))

# write.csv(dats_BLM, file = "C:/Users/DePrengm/Denver Botanic Gardens/Conservation - General/AllProjectsBySpecies/Eriogonum brandegeei/2020_Eriogonum-brandegeei_AprilGoebl_PVA/ErBr_scripts_May2021/erbr_BLM_demolag_20210603.csv")

# Add site names to BLM_erbr to match climate - but need dats_BLM to have climate
table(seasonAllsites$Group.3)
table(BLM_erbr$Site)
BLM_erbr$SiteName <- "CastleGarden_BLM"
BLM_erbr$SiteName[grep("DGBLM", BLM_erbr$Site)] <- "Droney_BLM"
BLM_erbr$SiteName[grep("GPQ", BLM_erbr$Site)] <- "GardenParkQuarry_BLM"
BLM_erbr$SiteName[grep("BB", BLM_erbr$Site)] <- "BigBend_BLM"
table(BLM_erbr$SiteName)

## Reshape climate now, not later
seasonAllsites_reshape <- reshape(seasonAllsites, idvar = c("Group.1","Group.3"), timevar = c("Group.2"), direction = "wide")

BLM_erbr.2 <- merge(dats_BLM, seasonAllsites_reshape, by.x = c("Year", "SiteName"), by.y = c("Group.1","Group.3"))
head(BLM_erbr.2)
BLM_erbr.2[is.na(BLM_erbr.2$Stage),]

write.csv(BLM_erbr.2, file = "C:/Users/DePrengm/Denver Botanic Gardens/Conservation - General/AllProjectsBySpecies/Eriogonum brandegeei/2020_Eriogonum-brandegeei_AprilGoebl_PVA/ErBr_scripts_May2021/erbr_BLM_demolag_20210603.csv")

# now merge BLM_erbr for other BLM_erbr.2 things
BLM_erbr.2 <- merge(BLM_erbr, seasonAllsites_reshape, by.x = c("Year", "SiteName"), by.y = c("Group.1","Group.3"))

## -----------------------------------------------------------------------------------


## MAKE TRANSECT UNIQUE (ie. transect, site combo) FOR USE AS PREDICTOR VARIABLE -----
BLM_erbr.2$TransectNew <- BLM_erbr.2$Site
## -----------------------------------------------------------------------------------


## ADD PROBABILITY OF REPRO RESPONSE VARIABLE ----------------------------------------
## Determine whether or not reproduction occurred 

BLM_erbr.2$InflYesNo <- BLM_erbr.2$Stage
## -----------------------------------------------------------------------------------
###############################################################################################


###############################################################################################
# Gross et al. 2002 
## START OF CHANGES NEEDED FOR JAGs --------------------------------------------------
dats_BLM <- BLM_erbr.2
## Make transect & year values numerical to use in jags as random effects 
dats_BLM$TransectNew.num <- as.factor(dats_BLM$Site)
dats_BLM$TransectNew.num <- as.numeric(dats_BLM$TransectNew.num) # unique number per transect across sites
TransectNew.num <- dats_BLM$TransectNew.num
table(TransectNew.num) # index for each transect; BigBend TransectNum.new == 1:53
table(dats_BLM$TransectNew.num, dats_BLM$Year, dats_BLM$SiteName)

dats_BLM$Year.num <- as.factor(dats_BLM$Year)
dats_BLM$Year.num <- as.numeric(dats_BLM$Year.num)   
Year.num <- dats_BLM$Year.num
table(Year.num) # index for each year survey; 1:5 == 2016:2020

## Make a linear index of transect-year combos
yrtranscombo <- 100*dats_BLM$TransectNew.num+dats_BLM$Year.num # Big bend first year == 3 == 2018 for <= 5300
table(yrtranscombo)

## Make dataframe w new plts (that are likely recent seedlings) for each transect & yr ------------
## Make df that will hold data containing new plants
years <- unique(dats_BLM$Year.num)
years <- years[order(years)]
dats_BLM.newPlts <- as.data.frame(rep(unique(dats_BLM$TransectNew.num), each=length(years)))
colnames(dats_BLM.newPlts) <- "TransectNew.num"
dats_BLM.newPlts$Year.num <- rep(years) # Now just a data frame with unique Transect number for each year

## Identify new plants 
newPlts <- dats_BLM[!is.na(dats_BLM$Stage),] %>% group_by(Tag) %>% slice(which.min(Year))   #Identify rows with 1st appearance for each plt
newPlts <- newPlts[newPlts$Year!=2016,]                           #Remove 2016 (first year of data collection) - 2018 first year for Big Bend
newPlts[which(newPlts$SiteName == "BigBend_BLM" & newPlts$Year == 2018),]
newPlts <- newPlts[-(which(newPlts$SiteName == "BigBend_BLM" & newPlts$Year == 2018)),]
# sz.cutoff <- 5                                                    #Sz cutoff, above which plt was likely not recently a seedling 
newPlts <- newPlts[newPlts$Stage == 0,]                   #Remove if reproductive, likely not reproductive if a first year plant
num.newPlts <- newPlts %>% group_by(TransectNew.num, Year.num) %>% summarise(num.newPlts=n())  #Count num new plts per yr & transect

### Edit 2021-09-04. To be able to bootstrap transitions without adding a fertility estimate to BLM's data, change newPlts to 'seedling'
dats_BLM_seedling <- dats_BLM
newPlts2 <- dats_BLM_seedling[!is.na(dats_BLM_seedling$Stage),] %>% group_by(Tag) %>% slice(which.min(Year))   #Identify rows with 1st appearance for each plt
newPlts2 <- newPlts2[newPlts2$Year!=2016,]                           #Remove 2016 (first year of data collection) - 2018 first year for Big Bend
newPlts2[which(newPlts2$SiteName == "BigBend_BLM" & newPlts2$Year == 2018),]
newPlts2 <- newPlts2[-(which(newPlts2$SiteName == "BigBend_BLM" & newPlts2$Year == 2018)),]
newPlts2 <- newPlts2[newPlts2$Stage == 0,]

write.csv(newPlts2, file = "C:/Users/DePrengm/Denver Botanic Gardens/Conservation - General/AllProjectsBySpecies/Eriogonum brandegeei/2020_Eriogonum-brandegeei_AprilGoebl_PVA/ErBr_scripts_May2021/erbr_newPlts2.csv")

## Add number of new plants to df of each transect & year
dats_BLM.newPlts <- left_join(dats_BLM.newPlts, num.newPlts, by=c("TransectNew.num", "Year.num"))
dats_BLM.newPlts$num.newPlts[is.na(dats_BLM.newPlts$num.newPlts)] <- 0   #Change NAs (no new plants) to zeros
dats_BLM.newPlts$num.newPlts[dats_BLM.newPlts$Year.num==1] <- NA         #Change new plts in 2016 (yr 1) to NA
# dats_BLM.newPlts[dats_BLM.newPlts$Year.num==3 & dats_BLM.newPlts$TransectNew.num < 54,] # change new in 2018 from BB to NA
dats_BLM.newPlts[dats_BLM.newPlts$Year.num==3 & dats_BLM.newPlts$TransectNew.num < 11,] # change since mistake, less than 11; 1:10 are BB
dats_BLM.newPlts$num.newPlts[dats_BLM.newPlts$Year.num < 4 & dats_BLM.newPlts$TransectNew.num < 11] <- NA # change new plts in 2016-2018 at Big Bend to NA

## Add column so new plts in t+1 match year t
dats_BLM.newPlts <- dats_BLM.newPlts %>% mutate(num.newPlts1=lead(num.newPlts)) 

BLM_newplts <- merge(dats_BLM.newPlts, dats_BLM[!duplicated(dats_BLM[,c("Year","SiteName","Site","TransectNew.num","Year.num")]),
                                                       c("Year","SiteName","Site","TransectNew.num","Year.num")],
                     by = c("Year.num","TransectNew.num"))

head(BLM_newplts)
write.csv(BLM_newplts, file = "C:/Users/DePrengm/Denver Botanic Gardens/Conservation - General/AllProjectsBySpecies/Eriogonum brandegeei/2020_Eriogonum-brandegeei_AprilGoebl_PVA/ErBr_scripts_May2021/erbr_BLM_newplts.csv",
     row.names = FALSE)


## Add climate variables to new plants data 
### add the data_id (row number column)
# dats_BLM <- dats_BLM_seedling
# dats_BLM_wide <- reshape(dats_BLM, # [,!(names(dats_BLM) %in% c("save","surv","lagsrtstg","lagforsurv"))], # take away these? Not needed when not doing JAGS model
#                          # idvar = names(dats_BLM)[c(1:8,12:15)], # This was before I added Data_id and changed other things!
#                          idvar = c("Year","SiteName","Transect","Meter","Tag..","Stage","Site","Data_id","Tag",
#                                    "save","surv","lagsrtstg","lagforsurv",
#                                    "TransectNew","InflYesNo","TransectNew.num","Year.num"), 
#                          timevar = "Group.2", 
#                          direction = "wide")  
# write.csv(dats_BLM, file = "C:/Users/DePrengm/Denver Botanic Gardens/Conservation - General/AllProjectsBySpecies/Eriogonum brandegeei/2020_Eriogonum-brandegeei_AprilGoebl_PVA/ErBr_scripts_May2021/erbr_BLM_demolag_20210603.csv")
# write.csv(dats_BLM_wide, file = "C:/Users/DePrengm/Denver Botanic Gardens/Conservation - General/AllProjectsBySpecies/Eriogonum brandegeei/2020_Eriogonum-brandegeei_AprilGoebl_PVA/ErBr_scripts_May2021/erbr_BLM_demolag_20210603.csv")

# Big bend and droney are the same
# ggplot(dats_BLM_wide, aes(Year, meanTmean.summer, colour = SiteName))+
#   geom_jitter(height = 0.05, width = 0.05)+
#   geom_line()
# 
# ggplot(dats_BLM_wide, aes(meanTmean.summer, Stage, colour = SiteName))+
#   geom_point()+
#   geom_smooth(method = "glm", method.args = list(family = "binomial"))
# 
# ggplot(dats_BLM_wide, aes(TotPpt.winter, Stage, colour = SiteName))+
#   geom_point()+
#   geom_smooth(method = "glm", method.args = list(family = "binomial"))
# 
# 
# ggplot(dats_BLM_wide, aes(Year, TotPpt.fall, colour = SiteName))+
#   geom_jitter()+
#   geom_line()

# dats_BLM.clim <- dats_BLM %>% dplyr::select(c(Year.num, PptFall, PptWinter, PptSummer, TempFall, TempWinter, TempSummer)) 
# clim <- unique(dats_BLM.clim)
# dats_BLM.newPlts <- dplyr::left_join(dats_BLM.newPlts, clim, by="Year.num")
# dats_BLM.clim <- dats_BLM_wide[,names(dats_BLM_wide)[c(11:18)]] ## ** HERE? do these numbers not line up with Data_id added? YES!
# dats_BLM.clim <- dats_BLM_wide[,c("TransectNew.num","Year.num","meanTmean.summer","TotPpt.summer","meanTmean.fall",
dats_BLM.clim <- dats_BLM[,c("TransectNew.num","Year.num","meanTmean.summer","TotPpt.summer","meanTmean.fall",
                                  "TotPpt.fall","meanTmean.winter","TotPpt.winter")]  
dats_BLM.clim <- dats_BLM.clim[!duplicated(dats_BLM.clim$TransectNew.num),]
dats_BLM.newPlts <- merge(dats_BLM.newPlts, 
                          dats_BLM.clim, by = c("Year.num","TransectNew.num")) # makes .x and .y when didn't before; now doesn't...
dats_BLM.newPlts <- dats_BLM.newPlts %>% mutate(TotPpt.fall1=lead(TotPpt.fall), 
                                                TotPpt.winter1=lead(TotPpt.winter), 
                                                TotPpt.summer1=lead(TotPpt.summer),
                                                meanTmean.fall1=lead(meanTmean.fall), 
                                                meanTmean.winter1=lead(meanTmean.winter),
                                                meanTmean.summer1=lead(meanTmean.summer))  
## Note: variables with '1' at the end (e.g. num.newPlts1, PptFall1) represent t+1 data   ## Why would you want next year climate data?


newplts <- dats_BLM.newPlts$num.newPlts1
newplt.trans <- dats_BLM.newPlts$TransectNew.num
newplt.yr <- dats_BLM.newPlts$Year.num

newPltlines <- length(dats_BLM.newPlts$TransectNew.num)

## Make a linear index of transect-year combos for new plts
newplt.yrtranscombo=100*newplt.trans+newplt.yr 
## ------------------------------------------------------------------------------------------------
write.csv(dats_BLM.clim, file = "C:/Users/DePrengm/Denver Botanic Gardens/Conservation - General/AllProjectsBySpecies/Eriogonum brandegeei/2020_Eriogonum-brandegeei_AprilGoebl_PVA/ErBr_scripts_May2021/erbr_BLM_climate_20210601.csv",
          row.names = FALSE)
write.csv(dats_BLM.newPlts, file = "C:/Users/DePrengm/Denver Botanic Gardens/Conservation - General/AllProjectsBySpecies/Eriogonum brandegeei/2020_Eriogonum-brandegeei_AprilGoebl_PVA/ErBr_scripts_May2021/erbr_BLMnewplts_climate_20210601.csv",
          row.names = FALSE)



## RUN ASSOCIATED JAGS MODEL ----------------------------------------------------------------------
# jags.mod <- run.jags('erbr_BLMJAGSmod_20210601.R', n.chains=3, data=dats, burnin=5000, thin=10, sample=30000, adapt=500, method='parallel')

#save(jags.mod, file='erbr_JAGSmod_c3t10s20b5_210406.rdata')
# saveRDS(jags.mod, "C:/Users/DePrengm/Denver Botanic Gardens/Conservation - General/AllProjectsBySpecies/Eriogonum brandegeei/2021_Eriogonum-brandegeei_SSA/Eriogonum-brandegeei_JAGS_output/erbr_JAGSmod_c3t10s30b5_210509.rds")
## ------------------------------------------------------------------------------------------------


## LOOK AT MODEL OUTPUT ---------------------------------------------------------------------------
# summary(jags.mod)
# plot(jags.mod)
# summ.mod <- summary(jags.mod)
# tail(summ.mod, n=20)
# gelman.diag(jags.mod, confidence = 0.95, transform=FALSE)
# 
# chains <- jags.mod$mcmc
# chains <- bind_rows(lapply(chains, as.data.frame))
# colMeds <- apply(chains,2,median)
# colSDs <- apply(chains,2,sd)
# 
# 
# ## Look at correlation b/w params 
# chains.1 <- chains %>% dplyr::select(!contains(c("randomeffect", "precision")))
# chains.1 <- chains.1 %>% dplyr::select(!c(deviance, resid.sum.sq))
# 
# cor.chains <- cor(chains.1)
# corrplot(cor.chains, method="circle", type="lower")


## ** Make bar graph comparing median param ests & 80% CIs b/w diff datasets **
## ------------------------------------------------------------------------------------------------




## -----------------------------------------------------------------------------------





## ------------------------------------------------------------------------------------------------
## CNHP does not have individuals tagged. Growth rate and reproductive rates are what we can pull. 
#   Reproduction by size - not connected but maybe in aggregate by transect? 
DroneyGulch_CNHP <- read.csv("C:/Users/DePrengm/Denver Botanic Gardens/Conservation - General/AllProjectsBySpecies/Eriogonum brandegeei/2021_Eriogonum-brandegeei_SSA/Eriogonum-brandegeei_CNHP_CNAP_demographicdata/Droney_GulchSWA_reformat_CNHP_2016-2021.csv")

ggplot(DroneyGulch_CNHP, aes(Year, Total.Clumps, colour = as.factor(Meter.Mark)))+
  geom_point()+
  stat_smooth(method = "glm", se = FALSE)+
  theme_bw()+
  scale_color_discrete(guide = FALSE)+
  ggtitle("Droney Gulch - Colorado Natural Heritage Program")

total2plot <- aggregate(Total.Clumps ~ Year, sum,data = DroneyGulch_CNHP)
ggplot(total2plot[total2plot$Year>2017,], aes(Year, Total.Clumps))+
  geom_point()+
  stat_smooth(method = "glm", se = TRUE)+
  theme_bw()+
  ggtitle("Droney Gulch - Colorado Natural Heritage Program")

## Count-based extinction probabilities and bootstrap confidence intervals
# popbio::countCDFex

tDGCNHP <- table(DroneyGulch_CNHP$Year,DroneyGulch_CNHP$Plot..) # a subset are 2016-2021, the remaining 2018-2021
plots6yrs <- table(DroneyGulch_CNHP$Year,DroneyGulch_CNHP$Plot..)[,(colSums(table(DroneyGulch_CNHP$Year,DroneyGulch_CNHP$Plot..))==6)]
count6yrs <- as.numeric(as.character(dimnames(plots6yrs)[[2]]))

total6yrs <- aggregate(Total.Clumps ~ Year, sum, data =  DroneyGulch_CNHP[DroneyGulch_CNHP$Plot.. %in% count6yrs,])
logN <- log(total6yrs$Total.Clumps[-c(1,6)]/total6yrs$Total.Clumps[-c(5:6)])
countCDFxt(mu = mean(logN), sig2 = var(logN), nt = 5, 
           Nc = 49, Ne = 1, tmax = 100, plot = TRUE)


## all plots 2018-2021
logN <- log(total2plot$Total.Clumps[total2plot$Year>2017][-1]/total2plot$Total.Clumps[total2plot$Year>2017][-4])
countCDFxt(mu = mean(logN), sig2 = var(logN), nt = 3, 
           Nc = 49, Ne = 1, tmax = 100, plot = TRUE)



# Each plot on its own
## Test
# x <- split(DroneyGulch_CNHP, DroneyGulch_CNHP$Plot..)[[1]]
# rm(x)

## Some counts of a transect are zero 
growthrates <- unlist(lapply(split(DroneyGulch_CNHP, DroneyGulch_CNHP$Plot..), function(x){
  logN <- log(x$Total.Clumps[-1]/x$Total.Clumps[-length(x$Total.Clumps)]) ## for the log( N_t+1 / N_t )
  }))



median(exp(growthrates[!is.infinite(growthrates) & !is.nan(growthrates)]))
growthratesgood <- growthrates[!is.infinite(growthrates) & !is.nan(growthrates)]
sd(growthratesgood)
mu_allgood <- mean(growthratesgood)
exp(mu_allgood)
sig2_allgood <- var(growthratesgood)
exp(c(mu_allgood - qt(0.975, nt - 1) * sqrt(sig2_allgood)/length(growthratesgood),
  mu + qt(0.975, length(growthratesgood) - 1) * sqrt(sig2_allgood)/length(growthratesgood)))

countCDFxt(mu =mu_allgood,sig2 = sig2_allgood, nt = length(growthratesgood), 
           Nc = floor(mean(DroneyGulch_CNHP$Total.Clumps)), Ne = 1, tmax = 30, plot = TRUE)
## Risk of extinction greater than 0.05 by about 4 years into the future 
