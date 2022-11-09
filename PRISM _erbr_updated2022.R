## April Goebl
## Script modified 10-21-20
## Collaboration with Denver Botanic Gardens
## Extract PRISM climate data for Eriogonum brandegeei





rm(list=ls())
dev.off()



## LOAD PACKAGES AND FUNCTIONS ----------------------------------------------
library(raster)
library(rgdal)
library(dplyr)
library(stringr)
library(tidyr)
library(prism)
## --------------------------------------------------------------------------




## EXTRACT PRISM CLIMATE DATA -----------------------------------------------

erbrlatlong <- data.frame(Lat = 38.5434, Long = -105.2184)
##From Michelle's 2018 report: Garden Park East macroplot origin= -105.195433, 38.543469 and Garden Park West macroplot origin= -105.199197, 38.544052


## MAKE LIST OF DESIRED YEAR & MONTH TO EXTRACT CLIMATE DATA FROM -----------
# year.span <- as.character(c(2003:2018))
year.span <- as.character(c(2003:2022))
month.span <- str_pad(as.character(c(1:12)), width=2, side="left", pad="0")

yyyymm <- NULL

for (yy in 1:length(year.span)) {
  for (mm in 1:length(month.span)) {
    yyyymm <- c(yyyymm,paste(year.span[yy], month.span[mm], sep=""))
  }
}

## Because we measure in early August each year so data for current year are the 12 preceeding months
## Remove fist 7 months (of 2013) and final 5 months (of 2018)  ## Do you mean 2003 and last of 2018?
yyyymm <- yyyymm[8:(length(yyyymm)-5)]
## Morph year & month labels into matrix
yyyymm <- as.data.frame(matrix(yyyymm, length(month.span), (length(year.span)-1)), row.names = FALSE)
colnames(yyyymm) <- as.character(c(2004:2022))
## ------------------------------------------------------------------------------



## CALCULATE CLIMATE VARIABLE FOR 3 SEASONS
#fall = Previous Aug-Nov
#Winter = Previous Dec - March
#Summer =  April - July



prism_set_dl_dir("Q:/Research/All_Projects_by_Species/Astragalus SPECIES/Astragalus_microcymbus/PRISM_asmiclimate")
get_prism_monthlys(type="ppt", mon = 1:12, keepZip = FALSE, years = 2022:2023)
get_prism_monthlys(type="tmean", mon = 1:12, keepZip = FALSE, years = 2003:2022)


m.ppt <- prism_archive_subset("ppt", "monthly", years = 2003:2022)
m.ppt <- do.call(rbind,lapply(m.ppt, function(x){
  x_rast <- pd_to_file(x)
  rastertemps <- raster(x_rast)
  data.frame(data = raster::extract(rastertemps,
                                    erbrlatlong[,c("Long","Lat")]),
             date = x)
}))

erbr.climate.monthly <- m.ppt %>%
  mutate(Year = as.numeric(substr(date,  (nchar(date)+1)-10,nchar(date)-6)),
         Month = as.numeric(substr(date, (nchar(date)+1)-6, nchar(date)-4)),
         # Make previous 12 months match year (previous to August survey)
         Prev12 = ifelse(Month > 7, Year+1, Year)) %>%
  #extract the climate variable name from the prism long name
  separate(date, sep = "_", into = c(NA,"Variable",NA,NA,NA,NA)) %>%
  mutate(season = case_when(Month %in% c(1:3,12) ~ "winter",
                            Month>=8 & Month<=11 ~ "fall",
                            Month>=4 & Month<=7 ~ "summer")) %>%
  group_by(season, Prev12, Variable) %>%
  dplyr::summarise(Value = sum(data), .groups = "keep")  %>%
  pivot_wider(names_from = c(Variable,season), values_from = Value) %>%
  dplyr::rename(Year = Prev12)


m.tmean <- prism_archive_subset("tmean", "monthly", years = 2003:2022)
m.tmean <- do.call(rbind,lapply(m.tmean, function(x){
  x_rast <- pd_to_file(x)
  rastertemps <- raster(x_rast)
  data.frame(data = raster::extract(rastertemps,
                                    erbrlatlong[,c("Long","Lat")]),
             date = x)
}))

erbr.climate.monthly.tmean <- m.tmean %>%
  mutate(Year = as.numeric(substr(date,  (nchar(date)+1)-10,nchar(date)-6)),
         Month = as.numeric(substr(date, (nchar(date)+1)-6, nchar(date)-4)),
         # Make previous 12 months match year (previous to August survey)
         Prev12 = ifelse(Month > 7, Year+1, Year)) %>%
  #extract the climate variable name from the prism long name
  separate(date, sep = "_", into = c(NA,"Variable",NA,NA,NA,NA)) %>%
  # mutate(season = "winter") %>%
  mutate(season = case_when(Month %in% c(1:3,12) ~ "winter",
                            Month>=8 & Month<=11 ~ "fall",
                            Month>=4 & Month<=7 ~ "summer")) %>%
  group_by(season, Prev12, Variable) %>%
  dplyr::summarise(Value = mean(data), .groups = "keep") %>%
  pivot_wider(names_from = c(Variable,season), values_from = Value) %>%
  dplyr::rename(Year = Prev12)

## TOTAL SEASONAL PRECIP ----------------------------------------------------------
# path.working <- "C:/Users/april/Dropbox/CU_Boulder_PhD/DBG_Internship/Prism_data_ppt/"
path.working <- "Q:/Research/All_Projects_by_Species/Astragalus SPECIES/Astragalus_microcymbus/PRISM_asmiclimate"

## Make a list of directories containing data for each year (must be unzipped)
dir.list <- list.files(path = path.working, full.names = FALSE, recursive = FALSE, include.dirs=TRUE)

## Loop over directory list and combine list of all bils from all years
bils.ppt <- NULL

for (bb in 1:length(dir.list)) {
  bils.ppt <- c(bils.ppt, list.files(path = paste(path.working,dir.list[bb], sep=""),
                                     pattern=".bil$", full.names = TRUE))
}

## Extract precip data and sum over months of each 'season'
ppt <- as.data.frame(matrix(NA, ncol(yyyymm), 4))
colnames(ppt) <- c("Year", "Tot_fall_ppt", "Tot_winter_ppt", "Tot_summer_ppt")
ppt$Year <- colnames(yyyymm)

for (rr in 1:nrow(ppt)) {
  ppt.fall <- 0
  ppt.winter <- 0
  ppt.summer <- 0

  for (pp in 1:4) {
    file.pos <- grep(yyyymm[pp,rr], bils.ppt)
    raster_file <- raster(bils.ppt[file.pos])                                           #Load specified raster
    ppt.fall <- ppt.fall + raster::extract(x = raster_file, erbrlatlong[,c("Long","Lat")])      #Extract ppt for given gps coord & month, sum to prev month
  }

  for (ww in 5:8) {
    file.pos <- grep(yyyymm[ww,rr], bils.ppt)
    raster_file <- raster(bils.ppt[file.pos])                                               #Load specified raster
    ppt.winter <- ppt.winter + raster::extract(x = raster_file, erbrlatlong[,c("Long","Lat")])      #Extract ppt for given gps coord & month, sum to prev month
  }

  for (ss in 9:12) {
    file.pos <- grep(yyyymm[ss,rr], bils.ppt)
    raster_file <- raster(bils.ppt[file.pos])                                               #Load specified raster
    ppt.summer <- ppt.summer + raster::extract(x = raster_file, erbrlatlong[,c("Long","Lat")])      #Extract ppt for given gps coord & month, sum to prev month
  }

  ppt[rr,2] <-  ppt.fall
  ppt[rr,3] <-  ppt.winter
  ppt[rr,4] <-  ppt.summer

}
## ------------------------------------------------------------------------------




## MEAN SEASONAL DAILY MEAN TEMP ---------------------------------------------------
path.temp <- "C:/Users/april/Dropbox/CU_Boulder_PhD/DBG_Internship/Prism_data_temp/"

## Make a list of directories containing data for each year (must be unzipped)
dir.temp <- list.files(path=path.temp, full.names=FALSE, recursive=FALSE, include.dirs=TRUE)

## Loop over directory list and combine list of all bils from all years
bils.temp <- NULL

for (bb in 1:length(dir.temp)) {
  bils.temp <- c(bils.temp, list.files(path = paste(path.temp,dir.temp[bb], sep=""),
                                       pattern=".bil$", full.names = TRUE))
}

## Extract mean temp data & avg over 12 months preceding annual demog data collection
temp <- as.data.frame(matrix(NA, ncol(yyyymm), 4))
colnames(temp) <- c("Year", "Mean_fall_temp", "Mean_winter_temp", "Mean_summer_temp")
temp$Year <- colnames(yyyymm)

for (rr in 1:nrow(temp)) {
  temp.fall <- 0
  temp.winter <- 0
  temp.summer <- 0

  for (pp in 1:4) {
    file.pos <- grep(yyyymm[pp,rr], bils.temp)
    raster_file <- raster(bils.temp[file.pos])                                          #Load specified raster
    temp.fall <- temp.fall + raster::extract(x = raster_file, erbrlatlong[,c("Long","Lat")])    #Extract temp for given gps coord & month, sum to prev month
  }

  for (ww in 5:8) {
    file.pos <- grep(yyyymm[ww,rr], bils.temp)
    raster_file <- raster(bils.temp[file.pos])                                              #Load specified raster
    temp.winter <- temp.winter + raster::extract(x = raster_file, erbrlatlong[,c("Long","Lat")])    #Extract temp for given gps coord & month, sum to prev month
  }

  for (pp in 9:12) {
    file.pos <- grep(yyyymm[pp,rr], bils.temp)
    raster_file <- raster(bils.temp[file.pos])                                              #Load specified raster
    temp.summer <- temp.summer + raster::extract(x = raster_file, erbrlatlong[,c("Long","Lat")])    #Extract temp for given gps coord & month, sum to prev month
  }

  temp[rr,2] <-  temp.fall/4                                            #Calc mean seasonal temp
  temp[rr,3] <-  temp.winter/4                                          #Calc mean seasonal temp
  temp[rr,4] <-  temp.summer/4                                          #Calc mean seasonal temp
}
## ------------------------------------------------------------------------------




## MEAN MONTHLY SEASONAL TMAX ---------------------------------------------------
path.tmax <- "C:/Users/april/Dropbox/CU_Boulder_PhD/DBG_Internship/Prism_data_tmax/"

## Make a list of directories containing data for each year (must be unzipped)
dir.tmax <- list.files(path=path.tmax, full.names=FALSE, recursive=FALSE, include.dirs=TRUE)

## Loop over directory list and combine list of all bils from all years
bils.tmax <- NULL

for (bb in 1:length(dir.tmax)) {
  bils.tmax <- c(bils.tmax, list.files(path=paste(path.tmax, dir.tmax[bb], sep=""),
                                       pattern=".bil$", full.names = TRUE))
}

## Extract max temp data & avg over 12 months preceding annual demog data collection
tmax <- as.data.frame(matrix(NA, ncol(yyyymm), 4))
colnames(tmax) <- c("Year", "Mean_fall_tmax", "Mean_winter_tmax", "Mean_summer_tmax")
tmax$Year <- colnames(yyyymm)

for (rr in 1:nrow(tmax)) {
  tmax.fall <- 0
  tmax.winter <- 0
  tmax.summer <- 0

  for (pp in 1:4) {
    file.pos <- grep(yyyymm[pp,rr], bils.tmax)
    raster_file <- raster(bils.tmax[file.pos])                                          #Load specified raster
    tmax.fall <- tmax.fall + raster::extract(x = raster_file, erbrlatlong[,c("Long","Lat")])    #Extract tmax for given gps coord & month, sum to prev month
  }

  for (ww in 5:8) {
    file.pos <- grep(yyyymm[ww,rr], bils.tmax)
    raster_file <- raster(bils.tmax[file.pos])                                              #Load specified raster
    tmax.winter <- tmax.winter + raster::extract(x = raster_file, erbrlatlong[,c("Long","Lat")])    #Extract tmax for given gps coord & month, sum to prev month
  }

  for (ss in 9:12) {
    file.pos <- grep(yyyymm[ss,rr], bils.tmax)
    raster_file <- raster(bils.tmax[file.pos])                                              #Load specified raster
    tmax.summer <- tmax.summer + raster::extract(x = raster_file, erbrlatlong[,c("Long","Lat")])    #Extract tmax for given gps coord & month, sum to prev month
  }

  tmax[rr,2] <-  tmax.fall/4                                                          #Calc mean monthly tmax
  tmax[rr,3] <-  tmax.winter/4                                                        #Calc mean monthly tmax
  tmax[rr,4] <-  tmax.summer/4                                                        #Calc mean monthly tmax

}
## ------------------------------------------------------------------------------




## MEAN MONTHLY SEASONAL TMIN ---------------------------------------------------
path.tmin <- "C:/Users/april/Dropbox/CU_Boulder_PhD/DBG_Internship/Prism_data_tmin/"

## Make a list of directories containing data for each year (must be unzipped)
dir.tmin <- list.files(path=path.tmin, full.names=FALSE, recursive=FALSE, include.dirs=TRUE)

## Loop over directory list and combine list of all bils from all years
bils.tmin <- NULL

for (bb in 1:length(dir.tmin)) {
  bils.tmin <- c(bils.tmin, list.files(path=paste(path.tmin, dir.tmin[bb], sep=""),
                                       pattern=".bil$", full.names = TRUE))
}

## Extract min temp data & avg over 12 months preceding annual demog data collection
tmin <- as.data.frame(matrix(NA, ncol(yyyymm), 4))
colnames(tmin) <- c("Year", "Mean_fall_tmin", "Mean_winter_tmin", "Mean_summer_tmin")
tmin$Year <- colnames(yyyymm)

for (rr in 1:nrow(tmin)) {
  tmin.fall <- 0
  tmin.winter <- 0
  tmin.summer <- 0

  for (pp in 1:4) {
    file.pos <- grep(yyyymm[pp,rr], bils.tmin)
    raster_file <- raster(bils.tmin[file.pos])                                          #Load specified raster
    tmin.fall <- tmin.fall + raster::extract(x = raster_file, erbrlatlong[,c("Long","Lat")])    #Extract tmin for given gps coord & month, sum to prev month
  }

  for (ww in 5:8) {
    file.pos <- grep(yyyymm[ww,rr], bils.tmin)
    raster_file <- raster(bils.tmin[file.pos])                                              #Load specified raster
    tmin.winter <- tmin.winter + raster::extract(x = raster_file, erbrlatlong[,c("Long","Lat")])    #Extract tmin for given gps coord & month, sum to prev month
  }

  for (ss in 9:12) {
    file.pos <- grep(yyyymm[ss,rr], bils.tmin)
    raster_file <- raster(bils.tmin[file.pos])                                              #Load specified raster
    tmin.summer <- tmin.summer + raster::extract(x = raster_file, erbrlatlong[,c("Long","Lat")])    #Extract tmin for given gps coord & month, sum to prev month
  }

  tmin[rr,2] <-  tmin.fall/4                                                          #Calc mean monthly tmin
  tmin[rr,3] <-  tmin.winter/4                                                        #Calc mean monthly tmin
  tmin[rr,4] <-  tmin.summer/4                                                        #Calc mean monthly tmin
}
## -------------------------------------------------------------------------------



## SAVE CLIMATE DATA --------------------------------------------------------------
## Combine different climate variables into one data frame
clim <- cbind(ppt, temp[,2:4], tmax[,2:4], tmin[,2:4])

## Save output
write.csv(clim, "erbr_climData3seasons_201022.csv", row.names=FALSE)
## -------------------------------------------------------------------------------





## LOOK AT CORRELATIONS OF CLIMATE VARS ------------------------------------------
## ** See erbr_dredge_tagClust.R **
#cor.clim <- cor(erbr[,6:11])
#corrplot(cor.clim, method="pie", type="lower")
#corrplot.mixed(cor.clim, lower.col="black")






## DOWNLOAD AND EXTRACT SEASONAL VALUES FROM PAST (E.G. 30 YEARS) ---------------

## FROM PRISM PACKAGE TUTORIAL
## https://github.com/ropensci/prism

## Set download directroy
prism_set_dl_dir("C:/Users/april/Dropbox/CU_Boulder_PhD/DBG_Internship/Prism_data")

## Download data
#get_prism_monthlys(type=c("ppt"), year=2001:2020, mon=1:12, keepZip=FALSE)
#get_prism_monthlys(type=c("tmean"), year=2003:2020, mon=1:12, keepZip=FALSE)
## --------------------------------------------------------------------------



## ASSIGN TEMPORAL RANGE FOR DESIRED CLIMATE DATA ---------------------------
erbrlatlong <- data.frame(Lat = 38.5434, Long = -105.2184)


## MAKE LIST OF DESIRED YEAR & MONTH TO EXTRACT CLIMATE DATA FROM
year.span <- as.character(c(1990:2020))
month.span <- str_pad(as.character(c(1:12)), width=2, side="left", pad="0")

yyyymm <- NULL

for (yy in 1:length(year.span)) {
  for (mm in 1:length(month.span)) {
    yyyymm <- c(yyyymm,paste(year.span[yy], month.span[mm], sep=""))
  }
}

## Remove fist 7 months and final 5 months to match the 'census year' (i.e. August-July)
## ** Note 2020 only had data upto July at the time of download **
yyyymm <- yyyymm[8:(length(yyyymm))]
## Morph year & month labels into matrix
yyyymm <- as.data.frame(matrix(yyyymm, length(month.span), (length(year.span)-1)), row.names = FALSE)
colnames(yyyymm) <- as.character(c(1991:2020))
## ------------------------------------------------------------------------------



## TOTAL SEASONAL PRECIP ----------------------------------------------------------
dirs.ppt <- prism_archive_subset("ppt", "monthly", year=1990:2020)
bils.ppt <- pd_to_file(dirs.ppt)


## Extract precip data and sum over months of each 'season'
ppt <- as.data.frame(matrix(NA, ncol(yyyymm), 4))
colnames(ppt) <- c("Year", "Tot_fall_ppt", "Tot_winter_ppt", "Tot_summer_ppt")
ppt$Year <- colnames(yyyymm)

for (rr in 1:nrow(ppt)) {
  ppt.fall <- 0
  ppt.winter <- 0
  ppt.summer <- 0

  for (pp in 1:4) {
    file.pos <- grep(yyyymm[pp,rr], bils.ppt)
    raster_file <- raster(bils.ppt[file.pos])                                           #Load specified raster
    ppt.fall <- ppt.fall + raster::extract(x = raster_file, erbrlatlong[,c("Long","Lat")])      #Extract ppt for given gps coord & month, sum to prev month
  }

  for (ww in 5:8) {
    file.pos <- grep(yyyymm[ww,rr], bils.ppt)
    raster_file <- raster(bils.ppt[file.pos])                                               #Load specified raster
    ppt.winter <- ppt.winter + raster::extract(x = raster_file, erbrlatlong[,c("Long","Lat")])      #Extract ppt for given gps coord & month, sum to prev month
  }

  for (ss in 9:12) {
    file.pos <- grep(yyyymm[ss,rr], bils.ppt)
    raster_file <- raster(bils.ppt[file.pos])                                               #Load specified raster
    ppt.summer <- ppt.summer + raster::extract(x = raster_file, erbrlatlong[,c("Long","Lat")])      #Extract ppt for given gps coord & month, sum to prev month
  }

  ppt[rr,2] <-  ppt.fall
  ppt[rr,3] <-  ppt.winter
  ppt[rr,4] <-  ppt.summer
}
## ------------------------------------------------------------------------------



## MEAN SEASONAL DAILY MEAN TEMP ---------------------------------------------------
dirs.temp <- prism_archive_subset("tmean", "monthly", year=1990:2020)
bils.temp <- pd_to_file(dirs.temp)


## Extract mean temp data & avg over 4 months of each season
temp <- as.data.frame(matrix(NA, ncol(yyyymm), 4))
colnames(temp) <- c("Year", "Mean_fall_temp", "Mean_winter_temp", "Mean_summer_temp")
temp$Year <- colnames(yyyymm)

for (rr in 1:nrow(temp)) {
  temp.fall <- 0
  temp.winter <- 0
  temp.summer <- 0

  for (pp in 1:4) {        #Aug-Nov
    file.pos <- grep(yyyymm[pp,rr], bils.temp)
    raster_file <- raster(bils.temp[file.pos])                                          #Load specified raster
    temp.fall <- temp.fall + raster::extract(x = raster_file, erbrlatlong[,c("Long","Lat")])       #Extract temp for given gps coord & month, sum to prev month
  }

  for (ww in 5:8) {        #Dec-Mar
    file.pos <- grep(yyyymm[ww,rr], bils.temp)
    raster_file <- raster(bils.temp[file.pos])                                           #Load specified raster
    temp.winter <- temp.winter + raster::extract(x = raster_file, erbrlatlong[,c("Long","Lat")])    #Extract temp for given gps coord & month, sum to prev month
  }

  for (pp in 9:12) {       #Apr-Jul
    file.pos <- grep(yyyymm[pp,rr], bils.temp)
    raster_file <- raster(bils.temp[file.pos])                                            #Load specified raster
    temp.summer <- temp.summer + raster::extract(x = raster_file, erbrlatlong[,c("Long","Lat")])    #Extract temp for given gps coord & month, sum to prev month
  }

  temp[rr,2] <-  temp.fall/4                                            #Calc mean seasonal temp
  temp[rr,3] <-  temp.winter/4                                          #Calc mean seasonal temp
  temp[rr,4] <-  temp.summer/4                                          #Calc mean seasonal temp
}
## ------------------------------------------------------------------------------

## column for each season and ppt and mean temp


## Update 2022-11-09 Michelle DePrenger-Levin, corrected current and previous year to split at 1:7 as current, 8:12 previous
clim <- erbr.climate.monthly %>%
  left_join(erbr.climate.monthly.tmean, by = "Year") %>%
  filter(Year < 2023 & Year > 2003)
clim <- cbind(erbr.climate.monthly, erbr.climate.monthly.tmean[,c(2:4)])
save(clim, file = paste("C:/Users/deprengm/OneDrive - Denver Botanic Gardens/P drive/hackathon/ErBr/erbr_climData3seas",
                      Sys.Date(),".Rdata", sep=""))
write.csv(clim, paste("C:/Users/deprengm/OneDrive - Denver Botanic Gardens/P drive/hackathon/ErBr/erbr_climData3seas",
                      Sys.Date(),".csv", sep=""), row.names = FALSE)

## SAVE CLIMATE DATA --------------------------------------------------------------
## Combine different climate variables into one data frame
clim <- cbind(ppt, temp[,2:4])

## Save output
write.csv(clim, "erbr_climData3seas30yr_210209.csv", row.names=FALSE)
## -------------------------------------------------------------------------------

