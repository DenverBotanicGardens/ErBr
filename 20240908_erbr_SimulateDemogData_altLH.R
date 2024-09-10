## April Goebl & Dan Doak
## Script modified 2024-09-04
## Collaboration with Denver Botanic Gardens on Eriogonum brandegeii modeling 
## Script written in response to reviews on Goebl et al 2024 Biol Conservation 
## Generate simulated data on lives of individual plants with known demographic parameters 
## Modify 'life history' to get additional simulated datasets



rm(list=ls())
setwd("C:/Users/Dan Doak/Desktop/Students/April/eriogonum models/manu fall2024/life histories")


## SET WD -----------------------------------------------------------------------------------------
#setwd("C:/Users/Dan Doak/Desktop/Students/April/eriogonum models/manu fall2024")
## ------------------------------------------------------------------------------------------------
#setwd("C:/Users/Dan Doak/Desktop/Students/April/eriogonum models/manu fall2024/life histories")


## LOAD DATA --------------------------------------------------------------------------------------
clim32yr <- read.csv("erbr_climData3seas32yr_221114.csv", header=TRUE)
erbr <- read.csv("erbr_TagClust2022_20230408.csv", header=TRUE)
erbr$Year <- as.factor(erbr$Year)
JAGSmodSumm_realDat <- readRDS("erbr_JAGSmodBestSUMM_c3t10s10b5_noYRE_20240829.rds")
mx.mean <- readRDS("erbrMeanMatrix_noYRE_P1k_20240715") #Load mean matrix variable
medParams <- readRDS("erbrMedParams_noYRE_20240803")    #Load median parameters variable 
## ------------------------------------------------------------------------------------------------



## LOAD PACKAGES ----------------------------------------------------------------------------------
library(dplyr)
library(matrixStats)
library(stringr)
## ------------------------------------------------------------------------------------------------




## GET MODIFIED VITAL RATE PARAMS AND MX TO SIMULATE DATA WITH ALTERNATIVE LIFE HISTORIES --------

## Fix each of the climate drivers at their mean values
## Use the real data from the study 
rel.yrs <- 2002:2021 #Select relevant subset of yrs (this matches what was used for stoch lam to keep Methods consistent) 
climYrs.rel <- clim32yr[clim32yr$Year>=2002 & clim32yr$Year<=2021,]
mn.clim <- as.matrix(colMeans(climYrs.rel[,-1]))
mn.clim <- as.data.frame(t(mn.clim))
## ------------------------------------------------------------------------------------------------


## Use the SSD from a mean matrix ---------------------------------------------------------------
## Obtain stable stage structure 
N.startNum <- 181                 #Set starting pop sz as 2007 obs sz

SSD=eigen(mx.mean)$vectors[,1]
SSD=Re(SSD/sum(SSD))
N.vecStart=N.startNum * SSD

# Specify min and max plt sz
minsize <- 1
maxsize <- (max(erbr$RosNew, na.rm=TRUE)) 
#maxsize <- 300  #For Fast Growth LH datasets ******


## For median sz estimation
## new size density estimation for median size estimation
pdfsz=density(erbr$RosNew, n=1024, cut=0, na.rm=TRUE) 
pdfsz2=cbind(pdfsz$x,pdfsz$y)
## This is a set of smoothed values that can then be used w weightedMedian in the matrixStats package to get a 'good' median for each class

bin.num <- 50  #Define number of bins 
  
## Improved method of finding median size/ bin mids (code from Dan)
vec.bin = c(minsize, minsize+1:bin.num*(maxsize-minsize)*(1/bin.num)) 
## Do this block to make medians the focal estimated size for each category
binmids = rep(NA, length(vec.bin)-1)

for(jj in 1:(length(vec.bin)-1)) {
## Set limits for subset according to bin breaks
    bounds <- c(vec.bin[jj], vec.bin[jj+1])
## Subset data according to bounds
    subsetszs <- pdfsz2[which(pdfsz2[,1] >= bounds[1] & pdfsz2[,1] < bounds[2]),]
    binmids[jj] <- weightedMedian(subsetszs[,1],subsetszs[,2])
  }

n.bin = length(binmids)



## ------------ Assign param values for intercepts & sz to generate alternative life histories -----------
# ####################################################################
# parameters for a fast growth, med survival, higher repro LH:
  medParams$grwth_RosCoef = medParams$grwth_RosCoef*1.1  
  medParams$grwth_intercept = 1.1   
# #
  medParams$surv_intercept=medParams$surv_intercept*1.2     
  medParams$surv_RosCoef=medParams$surv_RosCoef*0.8         
# #
  medParams$repro_intercept=medParams$repro_intercept+2 
  medParams$reproyesno_intercept = medParams$reproyesno_intercept+2
  medParams$reproyesno_RosCoef=medParams$reproyesno_RosCoef*1.2
# #######################################################

####################################################################
# parameters for a med growth, lower survival LH, high repro:
medParams$grwth_RosCoef = medParams$grwth_RosCoef*1.025
medParams$grwth_intercept = 1.025

medParams$surv_intercept=medParams$surv_intercept*1.3
medParams$surv_RosCoef=medParams$surv_RosCoef*0.5

medParams$repro_intercept=medParams$repro_intercept+5.5 #5.5
medParams$reproyesno_intercept = medParams$reproyesno_intercept+2
medParams$reproyesno_RosCoef=medParams$reproyesno_RosCoef*1.2
#######################################################

 



## Plug in median param vals, mean climate and starting plant size values into model formulas

    ## Survival (binomial) --
    pred.surv <- 1/(1+exp(-(medParams$surv_intercept + medParams$surv_RosCoef*log(binmids) + 
                              medParams$surv_PptWinterCoef*mn.clim$Tot_winter_ppt + 
                              medParams$surv_TempWinterCoef*mn.clim$Mean_winter_temp +
                              medParams$surv_TempSummerCoef*mn.clim$Mean_summer_temp +
                              medParams$surv_TempFallCoef*mn.clim$Mean_fall_temp)))
    
    ## Growth (negative binomial)--   
    pred.grwth <- exp(medParams$grwth_intercept + medParams$grwth_RosCoef*log(binmids) 
                          + medParams$grwth_TempFallCoef*mn.clim$Mean_fall_temp
                          + medParams$grwth_TempSummerCoef*mn.clim$Mean_summer_temp
                          + medParams$grwth_TempWinterCoef*mn.clim$Mean_winter_temp
                          + medParams$grwth_PptFallCoef*mn.clim$Tot_fall_ppt
                          + medParams$grwth_PptSummerCoef*mn.clim$Tot_summer_ppt
                          + medParams$grwth_PptWinterCoef*mn.clim$Tot_winter_ppt)
        
    pred.grwthVar <- exp(medParams$grwthvar_intercept + medParams$grwthvar_RosCoef*log(binmids)) 
        
    ## Probability of reproduction (binomial) --   
    pred.reproYesNo <- 1/(1+exp(-(medParams$reproyesno_intercept + medParams$reproyesno_RosCoef*log(binmids) +
                                      medParams$reproyesno_TempFallCoef*mn.clim$Mean_fall_temp +
                                      medParams$reproyesno_PptFallCoef*mn.clim$Tot_fall_ppt +
                                      medParams$reproyesno_PptSummerCoef*mn.clim$Tot_summer_ppt +
                                      medParams$reproyesno_TempSummerCoef*mn.clim$Mean_summer_temp +
                                      medParams$reproyesno_TempWinterCoef*mn.clim$Mean_winter_temp))) 
    
    ## Reproduction (negative binomial)
    pred.repro <- exp(medParams$repro_intercept + medParams$repro_RosCoef*log(binmids) + 
                                medParams$repro_TempFallCoef*mn.clim$Mean_fall_temp +
                                medParams$repro_PptFallCoef*mn.clim$Tot_fall_ppt +
                                medParams$repro_PptSummerCoef*mn.clim$Tot_summer_ppt +
                                medParams$repro_TempSummerCoef*mn.clim$Mean_summer_temp +
                                medParams$repro_TempWinterCoef*mn.clim$Mean_winter_temp)
          
      ## Seedling survival (binom)  
      sz.sdlg <- 1
      pred.survSdlg <- 1/(1+exp(-(medParams$surv_intercept + medParams$surv_RosCoef*log(sz.sdlg) + 
                                        medParams$surv_PptWinterCoef*mn.clim$Tot_winter_ppt + 
                                        medParams$surv_TempWinterCoef*mn.clim$Mean_winter_temp +
                                        medParams$surv_TempSummerCoef*mn.clim$Mean_summer_temp +
                                        medParams$surv_TempFallCoef*mn.clim$Mean_fall_temp)))
          
      ## Seedlings per inflor (neg binom)
      pred.numSdlg <- exp(medParams$newplt_intercept + log(pred.repro))
          
          
          ## Constructing matrix models
          grwth.mx <- matrix(data=NA, nrow=n.bin, ncol=n.bin)
          
          ## growth probs using cdf fn
          for (ss in 1:(n.bin)) {
            grwth.cdf <- pnorm(vec.bin, pred.grwth[ss], sqrt(pred.grwthVar[ss]))
            grwth <- grwth.cdf[2:length(vec.bin)] - grwth.cdf[1:(length(vec.bin)-1)]
            if (sum(grwth)>0) {   #If statement breaks code (puts NA's into mx) if sum of PDF is 0 (which happens if all prob is outside of sz bounds)
              grwth <- grwth/sum(grwth)
              grwth.mx[,ss] <- grwth
            } else {
              grwth.mx[,ss] <- NA
            } 
          } #End ss loop
          
          ## Make survival * growth matrix
          surv.grwth.mx <- grwth.mx * t(matrix(rep(pred.surv,(n.bin)),(n.bin)))
          mx1 <- surv.grwth.mx #growth and survival, without repro
          ## Add reproduction and recruitment
          mx <- matrix(0, (n.bin+1), (n.bin+1))
          mx[2:(n.bin+1), 2:(n.bin+1)] <- mx1
          mx[2,1] <- pred.survSdlg                                  #First column (seedling survival in element 2,1)
          mx[1,2:(n.bin+1)] <- pred.reproYesNo * pred.numSdlg       #First row (new seedlings)
          
          lam <- Re(eigen(mx)$values[1])              #Calculate lambda & store for each transect
          #mx.out.list[[length(mx.out.list) + 1]] <- mx
          


#plot(binmids,pred.grwth)
#          abline(a=0,b=1)
#          plot(binmids,pred.grwth, log='xy')
#          abline(a=0,b=1)
plot(binmids,pred.surv)      
#plot(binmids,pred.repro)

#print(lam) 
SSD=eigen(mx)$vectors[,1] 
plot(binmids,SSD[2:51])
## -------------------------------------------------------------------------------------------------

## Save mx for future runs 
#saveRDS(mx, file=paste("20240908", "_erbrMatrix_medGrLH", ".rds", sep=""))
mx.newLH <- readRDS("20240908_erbrMatrix_medGrLH.rds") #Load mean matrix variable
#saveRDS(mx, file=paste("20240908", "_erbrMatrix_hiLH", ".rds", sep=""))
#mx.newLH <- readRDS("20240908_erbrMatrix_hiLH.rds") #Load mean matrix variable

mx.newLH <- mx       
## -------------------------------------------------------------------------------------------------
          











## GENERATE SIMULATED DATA -------------------------------------------------------------------------         
## Start data set loop 
name <- as.character("SimDat50yrMedGrLH.") # For naming saved files below
n.datset <- 10
for (dd in 1:n.datset) {

print(dd)
  
  ## For a set number of years (e.g. 20), simulate climate variables for each year. Use the real data
  n.yrs <- 51 #Assign number of years (plus 1) to simulate climate data for 
  
  #Create empty variable to hold simulated climate data
  column.names <- colnames(clim32yr)
  sim.clim <- as.data.frame(matrix(NA, nrow=length(1:n.yrs), ncol=length(column.names)))
  colnames(sim.clim) <- column.names
  sim.clim$Year <- 1:n.yrs
   
  #List of random numbers that corresponds to a set of climate values 
  rel.yrs <- 2002:2021 #Select relevant subset of yrs (this matches what was used for stoch lam to keep Methods consistent) 
  climYrs.rel <- clim32yr[clim32yr$Year>=2002 & clim32yr$Year<=2021,]
  randVals.yr <- sample(1:nrow(climYrs.rel), size=n.yrs, prob=NULL, replace=TRUE)
  
  for (cc in 1:length(1:n.yrs)) {
    sim.clim[cc,2:7] <- climYrs.rel[randVals.yr[cc],2:7] 
    sim.clim$Clim_yr[cc] <- climYrs.rel[randVals.yr[cc],1]
  }
  ## ------------------------------------------------------------------------------------------------



  ## Use the SSD from a modified LH matrix (from above) to pick a set of plant starting sizes: 
  ## Use these as a set of multinomial probabilities for the starting size of each plant that is then simulated
  
  ## Obtain stable stage structure
  popSz.start <- 181                 #Set starting pop sz as 2007 obs sz
  
  SSD=eigen(mx.newLH)$vectors[,1]  #USE MATRIX GENERATED ABOVE
  SSD=Re(SSD/sum(SSD))
  N.vecStart=popSz.start * SSD
  
  # Specify min and max plt sz
  minsize <- 1
  maxsize <- (max(erbr$RosNew, na.rm=TRUE)) 
  
  ## For median sz estimation
  ##new size density estimation for median size estimation
  pdfsz=density(erbr$RosNew, n=1024, cut=0, na.rm=TRUE) 
  pdfsz2=cbind(pdfsz$x,pdfsz$y)
  ## This is a set of smoothed values that can then be used w weightedMedian in the matrixStats package to get a 'good' median for each class
  
  n.bin <- 50  #Define number of bins 
    
  ## Improved method of finding median size/ bin mids (code from Dan)
  vec.bin = c(minsize, minsize+1:n.bin*(maxsize-minsize)*(1/n.bin)) 
  ## Do this block to make medians the focal estimated size for each category
  binmids = rep(NA, length(vec.bin)-1)
  
  for(jj in 1:(length(vec.bin)-1)) {
  ## Set limits for subset according to bin breaks
      bounds <- c(vec.bin[jj], vec.bin[jj+1])
  ## Subset data according to bounds
      subsetszs <- pdfsz2[which(pdfsz2[,1] >= bounds[1] & pdfsz2[,1] < bounds[2]),]
      binmids[jj] <- weightedMedian(subsetszs[,1],subsetszs[,2])
    }
  
  binmids <- c(1, binmids)  
    
  ## Select starting sizes of plants for simulated data using SSD and median sz classes
  n.startPlts <- 200 #Num of starting plts
  N.startProbs <- (N.vecStart*1) / popSz.start
  sum(N.startProbs) #Should equal 1
  sz.startPlts <- sample(x=binmids, size=n.startPlts, replace=TRUE, prob=N.startProbs)  
  sz.startPlts <- round(sz.startPlts, digits=0)   #Round to nearest integer 
  ## ------------------------------------------------------------------------------------------------

  
  
  
  #Set dispersion param variables from runs of real data
  JAGSmodSumm_realDat   #Median values from JAGS run of real dat
  r.inf <- JAGSmodSumm_realDat[1,2]      
  r.sdlg <- JAGSmodSumm_realDat[2,2]
  
  
  #Matrices to hold sz & repro where rows are yrs and columns are plants
  #We will be getting and saving each plt's own fate in each yr w the vital rate functions
  mx.sz <- as.data.frame(matrix(NA, nrow=length(1:(n.yrs-1)), ncol=length(1:n.startPlts)))
  mx.reproInf <- as.data.frame(matrix(NA, nrow=length(1:(n.yrs-1)), ncol=length(1:n.startPlts)))
  mx.reproSdlg <- as.data.frame(matrix(NA, nrow=length(1:(n.yrs-1)), ncol=length(1:n.startPlts))) 

  
  #for 1 to the num of starting plants 
  for (pp in 1:n.startPlts) {  #Loop over starting plants
  
    sel.plt <- sz.startPlts[pp]
    
  
      for (yy in 2:(n.yrs-1)) {  #Loop over years
    
        ## Survival (binomial) --
        pred.surv <- 1/(1+exp(-(medParams$surv_intercept + medParams$surv_RosCoef*log(sel.plt) + 
                                  medParams$surv_PptWinterCoef*sim.clim[yy+1,]$Tot_winter_ppt + 
                                  medParams$surv_TempWinterCoef*sim.clim[yy+1,]$Mean_winter_temp +
                                  medParams$surv_TempSummerCoef*sim.clim[yy+1,]$Mean_summer_temp +
                                  medParams$surv_TempFallCoef*sim.clim[yy+1,]$Mean_fall_temp)))
        
        realzd.surv <- rbinom(n=1, size=1, prob=pred.surv)
        
        ##If survived, keep going, if realzd.surv=0, add zero to matrices, and end current loop (using break statement)
          if (realzd.surv==0) {
            mx.sz[yy,pp] <- 0            
            mx.sz[yy+1,pp] <- 0            #Add 2nd yr w zero following death so that if dead in missing data yr, death is recorded
            mx.reproInf[yy,pp] <- 0 
            mx.reproInf[yy+1,pp] <- 0 
            mx.reproSdlg[yy+1,pp] <- 0     
            break } 
        ## --
    
    
          ## Growth (negative binomial)--   
          pred.grwth <- exp(medParams$grwth_intercept + medParams$grwth_RosCoef*log(sel.plt) 
                            + medParams$grwth_TempFallCoef*sim.clim[yy+1,]$Mean_fall_temp
                            + medParams$grwth_TempSummerCoef*sim.clim[yy+1,]$Mean_summer_temp
                            + medParams$grwth_TempWinterCoef*sim.clim[yy+1,]$Mean_winter_temp
                            + medParams$grwth_PptFallCoef*sim.clim[yy+1,]$Tot_fall_ppt
                            + medParams$grwth_PptSummerCoef*sim.clim[yy+1,]$Tot_summer_ppt
                            + medParams$grwth_PptWinterCoef*sim.clim[yy+1,]$Tot_winter_ppt)
          
          pred.grwthVar <- exp(medParams$grwthvar_intercept + medParams$grwthvar_RosCoef*log(sel.plt)) 
          
          realzd.grwth <- rnorm(n=1, mean=pred.grwth, sd=sqrt(pred.grwthVar)) 
          realzd.grwth <- round(realzd.grwth, digits=0)   #Round to nearest integer since unit is rosettes
          ## Bound by largest and smallest sz class
          if (realzd.grwth < minsize) {   #If sz is equal to or less than 0, change to 1 (smallest sz) 
            realzd.grwth <- minsize
          }
          if (realzd.grwth > maxsize) {   #If sz is over upper sz bound, change to largest sz
            realzd.grwth <- maxsize
          }
          
          
          ## Enter realized size into size matrix
          mx.sz[yy,pp] <- realzd.grwth 
          ## --
  
          
          ## Probability of reproduction (binomial) --   
          pred.reproYesNo <- 1/(1+exp(-(medParams$reproyesno_intercept + medParams$reproyesno_RosCoef*log(sel.plt) +
                                        medParams$reproyesno_TempFallCoef*sim.clim[yy,]$Mean_fall_temp +
                                        medParams$reproyesno_PptFallCoef*sim.clim[yy,]$Tot_fall_ppt +
                                        medParams$reproyesno_PptSummerCoef*sim.clim[yy,]$Tot_summer_ppt +
                                        medParams$reproyesno_TempSummerCoef*sim.clim[yy,]$Mean_summer_temp +
                                        medParams$reproyesno_TempWinterCoef*sim.clim[yy,]$Mean_winter_temp))) 
      
            realzd.reproYesNo <- rbinom(n=1, size=1, prob=pred.reproYesNo)
           
            #If realized reproYesNo is 1, then continue. If 0 then enter 0 in repro matrix
            if (realzd.reproYesNo==0) {
              mx.reproInf[yy,pp] <- 0 
              mx.reproSdlg[yy+1,pp] <- 0 
              } else {
      
                
                
                ## Reproduction (negative binomial)
                pred.repro <- exp(medParams$repro_intercept + medParams$repro_RosCoef*log(sel.plt) + 
                                  medParams$repro_TempFallCoef*sim.clim[yy,]$Mean_fall_temp +
                                  medParams$repro_PptFallCoef*sim.clim[yy,]$Tot_fall_ppt +
                                  medParams$repro_PptSummerCoef*sim.clim[yy,]$Tot_summer_ppt +
                                  medParams$repro_TempSummerCoef*sim.clim[yy,]$Mean_summer_temp +
                                  medParams$repro_TempWinterCoef*sim.clim[yy,]$Mean_winter_temp)
                
                #This should be a random from a neg binomial, as that is how it is fit.  
                #Pull out the fitted r.inf variable and fit this way
                realzd.repro <- rnbinom(n=1, size=r.inf, mu=pred.repro)  
                if (realzd.repro < 1) {   #If num inflors less than 1, change to 1 (min number of infs when reproducing)
                  realzd.repro <- 1
                } 
            
                #Enter inf data into repro matrix
                mx.reproInf[yy,pp] <- realzd.repro
                
                # Seedlings (negative binomial)
                pred.numSdlg <- exp(medParams$newplt_intercept + log(realzd.repro))  
                realzd.numSdlg <- rnbinom(n=1, size=r.sdlg, mu=pred.numSdlg)       
                
                #Enter seedling data into seedling matrix
                mx.reproSdlg[yy+1,pp] <- realzd.numSdlg  } 
          ## --
        
            #Then, if the plt is still alive, you go to the next year and do the same.
          
  } ##End year loop
  
}   ##End starting number of plants loop


  ## Remove final rows with NAs, that get added if death occurs in last yr
  mx.sz <- mx.sz[1:(n.yrs-1),] 
  mx.reproInf <- mx.reproInf[1:(n.yrs-1),] 
  mx.reproSdlg <- mx.reproSdlg[1:(n.yrs-1),] 
  
  ## Add starting plant sizes to row 1
  mx.sz[1,] <- sz.startPlts
  
  
  ## Year 1 is starting plant sizes, there is no repro data 
  ## Year 2 is the 'realized' size and inflor number based on the starting size and climate 
  ## New seedlings can appear in year 3 or later based on the inflor numbers from year 2
  
  ## In size matrix, once dead, size is 0 for year died plus next year, and then NAs afterwards 
  ## In repro (inflors) matrix, once dead, repro is 0 for yr died plus next year, and then NAs afterwards
  ## 0 in this repro matrix can also mean there was no repro in that yr, in which case, non-NA and non-zero numbers follow 
  ## In repro (seedling) matrix, once dead, seedlings are 0 for yr died and NAs afterwards
  ## 0 in this seedling matrix can also mean no seedlings in that yr, in which case, non-NA values (0s or integers) follow 
  ## -------------------------------------------------------------------------------------------------

 #Lines 427-431:this makes the counter for needed seedlings each year:
  mx.survivors=mx.sz
  mx.survivors[is.na(mx.survivors)==FALSE]=1
  n.survivors=rowSums(mx.survivors,na.rm=TRUE)
  n.sdlgs.needed=200-n.survivors


  ## 4. After doing this, go back to the number of new plants in each year: for these, do the same approach as with the starting plts, but starting in the yr they are 'born'
  ## and then simulate them going forward. Use NA for the sz in the years before they are 'born'
  
  # Subset repro matrix to only contain cols (plts) with non-zero number of seedlings. Then loop thru that.
 # mx.sdlgYes <- mx.reproSdlg %>% select_if(~ sum(., na.rm=TRUE) > 0)  

  sdlg <- 1            #set size to be 1 ros for seedlings
  


    #for (rr in 1:ncol(mx.sdlgYes)) { #Loop over parent plants
    for (rr in 1:length(n.sdlgs.needed)) { #Loop over num seedlings needed
  
        realzd.grwth <- 1
      
        #if (n.sdlgs.needed[rr]==0) {
        #  break } 
        
        #sdlg.strtYrTemp <- which(mx.sdlgYes[,rr]>0)        #Find years when new seedlings appeared 
        #sdlg.perYr <- mx.sdlgYes[,rr][mx.sdlgYes[,rr]>0 & !is.na(mx.sdlgYes[,rr])]  #Find new sdlgs per year, if non-zero
        #sdlg.strtYr <- rep(sdlg.strtYrTemp, sdlg.perYr)     #List starting years for all new sdlgs, including if more than 1 in given yr
        sdlg.perYr <- n.sdlgs.needed[rr]
        colCount <- ncol(mx.sz)                             #Re-set column counter each loop based on new cols added in year loop
    
        
        ##Loop over number of new seedlings that a given parent plt had in 1 or more years
        #for (ss in 1:length(sdlg.strtYr)) { 
        for (ss in 1:sdlg.perYr) { 
            
          
          #startYr <- sdlg.strtYr[ss]
          startYr <- rr
          mx.sz[startYr,colCount+ss] <- sdlg  
          
    
          for (yy in startYr:(n.yrs-2)) {     #For each new seedling, loop over years
            
              plt.sz <- realzd.grwth * sdlg
            
              
              ## Survival (binomial) --  
              pred.surv <- 1/(1+exp(-(medParams$surv_intercept + medParams$surv_RosCoef*log(plt.sz) + 
                                            medParams$surv_PptWinterCoef*sim.clim[yy+1,]$Tot_winter_ppt + 
                                            medParams$surv_TempWinterCoef*sim.clim[yy+1,]$Mean_winter_temp +
                                            medParams$surv_TempSummerCoef*sim.clim[yy+1,]$Mean_summer_temp +
                                            medParams$surv_TempFallCoef*sim.clim[yy+1,]$Mean_fall_temp)))
              
              realzd.surv <- rbinom(n=1, size=1, prob=pred.surv)
            
              ##If survived, keep going, if realzd.survSdlg=0, add zero to matrix, and end current loop (using break statement)
              if (realzd.surv==0) {
                mx.sz[yy+1,colCount+ss] <- 0   
                mx.sz[yy+2,pp] <- 0            #Add 2nd yr w zero following death so that if dead in missing data yr, death is recorded
                break }  
              ## --
              
            
          
              ## Growth (negative binomial) --   
              pred.grwth <- exp(medParams$grwth_intercept + medParams$grwth_RosCoef*log(plt.sz)  
                              + medParams$grwth_TempFallCoef*sim.clim[yy+1,]$Mean_fall_temp
                              + medParams$grwth_TempSummerCoef*sim.clim[yy+1,]$Mean_summer_temp
                              + medParams$grwth_TempWinterCoef*sim.clim[yy+1,]$Mean_winter_temp
                              + medParams$grwth_PptFallCoef*sim.clim[yy+1,]$Tot_fall_ppt
                              + medParams$grwth_PptSummerCoef*sim.clim[yy+1,]$Tot_summer_ppt
                              + medParams$grwth_PptWinterCoef*sim.clim[yy+1,]$Tot_winter_ppt)
            
              pred.grwthVar <- exp(medParams$grwthvar_intercept + medParams$grwthvar_RosCoef*log(plt.sz)) 
            
              realzd.grwth <- rnorm(n=1, mean=pred.grwth, sd=sqrt(pred.grwthVar)) 
              realzd.grwth <- round(realzd.grwth, digits=0)   #Round to nearest integer since unit is rosettes
              
              ## Bound by largest and smallest sz class
              if (realzd.grwth < minsize) {   #If sz is equal to or less than 0, change to 1 (smallest sz) 
                realzd.grwth <- minsize
              }
              if (realzd.grwth > maxsize) {   #If sz is over upper sz bound, change to largest sz
                realzd.grwth <- maxsize
              }
            
          
            ##Enter realized size into size matrix
            mx.sz[yy+1,colCount+ss] <- realzd.grwth 
            ## --
            
          }  ## End year loop
        }    ## End new seedlings loop
    }        ## End parent plant loop

  ##Reproduction of new seedlings/ plants is not included here 
  ##We followed seedlings resulting from original plants but not subsequent seedlings that may have been 'born'

  ## Remove final rows with NAs, that get added if death occurs in last yr
  mx.sz <- mx.sz[1:(n.yrs-1),] 
  
  
  
  ## Add a year column
  mx.szWyr <- cbind.data.frame(1:(n.yrs-1), sim.clim$Clim_yr[(2:n.yrs)], mx.sz[1:(n.yrs-1),])
  colnames(mx.szWyr) <- c("Year", "Clim_yr", colnames(mx.sz))
  
  mx.reproWyr <- cbind.data.frame(1:(n.yrs-1), sim.clim$Clim_yr[(2:n.yrs)], mx.reproInf[1:(n.yrs-1),])
  colnames(mx.reproWyr) <- c("Year", "Clim_yr", colnames(mx.reproInf))
  
  ## Add columns to repro matrix so number matches that of size mx that represent new seedlings. Values should be 0 for no repro 
  num.addnCols <- ncol(mx.szWyr) - ncol(mx.reproWyr)
  addnCols <- as.data.frame(matrix(0, nrow=nrow(mx.reproWyr), ncol=num.addnCols))
  mx.reproWyrAdCol <- cbind(mx.reproWyr, addnCols)
  
  
  
  ## Modify output to match format of raw data so can be appropriately modified for use with JAGs code **
  ## 'Year' column can be e.g. 1-20 
  ## 'Clim_yr' column is e.g. 2004-2022
  ## 'Transect' column should be entered as 'TransectNew' where individuals are assigned to E.1-E.7 and W.1-W.5
  ## 'TagNew' column should be added where indivs have a unique tag ID that is E.### or W.###
  ## Therefore columns to have: 'Year', 'TransectNew', 'TagNew', 'RosNew', 'InflNew'
  ## Change output from above so there are rows for all individuals across all year (rows are years and plants)
  
  ## Combine all plants into 1 column, where all years for a given plt appear sequentially
  datComb <- NULL
  
  for (ll in 3:ncol(mx.szWyr)) {   
    temp <- cbind(mx.szWyr$Year, mx.szWyr$Clim_yr, mx.szWyr[,ll], mx.reproWyrAdCol[,ll])
    datComb <- rbind(datComb, temp)
  }
  
  datComb <- as.data.frame(datComb)
  colnames(datComb) <- c("Year", "ClimYr", "RosNew", "InflNew") 
  
  ## Assign plants to one of the following transects (random but equal distribution)
  tran <- c("E.1","E.2","E.3","E.4","E.5","E.6","E.7","W.1","W.2","W.3","W.4","W.5")
  num.plt <- ncol(mx.sz)
  tran.rep <- rep(tran, (ceiling(num.plt/length(tran))))
  tran.rep <- tran.rep[1:num.plt]
  tran.rand <- sample(tran.rep, replace=FALSE)
  datComb$TransectNew <- NULL
  datComb$TransectNew <- rep(tran.rand, each=(n.yrs-1), replace=FALSE)
  
  ## Assign unique tag ID to each plant 
  tag <- c(1:num.plt)
  tag.rep <- rep(tag, each=(n.yrs-1), replace=FALSE)
  datComb$TagNew <- paste(datComb$TransectNew, tag.rep, sep='.') #Tag new included site and transect; same format as real data
  
  
## --------------------------------------------------------------------------------------------
  ## save as csv 
  date <- Sys.Date()                                        #Enter date to be added to file name
  date <- str_replace_all(date, "-", "")
 # name <- as.character("SimDat20yrNoMissMedGrLH.")           #Enter name of file, e.g. Tagclust, 4to13, simulated data...
  
  write.csv(datComb, file=paste(date, "_erbr_", name, "NoMiss.", dd, ".csv", sep=""), row.names=FALSE)
  print(paste(date, name, dd, ".csv", sep=""))
  ## -----------------------------------------------------------------------------------------------------
  
  
  
  
  
  ## Modify output to contain missing years of data 
  ## Get the records for each individual, but then if a year is one of the missing data years, set that years data to NA 
  datComb1 <- datComb
  
  ## Assign years to be missing
  startConsecYrs <- 10
  yrs.missing <- as.integer(seq((startConsecYrs+1),(n.yrs-1), by=2))
  
  ## Change missing years to NA
  datComb1$RosNew[datComb1$Year %in% (yrs.missing)] <- NA 
  datComb1$InflNew[datComb1$Year %in% (yrs.missing)] <- NA 
  ## -----------------------------------------------------------------------------------------------------
  
  
  
  
  ## ADD CODE FROM erbr_1ReformatSIMdata_forJAGS SCRIPT -------------------------------------------------------
  
  
  ## LOAD PACKAGES AND FUNCTIONS --------------------------------------------------------------------
  library(dplyr)
  library(bazar)
  library(stringr) 
  library(tidyr)   
  ## ------------------------------------------------------------------------------------------------
  
  
  
  ## MODIFY FORM OF DATA ----------------------------------------------------------------------------
  ## Change zeros in Rosettes & Infls to NAs in RosNew & InflNew columns to indicate dead
  datComb1$RosNew[datComb1$RosNew==0] <- NA
  datComb1[datComb1$RosNew==0 & !is.na(datComb1$RosNew),] #Confirm that no RosNew=0
  datComb1$InflNew[datComb1$InflNew==0] <- NA
  datComb1$InflNew[is.na(datComb1$RosNew) & !is.na(datComb1$InflNew)] #Confirm
  datComb1$InflNew[is.na(datComb1$RosNew)] <- NA
  ## Change Infl to zero from NA if Rosettes has data 
  datComb1$InflNew[datComb1$RosNew>0 & !is.na(datComb1$RosNew) & is.na(datComb1$InflNew)] <- 0
  ## ------------------------------------------------------------------------------------------------
  
  
  
  ## Confirm no rows are duplicates in terms of TagNew and Year values
  datComb1[duplicated(datComb1[,c("TagNew","Year")]),]
  
  
  ## How many observed individuals (tag clusters) were there each year
  indivXyear <- datComb1 %>%
    dplyr::group_by(Year) %>%
    dplyr::summarise(Indivs = n_distinct(TagNew[RosNew > 0]))
  ## -----------------------------------------------------------------------------------
  
  
  
  
  ## START OF CHANGES NEEDED FOR JAGS --------------------------------------------------
  
  ## Make a column that indicates if a line should be kept: removes lines w NAs before a plt appeared or after died. Keeps lines that are NAs but bracketed by sz data yrs
  dats <- datComb1      #Placeholder for the new data
  dats$save <- 0        #Start with save = 0, or no
  dats$surv <- NA       #Column to show plt survival/ if plant is alive in current year
  
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
    left_join(clim32yr, by = c("ClimYr" = "Year"))
  ## -----------------------------------------------------------------------------------
  
  
  
  ## ADD PROBABILITY OF REPRO RESPONSE VARIABLE ----------------------------------------
  ## Determine whether or not reproduction occurred
  erbr.1$InflYesNo <- NA
  erbr.1$InflYesNo[erbr.1$InflNew > 0] <- 1
  erbr.1$InflYesNo[erbr.1$InflNew == 0] <- 0
  ## -----------------------------------------------------------------------------------
  
  
  
  ## SAVE FORMATTED DATA ---------------------------------------------------------------
  date <- Sys.Date()                                #Enter date to be added to file name
  date <- str_replace_all(date, "-", "")
#  name <- as.character("SimDat20yrMissMedGrLH.")           #Enter name of file, e.g. Tagclust, 4to13, simulated data...
  
  write.csv(erbr.1, file=paste(date, "_erbr_", name, "Miss.", dd, ".Format4JAGS", ".csv", sep=""), row.names=FALSE)
  print(paste(date, name, dd, ".Format4JAGS", ".csv", sep=""))
  ## -----------------------------------------------------------------------------------
  
  
  
}  ##End dataset loop

plot(erbr.1$Year,erbr.1$RosNew)

boxplot(erbr.1$RosNew~erbr.1$Year)
rowSums(mx.sdlgYes,na.rm=TRUE)

table(erbr.1$Year)
