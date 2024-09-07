## April Goebl & Dan Doak
## Script modified 2024-07-31
## Collaboration with Denver Botanic Gardens on Eriogonum brandegeii modeling 
## Script written in response to reviews on Goebl et al 2024 Biol Conservation 
## Generate simulated data on lives of individual plants with known demographic parameters 



rm(list=ls())



## SET WD -----------------------------------------------------------------------------------------
#setwd("C:/Users/Dan Doak/Desktop/Students/April/eriogonum models/manu fall2024")
## ------------------------------------------------------------------------------------------------



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



## Start data set loop here
n.datset <- 10
for (dd in 1:n.datset) {


  ## 1. First, for a set number of years (lets say 20), simulate climate variables for each year. 
  ## Use the real data from the study to choose sets of annual data values for the set of climate variables.
  n.yrs <- 21 #Assign number of years (plus 1) to simulate climate data for 
  
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




  ## 2. Then, use the SSD from a mean matrix from real data to pick a set of plant starting sizes: 
  ## you can actually just use these as a set of multinomial probabilities for the starting size of each plant that is then simulated.
  
  ## Obtain stable stage structure; modified from some of Dan's code 
  popSz.start <- 181                 #Set starting pop sz as 2007 obs sz
  
  SSD=eigen(mx.mean)$vectors[,1]
  SSD=Re(SSD/sum(SSD))
  N.vecStart=popSz.start * SSD
  
  
  # Specify min and max plt sz
  minsize <- 1
  maxsize <- (max(erbr$RosNew, na.rm=TRUE)) 
  
  ## For median sz estimation
  ## From Dan's bistorts to guppies code 
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
  n.startPlts <- 200 #Num of starting plts. Set much larger later (e.g. 500?) (make this considerably larger than you think you want to have)
  N.startProbs <- (N.vecStart*1) / popSz.start
  sum(N.startProbs) #Should equal 1
  sz.startPlts <- sample(x=binmids, size=n.startPlts, replace=TRUE, prob=N.startProbs)  
  sz.startPlts <- round(sz.startPlts, digits=0)   #Round to nearest integer 
  ## ------------------------------------------------------------------------------------------------

  


  ## 3. 
  
  #At the start of this code, make vectors of the parameter values for each of the vital rate functions 
  #these are the 'rules' that will determine each vital rate's predicted value for each plant in each year. I'd start off with the ones estimated by the analyses of the real data.
  
  #Load median parameters variable 
  #medParams <- readRDS("erbrMedParams_noYRE_20240803")
  
  #Make vector of param values for growth
  params.grwth <- c(medParams$grwth_intercept, medParams$grwth_RosCoef, medParams$grwth_TempFallCoef)
  #Make vector of param values for variance in growth
  params.grwthVar <- c(medParams$grwthvar_intercept, medParams$grwthvar_RosCoef)
  #Make vector of param values for survival
  params.surv <- c(medParams$surv_intercept, medParams$surv_RosCoef, medParams$surv_TempFallCoef)
  #Make vector of param values for probability of reproduction
  params.reproYesNo <- c(medParams$reproyesno_intercept, medParams$reproyesno_RosCoef, medParams$reproyesno_TempFallCoef)
  #Make vector of param values for amount of reproduction
  params.repro <- c(medParams$repro_intercept, medParams$repro_RosCoef, medParams$repro_TempFallCoef)
  #Make vector of param values for seedling survival
  params.survSdlg <- c(medParams$surv_intercept, medParams$surv_RosCoef, medParams$surv_TempFallCoef)
  #Make vector of param values for seedling amount
  params.numSdlg <- medParams$newplt_intercept  

  #Set dispersion param variables from runs of real data
  JAGSmodSumm_realDat   #Median values from JAGS run of real dat
  r.inf <- JAGSmodSumm_realDat[1,2]      
  r.sdlg <- JAGSmodSumm_realDat[2,2]
  

  ## Set number of desired years for particular simulation, e.g. 10, 20, 50
  #n.yrs <- 21 #Desired number of years of data plus 1 
  
  
  #Matrices to hold sz & repro where rows are yrs and columns are plants
  #We will be getting and saving each plt's own fate in each yr w the vital rate functions
  mx.sz <- as.data.frame(matrix(NA, nrow=length(1:(n.yrs-1)), ncol=length(1:n.startPlts)))
  mx.reproInf <- as.data.frame(matrix(NA, nrow=length(1:(n.yrs-1)), ncol=length(1:n.startPlts)))
  mx.reproSdlg <- as.data.frame(matrix(NA, nrow=length(1:(n.yrs-1)), ncol=length(1:n.startPlts))) 
  #mx.reproSdlg <- as.data.frame(matrix(NA, nrow=length(1:(n.yrs)), ncol=length(1:n.startPlts))) 
  
  

  #for 1 to the num of starting plants 
  for (pp in 1:n.startPlts) {  #Loop over starting plants
  
    sel.plt <- sz.startPlts[pp]
    
  
      for (yy in 2:(n.yrs-1)) {  #Loop over years
        ## ** Consider making the time loop start in yr 2, so that yr 1 is the starting values 
        ## ** This will make it easier to do this below, although you can also just loop over parent plts & yrs, as I think I am doing **
      
        #In the loop across years for one plant you do this:
        #In all cases, you'll want to make this a monte carlo: for e.g., you pick a random chance of surv from the prob of surv, choose 1 new sz from the distribution of possible new szs, etc. 
        #And, store the number of new seedlings produced that are predicted to be seen in the new year.
        #In this way, yr after yr, you get the data (sz, repro and surv) of each plt. So, populate a matrix w rows that are yrs, columns for each plt, that are the sz & repro & also 0 if the plt is dead.
        
        #In the first year, figure out probs of survival, growth, and reproduction for that year, given the fns for these vital rates, and the climate
        #Based on sz & that yr's clim variables (which you choose already) get all the predicted vital rates using these fns (e.g., pred.grwth)
        
        #Then, use the random variate functions to get the realized values for that plt in that yr. 
        #For e.g., you would get pred.surv and then use rbinom to get whether it lived (1) or died (0), based on the yr's clim & the plts sz. 
        #Similarly you would use pred.grwth and pred.grwthVar with rnorm to get a predicted size for the plant. This is the 'monte carlo' part.
        #Then, use these actual values to get the data for end of the yr: if surv=0, it is dead, but if not, it has the sz from rnorm. 
        #And reproduction is determined by the vital rates governing repro.
    
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
      
                
                #Both the number of inflors and number of seedlings are neg binomials. so, for these, you want to use the predicted mean numbers 
                #(e.g., the number of inflors given a plants size, climate, etc: repro_amount) and the dispersion parameter (e.g., r.inflors) 
                #to get the two parameters for a negbinomial and then use rnegbin to get a single value.
                
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




  ## 4. After doing this, go back to the number of new plants in each year: for these, do the same approach as with the starting plts, but starting in the yr they are 'born'
  ## and then simulate them going forward. Use NA for the sz in the years before they are 'born'
  
  # Subset repro matrix to only contain cols (plts) with non-zero number of seedlings. Then loop thru that.
  #mx.sdlgYes <- mx.reproSdlg[1:(n.yrs-1),] %>% select_if(funs(sum(., na.rm=TRUE) > 0)) #Exclude final yr of data to match sz matrix
  mx.sdlgYes <- mx.reproSdlg %>% select_if(~ sum(., na.rm=TRUE) > 0)  #Try this if error re. funs above
  #mx.sdlgYes <- mx.sdlgYes[1:(n.yrs-1),] #Remove final yr of data to match sz matrix
  
  sdlg <- 1            #set size to be 1 ros for seedlings
  


    for (rr in 1:ncol(mx.sdlgYes)) { #Loop over parent plants
    
        realzd.grwth <- 1
      
        sdlg.strtYrTemp <- which(mx.sdlgYes[,rr]>0)        #Find years when new seedlings appeared 
        sdlg.perYr <- mx.sdlgYes[,rr][mx.sdlgYes[,rr]>0 & !is.na(mx.sdlgYes[,rr])]  #Find new sdlgs per year, if non-zero
        sdlg.strtYr <- rep(sdlg.strtYrTemp, sdlg.perYr)     #List starting years for all new sdlgs, including if more than 1 in given yr
        colCount <- ncol(mx.sz)                             #Re-set column counter each loop based on new cols added in year loop
    
        
        ##Loop over number of new seedlings that a given parent plt had in 1 or more years
        for (ss in 1:length(sdlg.strtYr)) { 
          
          startYr <- sdlg.strtYr[ss]
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
  
  

  ## save as csv 
  date <- Sys.Date()                             #Enter date to be added to file name
  date <- str_replace_all(date, "-", "")
  name <- as.character("_erbr_SimDat20yrNoMiss.")           #Enter name of file, e.g. Tagclust, 4to13, simulated data...

  write.csv(datComb, file=paste(date, name, dd, ".csv", sep=""), row.names=FALSE)
  #print(paste(date, name, dd, ".csv", sep=""))
  ## -----------------------------------------------------------------------------------------------------
  
  
  
  ## MODIFY OUTPUT TO CONTAIN MISSING YEARS OF DATA ------------------------------------------------------ 
  ## Get the records for each individual, but then if a year is one of the missing data years, set that years data to NA 
  datComb1 <- datComb
  
  ## Assign years to be missing
  startConsecYrs <- 10
  yrs.missing <- as.integer(seq((startConsecYrs+1),(n.yrs-1), by=2))
  
  ## Change missing years to NA
  datComb1$RosNew[datComb1$Year %in% (yrs.missing)] <- NA 
  datComb1$InflNew[datComb1$Year %in% (yrs.missing)] <- NA 
  ## -----------------------------------------------------------------------------------------------------
  
  
  
  
 # rm(list=ls())
  
  ## LOAD DATA --------------------------------------------------------------------------------------
  #clim32yr <- read.csv("erbr_climData3seas32yr_221114.csv", header=TRUE)
  
  ## LOAD NO MISSING DATA TO PROCESS WITH THE FOLLOWING 
  #datComb1 <- read.csv("20240906_erbr_SimDat20yrNoMiss.1.csv", header = TRUE)
  #datComb1 <- read.csv("20240906_erbr_SimDat20yrNoMiss.2.csv", header = TRUE)
  #datComb1 <- read.csv("20240906_erbr_SimDat20yrNoMiss.3.csv", header = TRUE)
  #datComb1 <- read.csv("20240906_erbr_SimDat20yrNoMiss.4.csv", header = TRUE)
  #datComb1 <- read.csv("20240906_erbr_SimDat20yrNoMiss.5.csv", header = TRUE)
  #datComb1 <- read.csv("20240906_erbr_SimDat20yrNoMiss.6.csv", header = TRUE)
  #datComb1 <- read.csv("20240906_erbr_SimDat20yrNoMiss.7.csv", header = TRUE)
  #datComb1 <- read.csv("20240906_erbr_SimDat20yrNoMiss.8.csv", header = TRUE)
  #datComb1 <- read.csv("20240906_erbr_SimDat20yrNoMiss.9.csv", header = TRUE)
  #datComb1 <- read.csv("20240906_erbr_SimDat20yrNoMiss.10.csv", header = TRUE)
  
  ## ADD erbr_1ReformatSIMdata_forJAGS SCRIPT -------------------------------------------------------
  
  
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
    left_join(clim32yr, by = c("ClimYr" = "Year"))
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
  name <- as.character("SimDat20yrMiss.")        #Enter name of file, e.g. Tagclust, 4to13, simulated data...
  
  write.csv(erbr.1, file=paste(date, "_erbr_", name, dd, ".Format4JAGS", ".csv", sep=""), row.names=FALSE)
  ## -----------------------------------------------------------------------------------
  
  
  
}  ##End dataset loop



## FOR NO MISSING DATASETS, DONE MANUALLY 
#write.csv(erbr.1, file=paste(date, "_erbr_", "SimDat20yrNoMiss.10", ".Format4JAGS", ".csv", sep=""), row.names=FALSE)






  