## April Goebl & Dan Doak
## Script modified 2024-07-31
## Collaboration with Denver Botanic Gardens on Eriogonum brandegeii modeling 
## Script written in response to reviews on Goebl et al 2024 Biol Conservation 
## Generate simulated data on lives of individual plants with known demographic parameters 



rm(list=ls())



## SET WD -----------------------------------------------------------------------------------------
setwd("C:/Users/april/Dropbox/CU_Boulder_PhD/DBG_Internship")
## ------------------------------------------------------------------------------------------------



## LOAD DATA --------------------------------------------------------------------------------------
source('R_scripts/eigenall.r')
clim32yr <- read.csv("Manuscript/erbr_climData3seas32yr_221114.csv", header=TRUE)
erbr <- read.csv("Manuscript/erbr_TagClust2022_20230408.csv", header=TRUE)
erbr$Year <- as.factor(erbr$Year)
## ------------------------------------------------------------------------------------------------



## LOAD PACKAGES ----------------------------------------------------------------------------------
library(dplyr)
library(matrixStats)
## ------------------------------------------------------------------------------------------------




## 1. First, for a set number of years (lets say 30), simulate climate variables for each year. 
## What I would do is use the real data from the study to choose sets of annual data values for the set of climate variables.
num.yrs <- 51 #Assign number of years (plus 1) to simulate climate data for 

#Create empty variable to hold simulated climate data
column.names <- colnames(clim32yr)
sim.clim <- as.data.frame(matrix(NA, nrow=length(1:num.yrs), ncol=length(column.names)))
colnames(sim.clim) <- column.names
sim.clim$Year <- 1:num.yrs

#List of random numbers that corresponds to a set of climate values 
rel.yrs <- 2002:2021 #Select relevant subset of yrs (this matches what was used for stoch lam to keep Methods consistent) 
climYrs.rel <- clim32yr[clim32yr$Year>=2002 & clim32yr$Year<=2021,]
randVals.yr <- sample(1:nrow(climYrs.rel), size=num.yrs, prob=NULL, replace=TRUE)

for (cc in 1:length(1:num.yrs)) {
  sim.clim[cc,2:7] <- climYrs.rel[randVals.yr[cc],2:7] 
}
## ------------------------------------------------------------------------------------------------




## 2. Then, use the SSD from a mean matrix to pick a set of plant starting sizes: 
## you can actually just use these as a set of multinomial probabilities for the starting size of each plant that is then simulated.

## Obtain stable stage structure; modified from some of Dan's code 
N.startNum <- 181                 #Set starting pop sz as 2007 obs sz

#mx.num <- length(mx.list)
#mx.mean <- Reduce('+', mx.list)/mx.num  #Calculate mean of all matrices 

#date <- Sys.Date()        #Get date to be added to file name
#date <- str_replace_all(date, "-", "") 
#saveRDS(mx.mean, file=paste("C:/Users/april/Dropbox/CU_Boulder_PhD/DBG_Internship/erbrMeanMatrix_noYRE_P1k", date, sep="_"))

#Load mean matrix variable 
mx.mean <- readRDS("erbrMeanMatrix_noYRE_P1k_20240715")
## Check if first sz class is seedlings in the next year or not. **
## Currently model assumes plt is always in the smallest subsequent sz class if you survive. 

mnees=eigenall(mx.mean)
N.vecStart=N.startNum * mnees$stable.stage



# Specify min and max plt sz
minsize <- 1
maxsize <- (max(erbr$RosNew, na.rm=TRUE)) 

## For median sz estimation
## From Dan's bistorts to guppies code 
##new size density estimation for median size estimation
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

binmids <- c(1, binmids)  
  
  
## Select starting sizes of plants for simulated data using SSD and median sz classes
num.startPlts <- 100 #Num of starting plts. Set much larger later (e.g. 500?) (make this considerably larger than you think you want to have)
N.startProbs <- (N.vecStart*1) / N.startNum
sum(N.startProbs) #Should equal 1
sz.startPlts <- sample(x=binmids, size=num.startPlts, replace=TRUE, prob=N.startProbs)  
sz.startPlts <- round(sz.startPlts, digits=0)   #Round to nearest integer 
## ------------------------------------------------------------------------------------------------




## 3. 

#At the start of this code, you should make vectors of the parameter values for each of the vital rate functions 
#these are the 'rules' that will determine each vital rate's predicted value for each plant in each year. I'd start off with the ones estimated by the analyses of the real data.

## Calculate median param values from JAGS output from real data
#medParams <- as.data.frame(colMedians(as.matrix(chains)))
#medParams <-as.data.frame(t(as.data.frame(medParams)))
#colnames(medParams) <- colnames(chains)

#date <- Sys.Date()        #Get date to be added to file name
#date <- str_replace_all(date, "-", "") 
#saveRDS(medParams, file=paste("C:/Users/april/Dropbox/CU_Boulder_PhD/DBG_Internship/erbrMedParams_noYRE", date, sep="_"))

#Load median parameters variable 
medParams <- readRDS("erbrMedParams_noYRE_20240803")

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



## Set number of desired years for particular simulation, e.g. 10, 20, 50
num.yrs <- 21 #Desired number of years of data plus 1 


#Matrices to hold sz & repro where rows are yrs and columns are plants
mx.sz <- as.data.frame(matrix(NA, nrow=length(1:(num.yrs-1)), ncol=length(1:num.startPlts)))
mx.reproInf <- as.data.frame(matrix(NA, nrow=length(1:(num.yrs-1)), ncol=length(1:num.startPlts)))
mx.reproSdlg <- as.data.frame(matrix(NA, nrow=length(1:(num.yrs)), ncol=length(1:num.startPlts))) 

#One thing to note: you are never making a matrix for the vital rates, you are just getting each plt's own fate in each yr w the vital rate functions.



#for 1 to the num of starting plants 
for (pp in 1:num.startPlts) {  #Loop over starting plants

  sel.plt <- sz.startPlts[pp]
  
  
  for (yy in 1:(num.yrs-1)) {  #Loop over years
    ## ** AG: Consider making the time loop start in yr 2, so that yr 1 is the starting values 
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
    pred.surv <- 1/(1+exp(-(medParams$surv_intercept + medParams$surv_RosCoef*sel.plt + 
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
        mx.reproSdlg[yy+1,pp] <- 0     ##AG ** Does this need 2nd yr following death assigned to zero? ** 
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
              
              #* This should be a random from a neg binomial, as that is how it is fit. **** 
              #* So, you also need to pull out the fitted r.infls variable and fit this way
              #*erbr_JAGSmodbest_noYRE lines 92-99: Only loop over rows with >0 inflors, to compare estimated inflors to observed
              #for(i in 1:Nrows.w.inflors) {
                # p.infls[i] <- r.infls/(r.infls + repro_amount[rows.w.inflors[i]])
                #InflNew[rows.w.inflors[i]] ~ dnegbin(p.infls[i], r.infls) }
              #From rnbinom manual: size=dispersionParam, prob=size/(size+mu)
              #Then dispersionParam = (prob*mu) / (1-prob)
              #**Does this need to change to only include plts that are reproducing?
              r.inf <- (pred.reproYesNo*pred.repro) / (1 - pred.reproYesNo) #** Is prob=pred.reproYesNo correct here?
              realzd.repro <- rnbinom(n=1, size=r.inf, mu=pred.repro)       #** Is mu=pred.repro correct here?

              #Enter inf data into repro matrix
              mx.reproInf[yy,pp] <- realzd.repro
              
              # Seedlings per inflor (negative binomial)
              pred.numSdlg <- exp(medParams$newplt_intercept + log(pred.repro))  
              r.sdlg <- (pred.reproYesNo*pred.numSdlg) / (1-pred.reproYesNo)     
              realzd.numSdlg <- rnbinom(n=1, size=r.sdlg, mu=pred.numSdlg)       
              
              ##*What I would do here is take the product of the realized (result of a random draw) values for 
              ##did.the.plant.reproduce*number.inflors.if.reproducing*new.plts.per.inflor. I think that I have all the steps laid out correctly to get the realized # new plants. 
              ##so there is no need here for a new and different random draw, just the product of these other values. 
              ##*And this should be an integer already.
              realzd.numSdlg <- realzd.reproYesNo * realzd.repro * realzd.numSdlg  
              
              #Enter seedling data into seedling matrix
              mx.reproSdlg[yy+1,pp] <- realzd.numSdlg  } 
          ## ** Alternatively, consider adding up the number of seedlings each year, so that you just automatically get a total of new ones for each year. 
          ## or leave it and then add up so that you get a single vector of totals for each year, as I think it's doing below.**
          ## --
          
          ## ** Now the number of sdlgs, when not 0, seems too high *** 
            
            #Then, if the plt is still alive, you go to the next year and do the same.
          
  } ##End year loop
  
}   ##End starting number of plants loop

## ** Note that the starting sizes are not in the size output matrix. Does that seem right? 
## so year 1 is the 'realized' size and inflor number based on the starting size and climate. ** 
## New seedlings can appear in year 2 or later based on the inflor numbers from year 1. Does that seem correct? ** 

## ** In size matrix, once dead, size is 0 for year died and then NAs afterwards **
## ** In repro (inflors) matrix, once dead, repro is 0 for yr died and NAs afterwards
## ** 0 in this repro matrix can also mean there was no repro in that yr, in which case, non-NA and non-zero numbers follow **
## ** In repro (seedling) matrix, once dead, seedlings are 0 for yr died and NAs afterwards
## ** 0 in this seedling matrix can also mean no seedlings in that yr, in which case, non-NA values (0s or integers) follow ** 

## -------------------------------------------------------------------------------------------------






## 4. After doing this, you would go back to the number of new plants in each year: for these, do the same approach as with the starting plants, but starting in the year they are born
## and then simulate them going forward. I would use NA for the size in the years before they are 'born'.

## AG: Check if repro for a given year means num of seedlings in subsequent year or not, and if this is what it should be... **

# Subset repro matrix to only contain cols (plts) with non-zero number of seedlings. Then loop thru that.
#mx.sdlgYes <- mx.reproSdlg[1:(num.yrs-1),] %>% select_if(funs(sum(., na.rm=TRUE) > 0)) #Exclude final yr of data to match sz matrix
mx.sdlgYes <- mx.reproSdlg %>% select_if(~ sum(., na.rm=TRUE) > 0)  #Try this if error re. funs above
#mx.sdlgYes <- mx.sdlgYes[1:(num.yrs-1),] #Remove final year of data to match size matrix

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
      
    
      for (yy in startYr:(num.yrs-2)) {     #For each new seedling, loop over years
        
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
            mx.sz[yy+1,colCount+ss] <- 0   ## ** do we need two yrs of zero here following death? For missing data code below **
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
          if (realzd.grwth <= 0) {     #If size is equal to or less than 0, change to 1. Is there a better way to do this? **
            realzd.grwth <- 1
          }
        
          ##Enter realized size into size matrix
          mx.sz[yy+1,colCount+ss] <- realzd.grwth 
          ## --
          
      }  ## End year loop
    }    ## End new seedlings loop
}        ## End parent plant loop


## ** Reproduction of new seedlings/ plants is not included here yet ** 
## ** I reckon we can say that we followed seedlings resulting from original plants but not subsequent seedlings that may have been born **  





## Modify output to contain missing years of data 
## You'll get the records for each individual, but then say that if a year is one of the missing data years, just set that years data to NA. 
## Add a year column
mx.szWyr <- cbind.data.frame(1:(num.yrs-1), mx.sz[1:(num.yrs-1),])
colnames(mx.szWyr) <- c("Year", colnames(mx.sz))

mx.reproWyr <- cbind.data.frame(1:(num.yrs-1), mx.reproInf)
colnames(mx.reproWyr) <- c("Year", colnames(mx.reproInf))

## Assign years to be missing
startConsecYrs <- 10
yrs.missing <- as.integer(seq((startConsecYrs+1),(num.yrs-1), by=2))

## Change missing years to NA
mx.szWyr[mx.szWyr$Year %in% (yrs.missing),] <- NA 
mx.reproWyr[mx.reproWyr$Year %in% (yrs.missing),] <- NA 

## Add columns to repro matrix so number matches that of size mx. Values can be 0 for no repro ** change above in sdlg loop later **
num.addnCols <- ncol(mx.szWyr) - ncol(mx.reproWyr)
#mx.reproWyrAdCol 
addnCols <- as.data.frame(matrix(0, nrow=nrow(mx.reproWyr), ncol=num.addnCols))
mx.reproWyrAdCol <- cbind(mx.reproWyr, addnCols)

## Remove NA values from year column
mx.szWyr$Year <- 1:(num.yrs-1)
mx.reproWyrAdCol$Year <- 1:(num.yrs-1)


## In writing this, I am also realizing that there is one complication: if a plant dies one year, but it was not censused the next year, 
## then it should be censused the following year: in other words, usually a plant that is dead one year would have no records for subsequent years, 
## but in the case of one that died the missing year, there should be a census the next year. 
## The easiest way to deal with this might be to set two years of data to zero following death, for a plant that dies, 
## but then to say that if, after adding in the Nas for missing data, any record that is after a year recording zero is removed.
## -----------------------------------------------------------------------------------------------------




## Modify output to match format of raw data so can be appropriately modified for use with JAGs code **
## 'Site' column should be optional, but if adding, could split individuals between site 'E.' and 'W.'
## 'Year' column could be e.g. 1-50 or change to e.g. 2004-2022
## 'Transect' column could be entered as 'TransectNew' where individuals are assigned to E.1-E.7 and W.1-W.5
## 'TagNew' column should be added where indivs have a unique tag ID that is E.### or W.###
## Should be able to omit the following columns 'X', 'Y', 'DiameterX', 'DiameterY', 'Rust', 'BrType', 'InflBr', 'Comments'
## Therefore columns to have: 'Year', 'TransectNew', 'TagNew', 'RosNew', 'InflNew'
## Change output from above so there are rows for all individuals across all year (rows are years and plants)


## Combine all plants into 1 column, where all years for a given plt appear sequentially
datComb <- NULL

for (ll in 2:ncol(mx.szWyr)) {   
  temp <- cbind(mx.szWyr$Year, mx.szWyr[,ll], mx.reproWyrAdCol[,ll])
  datComb <- rbind(datComb, temp)
}

datComb <- as.data.frame(datComb)
colnames(datComb) <- c("Year", "RosNew", "InflNew") #Names cols 'Year', 'RosNew', 'InflNew'


## Assign plants to one of the following transects (random but equal distribution)
tran <- c("E.1","E.2","E.3","E.4","E.5","E.6","E.7","W.1","W.2","W.3","W.4","W.5")
num.plt <- ncol(mx.sz)
tran.rep <- rep(tran, (ceiling(num.plt/length(tran))))
tran.rep <- tran.rep[1:num.plt]
tran.rand <- sample(tran.rep, replace=FALSE)
datComb$TransectNew <- NULL
datComb$TransectNew <- rep(tran.rand, each=(num.yrs-1), replace=FALSE)


## Assign unique tag ID to each plant 
tag <- c(1:num.plt)
tag.rep <- rep(tag, each=(num.yrs-1), replace=FALSE)
datComb$TagNew <- paste(datComb$TransectNew, tag.rep, sep='.') #Tag new included site and transect; same format as real data


## ** save as csv ** 


## NOTE- from Dan: "What we want to do is make a data file of the simulated data that is exactly the same format as the real data, 
## so it can run through the same jags program. Then, we do this multiple times 
## and compare the distribution of fitted parameters to the ones that were used to create the simulated data." ** 

## **TO DO: Modify erbr_1ReformatData_forJAGS ** 
## Edit text to include methods, results, discussion on simulation and add figure(s) **

  