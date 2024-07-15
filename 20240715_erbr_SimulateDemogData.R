## April Goebl & Dan Doak
## Script started 2024-07-15
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
#mx.list <- readRDS("erbrIPMmxList_noYRE_P1k_20230503")  #Output generated from erbr_4IPMdeterm script 
## ------------------------------------------------------------------------------------------------



## 1. First, for a set number of years (lets say 30), simulate climate variables for each year. 
## What I would do is use the real data from the study to choose sets of annual data values for the set of climate variables.
num.yrs <- 30 #Assign number of years to simulate climate data for 

#Create empty variable to hold simulated climate data
column.names <- colnames(clim32yr)
sim.clim <- as.data.frame(matrix(NA, nrow=length(1:num.yrs), ncol=length(column.names)))
colnames(sim.clim) <- column.names
sim.clim$Year <- 1:num.yrs

#List of random numbers that corresponds to a set of climate values 
rel.yrs <- 2002:2021 #Select relevant subset of yrs (this matches what was used for stoch lam; will keep Methods consistent) 
climYrs.rel <- clim32yr[clim32yr$Year>=2002 & clim32yr$Year<=2021,]
randVals.yr <- sample(1:nrow(climYrs.rel), size=num.yrs, prob=NULL, replace=TRUE)

for (cc in 1:length(1:num.yrs)) {
  sim.clim[cc,2:7] <- climYrs.rel[randVals.yr[cc],2:7] 
}




## 2. Then, use the SSD from a mean matrix to pick a set of plant starting sizes: 
## you can actually just use these as a set of multinomial probabilities for the starting size of each plant that is then simulated.

## Obtain stable stage structure; modified from some of Dan's code 
N.startNum <- 181                  #Set starting pop sz as 2007 obs sz?

mx.num <- length(mx.list)
#mx.mean <- Reduce('+', mx.list)/mx.num  #Calculate mean of all matrices 

#date <- Sys.Date()        #Get date to be added to file name
#date <- str_replace_all(date, "-", "") 
#saveRDS(mx.mean, file=paste("C:/Users/april/Dropbox/CU_Boulder_PhD/DBG_Internship/erbrMeanMatrix_noYRE_P1k", date, sep="_"))

#Load mean matrix variable 
mx.mean <- readRDS("erbrMeanMatrix_noYRE_P1k_20240715")
#Does second row of matrix look funny (particularly column 1)? 

mnees=eigenall(mx.mean)
N.vecStart=N.startNum * mnees$stable.stage




## 3. for 1 to the number of starting plants (i would make this considerably larger than you think you want to have), you then do the following for each year of time:
  ## a. in the first year, figure out probs of survival, growth, and reproduction for that year, given the fns for these vital rates, and the climate. 
     ## In all cases, you'll want to make this a monte carlo: for example, you pick a random chance of survival from the prob of survival, choose one new size from the distribution of possible new sizes, etc. and, store the number of new seedlings produced that are predicted to be seen in the new year.

  ## b. in this way, year after year, you get the data (size, repro and survival) of each plant. So, I'd populate a matrix with rows that are years, columns for each plant, that are the size and repro and also zero if the plant is dead.

#Make matrix to hold size where rows are yrs and columns are plants
#Make matrix to hold repro where rows are yrs and columns are plants
#add zero if the plant is dead

num.startPlts <- 10 #Set much larger later (e.g. 500?)
for (pp in 1:num.startPlts) {
  
  for (yy in 1:num.yrs) {
    
  }
}

## DAN: I am stuck here. I am not sure how I figure out probs of surv, growth, repro given the fns for these vital rates and climate (where do the parameter values come from?).
## I don't quite understand: "make this a monte carlo: you pick a random chance of survival from the prob of survival, choose one new size from the distribution of possible new sizes, etc."
## I pasted some code here from our SSDM deterministic estimates. I feel like some of this code might be relevant, but not sure. Can you point me in the right direction? 

##### ***PASTED CODE FROM ONE SSDM SCRIPT********** ##################################################################
## Loop over different numbers of bins
for (bb in 1:length(bin.num)) {
  
  ## Improved method of finding median size/ bin mids (code from Dan)
  vec.bin = c(minsize, minsize+1:bin.num[bb]*(maxsize-minsize)*(1/bin.num[bb])) 
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
  truebinsizes = n.bin  
  
## Initialize variable with input data for predictions
in.data <- as.data.frame(binmids) 
colnames(in.data) <- "RosNew"


## Plug in median param vals, mean climate and selected size predictor values for sz classes/ loop into model formulas
## Exclude random transect effects here 

## Growth (neg binom)
pred.grwth <- exp(medParams$grwth_intercept + medParams$grwth_RosCoef*log(in.data$RosNew) 
                  + medParams$grwth_TempFallCoef*climMeans["Mean_fall_temp"]
                  + medParams$grwth_TempSummerCoef*climMeans["Mean_summer_temp"]
                  + medParams$grwth_TempWinterCoef*climMeans["Mean_winter_temp"]
                  + medParams$grwth_PptFallCoef*climMeans["Tot_fall_ppt"]
                  + medParams$grwth_PptSummerCoef*climMeans["Tot_summer_ppt"]
                  + medParams$grwth_PptWinterCoef*climMeans["Tot_winter_ppt"])

## Variance in growth (neg binom)
pred.grwthVar <- exp(medParams$grwthvar_intercept + medParams$grwthvar_RosCoef*log(in.data$RosNew)) 

## Survival (binom)  
pred.surv <- 1/(1+exp(-(medParams$surv_intercept + medParams$surv_RosCoef*in.data$RosNew + 
                          medParams$surv_PptWinterCoef*climMeans["Tot_winter_ppt"] + 
                          medParams$surv_TempWinterCoef*climMeans["Mean_winter_temp"] +
                          medParams$surv_TempSummerCoef*climMeans["Mean_summer_temp"] + 
                          medParams$surv_TempFallCoef*climMeans["Mean_fall_temp"])))

## Probability of reproducing (binom)  
pred.reproYesNo <- 1/(1+exp(-(medParams$reproyesno_intercept + medParams$reproyesno_RosCoef*log(in.data$RosNew) +
                                medParams$reproyesno_PptFallCoef*climMeans["Tot_fall_ppt"] +
                                medParams$reproyesno_PptSummerCoef*climMeans["Tot_summer_ppt"] +
                                medParams$reproyesno_TempFallCoef*climMeans["Mean_fall_temp"] +
                                medParams$reproyesno_TempSummerCoef*climMeans["Mean_summer_temp"] +
                                medParams$reproyesno_TempWinterCoef*climMeans["Mean_winter_temp"]))) 

## Reproduction (neg binom)
pred.repro <- exp(medParams$repro_intercept + medParams$repro_RosCoef*log(in.data$RosNew) + 
                    medParams$repro_PptFallCoef*climMeans["Tot_fall_ppt"] +
                    medParams$repro_PptSummerCoef*climMeans["Tot_summer_ppt"] +
                    medParams$repro_TempSummerCoef*climMeans["Mean_summer_temp"] +
                    medParams$repro_TempWinterCoef*climMeans["Mean_winter_temp"] +
                    medParams$repro_TempFallCoef*climMeans["Mean_fall_temp"])

## Seedling survival (binom)  
in.dataSdlg <- in.data[1,]  #Subset to keep only 1 row of in.data 
in.dataSdlg <- 1            #Change ros size to be 1 (sz of sdlg)
pred.survSdlg <- 1/(1+exp(-(medParams$surv_intercept + medParams$surv_RosCoef*log(in.dataSdlg) + 
                              medParams$surv_PptWinterCoef*climMeans["Tot_winter_ppt"] + 
                              medParams$surv_TempWinterCoef*climMeans["Mean_winter_temp"] +
                              medParams$surv_TempSummerCoef*climMeans["Mean_summer_temp"] + 
                              medParams$surv_TempFallCoef*climMeans["Mean_fall_temp"])))

## Seedlings per inflor (neg binom)
pred.numSdlg <- exp(medParams$newplt_intercept + log(pred.repro))


## From Dan's code
## Constructing matrix models
grwth.mx <- matrix(data=NA, nrow=n.bin, ncol=n.bin)
## Growth probs using cdf fn
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
mx1 <- surv.grwth.mx #Growth and survival, without repro
## Add reproduction and recruitment
mx <- matrix(0, (n.bin+1), (n.bin+1))
mx[2:(n.bin+1), 2:(n.bin+1)] <- mx1
mx[2,1] <- pred.survSdlg                                  #First column (seedling surv in element 2,1)
mx[1,2:(n.bin+1)] <- pred.numSdlg                         #First row (new seedlings)

lam.out.template[bb,2] <- Re(eigen(mx)$values[1])         #Calculate & store lambda

}  #End bb loop
###### ********END OF PASTED CODE******** ##################################################################################






## 4. After doing this, you would go back to the number of new plants in each year: for these, do the same approach as with the starting plants, but starting in the year they are born, and then simulate them going forward. I would use NA for the size in the years before they are 'born'.