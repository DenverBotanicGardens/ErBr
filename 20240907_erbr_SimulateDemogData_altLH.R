## April Goebl & Dan Doak
## Script modified 2024-09-04
## Collaboration with Denver Botanic Gardens on Eriogonum brandegeii modeling 
## Script written in response to reviews on Goebl et al 2024 Biol Conservation 
## Generate simulated data on lives of individual plants with known demographic parameters 
## Modify 'life history' to get additional simulated datasets



rm(list=ls())



## SET WD -----------------------------------------------------------------------------------------
#setwd("C:/Users/Dan Doak/Desktop/Students/April/eriogonum models/manu fall2024")
## ------------------------------------------------------------------------------------------------
setwd("C:/Users/Dan Doak/Desktop/Students/April/eriogonum models/manu fall2024/life histories")


## LOAD DATA --------------------------------------------------------------------------------------
clim32yr <- read.csv("erbr_climData3seas32yr_221114.csv", header=TRUE)
erbr <- read.csv("erbr_TagClust2022_20230408.csv", header=TRUE)
erbr$Year <- as.factor(erbr$Year)
#JAGSmodSumm_realDat <- readRDS("erbr_JAGSmodBestSUMM_c3t10s10b5_noYRE_20240829.rds")
mx.mean <- readRDS("erbrMeanMatrix_noYRE_P1k_20240715") #Load mean matrix variable
medParams <- readRDS("erbrMedParams_noYRE_20240803")    #Load median parameters variable 
## ------------------------------------------------------------------------------------------------



## LOAD PACKAGES ----------------------------------------------------------------------------------
library(dplyr)
library(matrixStats)
library(stringr)
## ------------------------------------------------------------------------------------------------




## 1. First, fix each of the climate drivers at their mean values
## Use the real data from the study 
rel.yrs <- 2002:2021 #Select relevant subset of yrs (this matches what was used for stoch lam to keep Methods consistent) 
climYrs.rel <- clim32yr[clim32yr$Year>=2002 & clim32yr$Year<=2021,]
mn.clim <- as.matrix(colMeans(climYrs.rel[,-1]))
mn.clim <- as.data.frame(t(mn.clim))
## ------------------------------------------------------------------------------------------------




## 2. Then, use the SSD from a mean matrix to pick a set of plant starting sizes: 
## you can actually just use these as a set of multinomial probabilities for the starting size of each plt that is then simulated

## Obtain stable stage structure; modified from some of Dan's code 
N.startNum <- 181                 #Set starting pop sz as 2007 obs sz

SSD=eigen(mx.mean)$vectors[,1]
SSD=Re(SSD/sum(SSD))
N.vecStart=N.startNum * SSD


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

#binmids <- c(1, binmids)  
n.bin = length(binmids)

  
## Select starting sizes of plants for simulated data using SSD and median sz classes
#num.startPlts <- 200 #Num of starting plts. Set much larger later (e.g. 500?) (make this considerably larger than you think you want to have)
#N.startProbs <- (N.vecStart*1) / N.startNum
#sum(N.startProbs) #Should equal 1
#sz.startPlts <- sample(x=binmids, size=num.startPlts, replace=TRUE, prob=N.startProbs)  
#sz.startPlts <- round(sz.startPlts, digits=0)   #Round to nearest integer 
## ------------------------------------------------------------------------------------------------




## 3. 

#At the start of this code, you should make vectors of the parameter values for each of the vital rate functions 
#these are the 'rules' that will determine each vital rate's predicted value for each plant in each year. I'd start off with the ones estimated by the analyses of the real data.

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
#JAGSmodSumm_realDat
#r.inf <- JAGSmodSumm_realDat[1,2]      
#r.sdlg <- JAGSmodSumm_realDat[2,2]


num.yrs <- 21 #Assign number of years (plus 1) to simulate data for 

#Matrices to hold sz & repro where rows are yrs and columns are plants
#mx.sz <- as.data.frame(matrix(NA, nrow=length(1:(num.yrs-1)), ncol=length(1:num.startPlts)))
#mx.reproInf <- as.data.frame(matrix(NA, nrow=length(1:(num.yrs-1)), ncol=length(1:num.startPlts)))
#mx.reproSdlg <- as.data.frame(matrix(NA, nrow=length(1:(num.yrs-1)), ncol=length(1:num.startPlts))) 


#Make variable to hold matrices?

#mx.out.list <- NULL        #List variable to store all matrices


# ####################################################################
# parameters for a fast growth, med survival, higher repro LH:
medParams$grwth_RosCoef = medParams$grwth_RosCoef*1.1
medParams$grwth_intercept = 1.1

medParams$surv_intercept=medParams$surv_intercept*1.2
medParams$surv_RosCoef=medParams$surv_RosCoef*0.8

medParams$repro_intercept=medParams$repro_intercept+2
medParams$reproyesno_intercept = medParams$reproyesno_intercept+2
medParams$reproyesno_RosCoef=medParams$reproyesno_RosCoef*1.2
# #######################################################

####################################################################
# # parameters for a med growth, lower survival LH, high repro:
# medParams$grwth_RosCoef = medParams$grwth_RosCoef*1.025
# medParams$grwth_intercept = 1.025
# 
# medParams$surv_intercept=medParams$surv_intercept*1.3
# medParams$surv_RosCoef=medParams$surv_RosCoef*0.5
# 
# medParams$repro_intercept=medParams$repro_intercept+5.5
# 
# medParams$reproyesno_intercept = medParams$reproyesno_intercept+2
# medParams$reproyesno_RosCoef=medParams$reproyesno_RosCoef*1.2
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
          
          
          
          
          ## From Dan's code
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

## -------------------------------------------------------------------------------------------------
plot(binmids,pred.grwth)
          abline(a=0,b=1)
          plot(binmids,pred.grwth, log='xy')
          abline(a=0,b=1)
plot(binmids,pred.surv)      
plot(binmids,pred.repro)
plot(binmids,pred.reproYesNo)

print(lam)




