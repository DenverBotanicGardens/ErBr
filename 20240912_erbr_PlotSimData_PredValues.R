## April Goebl & Dan Doak
## Script modified 2024-09-12
## Collaboration with CU Boulder and Denver Botanic Gardens on Eriogonum brandegeii modeling 
## Plot predicted size and growth rate results from modeling of simulated demographic data 
## to show JAGS model performance with missing data 


rm(list=ls())
graphics.off()


library(MASS)
library(dplyr)
library(matrixStats)
library(corrplot)
library(rgl)
library(viridis) 
library(car)
library(coda)
library(lme4)
library(resample)
library(gplots)
library(stringr)
library(plotrix)



## ASSIGN NAME VARIABLE FOR DESIRED DATASETS 
#date <- as.character("20240911")
#name <- as.character("SimDat20yr")
## ------------------------------------------------------------------------------------------------





## LOAD FILES -------------------------------------------------------------------------------------
clim32yr <- read.csv("erbr_climData3seas32yr_221114.csv", header=TRUE)

date <- "20240912"
name <- "SimDat20yr"
paramMMgrwth <- readRDS(file=paste(date, "_erbr_paramMMgrwthWint_", name, ".rds", sep=""))
paramMMsurv <- readRDS(file=paste(date, "_erbr_paramMMsurvWint_", name, ".rds", sep=""))
#glmm.grwth <- readRDS(file=paste(date, "_erbr_GLMMgrwth_", name, ".rds", sep=""))
#glmm.surv <- readRDS(file=paste(date, "_erbr_GLMMsurv_", name, ".rds", sep=""))


#summ.mod1 <- readRDS("20240906_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat20yr.1.rds")
summ.mod1 <- readRDS("20240906_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat20yr.1.rds")
summ.mod2 <- readRDS("20240906_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat20yr.2.rds")
summ.mod3 <- readRDS("20240906_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat20yr.3.rds")
summ.mod4 <- readRDS("20240906_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat20yr.4.rds")
summ.mod5 <- readRDS("20240906_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat20yr.5.rds")
summ.mod6 <- readRDS("20240907_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat20yr.6.rds")
summ.mod7 <- readRDS("20240907_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat20yr.7.rds")
summ.mod8 <- readRDS("20240907_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat20yr.8.rds")
summ.mod9 <- readRDS("20240907_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat20yr.9.rds")
summ.mod10 <- readRDS("20240907_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat20yr.10.rds")

simDat1 <- read.csv("C:/Users/april/Dropbox/CU_Boulder_PhD/DBG_Internship/Manuscript/MS_DataAndCode_archive/20240907_erbr_SimDat20yrNoMiss.1.Format4JAGS.csv", header=TRUE)
simDat2 <- read.csv("C:/Users/april/Dropbox/CU_Boulder_PhD/DBG_Internship/Manuscript/MS_DataAndCode_archive/20240907_erbr_SimDat20yrNoMiss.2.Format4JAGS.csv", header=TRUE)
simDat3 <- read.csv("C:/Users/april/Dropbox/CU_Boulder_PhD/DBG_Internship/Manuscript/MS_DataAndCode_archive/20240907_erbr_SimDat20yrNoMiss.3.Format4JAGS.csv", header=TRUE)
simDat4 <- read.csv("C:/Users/april/Dropbox/CU_Boulder_PhD/DBG_Internship/Manuscript/MS_DataAndCode_archive/20240907_erbr_SimDat20yrNoMiss.4.Format4JAGS.csv", header=TRUE)
simDat5 <- read.csv("C:/Users/april/Dropbox/CU_Boulder_PhD/DBG_Internship/Manuscript/MS_DataAndCode_archive/20240907_erbr_SimDat20yrNoMiss.5.Format4JAGS.csv", header=TRUE)
simDat6 <- read.csv("C:/Users/april/Dropbox/CU_Boulder_PhD/DBG_Internship/Manuscript/MS_DataAndCode_archive/20240907_erbr_SimDat20yrNoMiss.6.Format4JAGS.csv", header=TRUE)
simDat7 <- read.csv("C:/Users/april/Dropbox/CU_Boulder_PhD/DBG_Internship/Manuscript/MS_DataAndCode_archive/20240907_erbr_SimDat20yrNoMiss.7.Format4JAGS.csv", header=TRUE)
simDat8 <- read.csv("C:/Users/april/Dropbox/CU_Boulder_PhD/DBG_Internship/Manuscript/MS_DataAndCode_archive/20240907_erbr_SimDat20yrNoMiss.8.Format4JAGS.csv", header=TRUE)
simDat9 <- read.csv("C:/Users/april/Dropbox/CU_Boulder_PhD/DBG_Internship/Manuscript/MS_DataAndCode_archive/20240907_erbr_SimDat20yrNoMiss.9.Format4JAGS.csv", header=TRUE)
simDat10 <- read.csv("C:/Users/april/Dropbox/CU_Boulder_PhD/DBG_Internship/Manuscript/MS_DataAndCode_archive/20240907_erbr_SimDat20yrNoMiss.10.Format4JAGS.csv", header=TRUE)
## -----------------------------------------------------------------------------------------------





## GET CLIMATE VARIABLES ------------------------------------------------------------------------
#n.yrs <- 20 #Dataset length 
climYrs1 <- clim32yr[clim32yr$Year %in% unique(simDat1$ClimYr),]
climYrs2 <- clim32yr[clim32yr$Year %in% unique(simDat2$ClimYr),]
climYrs3 <- clim32yr[clim32yr$Year %in% unique(simDat3$ClimYr),]
climYrs4 <- clim32yr[clim32yr$Year %in% unique(simDat4$ClimYr),]
climYrs5 <- clim32yr[clim32yr$Year %in% unique(simDat5$ClimYr),]
climYrs6 <- clim32yr[clim32yr$Year %in% unique(simDat6$ClimYr),]
climYrs7 <- clim32yr[clim32yr$Year %in% unique(simDat7$ClimYr),]
climYrs8 <- clim32yr[clim32yr$Year %in% unique(simDat8$ClimYr),]
climYrs9 <- clim32yr[clim32yr$Year %in% unique(simDat9$ClimYr),]
climYrs10 <- clim32yr[clim32yr$Year %in% unique(simDat10$ClimYr),]
## ---------------------------------------------------------------------------------------------




## GET RANGE IF PLANT SIZES FOR MAKING PREDICTIONS ---------------------------------------------
## Initialize variable with size input data for predictions
#in.data <- as.data.frame(binmids) 
#in.data <- 1:151
#colnames(in.data) <- "RosNew"
#$RosNew)

#popSz.start <- 181                 #Set starting pop sz as 2007 obs sz

#SSD=eigen(mx.mean)$vectors[,1]
#SSD=Re(SSD/sum(SSD))
#N.vecStart=popSz.start * SSD


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
## ------------------------------------------------------------------------------------------




## --------------------------------------------------------------------------------------------
## OBTAIN TRUE VALUES FOR ADDING TO PLOTS
# Load data from previous JAGS runs of real data
medParams.realDat <- readRDS("erbrMedParams_noYRE_20240803")
medParams.realDatTr <- as.data.frame(t(medParams.realDat))
medParams.realDatTr <- as.data.frame(cbind(medParams.realDatTr$`colMedians(as.matrix(chains))`,colnames(medParams.realDat)))
colnames(medParams.realDatTr) <- c("realData", "Name")

#names.param <- rownames(summ.mod1)[26:41]
medParams.realDatTr <- medParams.realDatTr[55:70,]
medParams.realDatTr <- medParams.realDatTr[c(2:8,12:16),] #Remove intercept and GrwthVar  
medParams.realDatTr$ParamTitle <- c("Grwth Size","Grwth Fall Temp","Grwth Summer Temp","Grwth Winter Temp",
                                    "Grwth Fall Precip","Grwth Summer Precip","Grwth Winter Precip","Surv Size",
                                    "Surv Winter Precip","Surv Fall Temp","Surv Summer Temp","Surv Winter Temp")
medParams.realDatGrwth <- medParams.realDatTr[1:7,]
medParams.realDatSurv <- medParams.realDatTr[8:12,]
## ------------------------------------------------------------------------------------------





## Get Median param values from JAGS run ----------------------------------------------------
names.param <- rownames(summ.mod1)[c(26:33,36:41)]
names.paramTitles <- c("Grwth Intercept","Grwth Size","Grwth Fall Temp","Grwth Summer Temp","Grwth Winter Temp",
                       "Grwth Fall Precip","Grwth Summer Precip","Grwth Winter Precip","Surv Intercept","Surv Size",
                       "Surv Winter Precip","Surv Fall Temp","Surv Summer Temp","Surv Winter Temp")

medParams.1 <- summ.mod1[c(26:33,36:41),2]
medParams.2 <- summ.mod2[c(26:33,36:41),2]
medParams.3 <- summ.mod3[c(26:33,36:41),2]
medParams.4 <- summ.mod4[c(26:33,36:41),2]
medParams.5 <- summ.mod5[c(26:33,36:41),2]
medParams.6 <- summ.mod6[c(26:33,36:41),2]
medParams.7 <- summ.mod7[c(26:33,36:41),2]
medParams.8 <- summ.mod8[c(26:33,36:41),2]
medParams.9 <- summ.mod9[c(26:33,36:41),2]
medParams.10 <- summ.mod10[c(26:33,36:41),2]

medComb.grwth <- as.data.frame(cbind(medParams.1[1:8], medParams.2[1:8],medParams.3[1:8],
                                     medParams.4[1:8],medParams.5[1:8],medParams.6[1:8],
                                     medParams.7[1:8],medParams.8[1:8],medParams.9[1:8],medParams.10[1:8]))
## --------------------------------------------------------------------------------------




## Get GLMM param values ----------------------------------------------------------------
#paramsMM.grwthTmp <- as.data.frame(summary(glmm.grwth)$coefficients[1:nrow(summary(glmm.grwth)$coefficients),1])
#colnames(paramsMM.grwthTemp) <- paste("GLMM.",dd,sep="")
#paramsMM.grwth <- as.data.frame(c(paramsMM.grwth, paramsMM.grwthTmp))

#paramsMM.grwth <- as.data.frame(summary(glmm.grwth)$coefficients[1:nrow(summary(glmm.grwth)$coefficients),1])
#colnames(paramsMM.grwth) <- paste("GLMM.",1,sep="")
#paramsMM.grwth$ParamTitle <- c("Grwth Intercept","Grwth Size","Grwth Fall Precip","Grwth Winter Precip","Grwth Summer Precip",
#                               "Grwth Fall Temp","Grwth Winter Temp","Grwth Summer Temp")

#paramsMM.grwthOrd <- paramsMM.grwth[match(names.paramTitles[1:8], paramsMM.grwth$ParamTitle),]
#paramsMM.survOrd <- paramsMM.surv[match(names.paramTitles[8:12], paramsMM.surv$ParamTitle),]
## --------------------------------------





## OBTAIN VITAL RATE PREDICTIONS ----------------------------------------------

## Plug in median param vals, mean climate and selected size predictor values for sz classes/ loop into model formulas
## Exclude random transect effects here (alternatively, could try mean values)
## Growth (neg binom)
#pred.grwth <- exp(medParams$grwth_intercept + medParams$grwth_RosCoef*log(in.data$RosNew) 
#                  + medParams$grwth_TempFallCoef*climMeans["Mean_fall_temp"]
#                  + medParams$grwth_TempSummerCoef*climMeans["Mean_summer_temp"]
#                  + medParams$grwth_TempWinterCoef*climMeans["Mean_winter_temp"]
#                  + medParams$grwth_PptFallCoef*climMeans["Tot_fall_ppt"]
#                  + medParams$grwth_PptSummerCoef*climMeans["Tot_summer_ppt"]
#                  + medParams$grwth_PptWinterCoef*climMeans["Tot_winter_ppt"])

## Rep/ dataset 1
pred.grwth <- exp(medComb.grwth[1,10] + medComb.grwth[2,10]*log(in.data) 
                  + medComb.grwth[3,10]*clim$Mean_fall_temp
                  + medComb.grwth[4,10]*clim$Mean_summer_temp
                  + medComb.grwth[5,10]*clim$Mean_winter_temp
                  + medComb.grwth[6,10]*clim$Tot_fall_ppt
                  + medComb.grwth[7,10]*clim$Tot_summer_ppt
                  + medComb.grwth[8,10]*clim$Tot_winter_ppt)

pred.grwthGLM <- exp(paramsMM.grwthOrd[1,1] + paramsMM.grwthOrd[2,1]*log(in.data) 
                  + paramsMM.grwthOrd[3,1]*clim$Mean_fall_temp
                  + paramsMM.grwthOrd[4,1]*clim$Mean_summer_temp
                  + paramsMM.grwthOrd[5,1]*clim$Mean_winter_temp
                  + paramsMM.grwthOrd[6,1]*clim$Tot_fall_ppt
                  + paramsMM.grwthOrd[7,1]*clim$Tot_summer_ppt
                  + paramsMM.grwthOrd[8,1]*clim$Tot_winter_ppt)


## Survival (binom)  
pred.surv <- 1/(1+exp(-(medParams$surv_intercept + medParams$surv_RosCoef*in.data$RosNew + 
                          medParams$surv_PptWinterCoef*clim["Tot_winter_ppt"] + 
                          medParams$surv_TempWinterCoef*clim["Mean_winter_temp"] +
                          medParams$surv_TempSummerCoef*clim["Mean_summer_temp"] + 
                          medParams$surv_TempFallCoef*clim["Mean_fall_temp"])))



## PLOT
par(mfrow=c(1,1), mar=c(3.9,1.7,2.3,1.5))  #bottom, left, top and right 
#par(pty="s")

plot(pred.grwthGLM, pred.grwth)
plot(in.data, pred.grwth)
points(in.data, pred.grwthGLM, col="red")
