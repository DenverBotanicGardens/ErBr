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
name <- as.character("SimDat50yr")
#dats <- read.csv(file=paste(date,"_erbr_", name, dd, ".Format4JAGS", ".csv", sep=""), header=TRUE)
## ------------------------------------------------------------------------------------------------



## LOAD FILES -------------------------------------------------
clim32yr <- read.csv("erbr_climData3seas32yr_221114.csv", header=TRUE)

date <- "20240912"
name <- "SimDat50yr"
paramMMgrwth <- readRDS(file=paste(date, "_erbr_paramMMgrwth_", name, ".rds", sep=""))
paramMMsurv <- readRDS(file=paste(date, "_erbr_paramMMsurv_", name, ".rds", sep=""))
glmm.grwth <- readRDS(file=paste(date, "_erbr_GLMMgrwth_", name, ".rds", sep=""))
glmm.surv <- readRDS(file=paste(date, "_erbr_GLMMsurv_", name, ".rds", sep=""))


#summ.mod1 <- readRDS("20240906_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat20yr.1.rds")
#summ.mod1 <- readRDS("20240910_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat50yrMiss.1.rds")
#summ.mod2 <- readRDS("20240910_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat50yrMiss.2.rds")
#summ.mod3 <- readRDS("20240910_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat50yrMiss.3.rds")
#summ.mod4 <- readRDS("20240910_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat50yrMiss.4.rds")
#summ.mod5 <- readRDS("20240910_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat50yrMiss.5.rds")
#summ.mod6 <- readRDS("20240910_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat50yrMiss.6.rds")
#summ.mod7 <- readRDS("20240910_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat50yrMiss.7.rds")
#summ.mod8 <- readRDS("20240910_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat50yrMiss.8.rds")
#summ.mod9 <- readRDS("20240910_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat50yrMiss.9.rds")
summ.mod10 <- readRDS("20240910_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat50yrMiss.10.rds")
## --------------------------------



## GeT CLIMATE VARIABLES ----------
#n.yrs <- 51 #Assign number of years (plus 1) to simulate climate data for 

#Create empty variable to hold simulated climate data
#column.names <- colnames(clim32yr)
#sim.clim <- as.data.frame(matrix(NA, nrow=length(1:n.yrs), ncol=length(column.names)))
#colnames(sim.clim) <- column.names
#sim.clim$Year <- 1:n.yrs

#List of random numbers that corresponds to a set of climate values 
#rel.yrs <- 2002:2021 #Select relevant subset of yrs (this matches what was used for stoch lam to keep Methods consistent) 
#climYrs.rel <- clim32yr[clim32yr$Year>=2002 & clim32yr$Year<=2021,]
#randVals.yr <- sample(1:nrow(climYrs.rel), size=n.yrs, prob=NULL, replace=TRUE)

#for (cc in 1:length(1:n.yrs)) {
#  sim.clim[cc,2:7] <- climYrs.rel[randVals.yr[cc],2:7] 
#  sim.clim$Clim_yr[cc] <- climYrs.rel[randVals.yr[cc],1]
#}

## Example climate yrs 
SimDat <- read.csv("C:/Users/april/Dropbox/CU_Boulder_PhD/DBG_Internship/Manuscript/MS_DataAndCode_archive/20240907_erbr_SimDat20yrNoMiss.1.Format4JAGS.csv", header=TRUE)
SimDat$ClimYr[1:20]
clim <- clim32yr[clim32yr$Year==2015,]
## ----------------------------------------------------------------------------





## Get Median param values from JAGS run ---
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
## --------------------------------------




## Get GLMM param values -------------
#paramsMM.grwthTmp <- as.data.frame(summary(glmm.grwth)$coefficients[1:nrow(summary(glmm.grwth)$coefficients),1])
#colnames(paramsMM.grwthTemp) <- paste("GLMM.",dd,sep="")
#paramsMM.grwth <- as.data.frame(c(paramsMM.grwth, paramsMM.grwthTmp))

paramsMM.grwth <- as.data.frame(summary(glmm.grwth)$coefficients[1:nrow(summary(glmm.grwth)$coefficients),1])
colnames(paramsMM.grwth) <- paste("GLMM.",1,sep="")
paramsMM.grwth$ParamTitle <- c("Grwth Intercept","Grwth Size","Grwth Fall Precip","Grwth Winter Precip","Grwth Summer Precip",
                               "Grwth Fall Temp","Grwth Winter Temp","Grwth Summer Temp")

paramsMM.grwthOrd <- paramsMM.grwth[match(names.paramTitles[1:8], paramsMM.grwth$ParamTitle),]
paramsMM.survOrd <- paramsMM.surv[match(names.paramTitles[8:12], paramsMM.surv$ParamTitle),]
## --------------------------------------





## OBTAIN VITAL RATE PREDICTIONS ----------------------------------------------

## Initialize variable with size input data for predictions
#in.data <- as.data.frame(binmids) 
in.data <- 1:151
#colnames(in.data) <- "RosNew"
#$RosNew)

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
