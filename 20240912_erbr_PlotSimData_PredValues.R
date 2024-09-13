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
name <- as.character("SimDat20yr")
## ------------------------------------------------------------------------------------------------





## LOAD FILES -------------------------------------------------------------------------------------
## Real data
erbr <- read.csv("erbr_TagClust2022_20230408.csv", header=TRUE)
erbr$Year <- as.factor(erbr$Year)


## Climate data
clim32yr <- read.csv("erbr_climData3seas32yr_221114.csv", header=TRUE)


## Simulated data
dateSim <- "20240907"
pathSim <- "C:/Users/april/Dropbox/CU_Boulder_PhD/DBG_Internship/Manuscript/MS_DataAndCode_archive/"
simDat1 <- read.csv(file=paste(pathSim,dateSim,"_erbr_",name,"NoMiss.1.Format4JAGS.csv",sep=""), header=TRUE)
simDat2 <- read.csv(file=paste(pathSim,dateSim,"_erbr_",name,"NoMiss.2.Format4JAGS.csv",sep=""), header=TRUE)
simDat3 <- read.csv(file=paste(pathSim,dateSim,"_erbr_",name,"NoMiss.3.Format4JAGS.csv",sep=""), header=TRUE)
simDat4 <- read.csv(file=paste(pathSim,dateSim,"_erbr_",name,"NoMiss.4.Format4JAGS.csv",sep=""), header=TRUE)
simDat5 <- read.csv(file=paste(pathSim,dateSim,"_erbr_",name,"NoMiss.5.Format4JAGS.csv",sep=""), header=TRUE)
simDat6 <- read.csv(file=paste(pathSim,dateSim,"_erbr_",name,"NoMiss.6.Format4JAGS.csv",sep=""), header=TRUE)
simDat7 <- read.csv(file=paste(pathSim,dateSim,"_erbr_",name,"NoMiss.7.Format4JAGS.csv",sep=""), header=TRUE)
simDat8 <- read.csv(file=paste(pathSim,dateSim,"_erbr_",name,"NoMiss.8.Format4JAGS.csv",sep=""), header=TRUE)
simDat9 <- read.csv(file=paste(pathSim,dateSim,"_erbr_",name,"NoMiss.9.Format4JAGS.csv",sep=""), header=TRUE)
simDat10 <- read.csv(file=paste(pathSim,dateSim,"_erbr_",name,"NoMiss.10.Format4JAGS.csv",sep=""), header=TRUE)


## 'True' params from JAGS mod of real data 
medParams.realDat <- readRDS("erbrMedParams_noYRE_20240803")


## GLMM results
dateGLM <- "20240912"
#name <- "SimDat20yr"
paramsMM.grwth <- readRDS(file=paste(dateGLM, "_erbr_paramMMgrwthWint_", name, ".rds", sep=""))
paramsMM.surv <- readRDS(file=paste(dateGLM, "_erbr_paramMMsurvWint_", name, ".rds", sep=""))
#glmm.grwth <- readRDS(file=paste(date, "_erbr_GLMMgrwth_", name, ".rds", sep=""))
#glmm.surv <- readRDS(file=paste(date, "_erbr_GLMMsurv_", name, ".rds", sep=""))


## JAGS results
#summ.mod1 <- readRDS("20240906_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat20yr.1.rds")
dateSUMM <- "20240906"
summ.mod1 <- readRDS(file=paste(dateSUMM,"_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_",name,".1.rds", sep=""))
summ.mod2 <- readRDS(file=paste(dateSUMM,"_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_",name,".2.rds", sep=""))
summ.mod3 <- readRDS(file=paste(dateSUMM,"_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_",name,".3.rds", sep=""))
summ.mod4 <- readRDS(file=paste(dateSUMM,"_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_",name,".4.rds", sep=""))
summ.mod5 <- readRDS(file=paste(dateSUMM,"_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_",name,".5.rds", sep=""))
summ.mod6 <- readRDS(file=paste(dateSUMM,"_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_",name,".6.rds", sep=""))
summ.mod7 <- readRDS(file=paste(dateSUMM,"_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_",name,".7.rds", sep=""))
summ.mod8 <- readRDS(file=paste(dateSUMM,"_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_",name,".8.rds", sep=""))
summ.mod9 <- readRDS(file=paste(dateSUMM,"_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_",name,".9.rds", sep=""))
summ.mod10 <- readRDS(file=paste(dateSUMM,"_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_",name,".10.rds", sep=""))
#summ.mod10 <- readRDS("20240907_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat20yr.10.rds")
## -----------------------------------------------------------------------------------------------





## GET CLIMATE VARIABLES ------------------------------------------------------------------------
#n.yrs <- 20 #Dataset length 
climYrs1 <- clim32yr[clim32yr$Year %in% unique(simDat1$ClimYr),]
#climYrs1 <- cbind(as.factor(rep(1,nrow(climYrs1))), climYrs1)
climYrs2 <- clim32yr[clim32yr$Year %in% unique(simDat2$ClimYr),]
#climYrs2 <- cbind(as.factor(rep(2,nrow(climYrs2))), climYrs2)
climYrs3 <- clim32yr[clim32yr$Year %in% unique(simDat3$ClimYr),]
#climYrs3 <- cbind(as.factor(rep(3,nrow(climYrs3))), climYrs3)
climYrs4 <- clim32yr[clim32yr$Year %in% unique(simDat4$ClimYr),]
#climYrs4 <- cbind(as.factor(rep(4,nrow(climYrs4))), climYrs4)
climYrs5 <- clim32yr[clim32yr$Year %in% unique(simDat5$ClimYr),]
#climYrs5 <- cbind(as.factor(rep(5,nrow(climYrs5))), climYrs5)
climYrs6 <- clim32yr[clim32yr$Year %in% unique(simDat6$ClimYr),]
#climYrs6 <- cbind(as.factor(rep(6,nrow(climYrs6))), climYrs6)
climYrs7 <- clim32yr[clim32yr$Year %in% unique(simDat7$ClimYr),]
#climYrs7 <- cbind(as.factor(rep(7,nrow(climYrs7))), climYrs7)
climYrs8 <- clim32yr[clim32yr$Year %in% unique(simDat8$ClimYr),]
#climYrs8 <- cbind(as.factor(rep(8,nrow(climYrs8))), climYrs8)
climYrs9 <- clim32yr[clim32yr$Year %in% unique(simDat9$ClimYr),]
#climYrs9 <- cbind(as.factor(rep(9,nrow(climYrs9))), climYrs9)
climYrs10 <- clim32yr[clim32yr$Year %in% unique(simDat10$ClimYr),]
#climYrs10 <- cbind(as.factor(rep(10,nrow(climYrs10))), climYrs10)

climYrs <- rbind(climYrs1, climYrs2, climYrs3, climYrs4, climYrs5, climYrs6, climYrs7, climYrs8, climYrs9, climYrs10)
reps <- c(rep(1,nrow(climYrs1)), rep(2,nrow(climYrs2)), rep(3,nrow(climYrs3)), rep(4,nrow(climYrs4)),
          rep(5,nrow(climYrs5)), rep(6,nrow(climYrs6)), rep(7,nrow(climYrs7)), rep(8,nrow(climYrs8)),
          rep(9,nrow(climYrs9)), rep(10,nrow(climYrs10)))
climYrs$Rep <- as.factor(reps)
## ---------------------------------------------------------------------------------------------




## GET RANGE OF PLANT SIZES FOR MAKING PREDICTIONS ---------------------------------------------
#in.data <- as.data.frame(binmids) 
#in.data <- 1:151
#colnames(in.data) <- "RosNew"
#$RosNew)

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
## ------------------------------------------------------------------------------------------





## --------------------------------------------------------------------------------------------
## OBTAIN TRUE VALUES FOR ADDING TO PLOTS
# re-format data from previous JAGS runs of real data
medParams.realDatTr <- as.data.frame(t(medParams.realDat))
medParams.realDatTr <- as.data.frame(cbind(as.numeric(medParams.realDatTr$`colMedians(as.matrix(chains))`),colnames(medParams.realDat)))
colnames(medParams.realDatTr) <- c("realData", "Name")

medParams.realDatTr <- medParams.realDatTr[55:70,]
medParams.realDatTr <- medParams.realDatTr[c(1:8,11:16),] #Remove GrwthVar  
medParams.realDatTr$ParamTitle <- c("Grwth Intercept","Grwth Size","Grwth Fall Temp","Grwth Summer Temp","Grwth Winter Temp",
                                    "Grwth Fall Precip","Grwth Summer Precip","Grwth Winter Precip","Surv Intercept","Surv Size",
                                    "Surv Winter Precip","Surv Fall Temp","Surv Summer Temp","Surv Winter Temp")
medParams.realDatGrwth <- medParams.realDatTr[1:8,]
medParams.realDatGrwth$realData <- as.numeric(as.character(medParams.realDatGrwth$realData))
medParams.realDatSurv <- medParams.realDatTr[9:14,]
medParams.realDatSurv$realData <- as.numeric(as.character(medParams.realDatSurv$realData))
## ------------------------------------------------------------------------------------------






## GET MEDIAN PARAM VALUES FROM JAGS RUN ----------------------------------------------------
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

medComb.surv <- as.data.frame(cbind(medParams.1[9:14], medParams.2[9:14],medParams.3[9:14],
                                     medParams.4[9:14],medParams.5[9:14],medParams.6[9:14],
                                     medParams.7[9:14],medParams.8[9:14],medParams.9[9:14],medParams.10[9:14]))
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
## Plug in median param vals, annual climate and selected size predictor values 
## Exclude random transect effects here 

## Make matrix to hold output data
column.names <- c("DATASET","REP","CLIM YR","PLT SZ","SurvRate_JAGS","SurvRate_GLMM","SurvRate_GLMM","GrwthRate_JAGS","GrwthRate_GLMM","GrwthRate_True")
#output <- as.data.frame(matrix(NA,nrow=(10*20*51), ncol=length(column.names)))
colnames(output) <- column.names
output$DATASET <- rep(name, nrow(output))
#output$REP <- 

## ** make list of climate yr lengths for each yr, number of rows for each rep is therefore num yrs*51 ***

for (rr in 1:n.datset) {  #Dataset rep loop
#rr <- 1
  ## Select climate years for given rep
  climYrs.sel <- climYrs[climYrs$Rep==rr,]
  #assign("climYrs", unique(simDat1$ClimYr)
  
  #output$REP <- rep(rr, nrow(climYrs.sel))
  #output.rep <- subset(output, REP==rr) #Subset output dataframe by dataset replicate

  
  for (yy in 1:nrow(climYrs.sel))
  ## Climate year 
  yy <- 1
  #clim.sel <- climYrs.sel
  
  ## All sizes
  
  ## Growth (neg binom)
  ## JAGS
  pred.grwthJAGS <- exp(medComb.grwth[1,rr] + medComb.grwth[2,rr]*log(binmids) 
                    + medComb.grwth[3,rr]*climYrs.sel$Mean_fall_temp[yy]
                    + medComb.grwth[4,rr]*climYrs.sel$Mean_summer_temp[yy]
                    + medComb.grwth[5,rr]*climYrs.sel$Mean_winter_temp[yy]
                    + medComb.grwth[6,rr]*climYrs.sel$Tot_fall_ppt[yy]
                    + medComb.grwth[7,rr]*climYrs.sel$Tot_summer_ppt[yy]
                    + medComb.grwth[8,rr]*climYrs.sel$Tot_winter_ppt[yy])
  
  ## GLMM
  pred.grwthGLM <- exp(paramsMM.grwth[1,rr] + paramsMM.grwth[2,rr]*log(binmids) 
                    + paramsMM.grwth[3,rr]*climYrs.sel$Mean_fall_temp[yy]
                    + paramsMM.grwth[4,rr]*climYrs.sel$Mean_summer_temp[yy]
                    + paramsMM.grwth[5,rr]*climYrs.sel$Mean_winter_temp[yy]
                    + paramsMM.grwth[6,rr]*climYrs.sel$Tot_fall_ppt[yy]
                    + paramsMM.grwth[7,rr]*climYrs.sel$Tot_summer_ppt[yy]
                    + paramsMM.grwth[8,rr]*climYrs.sel$Tot_winter_ppt[yy])
  
  ## True params
  pred.grwthTrue <- exp(medParams.realDatGrwth$realData[1] + medParams.realDatGrwth$realData[2]*log(binmids) 
                        + medParams.realDatGrwth$realData[3]*climYrs.sel$Mean_fall_temp[yy]
                        + medParams.realDatGrwth$realData[4]*climYrs.sel$Mean_summer_temp[yy]
                        + medParams.realDatGrwth$realData[5]*climYrs.sel$Mean_winter_temp[yy]
                        + medParams.realDatGrwth$realData[6]*climYrs.sel$Tot_fall_ppt[yy]
                        + medParams.realDatGrwth$realData[7]*climYrs.sel$Tot_summer_ppt[yy]
                        + medParams.realDatGrwth$realData[8]*climYrs.sel$Tot_winter_ppt[yy])
  
  
  ## Survival (binom)
  ## JAGS
  pred.survJAGS <- 1/(1+exp(-(medComb.surv[1,rr] + medComb.surv[2,rr]*log(binmids) + 
                            medComb.surv[3,rr]*climYrs.sel$Tot_winter_ppt[yy] + 
                            medComb.surv[4,rr]*climYrs.sel$Mean_fall_temp[yy] +
                            medComb.surv[5,rr]*climYrs.sel$Mean_summer_temp[yy] +
                            medComb.surv[6,rr]*climYrs.sel$Mean_winter_temp[yy])))
  
  ## GLMM
  pred.survGLMM <- 1/(1+exp(-(paramsMM.surv[1,rr] + paramsMM.surv[2,rr]*log(binmids) + 
                                paramsMM.surv[3,rr]*climYrs.sel$Tot_winter_ppt[yy] + 
                                paramsMM.surv[4,rr]*climYrs.sel$Mean_fall_temp[yy] +
                                paramsMM.surv[5,rr]*climYrs.sel$Mean_summer_temp[yy] +
                                paramsMM.surv[6,rr]*climYrs.sel$Mean_winter_temp[yy])))
  
  ## True params
  pred.survTrue <- 1/(1+exp(-(medParams.realDatSurv$realData[1] + medParams.realDatSurv$realData[2]*log(binmids) + 
                                medParams.realDatSurv$realData[3]*climYrs.sel$Tot_winter_ppt[yy] + 
                                medParams.realDatSurv$realData[4]*climYrs.sel$Mean_fall_temp[yy] +
                                medParams.realDatSurv$realData[5]*climYrs.sel$Mean_summer_temp[yy] +
                                medParams.realDatSurv$realData[6]*climYrs.sel$Mean_winter_temp[yy])))
  
  
  output$REP <- rep(rr, )
  output$GrwthRate_JAGS[1:51] <- pred.grwthJAGS 
  output$GrwthRate_GLMM <- pred.grwthGLM 
  output$GrwthRate_True <- pred.grwthTrue 
  
  output$SurvRate_JAGS <- pred.survJAGS 
  output$SurvRate_GLMM <- pred.survGLM 
  output$SurvRate_True <- pred.survTrue 
  
  } ## End year loop
    
    
} ## End dataset replicates loop
  



## PLOT
par(mfrow=c(1,1), mar=c(3.9,1.7,2.3,1.5))  #bottom, left, top and right 
#par(pty="s")

plot(pred.grwthGLM, pred.grwth)
plot(in.data, pred.grwth)
points(in.data, pred.grwthGLM, col="red")
