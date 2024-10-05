## April Goebl & Dan Doak
## Script modified 2024-09-17
## Collaboration with CU Boulder and Denver Botanic Gardens on Eriogonum brandegeii modeling 
## Get predicted size and growth rate results from modeling of simulated demographic data 
## to show JAGS model performance with missing data versus GLMM estimates 


rm(list=ls())
graphics.off()


library(MASS)
library(dplyr)
library(matrixStats)
library(corrplot)
#library(rgl)
#library(viridis) 
library(car)
library(coda)
library(lme4)
library(resample)
library(gplots)
library(stringr)
library(plotrix)



## ASSIGN NAME VARIABLE FOR DESIRED DATASETS 
name <- as.character("SimDat40yr") #"SimDat20yr" "SimDat20yrHiGr" "SimDat20yrMedGr"
## ------------------------------------------------------------------------------------------------





## LOAD FILES -------------------------------------------------------------------------------------
## Real data
erbr <- read.csv("erbr_TagClust2022_20230408.csv", header=TRUE)
erbr$Year <- as.factor(erbr$Year)


## Climate data
clim32yr <- read.csv("erbr_climData3seas32yr_221114.csv", header=TRUE)
clim32yrMAXES <- read.csv("erbr_climData3seas32yr_MAXES.csv", header=TRUE)


## Simulated data for finding relevant climate years (could be miss or no-miss)
dateSim <- "20241005"
nameSim <-   "SimDat20yrMedGrNoMiss.srvCor.sdlgCor.grwthCor." #"SimDat20yrHiGrNoMiss.srvCor.sdlgCor.
pathSim <-  "C:/Users/april/Dropbox/CU_Boulder_PhD/DBG_Internship/Manuscript/MS_DataAndCode_archive/" #getwd() #
simDat1 <- read.csv(file=paste(pathSim,dateSim,"_erbr_",nameSim,"1.4JAGS.csv",sep=""), header=TRUE)
simDat2 <- read.csv(file=paste(pathSim,dateSim,"_erbr_",nameSim,"2.4JAGS.csv",sep=""), header=TRUE)
simDat3 <- read.csv(file=paste(pathSim,dateSim,"_erbr_",nameSim,"3.4JAGS.csv",sep=""), header=TRUE)
simDat4 <- read.csv(file=paste(pathSim,dateSim,"_erbr_",nameSim,"4.4JAGS.csv",sep=""), header=TRUE)
simDat5 <- read.csv(file=paste(pathSim,dateSim,"_erbr_",nameSim,"5.4JAGS.csv",sep=""), header=TRUE)
simDat6 <- read.csv(file=paste(pathSim,dateSim,"_erbr_",nameSim,"6.4JAGS.csv",sep=""), header=TRUE)
simDat7 <- read.csv(file=paste(pathSim,dateSim,"_erbr_",nameSim,"7.4JAGS.csv",sep=""), header=TRUE)
simDat8 <- read.csv(file=paste(pathSim,dateSim,"_erbr_",nameSim,"8.4JAGS.csv",sep=""), header=TRUE)
simDat9 <- read.csv(file=paste(pathSim,dateSim,"_erbr_",nameSim,"9.4JAGS.csv",sep=""), header=TRUE)
simDat10 <- read.csv(file=paste(pathSim,dateSim,"_erbr_",nameSim,"10.4JAGS.csv",sep=""), header=TRUE)



## 'True' params from JAGS mod of real data (or alternatative vales for alt LHs)
#medParams.realDat <- readRDS("erbrMedParams_noYRE_20240803")
medParams.realDat <- readRDS("erbrParams_MedGrAltLH_20240924.rds")
#medParams.realDat <- readRDS("erbrParams_HiGrAltLHcor_20240925.rds")





## GLMM results or files to run GLMM models below
name
#glmmSumm.grwth <- readRDS(file=paste("20241001_erbr_GLMMtmbSummGrwth_", name, ".rds", sep=""))
#glmmSumm.surv <- readRDS(file=paste("20241001_erbr_GLMMsummSurv_", name, ".rds", sep=""))


paramsMM.grwth <- readRDS(file=paste("20241005_erbr_paramMMtmbGrwth_", name, ".rds", sep=""))
seMM.grwth <- readRDS(file=paste("20241005_erbr_seMMtmbGrwth_", name,".rds", sep=""))
paramsMM.surv <- readRDS(file=paste("20241005_erbr_paramMMsurv_", name, ".rds", sep=""))
seMM.surv <- readRDS(file=paste("20241005_erbr_seMMsurv_", name, ".rds", sep=""))





## JAGS results
## MISSING
dateSUMM <- "20241002"
nameSUMM <-  "SimDat20yrMedGrMiss.srvCor.sdlgCor.grwthCor." #"SimDat20yrHiGrMiss.srvCor." #"SimDat20yrMedGrMiss.srvCor.sdlgCor."
summ.modMs1 <- readRDS(file=paste(dateSUMM,"_erbr_JAGSmodBestSUMM_",nameSUMM,"1.rds", sep=""))
summ.modMs2 <- readRDS(file=paste(dateSUMM,"_erbr_JAGSmodBestSUMM_",nameSUMM,"2.rds", sep=""))
summ.modMs3 <- readRDS(file=paste(dateSUMM,"_erbr_JAGSmodBestSUMM_",nameSUMM,"3.rds", sep=""))
summ.modMs4 <- readRDS(file=paste(dateSUMM,"_erbr_JAGSmodBestSUMM_",nameSUMM,"4.rds", sep=""))
summ.modMs5 <- readRDS(file=paste(dateSUMM,"_erbr_JAGSmodBestSUMM_",nameSUMM,"5.rds", sep=""))
summ.modMs6 <- readRDS(file=paste(dateSUMM,"_erbr_JAGSmodBestSUMM_",nameSUMM,"6.rds", sep=""))
summ.modMs7 <- readRDS(file=paste(dateSUMM,"_erbr_JAGSmodBestSUMM_",nameSUMM,"7.rds", sep=""))
summ.modMs8 <- readRDS(file=paste(dateSUMM,"_erbr_JAGSmodBestSUMM_",nameSUMM,"8.rds", sep=""))
summ.modMs9 <- readRDS(file=paste(dateSUMM,"_erbr_JAGSmodBestSUMM_",nameSUMM,"9.rds", sep=""))
summ.modMs10 <- readRDS(file=paste(dateSUMM,"_erbr_JAGSmodBestSUMM_",nameSUMM,"10.rds", sep=""))


## NO MISSING
dateSUMM <- "20240926"
nameSUMMno <- "SimDat20yrNoMiss.srvCor.sdlgCor."  #nameSim #"SimDat40yrMiss"
summ.modNo1 <- readRDS(file=paste(dateSUMM,"_erbr_JAGSmodBestSUMM_",nameSUMMno,"1.rds", sep=""))
summ.modNo2 <- readRDS(file=paste(dateSUMM,"_erbr_JAGSmodBestSUMM_",nameSUMMno,"2.rds", sep=""))
summ.modNo3 <- readRDS(file=paste(dateSUMM,"_erbr_JAGSmodBestSUMM_",nameSUMMno,"3.rds", sep=""))
summ.modNo4 <- readRDS(file=paste(dateSUMM,"_erbr_JAGSmodBestSUMM_",nameSUMMno,"4.rds", sep=""))
summ.modNo5 <- readRDS(file=paste(dateSUMM,"_erbr_JAGSmodBestSUMM_",nameSUMMno,"5.rds", sep=""))
summ.modNo6 <- readRDS(file=paste(dateSUMM,"_erbr_JAGSmodBestSUMM_",nameSUMMno,"6.rds", sep=""))
summ.modNo7 <- readRDS(file=paste(dateSUMM,"_erbr_JAGSmodBestSUMM_",nameSUMMno,"7.rds", sep=""))
summ.modNo8 <- readRDS(file=paste(dateSUMM,"_erbr_JAGSmodBestSUMM_",nameSUMMno,"8.rds", sep=""))
summ.modNo9 <- readRDS(file=paste(dateSUMM,"_erbr_JAGSmodBestSUMM_",nameSUMMno,"9.rds", sep=""))
summ.modNo10 <- readRDS(file=paste(dateSUMM,"_erbr_JAGSmodBestSUMM_",nameSUMMno,"10.rds", sep=""))
## -----------------------------------------------------------------------------------------------






## GET CLIMATE VARIABLES ------------------------------------------------------------------------
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

climYrs <- rbind(climYrs1, climYrs2, climYrs3, climYrs4, climYrs5, climYrs6, climYrs7, climYrs8, climYrs9, climYrs10)
reps <- c(rep(1,nrow(climYrs1)), rep(2,nrow(climYrs2)), rep(3,nrow(climYrs3)), rep(4,nrow(climYrs4)),
          rep(5,nrow(climYrs5)), rep(6,nrow(climYrs6)), rep(7,nrow(climYrs7)), rep(8,nrow(climYrs8)),
          rep(9,nrow(climYrs9)), rep(10,nrow(climYrs10)))
climYrs$Rep <- as.factor(reps)
## ---------------------------------------------------------------------------------------------

## SCALE CLIMATE VARIABLES ---------------------------------------------------------------------
climYrs$Tot_fall_ppt <- climYrs$Tot_fall_ppt / clim32yrMAXES$Tot_fall_ppt
climYrs$Tot_winter_ppt <- climYrs$Tot_winter_ppt / clim32yrMAXES$Tot_winter_ppt
climYrs$Tot_summer_ppt <- climYrs$Tot_summer_ppt / clim32yrMAXES$Tot_summer_ppt
climYrs$Mean_fall_temp <- climYrs$Mean_fall_temp / clim32yrMAXES$Mean_fall_temp
climYrs$Mean_winter_temp <- climYrs$Mean_winter_temp / clim32yrMAXES$Mean_winter_temp
climYrs$Mean_summer_temp <- climYrs$Mean_summer_temp / clim32yrMAXES$Mean_summer_temp
## ---------------------------------------------------------------------------------------------






## GET RANGE OF PLANT SIZES FOR MAKING PREDICTIONS ---------------------------------------------

# Specify min and max plt sz
minsize <- 1
#maxsize <- 151 #(max(erbr$RosNew, na.rm=TRUE)) 
#maxsize <- 600  ## ** Try this for Fast gr LH datasets **
maxsize <- 80   ## ** Try this for Med gr LH datasets **
 
## For median sz estimation
##new size density estimation for median size estimation
pdfsz=density(erbr$RosNew, n=1024, cut=0, na.rm=TRUE) 
pdfsz2=cbind(pdfsz$x,pdfsz$y)
## This is set of smoothed vals that can be used w weightedMedian in matrixStats package to get a 'good' median for each class

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



## For HI and MED GR
bin.num <- 50  #Define number of bins 

## Improved method of finding median size/ bin mids (code from Dan)
vec.bin = c(minsize, minsize+1:bin.num*(maxsize-minsize)*(1/bin.num)) 
## Do this block to make medians the focal estimated size for each category
binmids = rep(NA, length(vec.bin)-1)

for (bb in 1:length(vec.bin)-1) {
  bounds <- c(vec.bin[bb], vec.bin[bb+1])
  binmids[bb] <- median(bounds)
}

n.bin = length(binmids)
## ------------------------------------------------------------------------------------------





## --------------------------------------------------------------------------------------------
## OBTAIN TRUE VALUES FOR ADDING TO PLOTS
# re-format data from previous JAGS runs of real data
medParams.realDatTr <- as.data.frame(t(medParams.realDat))
medParams.realDatTr <- as.data.frame(cbind(as.numeric(medParams.realDatTr$`colMedians(as.matrix(chains))`),colnames(medParams.realDat)))
colnames(medParams.realDatTr) <- c("realData", "Name")

medParams.realDatTr <- medParams.realDatTr[55:70,]
medParams.realDatTr.grow.r = medParams.realDatTr[9:10,] #Save the growth variance (r) parameters
medParams.realDatTr <- medParams.realDatTr[c(1:8,11:16),] #Remove GrwthVar  

medParams.realDatTr$ParamTitle <- c("Grwth Intercept","Grwth Size","Grwth Fall Temp","Grwth Summer Temp","Grwth Winter Temp",
                                    "Grwth Fall Precip","Grwth Summer Precip","Grwth Winter Precip","Surv Intercept","Surv Size",
                                    "Surv Winter Precip","Surv Fall Temp","Surv Summer Temp","Surv Winter Temp")
medParams.realDatGrwth <- medParams.realDatTr[1:8,]
medParams.realDatGrwth$realData <- as.numeric(as.character(medParams.realDatGrwth$realData))
medParams.realDatSurv <- medParams.realDatTr[9:14,]
medParams.realDatSurv$realData <- as.numeric(as.character(medParams.realDatSurv$realData))


## RESCALE 'TRUE' PARAMS --------------------------------------------------------------------
tempFsc <- medParams.realDatGrwth$realData[medParams.realDatGrwth$ParamTitle=="Grwth Fall Temp"] * clim32yrMAXES$Mean_fall_temp
tempSsc <- medParams.realDatGrwth$realData[medParams.realDatGrwth$ParamTitle=="Grwth Summer Temp"] * clim32yrMAXES$Mean_summer_temp
tempWsc <- medParams.realDatGrwth$realData[medParams.realDatGrwth$ParamTitle=="Grwth Winter Temp"] * clim32yrMAXES$Mean_winter_temp
PptFsc <- medParams.realDatGrwth$realData[medParams.realDatGrwth$ParamTitle=="Grwth Fall Precip"] * clim32yrMAXES$Tot_fall_ppt
PptSsc <- medParams.realDatGrwth$realData[medParams.realDatGrwth$ParamTitle=="Grwth Summer Precip"] * clim32yrMAXES$Tot_summer_ppt
PptWsc <- medParams.realDatGrwth$realData[medParams.realDatGrwth$ParamTitle=="Grwth Winter Precip"] * clim32yrMAXES$Tot_winter_ppt

medParams.realDatGrwth$realData[medParams.realDatGrwth$ParamTitle=="Grwth Fall Temp"] <- tempFsc
medParams.realDatGrwth$realData[medParams.realDatGrwth$ParamTitle=="Grwth Summer Temp"] <- tempSsc
medParams.realDatGrwth$realData[medParams.realDatGrwth$ParamTitle=="Grwth Winter Temp"] <- tempWsc
medParams.realDatGrwth$realData[medParams.realDatGrwth$ParamTitle=="Grwth Fall Precip"] <- PptFsc
medParams.realDatGrwth$realData[medParams.realDatGrwth$ParamTitle=="Grwth Summer Precip"] <- PptSsc 
medParams.realDatGrwth$realData[medParams.realDatGrwth$ParamTitle=="Grwth Winter Precip"] <- PptWsc 


tempFsc.surv <- medParams.realDatSurv$realData[medParams.realDatSurv$ParamTitle=="Surv Fall Temp"] * clim32yrMAXES$Mean_fall_temp
tempSsc.surv <- medParams.realDatSurv$realData[medParams.realDatSurv$ParamTitle=="Surv Summer Temp"] * clim32yrMAXES$Mean_summer_temp
tempWsc.surv <- medParams.realDatSurv$realData[medParams.realDatSurv$ParamTitle=="Surv Winter Temp"] * clim32yrMAXES$Mean_winter_temp
PptFsc.surv <- medParams.realDatSurv$realData[medParams.realDatSurv$ParamTitle=="Surv Fall Precip"] * clim32yrMAXES$Tot_fall_ppt
PptSsc.surv <- medParams.realDatSurv$realData[medParams.realDatSurv$ParamTitle=="Surv Summer Precip"] * clim32yrMAXES$Tot_summer_ppt
PptWsc.surv <- medParams.realDatSurv$realData[medParams.realDatSurv$ParamTitle=="Surv Winter Precip"] * clim32yrMAXES$Tot_winter_ppt

medParams.realDatSurv$realData[medParams.realDatSurv$ParamTitle=="Surv Fall Temp"] <- tempFsc.surv
medParams.realDatSurv$realData[medParams.realDatSurv$ParamTitle=="Surv Summer Temp"] <- tempSsc.surv
medParams.realDatSurv$realData[medParams.realDatSurv$ParamTitle=="Surv Winter Temp"] <- tempWsc.surv
medParams.realDatSurv$realData[medParams.realDatSurv$ParamTitle=="Surv Fall Precip"] <- PptFsc.surv
medParams.realDatSurv$realData[medParams.realDatSurv$ParamTitle=="Surv Summer Precip"] <- PptSsc.surv 
medParams.realDatSurv$realData[medParams.realDatSurv$ParamTitle=="Surv Winter Precip"] <- PptWsc.surv 
## ------------------------------------------------------------------------------------------






## GET MEDIAN PARAM VALUES FROM JAGS RUN ----------------------------------------------------
names.param <- rownames(summ.modMs1)[c(26:33,36:41)]
names.paramTitles <- c("Grwth Intercept","Grwth Size","Grwth Fall Temp","Grwth Summer Temp","Grwth Winter Temp",
                       "Grwth Fall Precip","Grwth Summer Precip","Grwth Winter Precip","Surv Intercept","Surv Size",
                       "Surv Winter Precip","Surv Fall Temp","Surv Summer Temp","Surv Winter Temp")

## Missing data
medParamsMs.1 <- summ.modMs1[c(26:33,36:41),2]
medParamsMs.2 <- summ.modMs2[c(26:33,36:41),2]
medParamsMs.3 <- summ.modMs3[c(26:33,36:41),2]
medParamsMs.4 <- summ.modMs4[c(26:33,36:41),2]
medParamsMs.5 <- summ.modMs5[c(26:33,36:41),2]
medParamsMs.6 <- summ.modMs6[c(26:33,36:41),2]
medParamsMs.7 <- summ.modMs7[c(26:33,36:41),2]
medParamsMs.8 <- summ.modMs8[c(26:33,36:41),2]
medParamsMs.9 <- summ.modMs9[c(26:33,36:41),2]
medParamsMs.10 <- summ.modMs10[c(26:33,36:41),2]

medCombMs.grwth <- as.data.frame(cbind(medParamsMs.1[1:8], medParamsMs.2[1:8],medParamsMs.3[1:8],
                                     medParamsMs.4[1:8],medParamsMs.5[1:8],medParamsMs.6[1:8],
                                     medParamsMs.7[1:8],medParamsMs.8[1:8],medParamsMs.9[1:8],medParamsMs.10[1:8]))

medCombMs.surv <- as.data.frame(cbind(medParamsMs.1[9:14], medParamsMs.2[9:14],medParamsMs.3[9:14],
                                     medParamsMs.4[9:14],medParamsMs.5[9:14],medParamsMs.6[9:14],
                                     medParamsMs.7[9:14],medParamsMs.8[9:14],medParamsMs.9[9:14],medParamsMs.10[9:14]))



## Nomissing data
medParamsNo.1 <- summ.modNo1[c(26:33,36:41),2]
medParamsNo.2 <- summ.modNo2[c(26:33,36:41),2]
medParamsNo.3 <- summ.modNo3[c(26:33,36:41),2]
medParamsNo.4 <- summ.modNo4[c(26:33,36:41),2]
medParamsNo.5 <- summ.modNo5[c(26:33,36:41),2]
medParamsNo.6 <- summ.modNo6[c(26:33,36:41),2]
medParamsNo.7 <- summ.modNo7[c(26:33,36:41),2]
medParamsNo.8 <- summ.modNo8[c(26:33,36:41),2]
medParamsNo.9 <- summ.modNo9[c(26:33,36:41),2]
medParamsNo.10 <- summ.modNo10[c(26:33,36:41),2]

medCombNo.grwth <- as.data.frame(cbind(medParamsNo.1[1:8], medParamsNo.2[1:8],medParamsNo.3[1:8],
                                     medParamsNo.4[1:8],medParamsNo.5[1:8],medParamsNo.6[1:8],
                                     medParamsNo.7[1:8],medParamsNo.8[1:8],medParamsNo.9[1:8],medParamsNo.10[1:8]))

medCombNo.surv <- as.data.frame(cbind(medParamsNo.1[9:14], medParamsNo.2[9:14],medParamsNo.3[9:14],
                                    medParamsNo.4[9:14],medParamsNo.5[9:14],medParamsNo.6[9:14],
                                    medParamsNo.7[9:14],medParamsNo.8[9:14],medParamsNo.9[9:14],medParamsNo.10[9:14]))
## --------------------------------------------------------------------------------------





## EXTRACT ESTIMATES FROM GLMM SUMMARIES -----------------------------------------------
## --------------------------------------------------------------------------------------





## OBTAIN VITAL RATE PREDICTIONS ----------------------------------------------
## Plug in median param vals, annual climate and selected size predictor values 
## Exclude random transect effect 

## Make matrix to hold output data
column.names <- c("DATASET","REP","CLIM_YR","PLT_SZ","SurvRate_JAGS","SurvRate_GLMM","SurvRate_True",
                  "GrwthRate_JAGS","GrwthRate_GLMM","GrwthRate_True")

column.names <- c("DATASET","REP","CLIM_YR","PLT_SZ","SurvRate_JAGS","SurvRate_JAGSnoMs","SurvRate_GLMM","SurvRate_True",
                  "GrwthRate_JAGS","GrwthRate_JAGSnoMs","GrwthRate_GLMM","GrwthRate_True")

output <- NULL
n.datset <- 10

for (rr in 1:n.datset) {  #Dataset rep loop
  
  ## Select climate years for given rep
  climYrs.sel <- climYrs[climYrs$Rep==rr,]
  output.rep <- NULL
  
  for (yy in 1:nrow(climYrs.sel))  { #Year loop
  ## All sizes
  
  ## Growth (neg binom)
  ## JAGS
  pred.grwthJAGS <- exp(medCombMs.grwth[1,rr] + medCombMs.grwth[2,rr]*log(binmids) 
                    + medCombMs.grwth[3,rr]*climYrs.sel$Mean_fall_temp[yy]
                    + medCombMs.grwth[4,rr]*climYrs.sel$Mean_summer_temp[yy]
                    + medCombMs.grwth[5,rr]*climYrs.sel$Mean_winter_temp[yy]
                    + medCombMs.grwth[6,rr]*climYrs.sel$Tot_fall_ppt[yy]
                    + medCombMs.grwth[7,rr]*climYrs.sel$Tot_summer_ppt[yy]
                    + medCombMs.grwth[8,rr]*climYrs.sel$Tot_winter_ppt[yy])
  
  ## JAGS - no missing data
  #pred.grwthJAGSno <- exp(medCombNo.grwth[1,rr] + medCombNo.grwth[2,rr]*log(binmids) 
   #                     + medCombNo.grwth[3,rr]*climYrs.sel$Mean_fall_temp[yy]
    #                    + medCombNo.grwth[4,rr]*climYrs.sel$Mean_summer_temp[yy]
     #                   + medCombNo.grwth[5,rr]*climYrs.sel$Mean_winter_temp[yy]
      #                  + medCombNo.grwth[6,rr]*climYrs.sel$Tot_fall_ppt[yy]
       #                 + medCombNo.grwth[7,rr]*climYrs.sel$Tot_summer_ppt[yy]
        #                + medCombNo.grwth[8,rr]*climYrs.sel$Tot_winter_ppt[yy])
  
  ## GLMM
  pred.grwthGLMM <- exp(paramsMM.grwth[1,rr] + paramsMM.grwth[2,rr]*log(binmids) 
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
  
  pred.grow.r <- exp(as.numeric(as.character(medParams.realDatTr.grow.r[1,1])) + 
                     as.numeric(as.character(medParams.realDatTr.grow.r[2,1]))*log(binmids)) 
  dense.zero <- dnbinom(0, size=pred.grow.r, mu=pred.grwthTrue, log = FALSE)
  pred.grwthTrue <- dense.zero+pred.grwthTrue*(1- dense.zero)
  
  
  
  ## Survival (binom)
  ## JAGS
  pred.survJAGS <- 1/(1+exp(-(medCombMs.surv[1,rr] + medCombMs.surv[2,rr]*log(binmids) + 
                            medCombMs.surv[3,rr]*climYrs.sel$Tot_winter_ppt[yy] + 
                            medCombMs.surv[4,rr]*climYrs.sel$Mean_fall_temp[yy] +
                            medCombMs.surv[5,rr]*climYrs.sel$Mean_summer_temp[yy] +
                            medCombMs.surv[6,rr]*climYrs.sel$Mean_winter_temp[yy])))
  
  ## JAGS - no missing data
  #pred.survJAGSno <- 1/(1+exp(-(medCombNo.surv[1,rr] + medCombNo.surv[2,rr]*log(binmids) + 
   #                             medCombNo.surv[3,rr]*climYrs.sel$Tot_winter_ppt[yy] + 
    #                            medCombNo.surv[4,rr]*climYrs.sel$Mean_fall_temp[yy] +
     #                           medCombNo.surv[5,rr]*climYrs.sel$Mean_summer_temp[yy] +
      #                          medCombNo.surv[6,rr]*climYrs.sel$Mean_winter_temp[yy])))
  
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
  
  output.yr <- as.data.frame(matrix(NA,nrow=(length(binmids)), ncol=length(column.names)))
  colnames(output.yr) <- column.names
  output.yr$REP <- rr
  output.yr$CLIM_YR <- climYrs.sel$Year[yy]
  output.yr$PLT_SZ <- binmids
  
  output.yr$GrwthRate_JAGS <- pred.grwthJAGS 
  #output.yr$GrwthRate_JAGSnoMs <- pred.grwthJAGSno
  output.yr$GrwthRate_GLMM <- pred.grwthGLMM 
  output.yr$GrwthRate_True <- pred.grwthTrue 
  
  output.yr$SurvRate_JAGS <- pred.survJAGS 
  #output.yr$SurvRate_JAGSnoMs <- pred.survJAGSno 
  output.yr$SurvRate_GLMM <- pred.survGLMM 
  output.yr$SurvRate_True <- pred.survTrue 
  
  output.rep <- rbind(output.rep, output.yr)
  
  } ## End year loop
  
  output <- rbind(output, output.rep)

} ## End dataset replicates loop
 
 
output$DATASET <- name



## Save output 
write.csv(output, file="20241005_erbrSimDat20yrMedGr_predVals.csv", row.names=FALSE)
#file=paste(date, name, dd, ".csv", sep=""), row.names=FALSE)

## LOAD PREDICTIONS
#output <- read.csv("20240924_erbrSimDat20yrHiGr_predVals.csv", header=TRUE) ## **Try increasing max sz to 300 **
#output <- read.csv("20240924_erbrSimDat20yrMedGr_predVals.csv", header=TRUE)
#output <- read.csv("20240922_erbrSimDat40yr_predVals.csv", header=TRUE)
## -----------------------------------------------------------------------------------






## SHOW GRWTH AND SURV VS SZ FOR EACH LH --------------------------------------------- 
## This would be a good way to show the basic differences in the LHs 
#simDatComb <- simDatComb %>% mutate(RosNew1=lead(RosNew))  
#plot(simDatComb$RosNew, simDatComb$RosNew1, xlab="Plant size in year t", ylab="Plant size in year t+1",
#     main=name, pch=19,cex=0.6) 

## Do for mean climate year ***
plot(output$PLT_SZ, output$SurvRate_True, xlab="Plant size", ylab="Survival rate based on input parameters",
     pch=19, cex=0.6, main=name, ylim=c(0,1))






## PLOT VR PREDICTIONS ---------------------------------------------------------
par(pty="s")
par(mfrow=c(1,1), mar=c(4,4,2,1))  #bottom, left, top and right 

## Plot most recent rep from above
plot(binmids, pred.grwthJAGS)
     points(binmids, pred.grwthGLMM, col="red")
     points(binmids, pred.grwthTrue, col="blue")
     
plot(binmids, pred.survJAGS)
     points(binmids, pred.survGLMM, col="red")
     points(binmids, pred.survTrue, col="blue")

## Plot predicted VRS against plant size 
#plot(log(output$PLT_SZ),output$GrwthRate_True, main=name)  
#plot(output$PLT_SZ,output$GrwthRate_GLMM, main="20 years - Med Growth LH")  
#plot(log(output$PLT_SZ),output$GrwthRate_JAGS, main=name)  

#plot(log(output$PLT_SZ),output$SurvRate_True, main=name)  
#plot(output$PLT_SZ,output$SurvRate_GLMM, main="20 years - Med Growth LH")  
#plot(log(output$PLT_SZ),output$SurvRate_JAGS, main=name)  



## Plot reps separately 
## GROWTH
#par(mfrow=c(4,3), mar=c(4,4,2,1))  #bottom, left, top and right 
#for (rr in 1:10) {
#plot(output$PLT_SZ[output$REP==rr], output$GrwthRate_True[output$REP==rr], xlab="plant size", ylab="True predicted growth",
#     cex=0.65, pch=16)
#}
#plot.new()
#legend("center", "20 YEAR\nMEDIUM GRWTH LH", bty="n",cex=1.15)

#par(mfrow=c(4,3), mar=c(4,4,2,1))  #bottom, left, top and right 
#for (rr in 1:10) {
#  plot(output$PLT_SZ[output$REP==rr], output$GrwthRate_JAGS[output$REP==rr], xlab="plant size", ylab="MCMC predicted growth",
#       cex=0.65, pch=16, ylim=c(0,350))
#}
#plot.new()
#legend("center", "20 YEAR\nMEDIUM GRWTH LH", bty="n",cex=1.15)

#par(mfrow=c(4,3), mar=c(4,4,2,1))  #bottom, left, top and right 
#for (rr in 1:10) {
#  plot(output$PLT_SZ[output$REP==rr], output$GrwthRate_GLMM[output$REP==rr], xlab="plant size", ylab="MCMC predicted growth",
#       cex=0.65, pch=16, ylim=c(0,175))
#}
#plot.new()
#legend("center", "20 YEAR\nMEDIUM GRWTH LH", bty="n",cex=1.15)
#plot(output$PLT_SZ[output$REP==rr & output$CLIM_YR==2002], output$GrwthRate_True[output$REP==rr & output$CLIM_YR==2002])


## SURVIVAL
## Plot reps separately 
#par(mfrow=c(4,3), mar=c(4,4,2,1))  #bottom, left, top and right 
#for (rr in 1:10) {
#  plot(output$PLT_SZ[output$REP==rr], output$SurvRate_True[output$REP==rr], xlab="plant size", ylab="True predicted survival",
#       cex=0.65, pch=16, ylim=c(0.2,1))
#}
#plot.new()
#legend("center", "20 YEAR\nMEDIUM GRWTH LH", bty="n",cex=1.15)

#par(mfrow=c(4,3), mar=c(4,4,2,1))  #bottom, left, top and right 
#for (rr in 1:10) {
#  plot(output$PLT_SZ[output$REP==rr], output$SurvRate_JAGS[output$REP==rr], xlab="plant size", ylab="MCMC predicted survival",
#       cex=0.65, pch=16, ylim=c(0.2,1))
#}
#plot.new()
#legend("center", "20 YEAR\nMEDIUM GRWTH LH", bty="n",cex=1.15)

#par(mfrow=c(4,3), mar=c(4,4,2,1))  #bottom, left, top and right 
#for (rr in 1:10) {
#  plot(output$PLT_SZ[output$REP==rr], output$SurvRate_GLMM[output$REP==rr], xlab="plant size", ylab="MCMC predicted survival",
#       cex=0.65, pch=16, ylim=c(0.2,1))
#}
#plot.new()
#legend("center", "20 YEAR\nMEDIUM GRWTH LH", bty="n",cex=1.15)
## ------------------------------



## Histograms of size ------------
par(mfrow=c(1,1), mar=c(4,4,2,1))  #bottom, left, top and right 
simDatComb<-rbind(simDat1,simDat2,simDat3,simDat4,simDat5,simDat6,simDat7,simDat8,simDat9,simDat10)
hist(simDatComb$RosNew, xlim=c(0,100), breaks=40, xlab="Plant size", main=name)   
## ------------------------------



## Plot all reps, years, sizes together           
par(mfrow=c(2,2), mar=c(4,4,2,1))  #bottom, left, top and right 
par(pty="s")

## Growth   
mainTitle <- name
minAx <- 0
maxAx <- 50
legpos <- "bottomright"

plot(output$GrwthRate_True, output$GrwthRate_GLMM,col=alpha("grey40",0.5),main=mainTitle,
     xlab="True vital rate - GROWTH",ylab="Estimated vital rate - GROWTH", ylim=c(minAx,maxAx), xlim=c(minAx,maxAx))
abline(a=0, b=1)
legend(legpos, "GLMM no-missing data", col=alpha("grey40",0.5),pch=19, cex=1,horiz=FALSE, bty="y")

plot(output$GrwthRate_True, output$GrwthRate_JAGS,col=alpha("purple",0.5),main=mainTitle,
     xlab="True vital rate - GROWTH",ylab="Estimated vital rate - GROWTH",ylim=c(minAx,maxAx), xlim=c(minAx,maxAx))#,
abline(a=0, b=1)
legend(legpos, "MCMC missing data", col=alpha("purple",0.5),pch=19, cex=1,horiz=FALSE, bty="y")


plot(output$GrwthRate_True, output$GrwthRate_JAGSnoMs,col=alpha("pink",0.5),main=mainTitle,
     xlab="True vital rate - GROWTH",ylab="Estimated vital rate - GROWTH", ylim=c(minAx,maxAx), xlim=c(minAx,maxAx))
abline(a=0, b=1)
legend(legpos, "MCMC no-missing data", col=alpha("pink",0.5),pch=19, cex=1,horiz=FALSE, bty="y")
#legend("bottomright", c("GLMM","MCMC"), col=c(alpha("grey40",0.5),alpha("purple",0.3)), pch=19, cex=1,
#       horiz=FALSE, bty="y",seg.len=1)

 
## Survival 
minAx <- 0
maxAx <- 1
legPos <- "bottomright"
mainTitle <- name

plot(output$SurvRate_True, output$SurvRate_GLMM,col=alpha("grey40",0.5),main=mainTitle,
     xlab="True vital rate - SURVIVAL",ylab="Estimated vital rate - SURVIVAL",ylim=c(minAx,maxAx), xlim=c(minAx,maxAx))
abline(a=0, b=1)
legend(legPos, "GLMM no-missing data", col=alpha("grey40",0.5),pch=19, cex=1,horiz=FALSE, bty="y")

plot(output$SurvRate_True, output$SurvRate_JAGS,col=alpha("purple",0.5),main=mainTitle,
     xlab="True vital rate - SURVIVAL",ylab="Estimated vital rate - SURVIVAL",ylim=c(minAx,maxAx), xlim=c(minAx,maxAx))
abline(a=0, b=1)
legend(legPos, "MCMC missing data", col=alpha("purple",0.5),pch=19, cex=1,horiz=FALSE, bty="y")


plot(output$SurvRate_True, output$SurvRate_JAGSnoMs,col=alpha("pink",0.5),main=mainTitle,
     xlab="True vital rate - SURVIVAL",ylab="Estimated vital rate - SURVIVAL",ylim=c(minAx,maxAx), xlim=c(minAx,maxAx))
abline(a=0, b=1)
legend(legPos, "MCMC no-missing data", col=alpha("pink",0.5),pch=19, cex=1,horiz=FALSE, bty="y")
#points(output$SurvRate_True, output$SurvRate_JAGS, col=alpha("purple",0.3))
legend("bottomright", c("GLMM","MCMC missing","MCMC no-missing"), col=c(alpha("grey40",0.5),alpha("purple",0.5),alpha("pink",0.5)), 
       pch=19, cex=1.25,horiz=FALSE, bty="y",seg.len=1)
## -----------------------------------------------------------------------------------
     


     
     
     
     
     
## PLOT GLMM VS JAGS PARAMETER ESTIMATES -----------------------------------------
rownames(summ.modMs1)[c(26:33,36:41)]
     
## Get upper and lower 95% limits on JAGS estimates
     
## Missing data
low95ms.1 <- summ.modMs1[c(26:33,36:41),1]
low95ms.2 <- summ.modMs2[c(26:33,36:41),1]
low95ms.3 <- summ.modMs3[c(26:33,36:41),1]
low95ms.4 <- summ.modMs4[c(26:33,36:41),1]
low95ms.5 <- summ.modMs5[c(26:33,36:41),1]
low95ms.6 <- summ.modMs6[c(26:33,36:41),1]
low95ms.7 <- summ.modMs7[c(26:33,36:41),1]
low95ms.8 <- summ.modMs8[c(26:33,36:41),1]
low95ms.9 <- summ.modMs9[c(26:33,36:41),1]
low95ms.10 <- summ.modMs10[c(26:33,36:41),1]
     
up95ms.1 <- summ.modMs1[c(26:33,36:41),3]
up95ms.2 <- summ.modMs2[c(26:33,36:41),3]
up95ms.3 <- summ.modMs3[c(26:33,36:41),3]
up95ms.4 <- summ.modMs4[c(26:33,36:41),3]
up95ms.5 <- summ.modMs5[c(26:33,36:41),3]
up95ms.6 <- summ.modMs6[c(26:33,36:41),3]
up95ms.7 <- summ.modMs7[c(26:33,36:41),3]
up95ms.8 <- summ.modMs8[c(26:33,36:41),3]
up95ms.9 <- summ.modMs9[c(26:33,36:41),3]
up95ms.10 <- summ.modMs10[c(26:33,36:41),3]

lowCombMs.grwth <- as.data.frame(cbind(low95ms.1[1:8], low95ms.2[1:8],low95ms.3[1:8],
                                      low95ms.4[1:8],low95ms.5[1:8],low95ms.6[1:8],
                                      low95ms.7[1:8],low95ms.8[1:8],low95ms.9[1:8],low95ms.10[1:8]))

upCombMs.grwth <- as.data.frame(cbind(up95ms.1[1:8], up95ms.2[1:8],up95ms.3[1:8],
                                       up95ms.4[1:8],up95ms.5[1:8],up95ms.6[1:8],
                                       up95ms.7[1:8],up95ms.8[1:8],up95ms.9[1:8],up95ms.10[1:8]))

lowCombMs.surv <- as.data.frame(cbind(low95ms.1[9:14], low95ms.2[9:14],low95ms.3[9:14],
                                     low95ms.4[9:14],low95ms.5[9:14],low95ms.6[9:14],
                                     low95ms.7[9:14],low95ms.8[9:14],low95ms.9[9:14],low95ms.10[9:14])) 

upCombMs.surv <- as.data.frame(cbind(up95ms.1[9:14], up95ms.2[9:14],up95ms.3[9:14],
                                     up95ms.4[9:14],up95ms.5[9:14],up95ms.6[9:14],
                                     up95ms.7[9:14],up95ms.8[9:14],up95ms.9[9:14],up95ms.10[9:14])) 




## No-missing data
low95no.1 <- summ.modNo1[c(26:33,36:41),1]
low95no.2 <- summ.modNo2[c(26:33,36:41),1]
low95no.3 <- summ.modNo3[c(26:33,36:41),1]
low95no.4 <- summ.modNo4[c(26:33,36:41),1]
low95no.5 <- summ.modNo5[c(26:33,36:41),1]
low95no.6 <- summ.modNo6[c(26:33,36:41),1]
low95no.7 <- summ.modNo7[c(26:33,36:41),1]
low95no.8 <- summ.modNo8[c(26:33,36:41),1]
low95no.9 <- summ.modNo9[c(26:33,36:41),1]
low95no.10 <- summ.modNo10[c(26:33,36:41),1]

up95no.1 <- summ.modNo1[c(26:33,36:41),3]
up95no.2 <- summ.modNo2[c(26:33,36:41),3]
up95no.3 <- summ.modNo3[c(26:33,36:41),3]
up95no.4 <- summ.modNo4[c(26:33,36:41),3]
up95no.5 <- summ.modNo5[c(26:33,36:41),3]
up95no.6 <- summ.modNo6[c(26:33,36:41),3]
up95no.7 <- summ.modNo7[c(26:33,36:41),3]
up95no.8 <- summ.modNo8[c(26:33,36:41),3]
up95no.9 <- summ.modNo9[c(26:33,36:41),3]
up95no.10 <- summ.modNo10[c(26:33,36:41),3]

lowCombNo.grwth <- as.data.frame(cbind(low95no.1[1:8], low95no.2[1:8],low95no.3[1:8],
                                       low95no.4[1:8],low95no.5[1:8],low95no.6[1:8],
                                       low95no.7[1:8],low95no.8[1:8],low95no.9[1:8],low95no.10[1:8]))

upCombNo.grwth <- as.data.frame(cbind(up95no.1[1:8], up95no.2[1:8],up95no.3[1:8],
                                      up95no.4[1:8],up95no.5[1:8],up95no.6[1:8],
                                      up95no.7[1:8],up95no.8[1:8],up95no.9[1:8],up95no.10[1:8]))

lowCombNo.surv <- as.data.frame(cbind(low95no.1[9:14], low95no.2[9:14],low95no.3[9:14],
                                      low95no.4[9:14],low95no.5[9:14],low95no.6[9:14],
                                      low95no.7[9:14],low95no.8[9:14],low95no.9[9:14],low95no.10[9:14])) 

upCombNo.surv <- as.data.frame(cbind(up95no.1[9:14], up95no.2[9:14],up95no.3[9:14],
                                     up95no.4[9:14],up95no.5[9:14],up95no.6[9:14],
                                     up95no.7[9:14],up95no.8[9:14],up95no.9[9:14],up95no.10[9:14])) 
## ------------------------------------------------------------------------------------------------
     
     
     




## GLMM vs JAGS missing parameter estimates ---------------------------------------------------------------
## Plot Growth GLMM vs JAGS missing with GLMM SE and JAGS 95% quantiles 
par(mfrow=c(4,4), mar=c(4,4,2.3,1))  #bottom, left, top and right 
par(pty="s")
     
     for (nn in 1:nrow(paramsMM.grwth)) {
       
       maxLim <- max(c((as.numeric(paramsMM.grwth[nn,1:10]) + as.numeric(seMM.grwth[nn,1:10])),
                       as.numeric(upCombMs.grwth[nn,1:10]),
                       as.numeric(as.character(medParams.realDatGrwth$realData[nn])))) 
       minLim <- min(c((as.numeric(paramsMM.grwth[nn,1:10]) - as.numeric(seMM.grwth[nn,1:10])),
                       as.numeric(lowCombMs.grwth[nn,1:10]),
                       as.numeric(as.character(medParams.realDatGrwth$realData[nn])))) 
       
       plot(as.numeric(paramsMM.grwth[nn,1:10]), as.numeric(medCombMs.grwth[nn,1:10]),
            ylim=c(minLim,maxLim), xlim=c(minLim,maxLim), xlab="GLMM no-missing data", 
            ylab="MCMC missing data", main=paramsMM.grwth[nn,11], pch=19)
       abline(h=as.numeric(as.character(medParams.realDatGrwth$realData[nn])),
             v=as.numeric(as.character(medParams.realDatGrwth$realData[nn])), col="grey60")
       #abline(a=0, b=1)
       plotCI(x=as.numeric(paramsMM.grwth[nn,1:10]), y=as.numeric(medCombMs.grwth[nn,1:10]), 
              uiw=as.numeric(seMM.grwth[nn,1:10]), err="x", add=T, sfrac=0)
       plotCI(x=as.numeric(paramsMM.grwth[nn,1:10]), y=as.numeric(medCombMs.grwth[nn,1:10]), 
              liw=as.numeric(medCombMs.grwth[nn,1:10])-as.numeric(lowCombMs.grwth[nn,1:10]),
              uiw=as.numeric(upCombMs.grwth[nn,1:10])-as.numeric(medCombMs.grwth[nn,1:10]), 
              err="y", add=T, sfrac=0)
      points(as.numeric(as.character(medParams.realDatGrwth$realData[nn])), 
              as.numeric(as.character(medParams.realDatGrwth$realData[nn])), pch=8,
              cex=1.75, col="red")
     }
     

## Plot Survival GLMM vs JAGS missing with GLMM SE and JAGS 95% quantiles 
     for (nn in 1:nrow(paramsMM.surv)) {
       
       maxLim <- max(c((as.numeric(paramsMM.surv[nn,1:10]) + as.numeric(seMM.surv[nn,1:10])),
                       as.numeric(upCombMs.surv[nn,1:10]), 
                       as.numeric(as.character(medParams.realDatSurv$realData[nn])))) 
       minLim <- min(c((as.numeric(paramsMM.surv[nn,1:10]) - as.numeric(seMM.surv[nn,1:10])),
                       as.numeric(lowCombMs.surv[nn,1:10]),
                       as.numeric(as.character(medParams.realDatSurv$realData[nn])))) 
       
       plot(as.numeric(paramsMM.surv[nn,1:10]), as.numeric(medCombMs.surv[nn,1:10]),
            ylim=c(minLim,maxLim), xlim=c(minLim,maxLim), xlab="GLMM no-missing data", 
            ylab="MCMC missing data", main=paramsMM.surv[nn,11], pch=19)
       abline(h=as.numeric(as.character(medParams.realDatSurv$realData[nn])),
             v=as.numeric(as.character(medParams.realDatSurv$realData[nn])), col="grey60")
       plotCI(x=as.numeric(paramsMM.surv[nn,1:10]), y=as.numeric(medCombMs.surv[nn,1:10]), 
              uiw=as.numeric(seMM.surv[nn,1:10]), err="x", add=T, sfrac=0)
       plotCI(x=as.numeric(paramsMM.surv[nn,1:10]), y=as.numeric(medCombMs.surv[nn,1:10]), 
              liw=as.numeric(medCombMs.surv[nn,1:10])-as.numeric(lowCombMs.surv[nn,1:10]),
              uiw=as.numeric(upCombMs.surv[nn,1:10])-as.numeric(medCombMs.surv[nn,1:10]), 
              err="y", add=T, sfrac=0)
       points(as.numeric(as.character(medParams.realDatSurv$realData[nn])), 
              as.numeric(as.character(medParams.realDatSurv$realData[nn])), pch=8,
             cex=1.75, col="red")
     }

plot.new()
legend("center", name, bty="n",cex=1.2)
## ------------------------------------------------------------------------------------------------
     
     
     
     
     
     

## FOR 20 YR OBSERVED LH DATA SETS WHERE WE RAN JAGS ON NO-MISSING DATA 
## JAGS parameter estimates MISSING data VS NO-MISSING data ---------------------------------------
## Plot Growth JAGS with 95% quantiles 
par(mfrow=c(4,4), mar=c(4,4,2.3,1))  #bottom, left, top and right 
par(pty="s")
     
     for (nn in 1:nrow(medCombNo.grwth)) {
       
       maxLim <- max(c(as.numeric(upCombNo.grwth[nn,1:10])),
                       as.numeric(upCombMs.grwth[nn,1:10]),
                       as.numeric(as.character(medParams.realDatGrwth$realData[nn]))) 
       minLim <- min(c(as.numeric(lowCombNo.grwth[nn,1:10])),
                       as.numeric(lowCombMs.grwth[nn,1:10]),
                       as.numeric(as.character(medParams.realDatGrwth$realData[nn])))
       
       plot(as.numeric(medCombNo.grwth[nn,1:10]), as.numeric(medCombMs.grwth[nn,1:10]),
            ylim=c(minLim,maxLim), xlim=c(minLim,maxLim), xlab="MCMC no-missing data", 
            ylab="MCMC missing data", main=names.paramTitles[nn], pch=19)
       abline(h=as.numeric(as.character(medParams.realDatGrwth$realData[nn])),
              v=as.numeric(as.character(medParams.realDatGrwth$realData[nn])), col="grey60")       
       plotCI(x=as.numeric(medCombNo.grwth[nn,1:10]), y=as.numeric(medCombMs.grwth[nn,1:10]), 
              liw=as.numeric(medCombNo.grwth[nn,1:10])-as.numeric(lowCombNo.grwth[nn,1:10]),
              uiw=as.numeric(upCombNo.grwth[nn,1:10])-as.numeric(medCombNo.grwth[nn,1:10]), 
              err="x", add=T, sfrac=0)
       plotCI(x=as.numeric(medCombNo.grwth[nn,1:10]), y=as.numeric(medCombMs.grwth[nn,1:10]), 
              liw=as.numeric(medCombMs.grwth[nn,1:10])-as.numeric(lowCombMs.grwth[nn,1:10]),
              uiw=as.numeric(upCombMs.grwth[nn,1:10])-as.numeric(medCombMs.grwth[nn,1:10]), 
              err="y", add=T, sfrac=0)
       points(as.numeric(as.character(medParams.realDatGrwth$realData[nn])), 
              as.numeric(as.character(medParams.realDatGrwth$realData[nn])), pch=8,
              cex=1.75, col="red")
     }
     
     
## Plot Survival JAGS with JAGS 95% quantiles 
     for (nn in 1:nrow(medCombNo.surv)) {
       
       maxLim <- max(c(as.numeric(upCombNo.surv[nn,1:10])),
                     as.numeric(upCombMs.surv[nn,1:10]),
                     as.numeric(as.character(medParams.realDatSurv$realData[nn]))) 
       minLim <- min(c(as.numeric(lowCombNo.surv[nn,1:10])),
                     as.numeric(lowCombMs.surv[nn,1:10]),
                     as.numeric(as.character(medParams.realDatSurv$realData[nn])))
       
       plot(as.numeric(medCombNo.surv[nn,1:10]), as.numeric(medCombMs.surv[nn,1:10]),
            ylim=c(minLim,maxLim), xlim=c(minLim,maxLim), xlab="MCMC no-missing data", 
            ylab="MCMC missing data", main=names.paramTitles[nn+8], pch=19)
       abline(h=as.numeric(as.character(medParams.realDatSurv$realData[nn])),
              v=as.numeric(as.character(medParams.realDatSurv$realData[nn])), col="grey60")       
       plotCI(x=as.numeric(medCombNo.surv[nn,1:10]), y=as.numeric(medCombMs.surv[nn,1:10]), 
              liw=as.numeric(medCombNo.surv[nn,1:10])-as.numeric(lowCombNo.surv[nn,1:10]),
              uiw=as.numeric(upCombNo.surv[nn,1:10])-as.numeric(medCombNo.surv[nn,1:10]), 
              err="x", add=T, sfrac=0)
       plotCI(x=as.numeric(medCombNo.surv[nn,1:10]), y=as.numeric(medCombMs.surv[nn,1:10]), 
              liw=as.numeric(medCombMs.surv[nn,1:10])-as.numeric(lowCombMs.surv[nn,1:10]),
              uiw=as.numeric(upCombMs.surv[nn,1:10])-as.numeric(medCombMs.surv[nn,1:10]), 
              err="y", add=T, sfrac=0)
       points(as.numeric(as.character(medParams.realDatSurv$realData[nn])), 
              as.numeric(as.character(medParams.realDatSurv$realData[nn])), pch=8,
              cex=1.75, col="red")
     }
plot.new()
legend("center", name, bty="n",cex=1.2)
## ------------------------------------------------------------------------------------------------
     
     

## GLMM vs JAGS no-missing parameter estimates ------------------------------------------------------
## Plot Growth GLMM vs JAGS no-missing with GLMM SE and JAGS 95% quantiles 
par(mfrow=c(4,4), mar=c(4,4,2.3,1))  #bottom, left, top and right 
par(pty="s")
     
     for (nn in 1:nrow(paramsMM.grwth)) {
       
       maxLim <- max(c((as.numeric(paramsMM.grwth[nn,1:10]) + as.numeric(seMM.grwth[nn,1:10])),
                       as.numeric(upCombNo.grwth[nn,1:10]))) 
       minLim <- min(c((as.numeric(paramsMM.grwth[nn,1:10]) - as.numeric(seMM.grwth[nn,1:10])),
                       as.numeric(lowCombNo.grwth[nn,1:10]))) 
       
       plot(as.numeric(paramsMM.grwth[nn,1:10]), as.numeric(medCombNo.grwth[nn,1:10]),
            ylim=c(minLim,maxLim), xlim=c(minLim,maxLim), xlab="GLMM no-missing data", 
            ylab="MCMC no-missing data", main=paramsMM.grwth[nn,11], pch=19)
       abline(h=as.numeric(as.character(medParams.realDatGrwth$realData[nn])),
              v=as.numeric(as.character(medParams.realDatGrwth$realData[nn])), col="grey60")       
       plotCI(x=as.numeric(paramsMM.grwth[nn,1:10]), y=as.numeric(medCombNo.grwth[nn,1:10]), 
              uiw=as.numeric(seMM.grwth[nn,1:10]), err="x", add=T, sfrac=0)
       plotCI(x=as.numeric(paramsMM.grwth[nn,1:10]), y=as.numeric(medCombNo.grwth[nn,1:10]), 
              liw=as.numeric(medCombNo.grwth[nn,1:10])-as.numeric(lowCombNo.grwth[nn,1:10]),
              uiw=as.numeric(upCombNo.grwth[nn,1:10])-as.numeric(medCombNo.grwth[nn,1:10]), 
              err="y", add=T, sfrac=0)
       points(as.numeric(as.character(medParams.realDatGrwth$realData[nn])), 
              as.numeric(as.character(medParams.realDatGrwth$realData[nn])), pch=8,
              cex=1.75, col="red")
     }
     
     
## Plot Survival GLMM vs JAGS no-missing with GLMM SE and JAGS 95% quantiles 
     for (nn in 1:nrow(paramsMM.surv)) {
       
       maxLim <- max(c((as.numeric(paramsMM.surv[nn,1:10]) + as.numeric(seMM.surv[nn,1:10])),
                       as.numeric(upCombNo.surv[nn,1:10]))) 
       minLim <- min(c((as.numeric(paramsMM.surv[nn,1:10]) - as.numeric(seMM.surv[nn,1:10])),
                       as.numeric(lowCombNo.surv[nn,1:10]))) 
       
       plot(as.numeric(paramsMM.surv[nn,1:10]), as.numeric(medCombNo.surv[nn,1:10]),
            ylim=c(minLim,maxLim), xlim=c(minLim,maxLim), xlab="GLMM no-missing data", 
            ylab="MCMC no-missing data", main=paramsMM.surv[nn,11], pch=19)
       abline(h=as.numeric(as.character(medParams.realDatSurv$realData[nn])),
              v=as.numeric(as.character(medParams.realDatSurv$realData[nn])), col="grey60")            
       plotCI(x=as.numeric(paramsMM.surv[nn,1:10]), y=as.numeric(medCombNo.surv[nn,1:10]), 
              uiw=as.numeric(seMM.surv[nn,1:10]), err="x", add=T, sfrac=0)
       plotCI(x=as.numeric(paramsMM.surv[nn,1:10]), y=as.numeric(medCombNo.surv[nn,1:10]), 
              liw=as.numeric(medCombNo.surv[nn,1:10])-as.numeric(lowCombNo.surv[nn,1:10]),
              uiw=as.numeric(upCombNo.surv[nn,1:10])-as.numeric(medCombNo.surv[nn,1:10]), 
              err="y", add=T, sfrac=0)
       points(as.numeric(as.character(medParams.realDatSurv$realData[nn])), 
              as.numeric(as.character(medParams.realDatSurv$realData[nn])), pch=8,
              cex=1.75, col="red")
     }
## --------------------------------------------------------------------------------------------     
      
     