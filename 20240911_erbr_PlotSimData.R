## April Goebl & Dan Doak
## Script modified 2024-09-11
## Collaboration with CU Boulder and Denver Botanic Gardens on Eriogonum brandegeii modeling 
## Plot parameter results from modeling of simulated demographic data 
## to show JAGS model performance with missing data 


rm(list=ls())
graphics.off()


## LOAD PACKAGES AND FUNCTIONS --------------------------------------------------------------------
library(lme4)
library(ggplot2)
library(dplyr)
library(resample)
library(gplots)
library(matrixStats)
library(stringr)
library(plotrix)
## ------------------------------------------------------------------------------------------------




## ASSIGN NAME VARIABLE FOR DESIRED DATASETS 
#date <- as.character("20240911")
name <- as.character("SimDat20yr.")
#dats <- read.csv(file=paste(date,"_erbr_", name, dd, ".Format4JAGS", ".csv", sep=""), header=TRUE)
## ------------------------------------------------------------------------------------------------



## LOAD MODEL SUMMARY OUTPUT -------------------------------------------------
summ.mod1 <- readRDS("20240910_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat50yrMiss.1.rds")
summ.mod2 <- readRDS("20240910_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat50yrMiss.2.rds")
summ.mod3 <- readRDS("20240910_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat50yrMiss.3.rds")
summ.mod4 <- readRDS("20240910_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat50yrMiss.4.rds")
summ.mod5 <- readRDS("20240910_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat50yrMiss.5.rds")
summ.mod6 <- readRDS("20240910_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat50yrMiss.6.rds")
summ.mod7 <- readRDS("20240910_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat50yrMiss.7.rds")
summ.mod8 <- readRDS("20240910_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat50yrMiss.8.rds")
summ.mod9 <- readRDS("20240910_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat50yrMiss.9.rds")
summ.mod10 <- readRDS("20240910_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat50yrMiss.10.rds")

date <- as.character("20240906")
name
summ.mod1 <- readRDS(file=paste(date,"_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_", name, "1.rds", sep=""))
summ.mod2 <- readRDS(file=paste(date,"_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_", name, "2.rds", sep=""))
summ.mod3 <- readRDS(file=paste(date,"_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_", name, "3.rds", sep=""))
summ.mod4 <- readRDS(file=paste(date,"_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_", name, "4.rds", sep=""))
summ.mod5 <- readRDS(file=paste(date,"_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_", name, "5.rds", sep=""))
summ.mod6 <- readRDS(file=paste(date,"_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_", name, "6.rds", sep=""))
summ.mod7 <- readRDS(file=paste(date,"_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_", name, "7.rds", sep=""))
summ.mod8 <- readRDS(file=paste(date,"_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_", name, "8.rds", sep=""))
summ.mod9 <- readRDS(file=paste(date,"_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_", name, "9.rds", sep=""))
summ.mod10 <- readRDS(file=paste(date,"_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_", name, "10.rds", sep=""))
## -----------------------------------------------------------------------------------------------




## Make plot comparing median param ests of rep simulated datasets & compare GLMM estimate w/out missing data
## Show median param ests of diff datasets as points and upper and lower 95% limits (from JAGS summary)

## SUBSET OUTPUT TO JUST KEEP PARAMS OF INTEREST
#names.param <- colnames(as.matrix(summ.mod1$mcmc[1]))[26:41]
#names.param <- rownames(summ.mod1)[26:41]
#names.param <- names.param[c(2:8,12:16)] #Remove intercept and GrwthVar  
names.paramTitles <- c("Grwth Intercept","Grwth Size","Grwth Fall Temp","Grwth Summer Temp","Grwth Winter Temp",
                       "Grwth Fall Precip","Grwth Summer Precip","Grwth Winter Precip","Surv Intercept","Surv Size",
                       "Surv Winter Precip","Surv Fall Temp","Surv Summer Temp","Surv Winter Temp")

## Re-order parameter names for plotting 
#index<-c(1,5,3,7,11,9,13,2,14,6,4,8)
#names.paramOrd <- names.param[order(index)]
#names.paramTitlesOrd <- names.paramTitles[order(index)]
## -----------------------------------------------------------------------------------------------



## OBTAIN RELEVANT AND COMBINE MEDIAN PARAMETER VALUES FROM DIFF SIM DATASETS --------------------
#summ.mod1$summaries
medParams.1 <- summ.mod1[c(27:33,37:41),2]
medParams.2 <- summ.mod2[c(27:33,37:41),2]
medParams.3 <- summ.mod3[c(27:33,37:41),2]
medParams.4 <- summ.mod4[c(27:33,37:41),2]
medParams.5 <- summ.mod5[c(27:33,37:41),2]
medParams.6 <- summ.mod6[c(27:33,37:41),2]
medParams.7 <- summ.mod7[c(27:33,37:41),2]
medParams.8 <- summ.mod8[c(27:33,37:41),2]
medParams.9 <- summ.mod9[c(27:33,37:41),2]
medParams.10 <- summ.mod10[c(27:33,37:41),2]

mnParams.1 <- summ.mod1[c(27:33,37:41),4]
mnParams.2 <- summ.mod2[c(27:33,37:41),4]
mnParams.3 <- summ.mod3[c(27:33,37:41),4]
mnParams.4 <- summ.mod4[c(27:33,37:41),4]
mnParams.5 <- summ.mod5[c(27:33,37:41),4]
mnParams.6 <- summ.mod6[c(27:33,37:41),4]
mnParams.7 <- summ.mod7[c(27:33,37:41),4]
mnParams.8 <- summ.mod8[c(27:33,37:41),4]
mnParams.9 <- summ.mod9[c(27:33,37:41),4]
mnParams.10 <- summ.mod10[c(27:33,37:41),4]

low95.1 <- summ.mod1[c(27:33,37:41),1]
low95.2 <- summ.mod2[c(27:33,37:41),1]
low95.3 <- summ.mod3[c(27:33,37:41),1]
low95.4 <- summ.mod4[c(27:33,37:41),1]
low95.5 <- summ.mod5[c(27:33,37:41),1]
low95.6 <- summ.mod6[c(27:33,37:41),1]
low95.7 <- summ.mod7[c(27:33,37:41),1]
low95.8 <- summ.mod8[c(27:33,37:41),1]
low95.9 <- summ.mod9[c(27:33,37:41),1]
low95.10 <- summ.mod10[c(27:33,37:41),1]

up95.1 <- summ.mod1[c(27:33,37:41),3]
up95.2 <- summ.mod2[c(27:33,37:41),3]
up95.3 <- summ.mod3[c(27:33,37:41),3]
up95.4 <- summ.mod4[c(27:33,37:41),3]
up95.5 <- summ.mod5[c(27:33,37:41),3]
up95.6 <- summ.mod6[c(27:33,37:41),3]
up95.7 <- summ.mod7[c(27:33,37:41),3]
up95.8 <- summ.mod8[c(27:33,37:41),3]
up95.9 <- summ.mod9[c(27:33,37:41),3]
up95.10 <- summ.mod10[c(27:33,37:41),3]

#sd.1 <- summ.mod1[c(27:33,37:41),5]
#sd.2 <- summ.mod2[c(27:33,37:41),5]
#sd.3 <- summ.mod3[c(27:33,37:41),5]
#sd.4 <- summ.mod4[c(27:33,37:41),5]
#sd.5 <- summ.mod5[c(27:33,37:41),5]
#sd.6 <- summ.mod6[c(27:33,37:41),5]
#sd.7 <- summ.mod7[c(27:33,37:41),5]
#sd.8 <- summ.mod8[c(27:33,37:41),5]
#sd.9 <- summ.mod9[c(27:33,37:41),5]
#sd.10 <- summ.mod10[c(27:33,37:41),5]

## COMBINE
medComb.grwth <- as.data.frame(cbind(medParams.1[1:7], medParams.2[1:7],medParams.3[1:7],
                                     medParams.4[1:7],medParams.5[1:7],medParams.6[1:7],
                                     medParams.7[1:7],medParams.8[1:7],medParams.9[1:7],medParams.10[1:7]))
lowComb.grwth <- as.data.frame(cbind(low95.1[1:7], low95.2[1:7],low95.3[1:7],
                                     low95.4[1:7],low95.5[1:7],low95.6[1:7],
                                     low95.7[1:7],low95.8[1:7],low95.9[1:7],low95.10[1:7]))
upComb.grwth <- as.data.frame(cbind(up95.1[1:7], up95.2[1:7],up95.3[1:7],
                                     up95.4[1:7],up95.5[1:7],up95.6[1:7],
                                     up95.7[1:7],up95.8[1:7],up95.9[1:7],up95.10[1:7]))
#sdComb.grwth <- as.data.frame(cbind(sd.1[1:7], sd.2[1:7],sd.3[1:7],
#                                    sd.4[1:7],sd.5[1:7],sd.6[1:7],
#                                    sd.7[1:7],sd.8[1:7],sd.9[1:7],sd.10[1:7]))
mnComb.grwth <- as.data.frame(cbind(mnParams.1[1:7], mnParams.2[1:7],mnParams.3[1:7],
                                     mnParams.4[1:7],mnParams.5[1:7],mnParams.6[1:7],
                                     mnParams.7[1:7],mnParams.8[1:7],mnParams.9[1:7],mnParams.10[1:7]))

medComb.surv <- as.data.frame(cbind(medParams.1[8:12], medParams.2[8:12],medParams.3[8:12],
                                     medParams.4[8:12],medParams.5[8:12],medParams.6[8:12],
                                     medParams.7[8:12],medParams.8[8:12],medParams.9[8:12],medParams.10[8:12]))
lowComb.surv <- as.data.frame(cbind(low95.1[8:12], low95.2[8:12],low95.3[8:12],
                                     low95.4[8:12],low95.5[8:12],low95.6[8:12],
                                     low95.7[8:12],low95.8[8:12],low95.9[8:12],low95.10[8:12]))
upComb.surv <- as.data.frame(cbind(up95.1[8:12], up95.2[8:12],up95.3[8:12],
                                    up95.4[8:12],up95.5[8:12],up95.6[8:12],
                                    up95.7[8:12],up95.8[8:12],up95.9[8:12],up95.10[8:12]))
mnComb.surv <- as.data.frame(cbind(mnParams.1[8:12], mnParams.2[8:12],mnParams.3[8:12],
                                    mnParams.4[8:12],mnParams.5[8:12],mnParams.6[8:12],
                                    mnParams.7[8:12],mnParams.8[8:12],mnParams.9[8:12],mnParams.10[8:12]))
## ------------------------------------------------------------------------------------------------






## OBTAIN GLMM ESTIMATES -------------------------------------------------------------------------
## Loop over datasets 
date <- as.character("20240911")
#name
#name <- str_replace_all(name, ".", "")
name <- "SimDat20yrHiGrLH."

#modList.grwth <- NULL        #List variable to store all models
#modList.surv <- NULL        #List variable to store all models
paramsMM.grwth <- NULL
seMM.grwth <- NULL
paramsMM.surv <- NULL
seMM.surv <- NULL
n.datset <- 10

for (dd in 1:n.datset) {
  
  noMiss <- readRDS(file=paste(date, "_erbr_", name, "NoMiss.",dd,".4GLM",".rds", sep=""))

  ## Add t+1 climate, sz, & tag into erbr data 
  noMiss <- noMiss %>% mutate(TagNew1=lead(TagNew), RosNew1=lead(RosNew), Surv1=lead(surv))  
  
  noMiss <- noMiss[which(noMiss$TagNew == noMiss$TagNew1),]  #Remove lines with mis-matched individuals 
  
  ## Log size in all models to match JAGS; this makes ending size a linear function of starting size
  
  ## Growth (param order is same as for JAGS output)
  glmm.grwth <- glmer.nb(RosNew1 ~ log(RosNew) + TempFall + TempSummer + TempWinter + 
                           PptFall + PptSummer + PptWinter + (1|TransectNew), data=noMiss)
  
  #glmm.grwth <- glmer.nb(RosNew1 ~ log(RosNew) + PptFall + PptWinter + PptSummer + 
  #                         TempFall + TempWinter + TempSummer + (1|TransectNew), data=noMiss)
  
  ## Survival  (param order is same as for JAGS output)
  glmm.surv <- glmer(Surv1 ~ log(RosNew) + PptWinter + TempFall + TempSummer + 
                       TempWinter + (1|TransectNew), family=binomial(link='logit'), data=noMiss)
  #glmm.surv <- glmer(Surv1 ~ log(RosNew) + PptWinter + TempFall + (TempWinter) + 
  #                     (TempSummer) + (1|TransectNew), family=binomial(link='logit'), data=noMiss)
  
  ## Save model object in list
  #modList.grwth[[length(modList.grwth) + 1]] <- glmm.grwth
  #modList.surv[[length(modList.surv) + 1]] <- glmm.surv
  
  
  ## Extract parameter estimates, including intercept, and SEs from GLMMs
  paramsMM.grwthTmp <- as.data.frame(summary(glmm.grwth)$coefficients[1:nrow(summary(glmm.grwth)$coefficients),1])
  seMM.grwthTmp <- as.data.frame(summary(glmm.grwth)$coefficients[1:nrow(summary(glmm.grwth)$coefficients),2])
  
  colnames(paramsMM.grwthTmp) <- paste("GLMM.",dd,sep="")
  colnames(seMM.grwthTmp) <- paste("SE.",dd,sep="")#),"ParamTitle")
  
  paramsMM.grwth <- as.data.frame(c(paramsMM.grwth, paramsMM.grwthTmp))
  seMM.grwth <- as.data.frame(c(seMM.grwth, seMM.grwthTmp))
  
  
  paramsMM.survTmp <- as.data.frame(summary(glmm.surv)$coefficients[1:nrow(summary(glmm.surv)$coefficients),1])
  seMM.survTmp <- as.data.frame(summary(glmm.surv)$coefficients[1:nrow(summary(glmm.surv)$coefficients),2])
  
  colnames(paramsMM.survTmp) <- paste("GLMM.",dd,sep="")
  colnames(seMM.survTmp) <- paste("SE.",dd,sep="")
  paramsMM.surv <- as.data.frame(c(paramsMM.surv, paramsMM.survTmp))
  seMM.surv <- as.data.frame(c(seMM.surv, seMM.survTmp))
  
}

paramsMM.grwthTmp
paramsMM.grwth$ParamTitle <- names.paramTitles[1:8]
seMM.grwth$ParamTitle <- names.paramTitles[1:8]
paramsMM.survTmp
paramsMM.surv$ParamTitle <- names.paramTitles[9:14]
seMM.surv$ParamTitle <- names.paramTitles[9:14]

#paramsMM.grwth$ParamTitle <- c("Grwth Intercept","Grwth Size","Grwth Fall Precip","Grwth Winter Precip","Grwth Summer Precip",
#                               "Grwth Fall Temp","Grwth Winter Temp","Grwth Summer Temp")
#seMM.grwth$ParamTitle <- c("Grwth Intercept","Grwth Size","Grwth Fall Precip","Grwth Winter Precip","Grwth Summer Precip",
#                           "Grwth Fall Temp","Grwth Winter Temp","Grwth Summer Temp")
#paramsMM.surv$ParamTitle <- c("Surv Size","Surv Winter Precip",
#                              "Surv Fall Temp","Surv Winter Temp","Surv Summer Temp")
#seMM.surv$ParamTitle <- c("Surv Size","Surv Winter Precip",
#                          "Surv Fall Temp","Surv Winter Temp","Surv Summer Temp")


## Load saved GLMM results
#paramsMM.grwth <- readRDS(file=paste("20240908_erbr_paramMMgrwth_SimDat20yr",".rds", sep=""))
#seMM.grwth <- readRDS(file=paste("20240908_erbr_seMMgrwth_SimDat20yr",".rds", sep=""))
#paramsMM.surv <- readRDS(file=paste("20240908_erbr_paramMMsurv_SimDat20yr",".rds", sep=""))
#seMM.surv <- readRDS(file=paste("20240908_erbr_seMMsurv_SimDat20yr",".rds", sep=""))

## Re-order parameter names for plotting 
#paramsMM.grwthOrd <- paramsMM.grwth[match(names.paramTitles[1:7], paramsMM.grwth$ParamTitle),]
#seMM.grwthOrd <- seMM.grwth[match(names.paramTitles[1:7], seMM.grwth$ParamTitle),]
#paramsMM.survOrd <- paramsMM.surv[match(names.paramTitles[8:12], paramsMM.surv$ParamTitle),]
#seMM.survOrd <- seMM.surv[match(names.paramTitles[8:12], seMM.surv$ParamTitle),]


## Save GLMM param and se results 
date <- as.character("20240912")
name
name <- "SimDat20yrHiGr"
saveRDS(paramsMM.grwth, file=paste(date, "_erbr_paramMMgrwthWint_", name, ".rds", sep=""))
saveRDS(seMM.grwth, file=paste(date, "_erbr_seMMgrwthWint_", name, ".rds", sep=""))
saveRDS(paramsMM.surv, file=paste(date, "_erbr_paramMMsurvWint_", name, ".rds", sep=""))
saveRDS(seMM.surv, file=paste(date, "_erbr_seMMsurvWint_", name, ".rds", sep=""))

## Save GLMM model objects
#saveRDS(glmm.grwth, file=paste(date, "_erbr_GLMMgrwth_", name, ".rds", sep=""))
#saveRDS(glmm.surv, file=paste(date, "_erbr_GLMMsurv_", name, ".rds", sep=""))


#name <- "SimDat20yrMedGrLH."
## Load saved GLMM results
#paramsMM.grwthOrd <- readRDS(file=paste("20240911_erbr_paramMMgrwthOrd_", name, ".rds", sep=""))
#seMM.grwthOrd <- readRDS(file=paste("20240911_erbr_seMMgrwthOrd_", name,".rds", sep=""))
#paramsMM.survOrd <- readRDS(file=paste("20240911_erbr_paramMMsurvOrd_", name, ".rds", sep=""))
#seMM.survOrd <- readRDS(file=paste("20240911_erbr_seMMsurvOrd_", name, ".rds", sep=""))
## ------------------------------------------------------------------------------------------------------------






## ------------------------------------------------------------------
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
## -----------------------------------------------------




## PLOT RESULTS -----------------------------------------------------------------------------------------------

## Growth with JAGS SD
#par(mfrow=c(3,3), mar=c(3.9,1.7,2.3,1.5))  #bottom, left, top and right 
#par(pty="s")
#yMin <- min(as.numeric(medComb.grwth[nn,1:10])-as.numeric(lowComb.grwth[nn,1:10])) * 0.1
#yMax <- max(as.numeric(upComb.grwth[nn,1:10])-as.numeric(medComb.grwth[nn,1:10])) * 0.1
 
#for (nn in 1:nrow(paramsMM.grwthOrd)) {
  
  #maxLim <- max(c((as.numeric(paramsMM.grwthOrd[nn,1:10]) + as.numeric(seMM.grwth[nn,1:10])),
  #                as.numeric(upComb.grwth[nn,1:10]))) 
  #minLim <- min(c((as.numeric(paramsMM.grwthOrd[nn,1:10]) - as.numeric(seMM.grwth[nn,1:10])),
  #                as.numeric(lowComb.grwth[nn,1:10]))) 
 # maxLim <- max(c((as.numeric(paramsMM.grwthOrd[nn,1:10]) + as.numeric(seMM.grwth[nn,1:10])),
  #                (as.numeric(medComb.grwth[nn,1:10]) + as.numeric(sdComb.grwth[nn,1:10])))) 
  #minLim <- min(c((as.numeric(paramsMM.grwthOrd[nn,1:10]) - as.numeric(seMM.grwth[nn,1:10])),
   #               (as.numeric(medComb.grwth[nn,1:10]) - as.numeric(sdComb.grwth[nn,1:10])))) 
  
  #plot(as.numeric(paramsMM.grwthOrd[nn,1:10]), as.numeric(medComb.grwth[nn,1:10]),
  #     ylim=c(minLim,maxLim), xlim=c(minLim,maxLim), xlab="GLMM estimate", 
  #     ylab="JAGS estimate", main=paramsMM.grwthOrd[nn,11], pch=19)
  #abline(a=0, b=1)
  #plotCI(x=as.numeric(paramsMM.grwthOrd[nn,1:10]), y=as.numeric(medComb.grwth[nn,1:10]), 
  #       uiw=as.numeric(seMM.grwth[nn,1:10]), err="x", add=T, sfrac=0)
  #plotCI(x=as.numeric(paramsMM.grwthOrd[nn,1:10]), y=as.numeric(medComb.grwth[nn,1:10]), 
  #       liw=as.numeric(medComb.grwth[nn,1:10])-as.numeric(lowComb.grwth[nn,1:10]),
  #       uiw=as.numeric(upComb.grwth[nn,1:10])-as.numeric(medComb.grwth[nn,1:10]), 
  #       err="y", add=T, sfrac=0)
  #plotCI(x=as.numeric(paramsMM.grwthOrd[nn,1:10]), y=as.numeric(medComb.grwth[nn,1:10]), 
  #       uiw=as.numeric(seMM.grwthOrd[nn,1:10]), err="x", add=T, sfrac=0)
  #plotCI(x=as.numeric(paramsMM.grwthOrd[nn,1:10]), y=as.numeric(medComb.grwth[nn,1:10]), 
  #       uiw=as.numeric(sdComb.grwthOrd[nn,1:10]), err="y", add=T, sfrac=0)
#}


## Plot Growth with 95% quantiles 
par(mfrow=c(3,3), mar=c(3.9,2,2.3,1.5))  #bottom, left, top and right 
par(pty="s")

for (nn in 1:nrow(paramsMM.grwthOrd)) {
  
  maxLim <- max(c((as.numeric(paramsMM.grwthOrd[nn,1:10]) + as.numeric(seMM.grwthOrd[nn,1:10])),
                  as.numeric(upComb.grwth[nn,1:10]))) 
  minLim <- min(c((as.numeric(paramsMM.grwthOrd[nn,1:10]) - as.numeric(seMM.grwthOrd[nn,1:10])),
                  as.numeric(lowComb.grwth[nn,1:10]))) 
  
  plot(as.numeric(paramsMM.grwthOrd[nn,1:10]), as.numeric(medComb.grwth[nn,1:10]),
       ylim=c(minLim,maxLim), xlim=c(minLim,maxLim), xlab="GLMM estimate", 
       ylab="JAGS estimate", main=paramsMM.grwthOrd[nn,11], pch=19)
  abline(a=0, b=1)
  plotCI(x=as.numeric(paramsMM.grwthOrd[nn,1:10]), y=as.numeric(medComb.grwth[nn,1:10]), 
         uiw=as.numeric(seMM.grwthOrd[nn,1:10]), err="x", add=T, sfrac=0)
  plotCI(x=as.numeric(paramsMM.grwthOrd[nn,1:10]), y=as.numeric(medComb.grwth[nn,1:10]), 
         liw=as.numeric(medComb.grwth[nn,1:10])-as.numeric(lowComb.grwth[nn,1:10]),
         uiw=as.numeric(upComb.grwth[nn,1:10])-as.numeric(medComb.grwth[nn,1:10]), 
         err="y", add=T, sfrac=0)
  points(as.numeric(as.character(medParams.realDatGrwth$realData[nn])), 
         as.numeric(as.character(medParams.realDatGrwth$realData[nn])), pch=8,
         cex=1.75, col="red")
}



## Plot MEAN Growth with 95% quantiles 
par(mfrow=c(3,3), mar=c(3.9,2,2.3,1.5))  #bottom, left, top and right 
par(pty="s")

for (nn in 1:nrow(paramsMM.grwthOrd)) {
  
  maxLim <- max(c((as.numeric(paramsMM.grwthOrd[nn,1:10]) + as.numeric(seMM.grwthOrd[nn,1:10])),
                  as.numeric(upComb.grwth[nn,1:10]))) 
  minLim <- min(c((as.numeric(paramsMM.grwthOrd[nn,1:10]) - as.numeric(seMM.grwthOrd[nn,1:10])),
                  as.numeric(lowComb.grwth[nn,1:10]))) 
  
  plot(as.numeric(paramsMM.grwthOrd[nn,1:10]), as.numeric(mnComb.grwth[nn,1:10]),
       ylim=c(minLim,maxLim), xlim=c(minLim,maxLim), xlab="GLMM estimate", 
       ylab="JAGS estimate", main=paramsMM.grwthOrd[nn,11], pch=19)
  abline(a=0, b=1)
  plotCI(x=as.numeric(paramsMM.grwthOrd[nn,1:10]), y=as.numeric(mnComb.grwth[nn,1:10]), 
         uiw=as.numeric(seMM.grwthOrd[nn,1:10]), err="x", add=T, sfrac=0)
  plotCI(x=as.numeric(paramsMM.grwthOrd[nn,1:10]), y=as.numeric(mnComb.grwth[nn,1:10]), 
         liw=as.numeric(mnComb.grwth[nn,1:10])-as.numeric(lowComb.grwth[nn,1:10]),
         uiw=as.numeric(upComb.grwth[nn,1:10])-as.numeric(mnComb.grwth[nn,1:10]), 
         err="y", add=T, sfrac=0)
  points(as.numeric(as.character(medParams.realDatGrwth$realData[nn])), 
         as.numeric(as.character(medParams.realDatGrwth$realData[nn])), pch=8,
         cex=1.75, col="red")
}




## Plot Survival with 95% quantiles from JAGS and SE from GLMM
par(mfrow=c(2,3), mar=c(3.9,3.9,2.3,2))  #bottom, left, top and right 
par(pty="s")

for (nn in 1:nrow(paramsMM.survOrd)) {
  
  maxLim <- max(c((as.numeric(paramsMM.survOrd[nn,1:10]) + as.numeric(seMM.survOrd[nn,1:10])),
                  as.numeric(upComb.surv[nn,1:10]))) 
  minLim <- min(c((as.numeric(paramsMM.survOrd[nn,1:10]) - as.numeric(seMM.survOrd[nn,1:10])),
                  as.numeric(lowComb.surv[nn,1:10]))) 
  
  plot(as.numeric(paramsMM.survOrd[nn,1:10]), as.numeric(medComb.surv[nn,1:10]),
       ylim=c(minLim,maxLim), xlim=c(minLim,maxLim), xlab="GLMM estimate", 
       ylab="JAGS estimate", main=paramsMM.survOrd[nn,11], pch=19)
  abline(a=0, b=1)
  plotCI(x=as.numeric(paramsMM.survOrd[nn,1:10]), y=as.numeric(medComb.surv[nn,1:10]), 
         uiw=as.numeric(seMM.survOrd[nn,1:10]), err="x", add=T, sfrac=0)
  plotCI(x=as.numeric(paramsMM.survOrd[nn,1:10]), y=as.numeric(medComb.surv[nn,1:10]), 
         liw=as.numeric(medComb.surv[nn,1:10])-as.numeric(lowComb.surv[nn,1:10]),
         uiw=as.numeric(upComb.surv[nn,1:10])-as.numeric(medComb.surv[nn,1:10]), 
         err="y", add=T, sfrac=0)
  points(as.numeric(as.character(medParams.realDatSurv$realData[nn])), 
         as.numeric(as.character(medParams.realDatSurv$realData[nn])), pch=8,
         cex=1.75, col="red")
}
## -----------------------------------------------------------------------------




#plot(as.numeric(paramsMM.grwthOrd[2,1:10]), as.numeric(medComb.grwth[2,1:10]),
#     ylim=c(-0.25,0.4), xlim=c(-0.25,0.4), xlab="GLMM estimate", 
#     ylab="JAGS estimate", main=paramsMM.grwthOrd[2,11], pch=19)
#abline(a=0, b=1)
#plotCI(x=as.numeric(paramsMM.grwthOrd[2,1:10]), y=as.numeric(medComb.grwth[2,1:10]), 
#       uiw=as.numeric(seMM.grwth[2,1:10]), err="x", add=T, sfrac=0.005)
#plotCI(x=as.numeric(paramsMM.grwthOrd[2,1:10]), y=as.numeric(medComb.grwth[2,1:10]), 
#       liw=as.numeric(medComb.grwth[2,1:10])-as.numeric(lowComb.grwth[2,1:10]),
#       uiw=as.numeric(upComb.grwth[2,1:10])-as.numeric(medComb.grwth[2,1:10]), 
#       err="y", add=T, sfrac=0.001)


