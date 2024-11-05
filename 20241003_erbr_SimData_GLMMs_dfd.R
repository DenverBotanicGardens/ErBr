## April Goebl & Dan Doak
## Script modified 2024-10-01
## Collaboration with CU Boulder and Denver Botanic Gardens on Eriogonum brandegeii modeling 
## Run GLMM models on simulated data 


rm(list=ls())
graphics.off()

#install.packages("Matrix") 
#install.packages("TMB", type = "source")
#install.packages("glmmTMB", type = "source")
library(MASS)
library(dplyr)
library(matrixStats)
library(corrplot)
library(car)
library(coda)
library(lme4)
library(resample)
library(gplots)
library(stringr)
library(plotrix)
library(glmmTMB) # installed from tar for not the last version, but 1.1.9
library(TMB)
## ------------------------------------------------------------------------------------------------



setwd("C:/Users/Dan Doak/Desktop/Students/April/eriogonum models/manu fall2024/Sept 27 debugging simulateddata/checkingglmms")



## OBTAIN GLMM ESTIMATES ----------------------------------------------------------------

## Loop over datasets 
dateGLM <- as.character("20241001")
## ASSIGN NAME VARIABLE FOR DESIRED DATASETS 
name <- as.character("SimDat20yr")


#modList.grwth <- NULL        #List variable to store all model summaries
#modList.surv <- NULL        #List variable to store all models summaries
paramsMM.grwth <- NULL
seMM.grwth <- NULL    #Use SE in error bars in GLMM vs JAGS plots 
paramsMM.surv <- NULL
seMM.surv <- NULL
n.datset <- 10


for (dd in 1:n.datset) {
  
  ## Read in data
  noMiss <- readRDS(file=paste(dateGLM, "_erbr_", name, "NoMiss.srvCor.sdlgCor.grwthCor.",dd,".4GLM",".rds", sep=""))
  
  #lag the climate vars:
  
  nn=length(noMiss$TempFall)
  
  noMiss$TempFall.lag	   =c(noMiss$TempFall[2:(nn)], NA)
  noMiss$TempSummer.lag	 =c(noMiss$TempSummer[2:(nn)], NA)
  noMiss$TempWinter.lag	 =c(noMiss$TempWinter[2:(nn)], NA)
  noMiss$PptFall.lag	   =c(noMiss$PptFall[2:(nn)], NA)
  noMiss$PptWinter.lag	 =c(noMiss$PptWinter[2:(nn)], NA)
  noMiss$PptSummer.lag	 =c(noMiss$PptSummer[2:(nn)], NA)

  
  
  ## Add t+1 climate, sz, & tag into erbr data 
  noMiss <- noMiss %>% mutate(TagNew1=lead(TagNew), RosNew1=lead(RosNew), Surv1=lead(surv))  
  #Remove lines with mis-matched individuals 
  noMiss <- noMiss[which(noMiss$TagNew == noMiss$TagNew1),]  
  
  ## Note: Log size in all models to match JAGS; this makes ending size a linear function of starting size
  
  plot(noMiss$RosNew,noMiss$RosNew1)
  print(t(table(noMiss$Year)))
  
#lagged models------------------------------------------------
  glmm.grwth <- glmmTMB(RosNew1 ~ log(RosNew) + TempFall.lag + TempSummer.lag + TempWinter.lag +     PptFall.lag + PptSummer.lag + PptWinter.lag + (1|TransectNew), family=truncated_nbinom2(link = "log"), data=noMiss) 
  
  glmm.surv <- glmer(Surv1 ~ log(RosNew) + PptWinter.lag + TempFall.lag + TempSummer.lag + 
     TempWinter.lag + (1|TransectNew), family=binomial(link='logit'), data=noMiss)
  
  summary(glmm.grwth)
  summary(glmm.surv)
  
  #lagged models------------------------------------------------
  
  
  
 #  
 #  
 # jagmod <- as.data.frame(readRDS(file="20241001_erbr_JAGSmodBestSUMM_SimDat20yrNoMiss.srvCor.sdlgCor.grwthCor.1.rds"))
 # jagmod <- as.data.frame(readRDS(file="20241002_erbr_JAGSmodBestSUMM_SimDat20yrNoMiss.srvCor.sdlgCor.grwthCor.2.rds"))
 #  jagsparams=jagmod$Median
 #  
 #  ## Growth (param order is same as for JAGS output)
 #  glmm.grwth <- glmmTMB(RosNew1 ~ log(RosNew) + TempFall + TempSummer + TempWinter + 
 #                PptFall + PptSummer + PptWinter + (1|TransectNew), family=truncated_nbinom2(link = "log"), data=noMiss)
 #  
 #  glmm.grwth2 <- glmmTMB(RosNew1 ~ log(RosNew) + PptWinter+ TempFall + TempSummer + TempWinter + 
 #                          PptFall + PptSummer  + (1|TransectNew), family=truncated_nbinom2(link = "log"), data=noMiss)
 #  
 #  summary(glmm.grwth)
 #  summary(glmm.grwth2)
 #  
 #  
 #  ## Survival  (param order is same as for JAGS output)
 #  glmm.surv <- glmer(Surv1 ~ log(RosNew) + PptWinter + TempFall + TempSummer + 
 #                       TempWinter + (1|TransectNew), family=binomial(link='logit'), data=noMiss)
  
  
  ## Save model summaries in list
  #modList.grwth[[length(modList.grwth) + 1]] <- summary(glmm.grwth)
  #modList.surv[[length(modList.surv) + 1]] <- summary(glmm.surv)
  
  
  ## Extract parameter estimates, including intercept, and SEs from GLMMs
  ## Growth
  paramsMM.grwthTmp <- fixef(glmm.grwth)
  seMM.grwthTmp <- as.data.frame(summary.sdreport(glmm.grwth$sdr))$`Std. Error`[1:8]
  
  paramsMM.grwth <- as.data.frame(rbind(paramsMM.grwth, paramsMM.grwthTmp$cond))
  seMM.grwth <- as.data.frame(rbind(seMM.grwth, seMM.grwthTmp))
  
  
  ## Survival
  paramsMM.survTmp <- as.data.frame(summary(glmm.surv)$coefficients[1:nrow(summary(glmm.surv)$coefficients),1])
  seMM.survTmp <- as.data.frame(summary(glmm.surv)$coefficients[1:nrow(summary(glmm.surv)$coefficients),2])
  
  colnames(paramsMM.survTmp) <- paste("GLMM.",dd,sep="")
  colnames(seMM.survTmp) <- paste("SE.",dd,sep="")
  paramsMM.surv <- as.data.frame(c(paramsMM.surv, paramsMM.survTmp))
  seMM.surv <- as.data.frame(c(seMM.surv, seMM.survTmp))

}


paramsMM.grwth
paramsMM.grwthT <- as.data.frame(t(paramsMM.grwth))
colnames(paramsMM.grwthT) <- c("GLMM.1", "GLMM.2","GLMM.3","GLMM.4","GLMM.5","GLMM.6","GLMM.7","GLMM.8","GLMM.9","GLMM.10")
seMM.grwthT <- as.data.frame(t(seMM.grwth))
colnames(seMM.grwthT) <- c("GLMM.1", "GLMM.2","GLMM.3","GLMM.4","GLMM.5","GLMM.6","GLMM.7","GLMM.8","GLMM.9","GLMM.10")

paramsMM.surv
seMM.surv
names.paramTitles <- c("Grwth Intercept","Grwth Size","Grwth Fall Temp","Grwth Summer Temp","Grwth Winter Temp",
                       "Grwth Fall Precip","Grwth Summer Precip","Grwth Winter Precip","Surv Intercept","Surv Size",
                       "Surv Winter Precip","Surv Fall Temp","Surv Summer Temp","Surv Winter Temp")


paramsMM.grwthT$ParamTitle <- names.paramTitles[1:8]
seMM.grwthT$ParamTitle <- names.paramTitles[1:8]
paramsMM.surv$ParamTitle <- names.paramTitles[9:14]
seMM.surv$ParamTitle <- names.paramTitles[9:14]


## Save GLMM param and se results 
date <- as.character("20241003")
name

## Save GLMM model summaries
#saveRDS(modList.grwth, file=paste(date, "_erbr_GLMMtmbSummGrwth_", name, ".rds", sep=""))
#saveRDS(modList.surv, file=paste(date, "_erbr_GLMMsummSurv_", name, ".rds", sep=""))


saveRDS(paramsMM.grwthT, file=paste(date, "_erbr_paramMMgrwth_", name, ".rds", sep=""))
saveRDS(seMM.grwthT, file=paste(date, "_erbr_seMMgrwth_", name, ".rds", sep=""))
saveRDS(paramsMM.surv, file=paste(date, "_erbr_paramMMsurv_", name, ".rds", sep=""))
saveRDS(seMM.surv, file=paste(date, "_erbr_seMMsurv_", name, ".rds", sep=""))
## --------------------------------------
