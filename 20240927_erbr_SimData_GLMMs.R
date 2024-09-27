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
#library(corrplot)
#library(rgl)
#library(viridis) 
library(car)
#library(coda)
library(lme4)
library(resample)
#library(gplots)
library(stringr)
library(plotrix)
library(GLMMadaptive)



## ASSIGN NAME VARIABLE FOR DESIRED DATASETS 
name <- as.character("SimDat40yr") #"SimDat20yrMedGr" "SimDat20yr" "SimDat20yrHiGr"
## ------------------------------------------------------------------------------------------------






## OBTAIN GLMM ESTIMATES ----------------------------------------------------------------

## Loop over datasets 
dateGLM <- as.character("20240926")

modList.grwth <- NULL        #List variable to store all models
modList.surv <- NULL        #List variable to store all models
paramsMM.grwth <- NULL
seMM.grwth <- NULL    #Use SE in error bars in GLMM vs JAGS plots 
paramsMM.surv <- NULL
seMM.surv <- NULL
n.datset <- 10


for (dd in 1:n.datset) {
  
  ## Read in data
  noMiss <- readRDS(file=paste(dateGLM, "_erbr_", name, "NoMiss.srvCor.sdlgCor.",dd,".4GLM",".rds", sep=""))
  
  ## Add t+1 climate, sz, & tag into erbr data 
  noMiss <- noMiss %>% mutate(TagNew1=lead(TagNew), RosNew1=lead(RosNew), Surv1=lead(surv))  
  #Remove lines with mis-matched individuals 
  noMiss <- noMiss[which(noMiss$TagNew == noMiss$TagNew1),]  
  
  ## Note: Log size in all models to match JAGS; this makes ending size a linear function of starting size
  
  ## Growth
  ## Make re-scaled climate variables
  noMiss$TempFallSc <- noMiss$TempFall / max(noMiss$TempFall)
  noMiss$TempSummerSc <- noMiss$TempSummer / max(noMiss$TempSummer)
  noMiss$TempWinterSc <- noMiss$TempWinter / max(noMiss$TempWinter)
  noMiss$PptFallSc <- noMiss$PptFall / max(noMiss$PptFall)
  noMiss$PptSummerSc <- noMiss$PptSummer / max(noMiss$PptSummer)
  noMiss$PptWinterSc <- noMiss$PptWinter / max(noMiss$PptWinter)
  
  glmm.grwth <- glmer.nb(RosNew1 ~ log(RosNew) + TempFallSc + TempSummerSc + 
                           TempWinterSc + PptFallSc + PptSummerSc + PptWinterSc + (1|TransectNew), data=noMiss,
                           control=glmerControl(optCtrl=list(maxfun=50000)))
  
  glmm.grwth <- mixed_model(fixed=RosNew1 ~ log(RosNew) + TempFallSc + TempSummerSc + TempWinterSc + 
                              PptFallSc + PptSummerSc + PptWinterSc, random=~1|TransectNew, data=noMiss,
                            family=zi.negative.binomial(), zi_fixed=~1)#, zi_random=~1)
  
  
  
  ## Survival  (param order is same as for JAGS output)
  glmm.surv <- glmer(Surv1 ~ log(RosNew) + PptWinter + TempFall + TempSummer + 
                       TempWinter + (1|TransectNew), family=binomial(link='logit'), data=noMiss)
  
  ## Save model summaries in list
  modList.grwth[[length(modList.grwth) + 1]] <- summary(glmm.grwth)
  modList.surv[[length(modList.surv) + 1]] <- summary(glmm.surv)
  
  
  
  ## Extract parameter estimates, including intercept, and SEs from GLMMs
  paramsMM.grwthTmp <- as.data.frame(summary(glmm.grwth)$coefficients[1:nrow(summary(glmm.grwth)$coefficients),1])
  seMM.grwthTmp <- as.data.frame(summary(glmm.grwth)$coefficients[1:nrow(summary(glmm.grwth)$coefficients),2])
  
 colnames(paramsMM.grwthTmp) <- paste("GLMM.",dd,sep="")
 colnames(seMM.grwthTmp) <- paste("SE.",dd,sep="")
  
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
paramsMM.grwth
paramsMM.survTmp
names.paramTitles <- c("Grwth Intercept","Grwth Size","Grwth Fall Temp","Grwth Summer Temp","Grwth Winter Temp",
                       "Grwth Fall Precip","Grwth Summer Precip","Grwth Winter Precip","Surv Intercept","Surv Size",
                       "Surv Winter Precip","Surv Fall Temp","Surv Summer Temp","Surv Winter Temp")

paramsMM.grwth$ParamTitle <- names.paramTitles[1:8]
seMM.grwth$ParamTitle <- names.paramTitles[1:8]
paramsMM.surv
paramsMM.surv$ParamTitle <- names.paramTitles[9:14]
seMM.surv$ParamTitle <- names.paramTitles[9:14]


## Save GLMM param and se results 
date <- as.character("20240927")
name

## Save GLMM model summaries
saveRDS(modList.grwth, file=paste(date, "_erbr_GLMMsummGrwth_", name, ".rds", sep=""))
saveRDS(modList.surv, file=paste(date, "_erbr_GLMMsummSurv_", name, ".rds", sep=""))

saveRDS(paramsMM.grwth, file=paste(date, "_erbr_paramMMgrwthWint_", name, ".rds", sep=""))
saveRDS(seMM.grwth, file=paste(date, "_erbr_seMMgrwthWint_", name, ".rds", sep=""))
saveRDS(paramsMM.surv, file=paste(date, "_erbr_paramMMsurvWint_", name, ".rds", sep=""))
saveRDS(seMM.surv, file=paste(date, "_erbr_seMMsurvWint_", name, ".rds", sep=""))
## --------------------------------------






